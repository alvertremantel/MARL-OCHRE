use marl_config::stoich::{
    StoichEventKind, StoichRecord, StoichReservoir, StoichStage, StoichTickLedger, external_delta,
};
use marl_config::*;
use rayon::prelude::*;

const EXPLICIT_DIFFUSION_CFL_LIMIT: f32 = 1.0 / 6.0;
const MAX_DIFFUSION_SUBSTEPS: usize = 100_000;

pub fn validate_diffusion_config(sim: &SimulationConfig) -> Result<(), String> {
    if !sim.dt.is_finite() || sim.dt < 0.0 {
        return Err(format!("dt must be finite and nonnegative, got {}", sim.dt));
    }
    if !sim.alpha_eps.is_finite() || sim.alpha_eps < 0.0 {
        return Err(format!(
            "alpha_eps must be finite and nonnegative, got {}",
            sim.alpha_eps
        ));
    }
    if !sim.k_eps.is_finite() || sim.k_eps <= 0.0 {
        return Err(format!(
            "k_eps must be finite and positive, got {}",
            sim.k_eps
        ));
    }
    for (index, &d) in sim.d_voxel.iter().enumerate() {
        if !d.is_finite() || d < 0.0 {
            return Err(format!(
                "d_voxel[{index}] must be finite and nonnegative, got {d}"
            ));
        }
    }
    for (index, &decay) in sim.lambda_decay.iter().enumerate() {
        if !decay.is_finite() || decay < 0.0 {
            return Err(format!(
                "lambda_decay[{index}] must be finite and nonnegative, got {decay}"
            ));
        }
    }

    let required = required_stable_substeps(sim);
    if required > MAX_DIFFUSION_SUBSTEPS {
        return Err(format!(
            "diffusion requires {required} substeps, above maximum {MAX_DIFFUSION_SUBSTEPS}"
        ));
    }

    Ok(())
}

/// Minimum number of explicit substeps needed for the configured diffusion and
/// decay rates.
///
/// Diffusion coefficients are already expressed in voxel units, so the 3D
/// forward-Euler stability bound is `dt_sub * max(D) <= 1/6`. Explicit decay
/// also needs `dt_sub * max(lambda_decay) <= 1` to avoid clamping through zero.
pub fn stable_diffusion_substeps(sim: &SimulationConfig) -> usize {
    if sim.diffusion_substeps == 0 {
        return 0;
    }

    let required = required_stable_substeps(sim).min(MAX_DIFFUSION_SUBSTEPS);
    sim.diffusion_substeps.max(required.max(1))
}

fn required_stable_substeps(sim: &SimulationConfig) -> usize {
    let max_d = sim
        .d_voxel
        .iter()
        .copied()
        .filter(|d| d.is_finite())
        .fold(0.0f32, f32::max)
        .max(0.0);
    let max_decay = sim
        .lambda_decay
        .iter()
        .copied()
        .filter(|decay| decay.is_finite())
        .fold(0.0f32, f32::max)
        .max(0.0);
    if (max_d == 0.0 && max_decay == 0.0) || sim.dt <= 0.0 || !sim.dt.is_finite() {
        return 1;
    }

    let diffusion_required =
        (sim.dt as f64 * max_d as f64 / EXPLICIT_DIFFUSION_CFL_LIMIT as f64).ceil();
    let decay_required = (sim.dt as f64 * max_decay as f64).ceil();
    let required = diffusion_required.max(decay_required);
    if !required.is_finite() || required > usize::MAX as f64 {
        usize::MAX
    } else {
        required as usize
    }
}

/// 3D chemical concentration field.
///
/// Stores concentrations of all S_EXT chemical species at every voxel in
/// a flat array laid out as [z][y][x][species]. This memory order gives
/// cache-friendly access when iterating z-columns (the dominant access
/// pattern for diffusion and light attenuation).
///
/// A pre-allocated scratch buffer (`scratch`) eliminates per-substep
/// allocation during diffusion. The two buffers are swapped each substep
/// so no copying is ever needed.
#[derive(Clone)]
pub struct Field {
    grid: GridDims,
    /// Concentration data: grid_z * grid_y * grid_x * S_EXT floats
    pub data: Vec<f32>,
    /// Double-buffer for diffusion solver — same size as `data`.
    /// Swapped with `data` each substep to avoid allocation in the hot loop.
    scratch: Vec<f32>,
}

impl Field {
    pub fn new(grid: GridDims) -> Self {
        grid.validate().expect("invalid grid dimensions");
        let n = grid
            .field_float_count()
            .expect("validated grid field length must fit usize");
        Self {
            grid,
            data: vec![0.0; n],
            scratch: vec![0.0; n],
        }
    }

    pub fn new_default() -> Self {
        Self::new(GridDims::default())
    }

    #[inline]
    pub fn grid(&self) -> GridDims {
        self.grid
    }

    #[inline]
    pub fn grid_x(&self) -> usize {
        self.grid.x
    }

    #[inline]
    pub fn grid_y(&self) -> usize {
        self.grid.y
    }

    #[inline]
    pub fn grid_z(&self) -> usize {
        self.grid.z
    }

    #[inline]
    pub fn voxel_count(&self) -> usize {
        self.grid
            .voxel_count()
            .expect("validated grid voxel count must fit usize")
    }

    #[inline]
    pub fn field_len(&self) -> usize {
        self.grid
            .field_float_count()
            .expect("validated grid field length must fit usize")
    }

    #[inline]
    fn idx(&self, x: usize, y: usize, z: usize, s: usize) -> usize {
        ((z * self.grid.y + y) * self.grid.x + x) * S_EXT + s
    }

    #[inline]
    pub fn get(&self, x: usize, y: usize, z: usize, s: usize) -> f32 {
        self.data[self.idx(x, y, z, s)]
    }

    #[inline]
    pub fn set(&mut self, x: usize, y: usize, z: usize, s: usize, val: f32) {
        let i = self.idx(x, y, z, s);
        self.data[i] = val;
    }

    /// Read all species at a voxel
    pub fn read_voxel(&self, x: usize, y: usize, z: usize) -> [f32; S_EXT] {
        let mut out = [0.0f32; S_EXT];
        let base = self.idx(x, y, z, 0);
        out.copy_from_slice(&self.data[base..base + S_EXT]);
        out
    }

    pub fn species_totals(&self) -> [f32; S_EXT] {
        let mut totals = [0.0f32; S_EXT];
        for voxel in self.data.chunks_exact(S_EXT) {
            for s in 0..S_EXT {
                totals[s] += voxel[s].max(0.0);
            }
        }
        totals
    }

    /// Apply cell secretion/consumption deltas to a voxel
    pub fn apply_deltas(&mut self, x: usize, y: usize, z: usize, deltas: &[f32; S_EXT]) {
        let base = self.idx(x, y, z, 0);
        for (s, delta) in deltas.iter().enumerate() {
            self.data[base + s] = (self.data[base + s] + *delta).max(0.0);
        }
    }

    /// Helper: compute flat index without &self (needed for parallel closures).
    #[inline]
    fn idx_static(grid: GridDims, x: usize, y: usize, z: usize, s: usize) -> usize {
        ((z * grid.y + y) * grid.x + x) * S_EXT + s
    }

    /// Read a concentration from a raw data slice (used in parallel diffusion).
    #[inline]
    fn get_from(src: &[f32], grid: GridDims, x: usize, y: usize, z: usize, s: usize) -> f32 {
        src[Self::idx_static(grid, x, y, z, s)]
    }

    #[inline]
    fn niche_factor(structural: f32, sim: &SimulationConfig) -> f32 {
        let structural = structural.max(0.0);
        let denom = (sim.k_eps + structural).max(f32::EPSILON);
        (1.0 - sim.alpha_eps.max(0.0) * structural / denom).clamp(0.0, 1.0)
    }

    /// Run one diffusion substep for all species, parallelized over z-layers.
    ///
    /// Uses forward Euler integration of the 3D discrete Laplacian with
    /// Neumann (zero-flux) boundary conditions. Each voxel's new
    /// concentration is:
    ///
    ///   c' = c + dt * (D_local * laplacian(c) - lambda * c)
    ///
    /// where D_local is reduced by:
    ///   1. Niche construction (structural EPS deposits slow diffusion,
    ///      mimicking biofilm matrix)
    ///   2. Cell body exclusion — occupied voxels are **completely skipped**.
    ///      In a real microbial column, the extracellular medium (water/gel
    ///      between cells) is where diffusion happens. Cells are physical
    ///      objects that exclude the liquid phase. A voxel occupied by a
    ///      cell has no free water for chemicals to diffuse through.
    ///
    /// Occupied neighbors are treated as Neumann boundaries (zero-flux),
    /// identical to wall boundaries. This means chemicals cannot diffuse
    /// into, out of, or through cell-occupied space. Cells access the
    /// external medium through their transport machinery (see cell.rs),
    /// reading from adjacent empty voxels.
    ///
    /// This is the key mechanism that prevents grid saturation: interior
    /// cells in a dense colony are cut off from nutrients because the
    /// surrounding occupied voxels block diffusion. Only cells on the
    /// colony surface have access to the medium.
    ///
    /// The z-loop is embarrassingly parallel: each layer reads from the
    /// source buffer (immutable) and writes to its own slice of the
    /// scratch buffer. Rayon splits this across all available CPU cores.
    fn diffusion_step_inner(
        &mut self,
        dt_sub: f32,
        occupancy: Option<&[bool]>,
        sim: &SimulationConfig,
    ) {
        let grid = self.grid;
        let src = &self.data;
        let layer_size = grid.y * grid.x * S_EXT;
        let voxel_layer_size = grid.y * grid.x;

        self.scratch
            .par_chunks_mut(layer_size)
            .enumerate()
            .for_each(|(z, dst_layer)| {
                for y in 0..grid.y {
                    for x in 0..grid.x {
                        let occ_here =
                            occupancy.is_some_and(|o| o[z * voxel_layer_size + y * grid.x + x]);

                        // Occupied voxels are excluded from diffusion entirely.
                        // Their field concentrations are meaningless (chemicals
                        // inside cells are tracked in cell.internal[], not here).
                        // Just copy unchanged to maintain buffer consistency.
                        if occ_here {
                            for s in 0..S_EXT {
                                let local_idx = (y * grid.x + x) * S_EXT + s;
                                dst_layer[local_idx] = Self::get_from(src, grid, x, y, z, s);
                            }
                            continue;
                        }

                        // --- Empty voxel: compute conservative face fluxes ---

                        let base = (y * grid.x + x) * S_EXT;
                        let center_idx = Self::idx_static(grid, x, y, z, 0);
                        let niche_here = Self::niche_factor(src[center_idx + 7], sim);

                        // Check which neighbors are occupied or walls once per voxel.
                        // Walls and occupied neighbors are zero-flux faces.
                        let occ_check = |nx: usize, ny: usize, nz: usize| -> bool {
                            occupancy.is_some_and(|o| o[nz * voxel_layer_size + ny * grid.x + nx])
                        };
                        let neighbor = |nx: usize, ny: usize, nz: usize| {
                            if occ_check(nx, ny, nz) {
                                None
                            } else {
                                let idx = Self::idx_static(grid, nx, ny, nz, 0);
                                Some((idx, Self::niche_factor(src[idx + 7], sim)))
                            }
                        };

                        let xm = (x > 0).then(|| neighbor(x - 1, y, z)).flatten();
                        let xp = (x + 1 < grid.x).then(|| neighbor(x + 1, y, z)).flatten();
                        let ym = (y > 0).then(|| neighbor(x, y - 1, z)).flatten();
                        let yp = (y + 1 < grid.y).then(|| neighbor(x, y + 1, z)).flatten();
                        let zm = (z > 0).then(|| neighbor(x, y, z - 1)).flatten();
                        let zp = (z + 1 < grid.z).then(|| neighbor(x, y, z + 1)).flatten();
                        let neighbors = [xm, xp, ym, yp, zm, zp];

                        for s in 0..S_EXT {
                            let c = src[center_idx + s];
                            let base_d = sim.d_voxel[s].max(0.0);
                            let diffusion = neighbors
                                .iter()
                                .flatten()
                                .map(|&(neighbor_idx, niche_neighbor)| {
                                    let d_face = 0.5 * base_d * (niche_here + niche_neighbor);
                                    d_face * (src[neighbor_idx + s] - c)
                                })
                                .sum::<f32>();
                            let decay = sim.lambda_decay[s] * c;
                            let new_c = c + dt_sub * (diffusion - decay);

                            dst_layer[base + s] = new_c.max(0.0);
                        }
                    }
                }
            });

        // Swap buffers — the scratch buffer becomes the live data,
        // and the old data buffer becomes scratch for the next substep.
        std::mem::swap(&mut self.data, &mut self.scratch);
    }

    /// Run a full tick of diffusion with cell-body exclusion.
    ///
    /// Occupied voxels are completely excluded from diffusion — chemicals
    /// only move through the extracellular medium (empty voxels). This is
    /// the physically correct model: cells are solid objects that displace
    /// the liquid phase. Interior cells in a dense colony are cut off
    /// from nutrients, creating natural carrying capacity.
    pub fn diffuse_tick_with_cells(&mut self, occupancy: &[bool], sim: &SimulationConfig) {
        validate_diffusion_config(sim).expect("invalid diffusion configuration");
        assert_eq!(
            occupancy.len(),
            self.voxel_count(),
            "occupancy length must match field voxel count"
        );
        let substeps = stable_diffusion_substeps(sim);
        if substeps == 0 {
            return;
        }

        let dt_sub = sim.dt / substeps as f32;
        for _ in 0..substeps {
            self.diffusion_step_inner(dt_sub, Some(occupancy), sim);
        }
    }

    #[allow(dead_code)] // TODO: occupancy-free version kept for testing/benchmarking
    /// Run a full tick of diffusion (multiple substeps for stability)
    pub fn diffuse_tick(&mut self, sim: &SimulationConfig) {
        validate_diffusion_config(sim).expect("invalid diffusion configuration");
        let substeps = stable_diffusion_substeps(sim);
        if substeps == 0 {
            return;
        }

        let dt_sub = sim.dt / substeps as f32;
        for _ in 0..substeps {
            self.diffusion_step_inner(dt_sub, None, sim);
        }
    }

    /// Set boundary source terms (called once per tick before diffusion).
    /// Oxidant + carbon sourced from top (z=0), reductant from bottom (z=max).
    /// These are the only external inputs to the system — everything else is recycled.
    pub fn apply_boundary_sources(&mut self, sim: &SimulationConfig) {
        self.apply_boundary_sources_with_stoich(sim, None, false);
    }

    pub fn apply_boundary_sources_with_stoich(
        &mut self,
        sim: &SimulationConfig,
        mut stoich: Option<&mut StoichTickLedger>,
        keep_events: bool,
    ) {
        // Top face: oxidant (species 1) and carbon (species 3)
        for y in 0..self.grid.y {
            for x in 0..self.grid.x {
                let ox = self.get(x, y, 0, 1);
                let new_ox = (ox + sim.source_rate_oxidant).min(sim.c_max);
                self.set(x, y, 0, 1, new_ox);
                if let Some(ledger) = stoich.as_deref_mut() {
                    record_boundary_source(ledger, 1, new_ox - ox, keep_events);
                }
                let ca = self.get(x, y, 0, 3);
                let new_ca = (ca + sim.source_rate_carbon).min(sim.c_max);
                self.set(x, y, 0, 3, new_ca);
                if let Some(ledger) = stoich.as_deref_mut() {
                    record_boundary_source(ledger, 3, new_ca - ca, keep_events);
                }
            }
        }

        // Bottom face: reductant (species 2)
        for y in 0..self.grid.y {
            for x in 0..self.grid.x {
                let z = self.grid.z - 1;
                let re = self.get(x, y, z, 2);
                let new_re = (re + sim.source_rate_reductant).min(sim.c_max);
                self.set(x, y, z, 2, new_re);
                if let Some(ledger) = stoich.as_deref_mut() {
                    record_boundary_source(ledger, 2, new_re - re, keep_events);
                }
            }
        }
    }
}

fn record_boundary_source(
    ledger: &mut StoichTickLedger,
    species: usize,
    amount: f32,
    keep_events: bool,
) {
    if amount <= 0.0 {
        return;
    }
    ledger.record(
        StoichRecord::new(
            StoichStage::BoundarySources,
            StoichEventKind::BoundarySource,
            amount,
        )
        .model_delta(external_delta(species, amount))
        .balancing_reservoir(StoichReservoir::BoundaryInput)
        .species(species),
        keep_events,
    );
}

impl Default for Field {
    fn default() -> Self {
        Self::new_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init_deterministic_field(field: &mut Field) {
        for z in 0..field.grid_z() {
            for y in 0..field.grid_y() {
                for x in 0..field.grid_x() {
                    for s in 0..S_EXT {
                        let value = ((x * 13 + y * 17 + z * 19 + s * 23) % 541) as f32 * 0.001;
                        field.set(x, y, z, s, value);
                    }
                }
            }
        }
    }

    fn conservative_sim(species: usize, d: f32) -> SimulationConfig {
        let mut d_voxel = [0.0; S_EXT];
        d_voxel[species] = d;
        SimulationConfig {
            diffusion_substeps: 1,
            d_voxel,
            lambda_decay: [0.0; S_EXT],
            ..Default::default()
        }
    }

    fn total_species(field: &Field, species: usize) -> f64 {
        field
            .data
            .chunks_exact(S_EXT)
            .map(|voxel| voxel[species] as f64)
            .sum()
    }

    fn test_grid() -> GridDims {
        GridDims { x: 5, y: 4, z: 3 }
    }

    #[test]
    fn default_field_uses_default_grid_dimensions() {
        let field = Field::default();

        assert_eq!(field.grid(), GridDims::default());
        assert_eq!(field.grid_x(), GRID_X);
        assert_eq!(field.grid_y(), GRID_Y);
        assert_eq!(field.grid_z(), GRID_Z);
        assert_eq!(field.voxel_count(), GRID_X * GRID_Y * GRID_Z);
        assert_eq!(field.field_len(), GRID_X * GRID_Y * GRID_Z * S_EXT);
        assert_eq!(field.data.len(), field.field_len());
    }

    #[test]
    fn field_allocates_from_runtime_grid_dimensions() {
        let grid = test_grid();
        let field = Field::new(grid);

        assert_eq!(field.grid(), grid);
        assert_eq!(field.grid_x(), 5);
        assert_eq!(field.grid_y(), 4);
        assert_eq!(field.grid_z(), 3);
        assert_eq!(field.voxel_count(), 60);
        assert_eq!(field.field_len(), 60 * S_EXT);
        assert_eq!(field.data.len(), 60 * S_EXT);
    }

    #[test]
    fn stable_substeps_enforce_explicit_diffusion_bound() {
        let mut d_voxel = [0.0; S_EXT];
        d_voxel[1] = 1.5;
        let sim = SimulationConfig {
            diffusion_substeps: 1,
            dt: 1.0,
            d_voxel,
            ..Default::default()
        };

        assert_eq!(stable_diffusion_substeps(&sim), 9);
    }

    #[test]
    fn stable_substeps_enforce_explicit_decay_bound() {
        let mut d_voxel = [0.0; S_EXT];
        d_voxel[1] = 0.1;
        let mut lambda_decay = [0.0; S_EXT];
        lambda_decay[1] = 12.0;
        let sim = SimulationConfig {
            diffusion_substeps: 1,
            dt: 1.0,
            d_voxel,
            lambda_decay,
            ..Default::default()
        };

        assert_eq!(stable_diffusion_substeps(&sim), 12);
    }

    #[test]
    fn invalid_diffusion_coefficients_are_rejected() {
        let mut d_voxel = SimulationConfig::default().d_voxel;
        d_voxel[1] = f32::NAN;
        let sim = SimulationConfig {
            d_voxel,
            ..Default::default()
        };
        assert!(validate_diffusion_config(&sim).is_err());

        let mut lambda_decay = SimulationConfig::default().lambda_decay;
        lambda_decay[1] = -0.1;
        let sim = SimulationConfig {
            lambda_decay,
            ..Default::default()
        };
        assert!(validate_diffusion_config(&sim).is_err());

        let sim = SimulationConfig {
            k_eps: 0.0,
            ..Default::default()
        };
        assert!(validate_diffusion_config(&sim).is_err());
    }

    #[test]
    fn zero_diffusion_substeps_is_noop() {
        let sim = SimulationConfig {
            diffusion_substeps: 0,
            ..Default::default()
        };

        let mut field = Field::new(test_grid());
        init_deterministic_field(&mut field);
        let before = field.data.clone();
        let occupancy = vec![false; field.voxel_count()];

        field.diffuse_tick_with_cells(&occupancy, &sim);

        assert_eq!(field.data, before);
    }

    #[test]
    #[should_panic(expected = "occupancy length must match field voxel count")]
    fn occupancy_length_must_match_runtime_grid() {
        let sim = SimulationConfig {
            diffusion_substeps: 0,
            ..Default::default()
        };
        let mut field = Field::new(test_grid());
        let occupancy = vec![false; field.voxel_count() - 1];

        field.diffuse_tick_with_cells(&occupancy, &sim);
    }

    #[test]
    fn occupied_voxel_is_copied_unchanged() {
        let sim = SimulationConfig {
            diffusion_substeps: 1,
            ..Default::default()
        };

        let mut field = Field::new(test_grid());
        init_deterministic_field(&mut field);

        let x = field.grid_x() / 2;
        let y = field.grid_y() / 2;
        let z = field.grid_z() / 2;
        let before = field.read_voxel(x, y, z);

        let mut occupancy = vec![false; field.voxel_count()];
        occupancy[z * field.grid_y() * field.grid_x() + y * field.grid_x() + x] = true;

        field.diffuse_tick_with_cells(&occupancy, &sim);

        assert_eq!(field.read_voxel(x, y, z), before);
    }

    #[test]
    fn occupied_neighbor_is_zero_flux_barrier() {
        let sim = conservative_sim(1, 1.0);
        let mut field = Field::new(test_grid());
        let x = field.grid_x() / 2;
        let y = field.grid_y() / 2;
        let z = field.grid_z() / 2;
        field.set(x + 1, y, z, 1, 5.0);

        let mut occupancy = vec![false; field.voxel_count()];
        occupancy[z * field.grid_y() * field.grid_x() + y * field.grid_x() + x] = true;

        let before = total_species(&field, 1);
        field.diffuse_tick_with_cells(&occupancy, &sim);
        let after = total_species(&field, 1);

        assert_eq!(field.get(x, y, z, 1), 0.0);
        assert!(
            (after - before).abs() <= 1e-5,
            "mass changed across occupied barrier: before={before}, after={after}"
        );
    }

    #[test]
    fn heterogeneous_eps_diffusion_conserves_mass_without_decay() {
        let sim = conservative_sim(1, 1.0);

        let mut field = Field::new(GridDims { x: 6, y: 4, z: 3 });
        for z in 0..field.grid_z() {
            for y in 0..field.grid_y() {
                for x in 0..field.grid_x() {
                    let structural = if x < field.grid_x() / 2 { 0.0 } else { 10.0 };
                    field.set(x, y, z, 7, structural);
                }
            }
        }
        field.set(
            field.grid_x() / 2 - 1,
            field.grid_y() / 2,
            field.grid_z() / 2,
            1,
            2.0,
        );
        field.set(
            field.grid_x() / 2,
            field.grid_y() / 2,
            field.grid_z() / 2,
            1,
            1.0,
        );

        let before = total_species(&field, 1);
        field.diffuse_tick(&sim);
        let after = total_species(&field, 1);

        assert!(
            (after - before).abs() <= 1e-5,
            "mass changed under conservative diffusion: before={before}, after={after}"
        );
    }

    #[test]
    fn deterministic_diffusion_stays_finite_and_nonnegative() {
        let sim = SimulationConfig {
            diffusion_substeps: 1,
            ..Default::default()
        };

        let mut field = Field::new(test_grid());
        init_deterministic_field(&mut field);
        let occupancy = vec![false; field.voxel_count()];

        field.diffuse_tick_with_cells(&occupancy, &sim);

        for (index, value) in field.data.iter().copied().enumerate() {
            assert!(
                value.is_finite(),
                "non-finite concentration at index {index}"
            );
            assert!(
                value >= 0.0,
                "negative concentration at index {index}: {value}"
            );
        }
    }
}
