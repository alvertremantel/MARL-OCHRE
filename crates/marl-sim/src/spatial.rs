use marl_config::stoich::{
    StoichEventKind, StoichRecord, StoichReservoir, StoichStage, StoichTickLedger, external_delta,
};
use marl_config::{GridDims, S_EXT, SimulationConfig};
use marl_field::field::Field;
use rand::Rng;
use std::collections::{HashMap, HashSet};

/// The 6 face-neighbor offsets in 3D (±x, ±y, ±z).
/// Shared by all neighbor-scanning functions.
const FACE_OFFSETS: [(i32, i32, i32); 6] = [
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
];

fn collect_empty_neighbors(
    grid: GridDims,
    pos: [u16; 3],
    cell_map: &HashMap<[u16; 3], usize>,
    reserved: Option<&HashSet<[u16; 3]>>,
) -> ([[u16; 3]; 6], usize) {
    let mut neighbors = [[0u16; 3]; 6];
    let mut count = 0;
    for &(dx, dy, dz) in &FACE_OFFSETS {
        let nx = pos[0] as i32 + dx;
        let ny = pos[1] as i32 + dy;
        let nz = pos[2] as i32 + dz;
        if nx >= 0
            && nx < grid.x as i32
            && ny >= 0
            && ny < grid.y as i32
            && nz >= 0
            && nz < grid.z as i32
        {
            let p = [nx as u16, ny as u16, nz as u16];
            if !cell_map.contains_key(&p) && reserved.is_none_or(|r| !r.contains(&p)) {
                neighbors[count] = p;
                count += 1;
            }
        }
    }
    (neighbors, count)
}

/// Collect the positions of all in-bounds, empty face-neighbors of a voxel.
pub fn empty_neighbors(
    grid: GridDims,
    pos: [u16; 3],
    cell_map: &HashMap<[u16; 3], usize>,
) -> Vec<[u16; 3]> {
    let (neighbors, count) = collect_empty_neighbors(grid, pos, cell_map, None);
    neighbors[..count].to_vec()
}

/// Read the average chemical environment available to a cell from
/// surrounding extracellular space.
///
/// In a real microbial column, cells access dissolved chemicals from the
/// liquid medium around them — not from inside their own body. A cell
/// surrounded by other cells has limited access to nutrients: only empty
/// (liquid-filled) neighboring voxels contribute.
///
/// Returns the mean concentration across all empty face-neighbors. If
/// the cell is completely enclosed by other cells, returns all zeros —
/// the cell is cut off from external resources and will starve.
pub fn read_neighbor_environment(
    pos: [u16; 3],
    field: &Field,
    cell_map: &HashMap<[u16; 3], usize>,
) -> [f32; S_EXT] {
    let (neighbors, count) = collect_empty_neighbors(field.grid(), pos, cell_map, None);
    if count == 0 {
        return [0.0; S_EXT];
    }
    let mut sum = [0.0f32; S_EXT];
    for npos in &neighbors[..count] {
        let voxel = field.read_voxel(npos[0] as usize, npos[1] as usize, npos[2] as usize);
        for s in 0..S_EXT {
            sum[s] += voxel[s];
        }
    }
    let n = count as f32;
    for value in &mut sum {
        *value /= n;
    }
    sum
}

/// Distribute a cell's secretion/consumption deltas to surrounding empty voxels.
///
/// In the real world, metabolic byproducts are released into the liquid
/// medium around the cell, and consumed substrates are depleted from that
/// same medium. The deltas are split equally among all empty face-neighbors.
///
/// If no empty neighbors exist, deltas are lost — the cell cannot exchange
/// chemicals with a fully packed environment (waste heat / trapped products).
pub fn apply_deltas_to_neighbors(
    pos: [u16; 3],
    field: &mut Field,
    cell_map: &HashMap<[u16; 3], usize>,
    deltas: &[f32; S_EXT],
) -> [f32; S_EXT] {
    apply_deltas_to_neighbors_with_stoich(pos, field, cell_map, deltas, None, false, 0)
}

pub fn apply_deltas_to_neighbors_with_stoich(
    pos: [u16; 3],
    field: &mut Field,
    cell_map: &HashMap<[u16; 3], usize>,
    deltas: &[f32; S_EXT],
    mut stoich: Option<&mut StoichTickLedger>,
    keep_events: bool,
    actor_id: u64,
) -> [f32; S_EXT] {
    let mut accepted = [0.0f32; S_EXT];
    let (neighbors, count) = collect_empty_neighbors(field.grid(), pos, cell_map, None);
    if count == 0 {
        for (species, delta) in deltas.iter().enumerate() {
            if *delta != 0.0
                && let Some(ledger) = stoich.as_deref_mut()
            {
                ledger.record(
                    StoichRecord::new(
                        StoichStage::SpatialExchange,
                        StoichEventKind::ClampLoss,
                        delta.abs(),
                    )
                    .model_delta(external_delta(species, -*delta))
                    .balancing_reservoir(StoichReservoir::ClampLoss)
                    .actor(actor_id)
                    .species(species),
                    keep_events,
                );
            }
        }
        return accepted; // enclosed cell — deltas are lost
    }
    for s in 0..S_EXT {
        let delta = deltas[s];
        if delta == 0.0 {
            continue;
        }

        if delta > 0.0 {
            let mut split_deltas = [0.0f32; S_EXT];
            split_deltas[s] = delta / count as f32;
            for npos in &neighbors[..count] {
                field.apply_deltas(
                    npos[0] as usize,
                    npos[1] as usize,
                    npos[2] as usize,
                    &split_deltas,
                );
            }
            accepted[s] += delta;
            continue;
        }

        // Consumption is distributed in proportion to available mass so
        // heterogeneous neighbor concentrations do not silently clamp away
        // requested uptake in low-concentration voxels.
        let mut available = [0.0f32; 6];
        let mut total_available = 0.0;
        for (i, npos) in neighbors[..count].iter().enumerate() {
            let value = field.get(npos[0] as usize, npos[1] as usize, npos[2] as usize, s);
            available[i] = value.max(0.0);
            total_available += available[i];
        }
        if total_available <= f32::EPSILON {
            continue;
        }

        let requested = (-delta).min(total_available);
        for (i, npos) in neighbors[..count].iter().enumerate() {
            if available[i] <= 0.0 {
                continue;
            }
            let mut split_deltas = [0.0f32; S_EXT];
            split_deltas[s] = -requested * (available[i] / total_available);
            field.apply_deltas(
                npos[0] as usize,
                npos[1] as usize,
                npos[2] as usize,
                &split_deltas,
            );
        }
        accepted[s] -= requested;
    }
    if let Some(ledger) = stoich {
        for s in 0..S_EXT {
            let actual = accepted[s];
            let lost = deltas[s] - actual;
            if lost.abs() > f32::EPSILON {
                ledger.record(
                    StoichRecord::new(
                        StoichStage::SpatialExchange,
                        StoichEventKind::ClampLoss,
                        lost.abs(),
                    )
                    .model_delta(external_delta(s, -lost))
                    .balancing_reservoir(StoichReservoir::ClampLoss)
                    .actor(actor_id)
                    .species(s),
                    keep_events,
                );
            }
        }
    }
    accepted
}

/// Find an empty voxel `division_neighbor_distance` steps away along a face axis.
/// Daughters bud off and drift — they don't accrete like plant tissue.
/// Falls back to adjacent placement if no long-range positions are available.
pub fn find_empty_neighbor(
    grid: GridDims,
    pos: [u16; 3],
    occupied: &HashMap<[u16; 3], usize>,
    rng: &mut impl Rng,
    sim: &SimulationConfig,
) -> Option<[u16; 3]> {
    find_empty_neighbor_avoiding(grid, pos, occupied, None, rng, sim)
}

pub fn find_empty_neighbor_avoiding(
    grid: GridDims,
    pos: [u16; 3],
    occupied: &HashMap<[u16; 3], usize>,
    reserved: Option<&HashSet<[u16; 3]>>,
    rng: &mut impl Rng,
    sim: &SimulationConfig,
) -> Option<[u16; 3]> {
    let dist = sim.division_neighbor_distance as i32;
    // Try distance-first (gap between parent and daughter)
    let mut candidates = [[0u16; 3]; 6];
    let mut count = 0;
    for &(dx, dy, dz) in &FACE_OFFSETS {
        let nx = pos[0] as i32 + dx * dist;
        let ny = pos[1] as i32 + dy * dist;
        let nz = pos[2] as i32 + dz * dist;
        if nx >= 0
            && nx < grid.x as i32
            && ny >= 0
            && ny < grid.y as i32
            && nz >= 0
            && nz < grid.z as i32
        {
            let npos = [nx as u16, ny as u16, nz as u16];
            if !occupied.contains_key(&npos) && reserved.is_none_or(|r| !r.contains(&npos)) {
                candidates[count] = npos;
                count += 1;
            }
        }
    }
    if count > 0 {
        let idx = rng.random_range(0..count);
        return Some(candidates[idx]);
    }

    // Fallback: adjacent placement if surrounded at distance
    let (adjacent, adjacent_count) = collect_empty_neighbors(grid, pos, occupied, reserved);
    if adjacent_count == 0 {
        None
    } else {
        let idx = rng.random_range(0..adjacent_count);
        Some(adjacent[idx])
    }
}

pub fn nearby_cell_indices(
    grid: GridDims,
    pos: [u16; 3],
    cell_map: &HashMap<[u16; 3], usize>,
    radius: u8,
) -> Vec<usize> {
    if radius == 0 {
        return Vec::new();
    }

    let radius = radius as i32;
    let mut found = Vec::new();
    for dz in -radius..=radius {
        for dy in -radius..=radius {
            for dx in -radius..=radius {
                if dx == 0 && dy == 0 && dz == 0 {
                    continue;
                }
                let nx = pos[0] as i32 + dx;
                let ny = pos[1] as i32 + dy;
                let nz = pos[2] as i32 + dz;
                if nx < 0
                    || nx >= grid.x as i32
                    || ny < 0
                    || ny >= grid.y as i32
                    || nz < 0
                    || nz >= grid.z as i32
                {
                    continue;
                }
                let npos = [nx as u16, ny as u16, nz as u16];
                if let Some(&idx) = cell_map.get(&npos) {
                    found.push(idx);
                }
            }
        }
    }
    found
}

pub fn find_cell_neighbor(
    grid: GridDims,
    pos: [u16; 3],
    cell_map: &HashMap<[u16; 3], usize>,
) -> Option<usize> {
    for &(dx, dy, dz) in &FACE_OFFSETS {
        let nx = pos[0] as i32 + dx;
        let ny = pos[1] as i32 + dy;
        let nz = pos[2] as i32 + dz;
        if nx >= 0
            && nx < grid.x as i32
            && ny >= 0
            && ny < grid.y as i32
            && nz >= 0
            && nz < grid.z as i32
        {
            let npos = [nx as u16, ny as u16, nz as u16];
            if let Some(&idx) = cell_map.get(&npos) {
                return Some(idx);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consumption_is_weighted_by_available_neighbor_mass() {
        let grid = GridDims { x: 5, y: 5, z: 5 };
        let mut field = Field::new(grid);
        let cell_map = HashMap::from([([2, 2, 2], 0)]);
        field.set(3, 2, 2, 1, 0.0);
        field.set(1, 2, 2, 1, 1.0);

        let mut deltas = [0.0f32; S_EXT];
        deltas[1] = -0.5;
        apply_deltas_to_neighbors([2, 2, 2], &mut field, &cell_map, &deltas);

        assert_eq!(field.get(3, 2, 2, 1), 0.0);
        assert!((field.get(1, 2, 2, 1) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn enclosed_exchange_records_balanced_clamp_loss_without_field_change() {
        let grid = GridDims { x: 5, y: 5, z: 5 };
        let mut field = Field::new(grid);
        let center = [2, 2, 2];
        let cell_map = HashMap::from([
            (center, 0),
            ([3, 2, 2], 1),
            ([1, 2, 2], 2),
            ([2, 3, 2], 3),
            ([2, 1, 2], 4),
            ([2, 2, 3], 5),
            ([2, 2, 1], 6),
        ]);
        let mut deltas = [0.0f32; S_EXT];
        deltas[4] = 1.25;
        let mut ledger = StoichTickLedger::default();

        let accepted = apply_deltas_to_neighbors_with_stoich(
            center,
            &mut field,
            &cell_map,
            &deltas,
            Some(&mut ledger),
            true,
            42,
        );

        assert_eq!(accepted[4], 0.0);
        assert_eq!(field.get(3, 2, 2, 4), 0.0);
        assert_eq!(ledger.events.len(), 1);
        let event = &ledger.events[0];
        assert_eq!(event.kind, StoichEventKind::ClampLoss);
        assert!(event.balanced);
        assert!(
            event
                .residual
                .is_near_zero(marl_config::stoich::STOICH_TOLERANCE)
        );
        assert!(event.reservoir_delta.total_abs_sum() > 0.0);
    }

    #[test]
    fn empty_neighbor_scan_uses_face_neighbors_only() {
        let grid = GridDims { x: 4, y: 4, z: 4 };
        let cell_map = HashMap::from([([1, 1, 1], 0), ([2, 1, 1], 1)]);
        let neighbors = empty_neighbors(grid, [1, 1, 1], &cell_map);

        assert_eq!(neighbors.len(), 5);
        assert!(!neighbors.contains(&[2, 1, 1]));
        assert!(!neighbors.contains(&[2, 2, 1]));
    }

    #[test]
    fn nearby_cell_scan_uses_bounded_local_radius() {
        let grid = GridDims { x: 5, y: 5, z: 5 };
        let cell_map = HashMap::from([
            ([2, 2, 2], 0),
            ([3, 2, 2], 1),
            ([3, 3, 2], 2),
            ([4, 2, 2], 3),
        ]);

        let mut near = nearby_cell_indices(grid, [2, 2, 2], &cell_map, 1);
        near.sort_unstable();
        assert_eq!(near, vec![1, 2]);
        assert!(nearby_cell_indices(grid, [2, 2, 2], &cell_map, 0).is_empty());
    }

    #[test]
    fn division_neighbor_search_avoids_reserved_birth_positions() {
        let grid = GridDims { x: 4, y: 4, z: 4 };
        let sim = SimulationConfig {
            division_neighbor_distance: 1,
            ..SimulationConfig::default()
        };
        let cell_map = HashMap::from([
            ([1, 1, 1], 0),
            ([0, 1, 1], 1),
            ([1, 0, 1], 2),
            ([1, 2, 1], 3),
            ([1, 1, 0], 4),
            ([1, 1, 2], 5),
        ]);
        let reserved = HashSet::from([[2, 1, 1]]);
        let mut rng = rand::rng();

        let found = find_empty_neighbor_avoiding(
            grid,
            [1, 1, 1],
            &cell_map,
            Some(&reserved),
            &mut rng,
            &sim,
        );

        assert_eq!(found, None);
    }

    #[test]
    fn large_grid_neighbor_math_does_not_overflow_near_boundary() {
        let grid = GridDims {
            x: i16::MAX as usize,
            y: 1,
            z: 1,
        };
        let cell_map = HashMap::new();
        let pos = [(i16::MAX - 1) as u16, 0, 0];

        assert_eq!(empty_neighbors(grid, pos, &cell_map), vec![[32765, 0, 0]]);

        let sim = SimulationConfig {
            division_neighbor_distance: 2,
            ..SimulationConfig::default()
        };
        let mut rng = rand::rng();

        assert_eq!(
            find_empty_neighbor_avoiding(grid, pos, &cell_map, None, &mut rng, &sim),
            Some([32764, 0, 0])
        );
    }
}
