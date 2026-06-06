//! Compile-time grid constants and runtime configuration for MARL.
//!
//! This crate has no simulation dependencies — it is a pure data crate
//! consumed by every other crate in the workspace.

pub mod stoich;
use stoich::StoichEnforcement;

// ============================================================================
// DEFAULT GRID DIMENSIONS
// ============================================================================
// Suggested sizes:
//   64x64x32   — quick debug runs (~7 ticks/sec)
//   128x128x64 — calibration runs (~1 tick/sec est.)
//   256x256x128 — production runs (needs rayon, ~0.1 tick/sec est.)
pub const DEFAULT_GRID_X: usize = 128;
pub const DEFAULT_GRID_Y: usize = 128;
pub const DEFAULT_GRID_Z: usize = 64;

// Backward-compatible default dimensions for tests and GPU fallback paths.
pub const GRID_X: usize = DEFAULT_GRID_X;
pub const GRID_Y: usize = DEFAULT_GRID_Y;
pub const GRID_Z: usize = DEFAULT_GRID_Z;

const _: () = assert!(DEFAULT_GRID_X <= i16::MAX as usize);
const _: () = assert!(DEFAULT_GRID_Y <= i16::MAX as usize);
const _: () = assert!(DEFAULT_GRID_Z <= i16::MAX as usize);

// Species counts
pub const S_EXT: usize = 12; // external chemical species
pub const M_INT: usize = 16; // internal chemical species

// Reaction network
pub const R_MAX: usize = 16; // max reactions per cell
pub const S_RECEPTORS: usize = 8;
pub const S_TRANSPORTERS: usize = 8;
pub const S_EFFECTORS: usize = 8;

pub const EXT_ENERGY: usize = 0;
pub const EXT_OXIDANT: usize = 1;
pub const EXT_REDUCTANT: usize = 2;
pub const EXT_CARBON: usize = 3;
pub const EXT_ORGANIC: usize = 4;
pub const EXT_SIGNAL_A: usize = 5;
pub const EXT_SIGNAL_B: usize = 6;
pub const EXT_STRUCTURAL: usize = 7;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChemicalComposition {
    pub carbon_backbone: f32,
    pub oxidizing_power: f32,
    pub reducing_power: f32,
    pub phosphate_like_activation: f32,
    pub lipid_like_tail: f32,
    pub signal_group: f32,
    pub structural_group: f32,
    pub toxin_group: f32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExternalSpeciesDescriptor {
    pub species: usize,
    pub name: &'static str,
    pub composition: ChemicalComposition,
    pub bond_energy: f32,
    pub extracellular_stability: f32,
    pub membrane_permeability: f32,
    pub storage_density: f32,
    pub work_coupling: f32,
    pub default_diffusion: f32,
    pub default_decay: f32,
}

impl ExternalSpeciesDescriptor {
    pub const fn inert(species: usize, name: &'static str) -> Self {
        Self {
            species,
            name,
            composition: ChemicalComposition {
                carbon_backbone: 0.0,
                oxidizing_power: 0.0,
                reducing_power: 0.0,
                phosphate_like_activation: 0.0,
                lipid_like_tail: 0.0,
                signal_group: 0.0,
                structural_group: 0.0,
                toxin_group: 0.0,
            },
            bond_energy: 0.0,
            extracellular_stability: 1.0,
            membrane_permeability: 0.2,
            storage_density: 0.0,
            work_coupling: 0.0,
            default_diffusion: 0.3,
            default_decay: 0.01,
        }
    }
}

pub const EXTERNAL_SPECIES: [ExternalSpeciesDescriptor; S_EXT] = [
    ExternalSpeciesDescriptor {
        species: EXT_ENERGY,
        name: "free_energy",
        composition: ChemicalComposition {
            carbon_backbone: 0.0,
            oxidizing_power: 0.0,
            reducing_power: 0.0,
            phosphate_like_activation: 1.0,
            lipid_like_tail: 0.0,
            signal_group: 0.0,
            structural_group: 0.0,
            toxin_group: 0.0,
        },
        bond_energy: 1.0,
        extracellular_stability: 0.15,
        membrane_permeability: 0.05,
        storage_density: 0.1,
        work_coupling: 1.0,
        default_diffusion: 0.0,
        default_decay: 0.2,
    },
    ExternalSpeciesDescriptor {
        species: EXT_OXIDANT,
        name: "oxidant",
        composition: ChemicalComposition {
            carbon_backbone: 0.0,
            oxidizing_power: 1.0,
            reducing_power: 0.0,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.0,
            signal_group: 0.0,
            structural_group: 0.0,
            toxin_group: 0.25,
        },
        bond_energy: 0.0,
        extracellular_stability: 0.8,
        membrane_permeability: 0.7,
        storage_density: 0.0,
        work_coupling: 0.15,
        default_diffusion: 1.5,
        default_decay: 0.01,
    },
    ExternalSpeciesDescriptor {
        species: EXT_REDUCTANT,
        name: "reductant",
        composition: ChemicalComposition {
            carbon_backbone: 0.0,
            oxidizing_power: 0.0,
            reducing_power: 1.0,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.0,
            signal_group: 0.0,
            structural_group: 0.0,
            toxin_group: 0.0,
        },
        bond_energy: 0.6,
        extracellular_stability: 0.75,
        membrane_permeability: 0.55,
        storage_density: 0.1,
        work_coupling: 0.4,
        default_diffusion: 1.0,
        default_decay: 0.01,
    },
    ExternalSpeciesDescriptor {
        species: EXT_CARBON,
        name: "carbon",
        composition: ChemicalComposition {
            carbon_backbone: 1.0,
            oxidizing_power: 0.0,
            reducing_power: 0.0,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.0,
            signal_group: 0.0,
            structural_group: 0.0,
            toxin_group: 0.0,
        },
        bond_energy: 0.35,
        extracellular_stability: 0.9,
        membrane_permeability: 0.45,
        storage_density: 0.45,
        work_coupling: 0.2,
        default_diffusion: 1.2,
        default_decay: 0.005,
    },
    ExternalSpeciesDescriptor {
        species: EXT_ORGANIC,
        name: "organic",
        composition: ChemicalComposition {
            carbon_backbone: 0.8,
            oxidizing_power: 0.0,
            reducing_power: 0.1,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.1,
            signal_group: 0.0,
            structural_group: 0.0,
            toxin_group: 0.2,
        },
        bond_energy: 0.25,
        extracellular_stability: 0.55,
        membrane_permeability: 0.25,
        storage_density: 0.35,
        work_coupling: 0.1,
        default_diffusion: 0.8,
        default_decay: 0.03,
    },
    ExternalSpeciesDescriptor {
        species: EXT_SIGNAL_A,
        name: "signal_a",
        composition: ChemicalComposition {
            carbon_backbone: 0.2,
            oxidizing_power: 0.0,
            reducing_power: 0.0,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.0,
            signal_group: 1.0,
            structural_group: 0.0,
            toxin_group: 0.0,
        },
        bond_energy: 0.05,
        extracellular_stability: 0.45,
        membrane_permeability: 0.4,
        storage_density: 0.0,
        work_coupling: 0.0,
        default_diffusion: 0.5,
        default_decay: 0.05,
    },
    ExternalSpeciesDescriptor {
        species: EXT_SIGNAL_B,
        name: "signal_b",
        composition: ChemicalComposition {
            carbon_backbone: 0.2,
            oxidizing_power: 0.0,
            reducing_power: 0.0,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.0,
            signal_group: 1.0,
            structural_group: 0.0,
            toxin_group: 0.0,
        },
        bond_energy: 0.05,
        extracellular_stability: 0.45,
        membrane_permeability: 0.4,
        storage_density: 0.0,
        work_coupling: 0.0,
        default_diffusion: 0.5,
        default_decay: 0.05,
    },
    ExternalSpeciesDescriptor {
        species: EXT_STRUCTURAL,
        name: "structural",
        composition: ChemicalComposition {
            carbon_backbone: 0.4,
            oxidizing_power: 0.0,
            reducing_power: 0.0,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.2,
            signal_group: 0.0,
            structural_group: 1.0,
            toxin_group: 0.0,
        },
        bond_energy: 0.1,
        extracellular_stability: 0.95,
        membrane_permeability: 0.02,
        storage_density: 0.2,
        work_coupling: 0.0,
        default_diffusion: 0.1,
        default_decay: 0.002,
    },
    ExternalSpeciesDescriptor::inert(8, "spare_0"),
    ExternalSpeciesDescriptor::inert(9, "spare_1"),
    ExternalSpeciesDescriptor::inert(10, "spare_2"),
    ExternalSpeciesDescriptor::inert(11, "spare_3"),
];

pub fn external_species_name(species: usize) -> &'static str {
    EXTERNAL_SPECIES
        .get(species)
        .map(|descriptor| descriptor.name)
        .unwrap_or("unknown")
}

pub fn default_external_diffusion() -> [f32; S_EXT] {
    let mut out = [0.0; S_EXT];
    let mut i = 0;
    while i < S_EXT {
        out[i] = EXTERNAL_SPECIES[i].default_diffusion;
        i += 1;
    }
    out
}

pub fn default_external_decay() -> [f32; S_EXT] {
    let mut out = [0.0; S_EXT];
    let mut i = 0;
    while i < S_EXT {
        out[i] = EXTERNAL_SPECIES[i].default_decay;
        i += 1;
    }
    out
}

// ============================================================================
// RUNTIME CONFIGURATION — SimulationConfig + OutputConfig
// ============================================================================
// Grid dimensions, physics, chemistry, biology, and output parameters are now
// runtime-configurable via an optional TOML file and CLI overrides. Species and
// ruleset array sizes remain compile-time because they determine fixed-size
// cell/ruleset arrays.

/// Runtime grid dimensions for the 3D simulation domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(default)]
pub struct GridDims {
    #[serde(default = "default_grid_x")]
    pub x: usize,
    #[serde(default = "default_grid_y")]
    pub y: usize,
    #[serde(default = "default_grid_z")]
    pub z: usize,
}

fn default_grid_x() -> usize {
    DEFAULT_GRID_X
}

fn default_grid_y() -> usize {
    DEFAULT_GRID_Y
}

fn default_grid_z() -> usize {
    DEFAULT_GRID_Z
}

impl GridDims {
    pub fn validate(self) -> Result<(), String> {
        if self.x == 0 || self.y == 0 || self.z == 0 {
            return Err(format!(
                "grid dimensions must be nonzero, got {}x{}x{}",
                self.x, self.y, self.z
            ));
        }
        if self.x > i16::MAX as usize || self.y > i16::MAX as usize || self.z > i16::MAX as usize {
            return Err(format!(
                "grid dimensions must fit i16/u16 coordinate math, got {}x{}x{}",
                self.x, self.y, self.z
            ));
        }
        self.voxel_count()
            .ok_or_else(|| "grid voxel count overflow".to_string())?;
        Ok(())
    }

    pub fn voxel_count(self) -> Option<usize> {
        self.x.checked_mul(self.y)?.checked_mul(self.z)
    }

    pub fn field_float_count(self) -> Option<usize> {
        self.voxel_count()?.checked_mul(S_EXT)
    }
}

impl Default for GridDims {
    fn default() -> Self {
        Self {
            x: DEFAULT_GRID_X,
            y: DEFAULT_GRID_Y,
            z: DEFAULT_GRID_Z,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundaryFace {
    Top,
    Bottom,
}

impl BoundaryFace {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Top => "top",
            Self::Bottom => "bottom",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
pub struct BoundarySourceConfig {
    pub species: usize,
    pub face: BoundaryFace,
    pub rate: f32,
}

impl BoundarySourceConfig {
    pub fn validate(self) -> Result<(), String> {
        if self.species >= S_EXT {
            return Err(format!(
                "boundary source species {} is out of range for {S_EXT} external species",
                self.species
            ));
        }
        if !self.rate.is_finite() || self.rate < 0.0 {
            return Err(format!(
                "boundary source species {} has invalid rate {}",
                self.species, self.rate
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
pub struct BoundaryPrimeConfig {
    pub species: usize,
    pub face: BoundaryFace,
    pub concentration: f32,
}

impl BoundaryPrimeConfig {
    pub fn validate(self) -> Result<(), String> {
        if self.species >= S_EXT {
            return Err(format!(
                "boundary prime species {} is out of range for {S_EXT} external species",
                self.species
            ));
        }
        if !self.concentration.is_finite() || self.concentration < 0.0 {
            return Err(format!(
                "boundary prime species {} has invalid concentration {}",
                self.species, self.concentration
            ));
        }
        Ok(())
    }
}

/// Physics, chemistry, biology, and seeding parameters.
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(default)]
pub struct SimulationConfig {
    // Spatiotemporal
    pub dx: f32,
    pub dt: f32,
    pub diffusion_substeps: usize,

    // Diffusion & decay (arrays, length S_EXT)
    pub d_voxel: [f32; S_EXT],
    pub lambda_decay: [f32; S_EXT],

    // Boundary sources
    pub source_rate_oxidant: f32,
    pub source_rate_carbon: f32,
    pub source_rate_reductant: f32,
    pub boundary_sources: Vec<BoundarySourceConfig>,

    // Cell metabolism
    pub epsilon: f32,
    pub c_max: f32,
    pub lambda_maintenance: f32,
    pub hard_death_floor: f32,
    pub reaction_maintenance: f32,

    // Cell cycle
    pub base_division_prep: f32,
    pub prep_maintenance_multiplier: f32,
    pub rush_penalty_rate: f32,

    // Niche construction
    pub alpha_eps: f32,
    pub k_eps: f32,

    // Light
    pub light_efficiency: f32,
    pub surface_intensity: f32,
    pub cell_absorption: f32,
    pub chemical_absorption: f32,
    pub light_floor: f32,

    // Mutation
    pub mutation_stddev: f32,
    pub structural_mutation_rate_mult: f32,
    pub meta_mutation_rate: f32,
    pub meta_mutation_clamp_low: f32,
    pub meta_mutation_clamp_high: f32,
    pub hill_exponent_clamp_low: f32,
    pub hill_exponent_clamp_high: f32,
    pub active_reaction_threshold: f32,

    // Horizontal gene transfer
    pub hgt_enabled: bool,
    pub hgt_interval: u32,
    pub hgt_radius: u8,
    pub hgt_base_rate: f32,
    pub hgt_max_events_per_tick: usize,

    // Seeding geometry (canonical 200-layer units)
    pub seed_margin: u16,
    pub phototroph_z_lo: f32,
    pub phototroph_z_hi: f32,
    pub chemolithotroph_z_lo: f32,
    pub chemolithotroph_z_hi: f32,
    pub anaerobe_z_lo: f32,
    pub anaerobe_z_hi: f32,

    // Division neighbor search
    pub division_neighbor_distance: u8,

    // Field initialization boundary priming
    pub boundary_prime_layers: usize,
    pub boundary_prime_oxidant: f32,
    pub boundary_prime_carbon: f32,
    pub boundary_prime_reductant: f32,
    pub boundary_primes: Vec<BoundaryPrimeConfig>,

    // Stoichiometry accounting / enforcement
    pub stoich_enforcement: StoichEnforcement,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            dx: 100.0e-6,
            dt: 1.0,
            diffusion_substeps: 10,

            d_voxel: default_external_diffusion(),
            lambda_decay: default_external_decay(),

            source_rate_oxidant: 0.4,
            source_rate_carbon: 0.15,
            source_rate_reductant: 0.5,
            boundary_sources: Vec::new(),

            epsilon: 0.001,
            c_max: 10.0,
            lambda_maintenance: 0.12,
            hard_death_floor: 0.01,
            reaction_maintenance: 0.003,

            base_division_prep: 20.0,
            prep_maintenance_multiplier: 2.0,
            rush_penalty_rate: 0.05,

            alpha_eps: 0.8,
            k_eps: 2.0,

            light_efficiency: 0.0,
            surface_intensity: 1.0,
            cell_absorption: 0.3,
            chemical_absorption: 0.05,
            light_floor: 1e-7,

            mutation_stddev: 0.1,
            structural_mutation_rate_mult: 0.1,
            meta_mutation_rate: 0.01,
            meta_mutation_clamp_low: 0.001,
            meta_mutation_clamp_high: 0.5,
            hill_exponent_clamp_low: 0.5,
            hill_exponent_clamp_high: 8.0,
            active_reaction_threshold: 1e-9,

            hgt_enabled: false,
            hgt_interval: 10,
            hgt_radius: 1,
            hgt_base_rate: 0.02,
            hgt_max_events_per_tick: 100,

            seed_margin: 5,
            phototroph_z_lo: 0.0,
            phototroph_z_hi: 3.0,
            chemolithotroph_z_lo: 80.0,
            chemolithotroph_z_hi: 130.0,
            anaerobe_z_lo: 120.0,
            anaerobe_z_hi: 180.0,

            division_neighbor_distance: 2,

            boundary_prime_layers: 2,
            boundary_prime_oxidant: 0.5,
            boundary_prime_carbon: 0.3,
            boundary_prime_reductant: 0.5,
            boundary_primes: Vec::new(),

            stoich_enforcement: StoichEnforcement::Off,
        }
    }
}

impl SimulationConfig {
    pub fn effective_boundary_sources(&self) -> Vec<BoundarySourceConfig> {
        if !self.boundary_sources.is_empty() {
            return self.boundary_sources.clone();
        }
        vec![
            BoundarySourceConfig {
                species: EXT_OXIDANT,
                face: BoundaryFace::Top,
                rate: self.source_rate_oxidant,
            },
            BoundarySourceConfig {
                species: EXT_CARBON,
                face: BoundaryFace::Top,
                rate: self.source_rate_carbon,
            },
            BoundarySourceConfig {
                species: EXT_REDUCTANT,
                face: BoundaryFace::Bottom,
                rate: self.source_rate_reductant,
            },
        ]
    }

    pub fn effective_boundary_primes(&self) -> Vec<BoundaryPrimeConfig> {
        if !self.boundary_primes.is_empty() {
            return self.boundary_primes.clone();
        }
        vec![
            BoundaryPrimeConfig {
                species: EXT_OXIDANT,
                face: BoundaryFace::Top,
                concentration: self.boundary_prime_oxidant,
            },
            BoundaryPrimeConfig {
                species: EXT_CARBON,
                face: BoundaryFace::Top,
                concentration: self.boundary_prime_carbon,
            },
            BoundaryPrimeConfig {
                species: EXT_REDUCTANT,
                face: BoundaryFace::Bottom,
                concentration: self.boundary_prime_reductant,
            },
        ]
    }

    pub fn validate_chemistry(&self) -> Result<(), String> {
        for source in self.effective_boundary_sources() {
            source.validate()?;
        }
        for prime in self.effective_boundary_primes() {
            prime.validate()?;
        }
        Ok(())
    }
}

/// Logging cadence, snapshot selection, image toggles, and output directory.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BinaryCompression {
    #[default]
    None,
    Zstd,
}

impl BinaryCompression {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Zstd => "zstd",
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RulesetOutputMode {
    #[default]
    Off,
    LayerAverages,
    Full,
    Both,
}

impl RulesetOutputMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::LayerAverages => "layer_averages",
            Self::Full => "full",
            Self::Both => "both",
        }
    }

    pub fn is_enabled(self) -> bool {
        self != Self::Off
    }

    pub fn writes_layer_averages(self) -> bool {
        matches!(self, Self::LayerAverages | Self::Both)
    }

    pub fn writes_full_dump(self) -> bool {
        matches!(self, Self::Full | Self::Both)
    }
}

/// Logging cadence, snapshot selection, image toggles, and output directory.
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(default)]
pub struct OutputConfig {
    pub max_ticks: u32,
    pub stats_interval: u32,
    pub snapshot_interval: u32,
    pub image_interval: u32,
    pub seed_count: usize,
    pub output_dir: String,

    // Snapshot species indices for XZ cross-sections
    pub xz_snapshot_species: Vec<usize>,

    // XY slice depths (as fractions, 0.0..1.0, resolved at runtime)
    pub xy_slice_depths_frac: Vec<f32>,

    // Toggle image types
    pub write_binary_field: bool,
    pub write_binary_cells: bool,
    pub binary_compression: BinaryCompression,
    pub binary_compression_level: i32,
    pub ruleset_interval: u32,
    pub ruleset_output_mode: RulesetOutputMode,
    pub write_tick_log: bool,
    pub write_csv_snapshots: bool,
    pub write_stoich_summary: bool,
    pub write_stoich_tick_log: bool,
    pub write_stoich_v2_summary: bool,
    pub write_stoich_v2_events: bool,
    pub write_ancestry_map: bool,
    pub write_density_map: bool,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            max_ticks: 5000,
            stats_interval: 100,
            snapshot_interval: 500,
            image_interval: 500,
            seed_count: 30,
            output_dir: format!(
                "output/run_{}x{}x{}",
                DEFAULT_GRID_X, DEFAULT_GRID_Y, DEFAULT_GRID_Z
            ),
            xz_snapshot_species: Vec::new(),
            xy_slice_depths_frac: Vec::new(),
            write_binary_field: true,
            write_binary_cells: true,
            binary_compression: BinaryCompression::Zstd,
            binary_compression_level: 3,
            ruleset_interval: 1000,
            ruleset_output_mode: RulesetOutputMode::Off,
            write_tick_log: false,
            write_csv_snapshots: false,
            write_stoich_summary: false,
            write_stoich_tick_log: false,
            write_stoich_v2_summary: false,
            write_stoich_v2_events: false,
            write_ancestry_map: false,
            write_density_map: false,
        }
    }
}

/// Unified configuration: simulation + output.
#[derive(Debug, Clone, Default, serde::Deserialize)]
pub struct Config {
    #[serde(default)]
    pub grid: GridDims,
    #[serde(default)]
    pub simulation: SimulationConfig,
    #[serde(default)]
    pub output: OutputConfig,
}

impl Config {
    /// Load configuration with hierarchy:
    ///   1. Built-in defaults
    ///   2. TOML file override (`marl.toml` or `--config <path>`)
    ///   3. CLI flags override (run-control only: --ticks, --stats, etc.)
    pub fn load() -> Self {
        let mut cfg = Self::default();

        // 2. TOML file override
        let args: Vec<String> = std::env::args().collect();
        let mut config_path: Option<String> = None;
        let mut i = 1;
        while i < args.len() {
            if args[i] == "--config" && i + 1 < args.len() {
                config_path = Some(args[i + 1].clone());
                i += 2;
            } else {
                i += 1;
            }
        }
        let toml_path = config_path.unwrap_or_else(|| "marl.toml".to_string());
        if let Ok(content) = std::fs::read_to_string(&toml_path) {
            match toml::from_str::<Config>(&content) {
                Ok(parsed) => cfg = parsed,
                Err(err) => eprintln!("Warning: failed to parse config {}: {err}", toml_path),
            }
        }

        // 3. CLI override (run-control flags only)
        let mut i = 1;
        while i < args.len() {
            match args[i].as_str() {
                "--config" => i += 2,
                "--ticks" if i + 1 < args.len() => {
                    if let Ok(v) = args[i + 1].parse() {
                        cfg.output.max_ticks = v;
                    }
                    i += 2;
                }
                "--stats" if i + 1 < args.len() => {
                    if let Ok(v) = args[i + 1].parse() {
                        cfg.output.stats_interval = v;
                    }
                    i += 2;
                }
                "--snapshot" if i + 1 < args.len() => {
                    if let Ok(v) = args[i + 1].parse() {
                        cfg.output.snapshot_interval = v;
                    }
                    i += 2;
                }
                "--ruleset-interval" if i + 1 < args.len() => {
                    if let Ok(v) = args[i + 1].parse() {
                        cfg.output.ruleset_interval = v;
                    }
                    i += 2;
                }
                "--images" if i + 1 < args.len() => {
                    if let Ok(v) = args[i + 1].parse() {
                        cfg.output.image_interval = v;
                    }
                    i += 2;
                }
                "--seed" if i + 1 < args.len() => {
                    if let Ok(v) = args[i + 1].parse() {
                        cfg.output.seed_count = v;
                    }
                    i += 2;
                }
                "--output" if i + 1 < args.len() => {
                    cfg.output.output_dir = args[i + 1].clone();
                    i += 2;
                }
                _ => {
                    i += 1;
                }
            }
        }

        cfg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_config_defaults_match_binary_plan() {
        let out = OutputConfig::default();
        assert!(out.write_binary_field);
        assert!(out.write_binary_cells);
        assert_eq!(out.binary_compression, BinaryCompression::Zstd);
        assert_eq!(out.binary_compression_level, 3);
        assert_eq!(out.ruleset_interval, 1000);
        assert_eq!(out.ruleset_output_mode, RulesetOutputMode::Off);
        assert!(!out.write_stoich_summary);
        assert!(!out.write_stoich_tick_log);
        assert!(!out.write_stoich_v2_summary);
        assert!(!out.write_stoich_v2_events);
    }

    #[test]
    fn external_species_catalog_matches_species_slots() {
        assert_eq!(EXTERNAL_SPECIES.len(), S_EXT);
        for (index, descriptor) in EXTERNAL_SPECIES.iter().enumerate() {
            assert_eq!(descriptor.species, index);
            assert_eq!(external_species_name(index), descriptor.name);
            assert!(descriptor.default_diffusion.is_finite());
            assert!(descriptor.default_decay.is_finite());
            assert!(descriptor.default_diffusion >= 0.0);
            assert!(descriptor.default_decay >= 0.0);
        }
    }

    #[test]
    fn default_physics_comes_from_external_species_catalog() {
        let sim = SimulationConfig::default();
        assert_eq!(sim.d_voxel, default_external_diffusion());
        assert_eq!(sim.lambda_decay, default_external_decay());
        assert_eq!(sim.lambda_decay[EXT_ENERGY], 0.2);
        assert!(
            EXTERNAL_SPECIES[EXT_ENERGY].work_coupling > EXTERNAL_SPECIES[EXT_CARBON].work_coupling
        );
        assert!(
            EXTERNAL_SPECIES[EXT_CARBON].storage_density
                > EXTERNAL_SPECIES[EXT_ENERGY].storage_density
        );
    }

    #[test]
    fn grid_defaults_and_toml_override_work() {
        let cfg = Config::default();
        assert_eq!(cfg.grid, GridDims::default());
        assert!(cfg.grid.validate().is_ok());

        let cfg: Config = toml::from_str(
            r#"
            [grid]
            x = 64
            y = 64
            z = 32
            "#,
        )
        .unwrap();

        assert_eq!(
            cfg.grid,
            GridDims {
                x: 64,
                y: 64,
                z: 32
            }
        );
        assert!(cfg.grid.validate().is_ok());

        let cfg: Config = toml::from_str(
            r#"
            [grid]
            x = 96
            "#,
        )
        .unwrap();

        assert_eq!(
            cfg.grid,
            GridDims {
                x: 96,
                y: DEFAULT_GRID_Y,
                z: DEFAULT_GRID_Z,
            }
        );
    }

    #[test]
    fn invalid_grid_dimensions_are_rejected() {
        assert!(GridDims { x: 0, y: 64, z: 32 }.validate().is_err());
        assert!(
            GridDims {
                x: i16::MAX as usize + 1,
                y: 64,
                z: 32,
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn enum_labels_match_toml_strings() {
        assert_eq!(BinaryCompression::None.as_str(), "none");
        assert_eq!(BinaryCompression::Zstd.as_str(), "zstd");
        assert_eq!(RulesetOutputMode::Off.as_str(), "off");
        assert_eq!(RulesetOutputMode::LayerAverages.as_str(), "layer_averages");
        assert_eq!(RulesetOutputMode::Full.as_str(), "full");
        assert_eq!(RulesetOutputMode::Both.as_str(), "both");
        assert!(!RulesetOutputMode::Off.is_enabled());
        assert!(RulesetOutputMode::LayerAverages.is_enabled());
        assert!(RulesetOutputMode::Full.is_enabled());
        assert!(RulesetOutputMode::Both.is_enabled());
    }

    #[test]
    fn ruleset_mode_writes_layer_and_full() {
        assert!(!RulesetOutputMode::Off.writes_layer_averages());
        assert!(RulesetOutputMode::LayerAverages.writes_layer_averages());
        assert!(!RulesetOutputMode::Full.writes_layer_averages());
        assert!(RulesetOutputMode::Both.writes_layer_averages());

        assert!(!RulesetOutputMode::Off.writes_full_dump());
        assert!(!RulesetOutputMode::LayerAverages.writes_full_dump());
        assert!(RulesetOutputMode::Full.writes_full_dump());
        assert!(RulesetOutputMode::Both.writes_full_dump());
    }

    #[test]
    fn stoich_output_flags_deserialize_from_toml() {
        let cfg: Config = toml::from_str(
            r#"
            [output]
            write_stoich_summary = true
            write_stoich_tick_log = true
            write_stoich_v2_summary = true
            write_stoich_v2_events = true
            "#,
        )
        .unwrap();

        assert!(cfg.output.write_stoich_summary);
        assert!(cfg.output.write_stoich_tick_log);
        assert!(cfg.output.write_stoich_v2_summary);
        assert!(cfg.output.write_stoich_v2_events);
    }

    #[test]
    fn stoich_enforcement_deserializes_from_toml() {
        let cfg: Config = toml::from_str(
            r#"
            [simulation]
            stoich_enforcement = "strict"
            "#,
        )
        .unwrap();

        assert_eq!(cfg.simulation.stoich_enforcement, StoichEnforcement::Strict);
    }

    #[test]
    fn explicit_boundary_chemistry_replaces_legacy_sources() {
        let cfg: Config = toml::from_str(
            r#"
            [simulation]
            source_rate_oxidant = 9.0
            boundary_sources = [
                { species = 4, face = "bottom", rate = 0.25 },
            ]
            boundary_primes = [
                { species = 5, face = "top", concentration = 0.75 },
            ]
            "#,
        )
        .unwrap();

        let sources = cfg.simulation.effective_boundary_sources();
        assert_eq!(sources.len(), 1);
        assert_eq!(sources[0].species, EXT_ORGANIC);
        assert_eq!(sources[0].face, BoundaryFace::Bottom);
        assert_eq!(sources[0].rate, 0.25);

        let primes = cfg.simulation.effective_boundary_primes();
        assert_eq!(primes.len(), 1);
        assert_eq!(primes[0].species, EXT_SIGNAL_A);
        assert_eq!(primes[0].face, BoundaryFace::Top);
        assert_eq!(primes[0].concentration, 0.75);
        assert!(cfg.simulation.validate_chemistry().is_ok());
    }

    #[test]
    fn chemistry_validation_rejects_bad_species_and_rates() {
        let mut sim = SimulationConfig {
            boundary_sources: vec![BoundarySourceConfig {
                species: S_EXT,
                face: BoundaryFace::Top,
                rate: 0.1,
            }],
            ..Default::default()
        };
        assert!(sim.validate_chemistry().is_err());

        sim.boundary_sources = vec![BoundarySourceConfig {
            species: EXT_CARBON,
            face: BoundaryFace::Top,
            rate: -0.1,
        }];
        assert!(sim.validate_chemistry().is_err());

        sim.boundary_sources.clear();
        sim.boundary_primes = vec![BoundaryPrimeConfig {
            species: EXT_CARBON,
            face: BoundaryFace::Top,
            concentration: f32::NAN,
        }];
        assert!(sim.validate_chemistry().is_err());
    }

    #[test]
    fn hgt_defaults_keep_runtime_disabled() {
        let sim = SimulationConfig::default();
        assert!(!sim.hgt_enabled);
        assert_eq!(sim.hgt_interval, 10);
        assert_eq!(sim.hgt_radius, 1);
        assert_eq!(sim.hgt_base_rate, 0.02);
        assert_eq!(sim.hgt_max_events_per_tick, 100);
    }

    #[test]
    fn hgt_config_deserializes_from_toml() {
        let cfg: Config = toml::from_str(
            r#"
            [simulation]
            hgt_enabled = true
            hgt_interval = 3
            hgt_radius = 2
            hgt_base_rate = 0.5
            hgt_max_events_per_tick = 7
            "#,
        )
        .unwrap();

        assert!(cfg.simulation.hgt_enabled);
        assert_eq!(cfg.simulation.hgt_interval, 3);
        assert_eq!(cfg.simulation.hgt_radius, 2);
        assert_eq!(cfg.simulation.hgt_base_rate, 0.5);
        assert_eq!(cfg.simulation.hgt_max_events_per_tick, 7);
    }
}
