//! Compile-time grid constants and runtime configuration for MARL.
//!
//! This crate has no simulation dependencies — it is a pure data crate
//! consumed by every other crate in the workspace.

pub mod stoich;

// ============================================================================
// GRID DIMENSIONS — compile-time constants
// ============================================================================
// These must be const because they determine array sizes throughout the code.
// To run at a different grid size, change these values and recompile:
//   cargo build --release    (~2 seconds incremental)
//
// Suggested sizes:
//   64x64x32   — quick debug runs (~7 ticks/sec)
//   128x128x64 — calibration runs (~1 tick/sec est.)
//   256x256x128 — production runs (needs rayon, ~0.1 tick/sec est.)
pub const GRID_X: usize = 128;
pub const GRID_Y: usize = 128;
pub const GRID_Z: usize = 64;

// Species counts
pub const S_EXT: usize = 12; // external chemical species
pub const M_INT: usize = 16; // internal chemical species

// Reaction network
pub const R_MAX: usize = 16; // max reactions per cell
pub const S_RECEPTORS: usize = 8;
pub const S_TRANSPORTERS: usize = 8;
pub const S_EFFECTORS: usize = 8;

// ============================================================================
// RUNTIME CONFIGURATION — SimulationConfig + OutputConfig
// ============================================================================
// All physics, chemistry, biology, and output parameters are now runtime-
// configurable via an optional TOML file and CLI overrides. Only the array-
// size constants above remain compile-time.

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
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            dx: 100.0e-6,
            dt: 1.0,
            diffusion_substeps: 10,

            d_voxel: [0.0, 1.5, 1.0, 1.2, 0.8, 0.5, 0.5, 0.1, 0.3, 0.3, 0.3, 0.3],
            lambda_decay: [
                0.0, 0.01, 0.01, 0.005, 0.03, 0.05, 0.05, 0.002, 0.01, 0.01, 0.01, 0.01,
            ],

            source_rate_oxidant: 0.4,
            source_rate_carbon: 0.15,
            source_rate_reductant: 0.5,

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
        }
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
            output_dir: format!("output/run_{}x{}x{}", GRID_X, GRID_Y, GRID_Z),
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
            write_ancestry_map: false,
            write_density_map: false,
        }
    }
}

/// Unified configuration: simulation + output.
#[derive(Debug, Clone, Default, serde::Deserialize)]
pub struct Config {
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
}
