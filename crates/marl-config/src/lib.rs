//! Compile-time grid constants and runtime configuration for MARL.
//!
//! This crate has no simulation dependencies — it is a pure data crate
//! consumed by every other crate in the workspace.

pub mod stoich;
use serde::Serialize;
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

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
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

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
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

fn enzyme_descriptor(species: usize, name: &'static str) -> ExternalSpeciesDescriptor {
    ExternalSpeciesDescriptor {
        species,
        name,
        composition: ChemicalComposition {
            carbon_backbone: 0.8,
            oxidizing_power: 0.0,
            reducing_power: 0.0,
            phosphate_like_activation: 0.0,
            lipid_like_tail: 0.0,
            signal_group: 0.0,
            structural_group: 0.8,
            toxin_group: 0.0,
        },
        bond_energy: 0.15,
        extracellular_stability: 0.4,
        membrane_permeability: 0.0,
        storage_density: 0.1,
        work_coupling: 0.05,
        default_diffusion: 0.0,
        default_decay: 0.0,
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
    external_species_descriptor(species).name
}

pub fn external_species_descriptor(species: usize) -> ExternalSpeciesDescriptor {
    EXTERNAL_SPECIES
        .get(species)
        .copied()
        .unwrap_or_else(|| ExternalSpeciesDescriptor::inert(species, "unknown"))
}

pub fn internal_species_descriptor(species: usize) -> ExternalSpeciesDescriptor {
    match species {
        EXT_ENERGY | EXT_OXIDANT | EXT_REDUCTANT | EXT_CARBON | EXT_ORGANIC => {
            external_species_descriptor(species)
        }
        5 => enzyme_descriptor(species, "enzyme_a"),
        6 => enzyme_descriptor(species, "enzyme_b"),
        7 => {
            let mut descriptor = external_species_descriptor(EXT_STRUCTURAL);
            descriptor.name = "carbon_reserve";
            descriptor
        }
        _ => ExternalSpeciesDescriptor::inert(species, "inactive"),
    }
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

fn descriptor_composition_load(composition: &ChemicalComposition) -> f32 {
    composition.carbon_backbone.max(0.0)
        + composition.oxidizing_power.max(0.0)
        + composition.reducing_power.max(0.0)
        + composition.phosphate_like_activation.max(0.0)
        + composition.lipid_like_tail.max(0.0)
        + composition.signal_group.max(0.0)
        + composition.structural_group.max(0.0)
        + composition.toxin_group.max(0.0)
}

fn descriptor_composition_overlap(a: &ChemicalComposition, b: &ChemicalComposition) -> f32 {
    let shared = a.carbon_backbone.min(b.carbon_backbone).max(0.0)
        + a.oxidizing_power.min(b.oxidizing_power).max(0.0)
        + a.reducing_power.min(b.reducing_power).max(0.0)
        + a.phosphate_like_activation
            .min(b.phosphate_like_activation)
            .max(0.0)
        + a.lipid_like_tail.min(b.lipid_like_tail).max(0.0)
        + a.signal_group.min(b.signal_group).max(0.0)
        + a.structural_group.min(b.structural_group).max(0.0)
        + a.toxin_group.min(b.toxin_group).max(0.0);
    let total = descriptor_composition_load(a).max(descriptor_composition_load(b));
    if total <= f32::EPSILON {
        1.0
    } else {
        (shared / total).clamp(0.0, 1.0)
    }
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
    pub transport_energy_cost_scale: f32,
    pub transport_permeability_cost_weight: f32,
    pub transport_composition_cost_weight: f32,
    pub reaction_descriptor_coupling_strength: f32,
    pub reaction_descriptor_min_factor: f32,
    pub reaction_descriptor_max_factor: f32,
    pub reaction_leakage_strength: f32,
    pub reaction_leakage_max_fraction: f32,
    pub reaction_byproduct_strength: f32,
    pub reaction_byproduct_max_fraction: f32,

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
            transport_energy_cost_scale: 0.02,
            transport_permeability_cost_weight: 1.0,
            transport_composition_cost_weight: 0.5,
            reaction_descriptor_coupling_strength: 0.25,
            reaction_descriptor_min_factor: 0.25,
            reaction_descriptor_max_factor: 1.5,
            reaction_leakage_strength: 0.15,
            reaction_leakage_max_fraction: 0.5,
            reaction_byproduct_strength: 0.08,
            reaction_byproduct_max_fraction: 0.25,

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

    pub fn transport_energy_cost_per_unit(&self, species: usize) -> f32 {
        let descriptor = external_species_descriptor(species);
        let permeability = descriptor.membrane_permeability.clamp(0.0, 1.0);
        let permeability_cost =
            (1.0 - permeability) * self.transport_permeability_cost_weight.max(0.0);
        let composition_cost = descriptor_composition_load(&descriptor.composition)
            * self.transport_composition_cost_weight.max(0.0);
        let cost =
            self.transport_energy_cost_scale.max(0.0) * (permeability_cost + composition_cost);
        if cost.is_finite() { cost } else { 0.0 }
    }

    pub fn reaction_descriptor_factor(
        &self,
        substrate: usize,
        product: usize,
        cofactor: u8,
    ) -> f32 {
        let strength = self.reaction_descriptor_coupling_strength.max(0.0);
        if strength <= f32::EPSILON {
            return 1.0;
        }

        let substrate = internal_species_descriptor(substrate);
        let product = internal_species_descriptor(product);
        let cofactor = (cofactor != stoich::NO_COFACTOR)
            .then(|| internal_species_descriptor(cofactor as usize));

        let mut composition_fit = 0.5
            + 0.5 * descriptor_composition_overlap(&substrate.composition, &product.composition);
        if product.work_coupling >= 0.75 {
            composition_fit = composition_fit.max(0.9);
        }
        let cofactor_bond = cofactor
            .map(|descriptor| descriptor.bond_energy)
            .unwrap_or(0.0);
        let cofactor_redox = cofactor
            .map(|descriptor| {
                descriptor
                    .composition
                    .oxidizing_power
                    .max(descriptor.composition.reducing_power)
            })
            .unwrap_or(0.0);
        let energy_supply = substrate.bond_energy + 0.5 * cofactor_bond + 0.5 * cofactor_redox;
        let energy_need = product.bond_energy * (0.5 + product.work_coupling);
        let energy_fit = if energy_need <= 0.01 {
            1.0
        } else {
            (0.5 + 0.5 * energy_supply / energy_need).clamp(0.25, 1.5)
        };
        let role_fit = if product.work_coupling >= 0.75 && energy_supply > 0.25 {
            1.3
        } else if product.storage_density >= 0.3 && substrate.composition.carbon_backbone > 0.2 {
            1.1
        } else if product.composition.structural_group >= 0.5
            && substrate.composition.carbon_backbone > 0.2
        {
            1.05
        } else if product.composition.signal_group >= 0.5
            && substrate.composition.signal_group < 0.2
        {
            0.75
        } else {
            1.0
        };

        let raw = composition_fit * energy_fit * role_fit;
        let factor = 1.0 + strength * (raw - 1.0);
        let min = self.reaction_descriptor_min_factor.max(0.0);
        let max = self.reaction_descriptor_max_factor.max(min);
        if factor.is_finite() {
            factor.clamp(min, max)
        } else {
            1.0
        }
    }

    pub fn reaction_descriptor_leak_fraction(
        &self,
        substrate: usize,
        product: usize,
        cofactor: u8,
    ) -> f32 {
        let strength = self.reaction_leakage_strength.max(0.0);
        if strength <= f32::EPSILON {
            return 0.0;
        }
        let factor = self.reaction_descriptor_factor(substrate, product, cofactor);
        let leak = strength * (1.0 - factor).max(0.0);
        if leak.is_finite() {
            leak.clamp(0.0, self.reaction_leakage_max_fraction.max(0.0))
        } else {
            0.0
        }
    }

    pub fn reaction_descriptor_byproduct_fraction(
        &self,
        substrate: usize,
        product: usize,
        cofactor: u8,
    ) -> f32 {
        let strength = self.reaction_byproduct_strength.max(0.0);
        if strength <= f32::EPSILON {
            return 0.0;
        }
        let factor = self.reaction_descriptor_factor(substrate, product, cofactor);
        let fraction = strength * (1.0 - factor).max(0.0);
        if fraction.is_finite() {
            fraction.clamp(0.0, self.reaction_byproduct_max_fraction.max(0.0))
        } else {
            0.0
        }
    }

    pub fn reaction_descriptor_byproduct_species(
        &self,
        substrate: usize,
        product: usize,
        cofactor: u8,
    ) -> usize {
        let substrate = internal_species_descriptor(substrate);
        let product = internal_species_descriptor(product);
        let cofactor = (cofactor != stoich::NO_COFACTOR)
            .then(|| internal_species_descriptor(cofactor as usize));

        let signal_load = substrate.composition.signal_group
            + product.composition.signal_group
            + cofactor
                .map(|descriptor| descriptor.composition.signal_group)
                .unwrap_or(0.0);
        if signal_load >= 0.5 {
            return if product.composition.signal_group >= substrate.composition.signal_group {
                EXT_SIGNAL_A
            } else {
                EXT_SIGNAL_B
            };
        }

        let structural_load = substrate.composition.structural_group
            + product.composition.structural_group
            + product.storage_density;
        if structural_load >= 0.7 {
            return EXT_STRUCTURAL;
        }

        EXT_ORGANIC
    }

    pub fn validate_chemistry(&self) -> Result<(), String> {
        for (name, value) in [
            (
                "transport_energy_cost_scale",
                self.transport_energy_cost_scale,
            ),
            (
                "transport_permeability_cost_weight",
                self.transport_permeability_cost_weight,
            ),
            (
                "transport_composition_cost_weight",
                self.transport_composition_cost_weight,
            ),
            (
                "reaction_descriptor_coupling_strength",
                self.reaction_descriptor_coupling_strength,
            ),
            (
                "reaction_descriptor_min_factor",
                self.reaction_descriptor_min_factor,
            ),
            (
                "reaction_descriptor_max_factor",
                self.reaction_descriptor_max_factor,
            ),
            ("reaction_leakage_strength", self.reaction_leakage_strength),
            (
                "reaction_leakage_max_fraction",
                self.reaction_leakage_max_fraction,
            ),
            (
                "reaction_byproduct_strength",
                self.reaction_byproduct_strength,
            ),
            (
                "reaction_byproduct_max_fraction",
                self.reaction_byproduct_max_fraction,
            ),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(format!(
                    "{name} must be finite and nonnegative, got {value}"
                ));
            }
        }
        if self.reaction_descriptor_max_factor < self.reaction_descriptor_min_factor {
            return Err(format!(
                "reaction_descriptor_max_factor ({}) must be >= reaction_descriptor_min_factor ({})",
                self.reaction_descriptor_max_factor, self.reaction_descriptor_min_factor
            ));
        }
        if self.reaction_leakage_max_fraction > 1.0 {
            return Err(format!(
                "reaction_leakage_max_fraction ({}) must be <= 1.0",
                self.reaction_leakage_max_fraction
            ));
        }
        if self.reaction_byproduct_max_fraction > 1.0 {
            return Err(format!(
                "reaction_byproduct_max_fraction ({}) must be <= 1.0",
                self.reaction_byproduct_max_fraction
            ));
        }
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
    fn transport_cost_comes_from_permeability_and_composition() {
        let sim = SimulationConfig::default();
        assert!(sim.transport_energy_cost_per_unit(EXT_ENERGY) > 0.0);
        assert!(
            sim.transport_energy_cost_per_unit(EXT_STRUCTURAL)
                > sim.transport_energy_cost_per_unit(EXT_CARBON)
        );

        let free_transport = SimulationConfig {
            transport_energy_cost_scale: 0.0,
            ..SimulationConfig::default()
        };
        assert_eq!(
            free_transport.transport_energy_cost_per_unit(EXT_STRUCTURAL),
            0.0
        );
    }

    #[test]
    fn reaction_descriptor_factor_uses_bond_energy_and_roles() {
        let sim = SimulationConfig::default();
        let energy_factor =
            sim.reaction_descriptor_factor(EXT_REDUCTANT, EXT_ENERGY, EXT_OXIDANT as u8);
        let toxicity_factor =
            sim.reaction_descriptor_factor(EXT_ENERGY, EXT_CARBON, EXT_OXIDANT as u8);
        let reserve_factor = sim.reaction_descriptor_factor(EXT_CARBON, 7, 0xFF);
        let inactive_factor = sim.reaction_descriptor_factor(EXT_CARBON, 8, 0xFF);
        assert!(energy_factor > toxicity_factor);
        assert!(energy_factor > 1.0);
        assert!(reserve_factor > inactive_factor);

        let disabled = SimulationConfig {
            reaction_descriptor_coupling_strength: 0.0,
            ..SimulationConfig::default()
        };
        assert_eq!(
            disabled.reaction_descriptor_factor(EXT_REDUCTANT, EXT_ENERGY, EXT_OXIDANT as u8),
            1.0
        );
        assert_eq!(
            disabled.reaction_descriptor_leak_fraction(EXT_CARBON, EXT_ENERGY, 0xFF),
            0.0
        );

        let poor_coupling = SimulationConfig {
            reaction_descriptor_coupling_strength: 1.0,
            reaction_leakage_strength: 1.0,
            ..SimulationConfig::default()
        };
        assert!(
            poor_coupling.reaction_descriptor_leak_fraction(EXT_CARBON, EXT_ENERGY, 0xFF) > 0.0
        );
        assert!(
            poor_coupling.reaction_descriptor_byproduct_fraction(EXT_CARBON, EXT_ENERGY, 0xFF)
                > 0.0
        );
        assert_eq!(
            poor_coupling.reaction_descriptor_byproduct_species(EXT_CARBON, EXT_ENERGY, 0xFF),
            EXT_ORGANIC
        );
        assert_eq!(
            poor_coupling.reaction_descriptor_byproduct_species(EXT_CARBON, 7, 0xFF),
            EXT_STRUCTURAL
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

        let mut sim = SimulationConfig {
            transport_energy_cost_scale: -1.0,
            ..Default::default()
        };
        assert!(sim.validate_chemistry().is_err());

        sim.transport_energy_cost_scale = 0.0;
        sim.transport_permeability_cost_weight = f32::NAN;
        assert!(sim.validate_chemistry().is_err());

        let sim = SimulationConfig {
            reaction_descriptor_min_factor: 2.0,
            reaction_descriptor_max_factor: 1.0,
            ..Default::default()
        };
        assert!(sim.validate_chemistry().is_err());

        let sim = SimulationConfig {
            reaction_leakage_strength: f32::NAN,
            ..Default::default()
        };
        assert!(sim.validate_chemistry().is_err());

        let sim = SimulationConfig {
            reaction_leakage_max_fraction: 1.1,
            ..Default::default()
        };
        assert!(sim.validate_chemistry().is_err());

        let sim = SimulationConfig {
            reaction_byproduct_strength: f32::NAN,
            ..Default::default()
        };
        assert!(sim.validate_chemistry().is_err());

        let sim = SimulationConfig {
            reaction_byproduct_max_fraction: 1.1,
            ..Default::default()
        };
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
