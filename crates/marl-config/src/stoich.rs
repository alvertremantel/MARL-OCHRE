use crate::M_INT;

pub const LIGHT_INTERNAL_SPECIES: usize = M_INT - 1;
pub const NO_COFACTOR: u8 = 0xFF;

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SpeciesRole {
    Material,
    Energy,
    Catalyst,
    Structural,
    Inactive,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize)]
pub struct MaterialVector {
    pub c: f32,
    pub h: f32,
    pub o: f32,
    pub s: f32,
}

impl MaterialVector {
    pub const ZERO: Self = Self {
        c: 0.0,
        h: 0.0,
        o: 0.0,
        s: 0.0,
    };

    pub const fn new(c: f32, h: f32, o: f32, s: f32) -> Self {
        Self { c, h, o, s }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct SpeciesMetadata {
    pub index: usize,
    pub name: &'static str,
    pub role: SpeciesRole,
    pub material: MaterialVector,
    pub redox_equiv: f32,
    pub energy_equiv: f32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize)]
pub struct StoichBudgetDelta {
    pub c: f32,
    pub h: f32,
    pub o: f32,
    pub s: f32,
    pub redox: f32,
    pub energy: f32,
}

impl StoichBudgetDelta {
    pub fn add_species(&mut self, species: SpeciesMetadata, amount: f32) {
        self.c += species.material.c * amount;
        self.h += species.material.h * amount;
        self.o += species.material.o * amount;
        self.s += species.material.s * amount;
        self.redox += species.redox_equiv * amount;
        self.energy += species.energy_equiv * amount;
    }

    pub fn add_delta(&mut self, other: Self) {
        self.c += other.c;
        self.h += other.h;
        self.o += other.o;
        self.s += other.s;
        self.redox += other.redox;
        self.energy += other.energy;
    }

    pub fn material_abs_sum(self) -> f32 {
        self.c.abs() + self.h.abs() + self.o.abs() + self.s.abs()
    }

    pub fn total_abs_sum(self) -> f32 {
        self.material_abs_sum() + self.redox.abs() + self.energy.abs()
    }
}

#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
pub struct StoichTickLedger {
    pub reaction_count: u64,
    pub active_flux: f32,
    pub imbalanced_reaction_count: u64,
    pub unknown_species_flux: f32,
    pub carbon_to_energy_flux: f32,
    pub reductant_to_energy_flux: f32,
    pub gross_material_abs: f32,
    pub gross_total_abs: f32,
    pub delta: StoichBudgetDelta,
}

impl StoichTickLedger {
    pub fn record_legacy_reaction(
        &mut self,
        substrate: u8,
        product: u8,
        catalyst: u8,
        cofactor: u8,
        flux: f32,
    ) -> StoichBudgetDelta {
        if flux <= 0.0 || !flux.is_finite() {
            return StoichBudgetDelta::default();
        }

        let mut delta = StoichBudgetDelta::default();
        let substrate_meta = internal_species(substrate as usize);
        let product_meta = internal_species(product as usize);
        let catalyst_meta = internal_species(catalyst as usize);

        delta.add_species(substrate_meta, -flux);
        delta.add_species(product_meta, flux);

        if cofactor != NO_COFACTOR {
            delta.add_species(internal_species(cofactor as usize), -0.5 * flux);
        }

        self.reaction_count += 1;
        self.active_flux += flux;
        self.gross_material_abs += delta.material_abs_sum();
        self.gross_total_abs += delta.total_abs_sum();
        self.delta.add_delta(delta);

        let has_unknown = substrate_meta.role == SpeciesRole::Inactive
            || product_meta.role == SpeciesRole::Inactive
            || catalyst_meta.role == SpeciesRole::Inactive
            || (cofactor != NO_COFACTOR
                && internal_species(cofactor as usize).role == SpeciesRole::Inactive);
        if has_unknown {
            self.unknown_species_flux += flux;
        }

        if delta.total_abs_sum() > 1e-6 {
            self.imbalanced_reaction_count += 1;
        }
        if product_meta.role == SpeciesRole::Energy && substrate_meta.material.c > 0.0 {
            self.carbon_to_energy_flux += flux;
        }
        if product_meta.role == SpeciesRole::Energy && substrate_meta.redox_equiv < 0.0 {
            self.reductant_to_energy_flux += flux;
        }

        delta
    }

    pub fn has_activity(self) -> bool {
        self.reaction_count > 0 || self.delta.total_abs_sum() > 0.0
    }
}

#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
pub struct StoichRunLedger {
    pub tick_count: u64,
    pub reaction_count: u64,
    pub active_flux: f32,
    pub imbalanced_reaction_count: u64,
    pub unknown_species_flux: f32,
    pub carbon_to_energy_flux: f32,
    pub reductant_to_energy_flux: f32,
    pub gross_material_abs: f32,
    pub gross_total_abs: f32,
    pub delta: StoichBudgetDelta,
}

impl StoichRunLedger {
    pub fn add_tick(&mut self, tick: &StoichTickLedger) {
        self.tick_count += 1;
        self.reaction_count += tick.reaction_count;
        self.active_flux += tick.active_flux;
        self.imbalanced_reaction_count += tick.imbalanced_reaction_count;
        self.unknown_species_flux += tick.unknown_species_flux;
        self.carbon_to_energy_flux += tick.carbon_to_energy_flux;
        self.reductant_to_energy_flux += tick.reductant_to_energy_flux;
        self.gross_material_abs += tick.gross_material_abs;
        self.gross_total_abs += tick.gross_total_abs;
        self.delta.add_delta(tick.delta);
    }

    pub fn material_abs_sum(self) -> f32 {
        self.gross_material_abs
    }

    pub fn net_material_abs_sum(self) -> f32 {
        self.delta.material_abs_sum()
    }

    pub fn total_abs_sum(self) -> f32 {
        self.gross_total_abs
    }
}

pub const fn external_species(index: usize) -> SpeciesMetadata {
    match index {
        1 => species(
            index,
            "oxidant",
            SpeciesRole::Material,
            MaterialVector::new(0.0, 0.0, 2.0, 0.0),
            4.0,
            0.0,
        ),
        2 => species(
            index,
            "reductant",
            SpeciesRole::Material,
            MaterialVector::new(0.0, 2.0, 0.0, 1.0),
            -4.0,
            0.0,
        ),
        3 => species(
            index,
            "carbon",
            SpeciesRole::Material,
            MaterialVector::new(1.0, 0.0, 2.0, 0.0),
            4.0,
            0.0,
        ),
        4 => species(
            index,
            "organic",
            SpeciesRole::Material,
            MaterialVector::new(1.0, 2.0, 1.0, 0.0),
            0.0,
            0.0,
        ),
        7 => species(
            index,
            "structural_eps",
            SpeciesRole::Structural,
            MaterialVector::new(1.0, 2.0, 1.0, 0.0),
            0.0,
            0.0,
        ),
        _ => inactive_species(index),
    }
}

pub const fn internal_species(index: usize) -> SpeciesMetadata {
    match index {
        0 => species(
            index,
            "energy",
            SpeciesRole::Energy,
            MaterialVector::ZERO,
            0.0,
            1.0,
        ),
        1 => species(
            index,
            "oxidant",
            SpeciesRole::Material,
            MaterialVector::new(0.0, 0.0, 2.0, 0.0),
            4.0,
            0.0,
        ),
        2 => species(
            index,
            "reductant",
            SpeciesRole::Material,
            MaterialVector::new(0.0, 2.0, 0.0, 1.0),
            -4.0,
            0.0,
        ),
        3 => species(
            index,
            "carbon",
            SpeciesRole::Material,
            MaterialVector::new(1.0, 0.0, 2.0, 0.0),
            4.0,
            0.0,
        ),
        4 => species(
            index,
            "organic",
            SpeciesRole::Material,
            MaterialVector::new(1.0, 2.0, 1.0, 0.0),
            0.0,
            0.0,
        ),
        5 => species(
            index,
            "enzyme_a",
            SpeciesRole::Structural,
            MaterialVector::new(5.0, 7.0, 2.0, 1.0),
            0.0,
            0.0,
        ),
        6 => species(
            index,
            "enzyme_b",
            SpeciesRole::Structural,
            MaterialVector::new(5.0, 7.0, 2.0, 1.0),
            0.0,
            0.0,
        ),
        7 => species(
            index,
            "carbon_reserve",
            SpeciesRole::Structural,
            MaterialVector::new(1.0, 2.0, 1.0, 0.0),
            0.0,
            0.0,
        ),
        LIGHT_INTERNAL_SPECIES => species(
            index,
            "light",
            SpeciesRole::Catalyst,
            MaterialVector::ZERO,
            0.0,
            0.0,
        ),
        _ => inactive_species(index),
    }
}

const fn species(
    index: usize,
    name: &'static str,
    role: SpeciesRole,
    material: MaterialVector,
    redox_equiv: f32,
    energy_equiv: f32,
) -> SpeciesMetadata {
    SpeciesMetadata {
        index,
        name,
        role,
        material,
        redox_equiv,
        energy_equiv,
    }
}

const fn inactive_species(index: usize) -> SpeciesMetadata {
    SpeciesMetadata {
        index,
        name: "inactive",
        role: SpeciesRole::Inactive,
        material: MaterialVector::ZERO,
        redox_equiv: 0.0,
        energy_equiv: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{M_INT, S_EXT};

    #[test]
    fn every_current_slot_has_metadata() {
        for index in 0..S_EXT {
            assert_eq!(external_species(index).index, index);
        }
        for index in 0..M_INT {
            assert_eq!(internal_species(index).index, index);
        }
        assert_eq!(internal_species(0).role, SpeciesRole::Energy);
        assert_eq!(
            internal_species(LIGHT_INTERNAL_SPECIES).role,
            SpeciesRole::Catalyst
        );
    }

    #[test]
    fn balanced_material_neutral_reaction_has_no_material_delta() {
        let mut ledger = StoichTickLedger::default();
        let delta =
            ledger.record_legacy_reaction(4, 7, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 1.0);
        assert_eq!(delta.material_abs_sum(), 0.0);
        assert_eq!(ledger.imbalanced_reaction_count, 0);
    }

    #[test]
    fn carbon_to_energy_reaction_is_reported_as_imbalance() {
        let mut ledger = StoichTickLedger::default();
        let delta =
            ledger.record_legacy_reaction(3, 0, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 2.0);
        assert_eq!(delta.c, -2.0);
        assert_eq!(delta.energy, 2.0);
        assert_eq!(ledger.carbon_to_energy_flux, 2.0);
        assert_eq!(ledger.imbalanced_reaction_count, 1);
    }

    #[test]
    fn cofactor_consumption_is_included_in_delta() {
        let mut ledger = StoichTickLedger::default();
        let delta = ledger.record_legacy_reaction(2, 0, 5, 1, 2.0);
        assert_eq!(delta.s, -2.0);
        assert_eq!(delta.o, -2.0);
        assert_eq!(delta.energy, 2.0);
        assert_eq!(ledger.reductant_to_energy_flux, 2.0);
    }

    #[test]
    fn light_catalyst_contributes_no_material_delta() {
        let mut ledger = StoichTickLedger::default();
        let delta =
            ledger.record_legacy_reaction(4, 7, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 3.0);
        assert_eq!(delta.total_abs_sum(), 0.0);
        assert_eq!(ledger.unknown_species_flux, 0.0);
    }

    #[test]
    fn inactive_species_is_reported_as_unknown_flux() {
        let mut ledger = StoichTickLedger::default();
        let delta = ledger.record_legacy_reaction(8, 0, 5, NO_COFACTOR, 1.0);
        assert_eq!(delta.energy, 1.0);
        assert_eq!(ledger.unknown_species_flux, 1.0);
    }

    #[test]
    fn gross_imbalance_accumulates_even_when_net_cancels() {
        let mut run = StoichRunLedger::default();

        let mut tick_a = StoichTickLedger::default();
        tick_a.record_legacy_reaction(3, 0, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 1.0);
        run.add_tick(&tick_a);

        let mut tick_b = StoichTickLedger::default();
        tick_b.delta.c = 1.0;
        tick_b.delta.o = 2.0;
        tick_b.gross_material_abs = 3.0;
        tick_b.gross_total_abs = 3.0;
        run.add_tick(&tick_b);

        assert_eq!(run.delta.c, 0.0);
        assert_eq!(run.delta.o, 0.0);
        assert_eq!(run.net_material_abs_sum(), 0.0);
        assert_eq!(run.material_abs_sum(), 6.0);
        assert_eq!(run.total_abs_sum(), 11.0);
    }
}
