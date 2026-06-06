use crate::M_INT;

pub const LIGHT_INTERNAL_SPECIES: usize = M_INT - 1;
pub const NO_COFACTOR: u8 = 0xFF;
pub const STRICT_TEMPLATE_NONE: u8 = 0xFF;
pub const STOICH_TOLERANCE: f32 = 1e-4;
pub const STOICH_STAGE_COUNT: usize = 11;
pub const STOICH_RESERVOIR_COUNT: usize = 7;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StoichEnforcement {
    #[default]
    Off,
    Audit,
    Strict,
}

impl StoichEnforcement {
    pub fn is_enabled(self) -> bool {
        self != Self::Off
    }

    pub fn is_strict(self) -> bool {
        self == Self::Strict
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Audit => "audit",
            Self::Strict => "strict",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SpeciesRole {
    Material,
    Energy,
    Catalyst,
    Structural,
    Inactive,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StoichReservoir {
    BoundaryInput,
    PhotonInput,
    HeatSink,
    AbioticDecaySink,
    RemovedBiomass,
    UnmodeledMaterialSink,
    ClampLoss,
}

impl StoichReservoir {
    pub const ALL: [Self; STOICH_RESERVOIR_COUNT] = [
        Self::BoundaryInput,
        Self::PhotonInput,
        Self::HeatSink,
        Self::AbioticDecaySink,
        Self::RemovedBiomass,
        Self::UnmodeledMaterialSink,
        Self::ClampLoss,
    ];

    pub const fn index(self) -> usize {
        match self {
            Self::BoundaryInput => 0,
            Self::PhotonInput => 1,
            Self::HeatSink => 2,
            Self::AbioticDecaySink => 3,
            Self::RemovedBiomass => 4,
            Self::UnmodeledMaterialSink => 5,
            Self::ClampLoss => 6,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StoichStage {
    BoundaryPriming,
    BoundarySources,
    DiffusionDecay,
    SpatialExchange,
    CellTransport,
    Reactions,
    Maintenance,
    Effectors,
    Division,
    Death,
    Light,
}

impl StoichStage {
    pub const ALL: [Self; STOICH_STAGE_COUNT] = [
        Self::BoundaryPriming,
        Self::BoundarySources,
        Self::DiffusionDecay,
        Self::SpatialExchange,
        Self::CellTransport,
        Self::Reactions,
        Self::Maintenance,
        Self::Effectors,
        Self::Division,
        Self::Death,
        Self::Light,
    ];

    pub const fn index(self) -> usize {
        match self {
            Self::BoundaryPriming => 0,
            Self::BoundarySources => 1,
            Self::DiffusionDecay => 2,
            Self::SpatialExchange => 3,
            Self::CellTransport => 4,
            Self::Reactions => 5,
            Self::Maintenance => 6,
            Self::Effectors => 7,
            Self::Division => 8,
            Self::Death => 9,
            Self::Light => 10,
        }
    }

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::BoundaryPriming => "boundary_priming",
            Self::BoundarySources => "boundary_sources",
            Self::DiffusionDecay => "diffusion_decay",
            Self::SpatialExchange => "spatial_exchange",
            Self::CellTransport => "cell_transport",
            Self::Reactions => "reactions",
            Self::Maintenance => "maintenance",
            Self::Effectors => "effectors",
            Self::Division => "division",
            Self::Death => "death",
            Self::Light => "light",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StoichEventKind {
    BoundaryPrime,
    BoundarySource,
    DiffusionDecay,
    SpatialExchange,
    Transport,
    TransportEnergyCost,
    Effector,
    LegacyReaction,
    BalancedReaction,
    Maintenance,
    ExpressionMaintenance,
    DivisionSplit,
    DeathRemoval,
    LightAvailability,
    StrictRejection,
    ClampLoss,
}

impl StoichEventKind {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::BoundaryPrime => "boundary_prime",
            Self::BoundarySource => "boundary_source",
            Self::DiffusionDecay => "diffusion_decay",
            Self::SpatialExchange => "spatial_exchange",
            Self::Transport => "transport",
            Self::TransportEnergyCost => "transport_energy_cost",
            Self::Effector => "effector",
            Self::LegacyReaction => "legacy_reaction",
            Self::BalancedReaction => "balanced_reaction",
            Self::Maintenance => "maintenance",
            Self::ExpressionMaintenance => "expression_maintenance",
            Self::DivisionSplit => "division_split",
            Self::DeathRemoval => "death_removal",
            Self::LightAvailability => "light_availability",
            Self::StrictRejection => "strict_rejection",
            Self::ClampLoss => "clamp_loss",
        }
    }
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

    pub fn negated(self) -> Self {
        Self {
            c: -self.c,
            h: -self.h,
            o: -self.o,
            s: -self.s,
            redox: -self.redox,
            energy: -self.energy,
        }
    }

    pub fn material_abs_sum(self) -> f32 {
        self.c.abs() + self.h.abs() + self.o.abs() + self.s.abs()
    }

    pub fn total_abs_sum(self) -> f32 {
        self.material_abs_sum() + self.redox.abs() + self.energy.abs()
    }

    pub fn is_near_zero(self, tolerance: f32) -> bool {
        self.total_abs_sum() <= tolerance
    }
}

#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
pub struct StoichStageSummary {
    pub event_count: u64,
    pub strict_rejection_count: u64,
    pub imbalanced_event_count: u64,
    pub gross_model_abs: f32,
    pub gross_reservoir_abs: f32,
    pub gross_residual_abs: f32,
    pub net_model_delta: StoichBudgetDelta,
    pub net_reservoir_delta: StoichBudgetDelta,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct StoichEvent {
    pub stage: StoichStage,
    pub kind: StoichEventKind,
    pub actor_id: u64,
    pub species_index: i16,
    pub template_id: u8,
    pub amount: f32,
    pub model_delta: StoichBudgetDelta,
    pub reservoir_delta: StoichBudgetDelta,
    pub residual: StoichBudgetDelta,
    pub residual_abs: f32,
    pub balanced: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct StoichRecord {
    pub stage: StoichStage,
    pub kind: StoichEventKind,
    pub actor_id: u64,
    pub species_index: i16,
    pub template_id: u8,
    pub amount: f32,
    pub model_delta: StoichBudgetDelta,
    pub reservoir: Option<StoichReservoir>,
    pub reservoir_delta: StoichBudgetDelta,
}

#[derive(Debug, Clone, Copy)]
pub struct LegacyReactionAudit {
    pub substrate: u8,
    pub product: u8,
    pub catalyst: u8,
    pub cofactor: u8,
    pub flux: f32,
    pub actor_id: u64,
}

impl LegacyReactionAudit {
    pub fn new(substrate: u8, product: u8, catalyst: u8, cofactor: u8, flux: f32) -> Self {
        Self {
            substrate,
            product,
            catalyst,
            cofactor,
            flux,
            actor_id: 0,
        }
    }

    pub fn actor(mut self, actor_id: u64) -> Self {
        self.actor_id = actor_id;
        self
    }
}

impl StoichRecord {
    pub fn new(stage: StoichStage, kind: StoichEventKind, amount: f32) -> Self {
        Self {
            stage,
            kind,
            actor_id: 0,
            species_index: -1,
            template_id: STRICT_TEMPLATE_NONE,
            amount,
            model_delta: StoichBudgetDelta::default(),
            reservoir: None,
            reservoir_delta: StoichBudgetDelta::default(),
        }
    }

    pub fn model_delta(mut self, delta: StoichBudgetDelta) -> Self {
        self.model_delta = delta;
        self
    }

    pub fn balancing_reservoir(mut self, reservoir: StoichReservoir) -> Self {
        self.reservoir = Some(reservoir);
        self.reservoir_delta = self.model_delta.negated();
        self
    }

    pub fn actor(mut self, actor_id: u64) -> Self {
        self.actor_id = actor_id;
        self
    }

    pub fn species(mut self, species_index: usize) -> Self {
        self.species_index = species_index as i16;
        self
    }

    pub fn template(mut self, template_id: u8) -> Self {
        self.template_id = template_id;
        self
    }
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct StoichTickLedger {
    pub reaction_count: u64,
    pub active_flux: f32,
    pub imbalanced_reaction_count: u64,
    pub unknown_species_flux: f32,
    pub carbon_to_energy_flux: f32,
    pub reductant_to_energy_flux: f32,
    pub gross_material_abs: f32,
    pub gross_total_abs: f32,
    pub strict_rejection_count: u64,
    pub delta: StoichBudgetDelta,
    pub stage_summaries: [StoichStageSummary; STOICH_STAGE_COUNT],
    pub reservoir_deltas: [StoichBudgetDelta; STOICH_RESERVOIR_COUNT],
    pub events: Vec<StoichEvent>,
}

impl Default for StoichTickLedger {
    fn default() -> Self {
        Self {
            reaction_count: 0,
            active_flux: 0.0,
            imbalanced_reaction_count: 0,
            unknown_species_flux: 0.0,
            carbon_to_energy_flux: 0.0,
            reductant_to_energy_flux: 0.0,
            gross_material_abs: 0.0,
            gross_total_abs: 0.0,
            strict_rejection_count: 0,
            delta: StoichBudgetDelta::default(),
            stage_summaries: [StoichStageSummary::default(); STOICH_STAGE_COUNT],
            reservoir_deltas: [StoichBudgetDelta::default(); STOICH_RESERVOIR_COUNT],
            events: Vec::new(),
        }
    }
}

impl StoichTickLedger {
    pub fn record(&mut self, record: StoichRecord, keep_event: bool) -> StoichBudgetDelta {
        if record.amount < 0.0 || !record.amount.is_finite() {
            return StoichBudgetDelta::default();
        }

        let mut residual = record.model_delta;
        residual.add_delta(record.reservoir_delta);
        let residual_abs = residual.total_abs_sum();
        let balanced = residual.is_near_zero(STOICH_TOLERANCE);
        let stage = &mut self.stage_summaries[record.stage.index()];
        stage.event_count += 1;
        stage.gross_model_abs += record.model_delta.total_abs_sum();
        stage.gross_reservoir_abs += record.reservoir_delta.total_abs_sum();
        stage.gross_residual_abs += residual_abs;
        stage.net_model_delta.add_delta(record.model_delta);
        stage.net_reservoir_delta.add_delta(record.reservoir_delta);
        if !balanced {
            stage.imbalanced_event_count += 1;
        }
        if record.kind == StoichEventKind::StrictRejection {
            stage.strict_rejection_count += 1;
            self.strict_rejection_count += 1;
        }
        if let Some(reservoir) = record.reservoir {
            self.reservoir_deltas[reservoir.index()].add_delta(record.reservoir_delta);
        }
        self.delta.add_delta(record.model_delta);
        self.gross_material_abs += record.model_delta.material_abs_sum();
        self.gross_total_abs +=
            record.model_delta.total_abs_sum() + record.reservoir_delta.total_abs_sum();

        if keep_event {
            self.events.push(StoichEvent {
                stage: record.stage,
                kind: record.kind,
                actor_id: record.actor_id,
                species_index: record.species_index,
                template_id: record.template_id,
                amount: record.amount,
                model_delta: record.model_delta,
                reservoir_delta: record.reservoir_delta,
                residual,
                residual_abs,
                balanced,
            });
        }

        residual
    }

    pub fn record_legacy_reaction(
        &mut self,
        audit: LegacyReactionAudit,
        keep_event: bool,
    ) -> StoichBudgetDelta {
        let LegacyReactionAudit {
            substrate,
            product,
            catalyst,
            cofactor,
            flux,
            actor_id,
        } = audit;
        if flux <= 0.0 || !flux.is_finite() {
            return StoichBudgetDelta::default();
        }

        let delta = legacy_reaction_delta(substrate, product, cofactor, flux);
        let substrate_meta = internal_species(substrate as usize);
        let product_meta = internal_species(product as usize);
        let catalyst_meta = internal_species(catalyst as usize);
        let cofactor_meta = (cofactor != NO_COFACTOR).then(|| internal_species(cofactor as usize));

        self.reaction_count += 1;
        self.active_flux += flux;

        let has_unknown = substrate_meta.role == SpeciesRole::Inactive
            || product_meta.role == SpeciesRole::Inactive
            || catalyst_meta.role == SpeciesRole::Inactive
            || cofactor_meta.is_some_and(|meta| meta.role == SpeciesRole::Inactive);
        if has_unknown {
            self.unknown_species_flux += flux;
        }

        if delta.total_abs_sum() > STOICH_TOLERANCE {
            self.imbalanced_reaction_count += 1;
        }
        if product_meta.role == SpeciesRole::Energy && substrate_meta.material.c > 0.0 {
            self.carbon_to_energy_flux += flux;
        }
        if product_meta.role == SpeciesRole::Energy && substrate_meta.redox_equiv < 0.0 {
            self.reductant_to_energy_flux += flux;
        }

        self.record(
            StoichRecord::new(
                StoichStage::Reactions,
                StoichEventKind::LegacyReaction,
                flux,
            )
            .model_delta(delta)
            .actor(actor_id),
            keep_event,
        );
        delta
    }

    pub fn record_balanced_reaction(
        &mut self,
        template: BalancedReactionTemplate,
        flux: f32,
        actor_id: u64,
        keep_event: bool,
    ) -> StoichBudgetDelta {
        if flux <= 0.0 || !flux.is_finite() {
            return StoichBudgetDelta::default();
        }

        let delta = template.legacy_delta(flux);
        self.reaction_count += 1;
        self.active_flux += flux;
        if template.product == 0 && internal_species(template.substrate as usize).material.c > 0.0 {
            self.carbon_to_energy_flux += flux;
        }
        if template.product == 0 && internal_species(template.substrate as usize).redox_equiv < 0.0
        {
            self.reductant_to_energy_flux += flux;
        }
        self.record(
            StoichRecord::new(
                StoichStage::Reactions,
                StoichEventKind::BalancedReaction,
                flux,
            )
            .model_delta(delta)
            .balancing_reservoir(template.reservoir)
            .actor(actor_id)
            .template(template.id),
            keep_event,
        );
        delta
    }

    pub fn record_strict_rejection(
        &mut self,
        stage: StoichStage,
        amount: f32,
        actor_id: u64,
        template_id: u8,
        keep_event: bool,
    ) {
        self.record(
            StoichRecord::new(stage, StoichEventKind::StrictRejection, amount.max(0.0))
                .actor(actor_id)
                .template(template_id),
            keep_event,
        );
    }

    pub fn has_activity(&self) -> bool {
        self.reaction_count > 0
            || self.delta.total_abs_sum() > 0.0
            || self.strict_rejection_count > 0
    }
}

#[derive(Debug, Clone, serde::Serialize)]
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
    pub strict_rejection_count: u64,
    pub delta: StoichBudgetDelta,
    pub stage_summaries: [StoichStageSummary; STOICH_STAGE_COUNT],
    pub reservoir_deltas: [StoichBudgetDelta; STOICH_RESERVOIR_COUNT],
}

impl Default for StoichRunLedger {
    fn default() -> Self {
        Self {
            tick_count: 0,
            reaction_count: 0,
            active_flux: 0.0,
            imbalanced_reaction_count: 0,
            unknown_species_flux: 0.0,
            carbon_to_energy_flux: 0.0,
            reductant_to_energy_flux: 0.0,
            gross_material_abs: 0.0,
            gross_total_abs: 0.0,
            strict_rejection_count: 0,
            delta: StoichBudgetDelta::default(),
            stage_summaries: [StoichStageSummary::default(); STOICH_STAGE_COUNT],
            reservoir_deltas: [StoichBudgetDelta::default(); STOICH_RESERVOIR_COUNT],
        }
    }
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
        self.strict_rejection_count += tick.strict_rejection_count;
        self.delta.add_delta(tick.delta);
        for (run, tick) in self.stage_summaries.iter_mut().zip(tick.stage_summaries) {
            run.event_count += tick.event_count;
            run.strict_rejection_count += tick.strict_rejection_count;
            run.imbalanced_event_count += tick.imbalanced_event_count;
            run.gross_model_abs += tick.gross_model_abs;
            run.gross_reservoir_abs += tick.gross_reservoir_abs;
            run.gross_residual_abs += tick.gross_residual_abs;
            run.net_model_delta.add_delta(tick.net_model_delta);
            run.net_reservoir_delta.add_delta(tick.net_reservoir_delta);
        }
        for (run, tick) in self.reservoir_deltas.iter_mut().zip(tick.reservoir_deltas) {
            run.add_delta(tick);
        }
    }

    pub fn material_abs_sum(&self) -> f32 {
        self.gross_material_abs
    }

    pub fn net_material_abs_sum(&self) -> f32 {
        self.delta.material_abs_sum()
    }

    pub fn total_abs_sum(&self) -> f32 {
        self.gross_total_abs
    }

    pub fn gross_residual_abs(&self) -> f32 {
        self.stage_summaries
            .iter()
            .map(|stage| stage.gross_residual_abs)
            .sum()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum BalancedTemplateKind {
    PhototrophyEnergy,
    PhototrophyOxidant,
    CarbonFixation,
    EnzymeSynthesisA,
    EnzymeSynthesisB,
    ChemolithotrophyEnergy,
    OxidantProcessing,
    ReserveBurn,
    AnaerobicEnergy,
    OxidantToxicity,
}

#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct BalancedReactionTemplate {
    pub id: u8,
    pub kind: BalancedTemplateKind,
    pub substrate: u8,
    pub product: u8,
    pub catalyst: u8,
    pub cofactor: u8,
    pub reservoir: StoichReservoir,
}

impl BalancedReactionTemplate {
    pub fn legacy_delta(self, flux: f32) -> StoichBudgetDelta {
        legacy_reaction_delta(self.substrate, self.product, self.cofactor, flux)
    }

    pub fn matches(self, substrate: u8, product: u8, catalyst: u8, cofactor: u8) -> bool {
        self.substrate == substrate
            && self.product == product
            && self.catalyst == catalyst
            && self.cofactor == cofactor
    }
}

pub const BALANCED_REACTION_TEMPLATES: [BalancedReactionTemplate; 13] = [
    template(
        0,
        BalancedTemplateKind::PhototrophyEnergy,
        3,
        0,
        15,
        NO_COFACTOR,
        StoichReservoir::PhotonInput,
    ),
    template(
        1,
        BalancedTemplateKind::PhototrophyOxidant,
        3,
        1,
        15,
        NO_COFACTOR,
        StoichReservoir::PhotonInput,
    ),
    template(
        2,
        BalancedTemplateKind::CarbonFixation,
        3,
        4,
        0,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        3,
        BalancedTemplateKind::EnzymeSynthesisA,
        3,
        5,
        6,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        4,
        BalancedTemplateKind::EnzymeSynthesisB,
        3,
        6,
        5,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        5,
        BalancedTemplateKind::ChemolithotrophyEnergy,
        2,
        0,
        5,
        1,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        6,
        BalancedTemplateKind::OxidantProcessing,
        1,
        4,
        0,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        7,
        BalancedTemplateKind::ReserveBurn,
        7,
        0,
        5,
        NO_COFACTOR,
        StoichReservoir::HeatSink,
    ),
    template(
        8,
        BalancedTemplateKind::AnaerobicEnergy,
        2,
        0,
        5,
        NO_COFACTOR,
        StoichReservoir::HeatSink,
    ),
    template(
        9,
        BalancedTemplateKind::OxidantToxicity,
        0,
        3,
        1,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        10,
        BalancedTemplateKind::CarbonFixation,
        3,
        4,
        0,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        11,
        BalancedTemplateKind::EnzymeSynthesisA,
        3,
        5,
        6,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
    template(
        12,
        BalancedTemplateKind::EnzymeSynthesisB,
        3,
        6,
        5,
        NO_COFACTOR,
        StoichReservoir::UnmodeledMaterialSink,
    ),
];

pub const fn template(
    id: u8,
    kind: BalancedTemplateKind,
    substrate: u8,
    product: u8,
    catalyst: u8,
    cofactor: u8,
    reservoir: StoichReservoir,
) -> BalancedReactionTemplate {
    BalancedReactionTemplate {
        id,
        kind,
        substrate,
        product,
        catalyst,
        cofactor,
        reservoir,
    }
}

pub fn balanced_template_for_reaction(
    substrate: u8,
    product: u8,
    catalyst: u8,
    cofactor: u8,
) -> Option<BalancedReactionTemplate> {
    BALANCED_REACTION_TEMPLATES
        .iter()
        .copied()
        .find(|template| template.matches(substrate, product, catalyst, cofactor))
}

pub fn template_by_id(id: u8) -> Option<BalancedReactionTemplate> {
    BALANCED_REACTION_TEMPLATES
        .iter()
        .copied()
        .find(|template| template.id == id)
}

pub fn legacy_reaction_delta(
    substrate: u8,
    product: u8,
    cofactor: u8,
    flux: f32,
) -> StoichBudgetDelta {
    let mut delta = StoichBudgetDelta::default();
    delta.add_species(internal_species(substrate as usize), -flux);
    delta.add_species(internal_species(product as usize), flux);

    if cofactor != NO_COFACTOR {
        delta.add_species(internal_species(cofactor as usize), -0.5 * flux);
    }
    delta
}

pub fn external_delta(species: usize, amount: f32) -> StoichBudgetDelta {
    let mut delta = StoichBudgetDelta::default();
    delta.add_species(external_species(species), amount);
    delta
}

pub fn internal_delta(species: usize, amount: f32) -> StoichBudgetDelta {
    let mut delta = StoichBudgetDelta::default();
    delta.add_species(internal_species(species), amount);
    delta
}

pub fn transfer_delta(
    ext_species: usize,
    int_species: usize,
    amount_into_internal: f32,
) -> StoichBudgetDelta {
    let mut delta = StoichBudgetDelta::default();
    delta.add_species(external_species(ext_species), -amount_into_internal);
    delta.add_species(internal_species(int_species), amount_into_internal);
    delta
}

pub fn internal_pool_delta(internal: &[f32; M_INT], sign: f32) -> StoichBudgetDelta {
    let mut delta = StoichBudgetDelta::default();
    for (index, amount) in internal.iter().enumerate() {
        if index == LIGHT_INTERNAL_SPECIES || *amount == 0.0 {
            continue;
        }
        delta.add_species(internal_species(index), sign * *amount);
    }
    delta
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
        let delta = ledger.record_legacy_reaction(
            LegacyReactionAudit::new(4, 7, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 1.0),
            false,
        );
        assert_eq!(delta.material_abs_sum(), 0.0);
        assert_eq!(ledger.imbalanced_reaction_count, 0);
    }

    #[test]
    fn carbon_to_energy_reaction_is_reported_as_imbalance() {
        let mut ledger = StoichTickLedger::default();
        let delta = ledger.record_legacy_reaction(
            LegacyReactionAudit::new(3, 0, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 2.0),
            false,
        );
        assert_eq!(delta.c, -2.0);
        assert_eq!(delta.energy, 2.0);
        assert_eq!(ledger.carbon_to_energy_flux, 2.0);
        assert_eq!(ledger.imbalanced_reaction_count, 1);
    }

    #[test]
    fn cofactor_consumption_is_included_in_delta() {
        let mut ledger = StoichTickLedger::default();
        let delta = ledger.record_legacy_reaction(LegacyReactionAudit::new(2, 0, 5, 1, 2.0), false);
        assert_eq!(delta.s, -2.0);
        assert_eq!(delta.o, -2.0);
        assert_eq!(delta.energy, 2.0);
        assert_eq!(ledger.reductant_to_energy_flux, 2.0);
    }

    #[test]
    fn light_catalyst_contributes_no_material_delta() {
        let mut ledger = StoichTickLedger::default();
        let delta = ledger.record_legacy_reaction(
            LegacyReactionAudit::new(4, 7, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 3.0),
            false,
        );
        assert_eq!(delta.total_abs_sum(), 0.0);
        assert_eq!(ledger.unknown_species_flux, 0.0);
    }

    #[test]
    fn inactive_species_is_reported_as_unknown_flux() {
        let mut ledger = StoichTickLedger::default();
        let delta = ledger
            .record_legacy_reaction(LegacyReactionAudit::new(8, 0, 5, NO_COFACTOR, 1.0), false);
        assert_eq!(delta.energy, 1.0);
        assert_eq!(ledger.unknown_species_flux, 1.0);
    }

    #[test]
    fn gross_imbalance_accumulates_even_when_net_cancels() {
        let mut run = StoichRunLedger::default();

        let mut tick_a = StoichTickLedger::default();
        tick_a.record_legacy_reaction(
            LegacyReactionAudit::new(3, 0, LIGHT_INTERNAL_SPECIES as u8, NO_COFACTOR, 1.0),
            false,
        );
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

    #[test]
    fn reservoir_balancing_closes_event_residual() {
        let mut ledger = StoichTickLedger::default();
        let delta = external_delta(3, 2.0);

        ledger.record(
            StoichRecord::new(
                StoichStage::BoundarySources,
                StoichEventKind::BoundarySource,
                2.0,
            )
            .model_delta(delta)
            .balancing_reservoir(StoichReservoir::BoundaryInput),
            true,
        );

        assert_eq!(ledger.events.len(), 1);
        assert!(ledger.events[0].balanced);
        assert_eq!(ledger.events[0].residual.total_abs_sum(), 0.0);
    }

    #[test]
    fn every_balanced_template_closes_with_its_reservoir() {
        for template in BALANCED_REACTION_TEMPLATES {
            let delta = template.legacy_delta(1.0);
            let mut residual = delta;
            residual.add_delta(delta.negated());
            assert!(
                residual.is_near_zero(STOICH_TOLERANCE),
                "template {:?} did not close",
                template.kind
            );
        }
    }

    #[test]
    fn strict_mode_rejects_unknown_legacy_template() {
        assert!(balanced_template_for_reaction(4, 0, 5, NO_COFACTOR).is_none());
        assert!(balanced_template_for_reaction(3, 0, 15, NO_COFACTOR).is_some());
    }
}
