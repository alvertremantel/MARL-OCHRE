//! Cell agent module — the evolvable catalytic core of the simulation.
//!
//! Each cell is a small biochemical computer: it senses its environment
//! (receptors), exchanges chemicals across its membrane (transporters),
//! runs internal catalytic reactions, secretes products (effectors), and
//! makes life/death/division decisions based on energy levels.
//!
//! The entire behavior of a cell is encoded in its **Ruleset** — a
//! collection of parameter arrays that define receptor sensitivities,
//! transport rates, reaction kinetics, etc. Rulesets mutate during
//! division, so evolution acts on these parameters directly. There is
//! no fitness function; selection emerges from the chemistry.
//!
//! ## The 5-Phase Cell Tick
//!
//! Each simulation tick, every cell runs through five phases in order:
//!
//! 1. **Receptor pass** — sense external chemical concentrations using
//!    Hill-function activations. Transporters can evolve receptor gates
//!    that amplify or suppress membrane flux from these activations.
//!
//! 2. **Transport pass** — move chemicals across the cell membrane.
//!    Uptake brings external species into the cell; secretion pushes
//!    internal species out. Both use Michaelis-Menten saturation.
//!
//! 3. **Intracellular reactions** — the catalytic network. Each reaction
//!    consumes a substrate and produces a product, catalyzed by a third
//!    species. Michaelis-Menten kinetics with an epsilon background rate
//!    (see `research-bootstrapping.md` for why epsilon matters).
//!
//! 4. **Effector pass** — secrete internal species into the field when
//!    concentrations exceed a threshold. Quiescent cells skip this.
//!
//! 5. **Fate decision** — based on energy (internal[0]):
//!    - Energy > division_energy → divide
//!    - Energy < death_energy → die
//!    - Energy < quiescence_energy → go dormant (skip effectors)
//!    - Otherwise → active

use marl_config::stoich::{
    BALANCED_REACTION_TEMPLATES, LegacyReactionAudit, StoichEventKind, StoichRecord,
    StoichReservoir, StoichStage, StoichTickLedger, balanced_template_for_reaction, internal_delta,
    transfer_delta,
};
use marl_config::*;

const LIGHT_SPECIES: usize = M_INT - 1;
const CHEMICAL_INT_SPECIES: usize = M_INT - 1;
const TRANSPORT_GATE_ACTIVATION_LIMIT: f32 = 1.0;
const TRANSPORT_GATE_WEIGHT_LIMIT: f32 = 3.0;
const TRANSPORT_GATE_FACTOR_MAX: f32 = 4.0;

/// Receptor parameters: Hill-function sensor for an external chemical.
///
/// The Hill equation models cooperative binding:
///   activation = gain * c^n / (k_half^n + c^n)
///
/// - `k_half`: the concentration at which activation reaches 50%
/// - `n_hill`: cooperativity coefficient (1 = hyperbolic, >1 = sigmoidal switch)
/// - `gain`: maximum output scaling
#[derive(Clone, Debug)]
pub struct ReceptorParams {
    pub k_half: f32,
    pub n_hill: f32,
    pub gain: f32,
}

/// Membrane transporter: moves one chemical species between the external
/// field and the cell's internal pool.
///
/// Each transporter is specific to one (ext_species, int_species) pair.
/// Uptake and secretion rates are independent — a transporter can do both,
/// creating a net flux direction based on concentration gradients. The
/// optional receptor gate multiplies both rates by
/// `clamp(1 + gate_weight * activation[gate_receptor], 0, max)`.
/// A zero gate weight preserves unconditional transport.
#[derive(Clone, Debug)]
pub struct TransportParams {
    pub uptake_rate: f32,
    pub secrete_rate: f32,
    pub ext_species: u8,
    pub int_species: u8,
    pub gate_receptor: u8,
    pub gate_weight: f32,
}

/// A single catalytic reaction in the cell's metabolic network.
///
/// Models enzyme kinetics: substrate → product, catalyzed by a third
/// internal species. Rate follows Michaelis-Menten:
///   rate = v_max * [S]/(k_m + [S]) * (epsilon + [catalyst]/(k_cat + [catalyst]))
///
/// The epsilon term ensures a tiny background rate even without catalyst,
/// preventing evolutionary dead ends where a useful reaction can never
/// start because the catalyst doesn't exist yet.
///
/// Optional cofactor: if cofactor != 0xFF, the reaction also requires
/// (and partially consumes) a second internal species.
#[derive(Clone, Debug)]
pub struct Reaction {
    pub substrate: u8, // consumed internal species
    pub product: u8,   // produced internal species
    pub catalyst: u8,  // enzyme/catalyst species
    pub cofactor: u8,  // second required species, 0xFF = none
    pub k_m: f32,      // Michaelis constant for substrate
    pub v_max: f32,    // maximum reaction rate
    pub k_cat: f32,    // half-saturation for catalyst
}

/// Effector: conditional secretion of an internal species into the field.
///
/// When internal[int_species] exceeds `threshold`, the cell secretes at
/// the given rate. This creates emergent signaling — cells that accumulate
/// waste products or metabolic byproducts automatically release them,
/// making those chemicals available to neighboring cells.
#[derive(Clone, Debug)]
pub struct EffectorParams {
    pub threshold: f32,
    pub rate: f32,
    pub int_species: u8,
    pub ext_species: u8,
}

/// Energy thresholds and cell cycle timing that determine cell fate.
///
/// There is no explicit fitness function — a cell that can accumulate
/// energy above `division_energy` AND sustain the costly preparation
/// phase will reproduce. One that can't maintain energy above
/// `death_energy` will die. Selection is entirely thermodynamic.
///
/// The division prep phase models DNA replication, organelle duplication,
/// and membrane synthesis — biologically expensive processes that take
/// time. Cells can evolve shorter prep times but pay exponentially
/// more energy to rush.
#[derive(Clone, Debug)]
pub struct FateParams {
    pub division_energy: f32,
    pub death_energy: f32,
    pub quiescence_energy: f32,
    /// How many ticks the cell spends in division prep. Evolvable.
    /// Lower = faster division but higher energy cost during prep.
    pub division_prep_ticks: f32,
}

/// The complete evolvable genotype of a cell.
///
/// Contains all parameters that define the cell's behavior: how it senses,
/// what it transports, which reactions it catalyzes, when it secretes,
/// and at what energy levels it divides/dies. All of these mutate during
/// reproduction, so the ruleset is the unit of evolution.
///
/// Size: ~314 bytes (8 receptors + 8 transporters + 16 reactions +
/// 8 effectors + fate params + mutation/HGT rates).
#[derive(Clone, Debug)]
pub struct Ruleset {
    pub receptors: [ReceptorParams; S_RECEPTORS],
    pub transport: [TransportParams; S_TRANSPORTERS],
    pub reactions: [Reaction; R_MAX],
    pub effectors: [EffectorParams; S_EFFECTORS],
    pub fate: FateParams,
    /// Probability of accepting a horizontal gene transfer event.
    /// Evolvable — evolution can select for or against gene sharing.
    pub hgt_propensity: f32,
    /// Per-parameter probability of point mutation during division.
    /// Also evolvable (meta-evolution: evolution of evolvability).
    pub mutation_rate: f32,
}

/// Events produced by a cell tick — the cell's "decision" for this timestep.
#[derive(Clone, Debug)]
pub enum CellEvent {
    None,
    Division,
    Death,
    Quiescence,
}

/// The complete state of a single cell at one point in time.
///
/// Each cell occupies exactly one voxel in the 3D grid. Its position,
/// internal chemical concentrations, and ruleset together define
/// everything about it. The lineage_id is a random tag assigned at
/// birth for tracking evolutionary lineages.
#[derive(Clone, Debug)]
pub struct CellState {
    pub pos: [u16; 3],
    pub lineage_id: u64,
    pub age: u32,
    pub internal: [f32; M_INT],
    pub ruleset: Ruleset,
    pub quiescent: bool,
    /// Which original starter metabolism this cell descends from.
    /// 0 = phototroph, 1 = chemolithotroph, 2 = anaerobe.
    /// Inherited at division, never mutated — permanent ancestral marker.
    pub starter_type: u8,
    /// Division prep countdown. 0 = not in prep. When > 0, the cell is
    /// preparing to divide and pays extra maintenance. Reaches 0 → divide.
    pub prep_remaining: u16,
}

fn finite_nonnegative(value: f32) -> f32 {
    if value.is_finite() {
        value.max(0.0)
    } else {
        0.0
    }
}

fn transport_gate_factor(tp: &TransportParams, activation: &[f32; S_RECEPTORS]) -> f32 {
    if !tp.gate_weight.is_finite() {
        return 1.0;
    }
    let weight = tp
        .gate_weight
        .clamp(-TRANSPORT_GATE_WEIGHT_LIMIT, TRANSPORT_GATE_WEIGHT_LIMIT);
    if weight.abs() <= f32::EPSILON {
        return 1.0;
    }
    let Some(signal) = activation.get(tp.gate_receptor as usize).copied() else {
        return 1.0;
    };
    let factor = 1.0 + weight * signal;
    if factor.is_finite() {
        factor.clamp(0.0, TRANSPORT_GATE_FACTOR_MAX)
    } else {
        1.0
    }
}

impl CellState {
    /// Run one complete cell update tick. Returns field deltas and an event.
    /// This is the 5-phase update:
    ///   Phase 1: Receptor pass (read external field, compute Hill activations)
    ///   Phase 2: Transport pass (membrane crossing, external <-> internal)
    ///   Phase 3: Intracellular reactions (catalytic network with epsilon background)
    ///   Phase 4: Effector pass (internal -> external secretion)
    ///   Phase 5: Fate decision (energy thresholds)
    pub fn tick(
        &mut self,
        ext_conc: &[f32; S_EXT],
        light: f32,
        sim: &SimulationConfig,
    ) -> ([f32; S_EXT], CellEvent) {
        self.tick_with_stoich(ext_conc, light, sim, None, false, false)
    }

    pub fn tick_with_stoich(
        &mut self,
        ext_conc: &[f32; S_EXT],
        light: f32,
        sim: &SimulationConfig,
        mut stoich: Option<&mut StoichTickLedger>,
        record_full_stoich: bool,
        keep_stoich_events: bool,
    ) -> ([f32; S_EXT], CellEvent) {
        let mut field_deltas = [0.0f32; S_EXT];
        // Cell tick runs once per full tick
        let dt = sim.dt;

        // === PHASE 1: RECEPTOR PASS ===
        // Activations modulate transport rates through each transporter's
        // evolvable receptor gate.
        let mut activation = [0.0f32; S_RECEPTORS];
        for i in 0..S_RECEPTORS {
            let r = &self.ruleset.receptors[i];
            if i < S_EXT
                && r.k_half.is_finite()
                && r.n_hill.is_finite()
                && r.gain.is_finite()
                && r.gain.abs() > sim.active_reaction_threshold
            {
                let c = ext_conc[i];
                if !c.is_finite() {
                    continue;
                }
                // Guard against negative/zero base in powf
                let k = r.k_half.max(1e-6);
                let n = r
                    .n_hill
                    .clamp(sim.hill_exponent_clamp_low, sim.hill_exponent_clamp_high);
                let kn = k.powf(n);
                let cn = c.max(0.0).powf(n);
                let value = r.gain * cn / (kn + cn + 1e-9);
                if value.is_finite() {
                    activation[i] = value.clamp(
                        -TRANSPORT_GATE_ACTIVATION_LIMIT,
                        TRANSPORT_GATE_ACTIVATION_LIMIT,
                    );
                }
            }
        }

        // === PHASE 2: TRANSPORT PASS ===
        let mut valid_transport = [false; S_TRANSPORTERS];
        let mut transport_ext = [0usize; S_TRANSPORTERS];
        let mut transport_int = [0usize; S_TRANSPORTERS];
        let mut uptake_request = [0.0f32; S_TRANSPORTERS];
        let mut secretion_request = [0.0f32; S_TRANSPORTERS];
        let mut uptake_requested_by_ext = [0.0f32; S_EXT];
        let mut secretion_requested_by_int = [0.0f32; M_INT];

        for i in 0..S_TRANSPORTERS {
            let tp = &self.ruleset.transport[i];
            let ext_idx = tp.ext_species as usize;
            let int_idx = tp.int_species as usize;
            if ext_idx >= S_EXT || int_idx >= CHEMICAL_INT_SPECIES {
                continue;
            }

            let ext_available = finite_nonnegative(ext_conc[ext_idx]);
            let internal_available = finite_nonnegative(self.internal[int_idx]);
            let gate_factor = transport_gate_factor(tp, &activation);
            let uptake_rate = finite_nonnegative(tp.uptake_rate) * gate_factor;
            let secretion_rate = finite_nonnegative(tp.secrete_rate) * gate_factor;
            let secretion = secretion_rate * internal_available / (1.0 + internal_available);
            let uptake = uptake_rate * ext_available / (1.0 + ext_available);
            let uptake_request_amount = uptake * dt;
            let secretion_request_amount = secretion * dt;
            let requested_amount = uptake_request_amount + secretion_request_amount;
            if sim.stoich_enforcement.is_strict()
                && requested_amount > sim.active_reaction_threshold
                && !transfer_delta(ext_idx, int_idx, 1.0).is_near_zero(1e-4)
            {
                if record_full_stoich && let Some(ledger) = stoich.as_deref_mut() {
                    ledger.record_strict_rejection(
                        StoichStage::CellTransport,
                        requested_amount,
                        self.lineage_id,
                        marl_config::stoich::STRICT_TEMPLATE_NONE,
                        keep_stoich_events,
                    );
                }
                continue;
            }

            valid_transport[i] = true;
            transport_ext[i] = ext_idx;
            transport_int[i] = int_idx;
            secretion_request[i] = secretion_request_amount;
            uptake_request[i] = uptake_request_amount;
            secretion_requested_by_int[int_idx] += secretion_request[i];
            uptake_requested_by_ext[ext_idx] += uptake_request[i];
        }

        let mut secretion_scale_by_int = [1.0f32; M_INT];
        for int_idx in 0..CHEMICAL_INT_SPECIES {
            let requested = secretion_requested_by_int[int_idx];
            let available = self.internal[int_idx].max(0.0);
            if requested > available && requested > 0.0 {
                secretion_scale_by_int[int_idx] = available / requested;
            }
        }

        let mut accepted_secretion_by_int = [0.0f32; M_INT];
        let mut accepted_secretion = [0.0f32; S_TRANSPORTERS];
        for i in 0..S_TRANSPORTERS {
            if !valid_transport[i] {
                continue;
            }
            let int_idx = transport_int[i];
            let amount = secretion_request[i] * secretion_scale_by_int[int_idx];
            accepted_secretion_by_int[int_idx] += amount;
            accepted_secretion[i] = amount;
        }

        let mut uptake_scale_by_ext = [1.0f32; S_EXT];
        for ext_idx in 0..S_EXT {
            let requested = uptake_requested_by_ext[ext_idx];
            let available = ext_conc[ext_idx].max(0.0);
            if requested > available && requested > 0.0 {
                uptake_scale_by_ext[ext_idx] = available / requested;
            }
        }

        let mut provisional_uptake = [0.0f32; S_TRANSPORTERS];
        let mut uptake_requested_by_int = [0.0f32; M_INT];
        for i in 0..S_TRANSPORTERS {
            if !valid_transport[i] {
                continue;
            }
            let ext_idx = transport_ext[i];
            let int_idx = transport_int[i];
            let amount = uptake_request[i] * uptake_scale_by_ext[ext_idx];
            provisional_uptake[i] = amount;
            uptake_requested_by_int[int_idx] += amount;
        }

        let mut uptake_scale_by_int = [1.0f32; M_INT];
        for int_idx in 0..CHEMICAL_INT_SPECIES {
            let requested = uptake_requested_by_int[int_idx];
            let post_secretion = self.internal[int_idx] - accepted_secretion_by_int[int_idx];
            let headroom = (sim.c_max - post_secretion).max(0.0);
            if requested > headroom && requested > 0.0 {
                uptake_scale_by_int[int_idx] = headroom / requested;
            }
        }

        let mut accepted_uptake_by_int = [0.0f32; M_INT];
        let mut accepted_uptake = [0.0f32; S_TRANSPORTERS];
        for i in 0..S_TRANSPORTERS {
            if !valid_transport[i] {
                continue;
            }
            let int_idx = transport_int[i];
            let amount = provisional_uptake[i] * uptake_scale_by_int[int_idx];
            accepted_uptake_by_int[int_idx] += amount;
            accepted_uptake[i] = amount;
        }

        let requested_transport_energy_cost = (0..S_TRANSPORTERS)
            .filter(|&i| valid_transport[i])
            .map(|i| {
                let ext_idx = transport_ext[i];
                (accepted_secretion[i] + accepted_uptake[i])
                    * sim.transport_energy_cost_per_unit(ext_idx)
            })
            .sum::<f32>();
        let net_energy_transport = accepted_uptake_by_int[0] - accepted_secretion_by_int[0];
        let energy_available = self.internal[0].max(0.0);
        let transport_scale = if requested_transport_energy_cost <= 0.0
            || !requested_transport_energy_cost.is_finite()
        {
            1.0
        } else {
            let uncovered_cost = requested_transport_energy_cost - net_energy_transport;
            if uncovered_cost <= 0.0 {
                1.0
            } else {
                (energy_available / uncovered_cost).clamp(0.0, 1.0)
            }
        };

        accepted_secretion_by_int = [0.0; M_INT];
        accepted_uptake_by_int = [0.0; M_INT];
        for i in 0..S_TRANSPORTERS {
            if !valid_transport[i] {
                continue;
            }
            let int_idx = transport_int[i];
            let ext_idx = transport_ext[i];
            let amount = accepted_secretion[i] * transport_scale;
            accepted_secretion_by_int[int_idx] += amount;
            field_deltas[ext_idx] += amount;
            if amount > 0.0
                && record_full_stoich
                && let Some(ledger) = stoich.as_deref_mut()
            {
                ledger.record(
                    StoichRecord::new(
                        StoichStage::CellTransport,
                        StoichEventKind::Transport,
                        amount,
                    )
                    .model_delta(transfer_delta(ext_idx, int_idx, -amount))
                    .actor(self.lineage_id)
                    .species(ext_idx),
                    keep_stoich_events,
                );
            }
        }

        for i in 0..S_TRANSPORTERS {
            if !valid_transport[i] {
                continue;
            }
            let int_idx = transport_int[i];
            let ext_idx = transport_ext[i];
            let amount = accepted_uptake[i] * transport_scale;
            accepted_uptake_by_int[int_idx] += amount;
            field_deltas[ext_idx] -= amount;
            if amount > 0.0
                && record_full_stoich
                && let Some(ledger) = stoich.as_deref_mut()
            {
                ledger.record(
                    StoichRecord::new(
                        StoichStage::CellTransport,
                        StoichEventKind::Transport,
                        amount,
                    )
                    .model_delta(transfer_delta(ext_idx, int_idx, amount))
                    .actor(self.lineage_id)
                    .species(ext_idx),
                    keep_stoich_events,
                );
            }
        }

        for int_idx in 0..CHEMICAL_INT_SPECIES {
            self.internal[int_idx] = (self.internal[int_idx] - accepted_secretion_by_int[int_idx]
                + accepted_uptake_by_int[int_idx])
                .clamp(0.0, sim.c_max);
        }

        let transport_energy_cost = requested_transport_energy_cost * transport_scale;
        self.internal[0] = (self.internal[0] - transport_energy_cost).max(0.0);
        if record_full_stoich
            && transport_energy_cost > 0.0
            && let Some(ledger) = stoich.as_deref_mut()
        {
            ledger.record(
                StoichRecord::new(
                    StoichStage::CellTransport,
                    StoichEventKind::TransportEnergyCost,
                    transport_energy_cost,
                )
                .model_delta(internal_delta(0, -transport_energy_cost))
                .balancing_reservoir(StoichReservoir::HeatSink)
                .actor(self.lineage_id)
                .species(0),
                keep_stoich_events,
            );
        }

        // Light is stored as a pseudo-internal concentration so reactions can use it as catalyst.
        // Internal species 15 (last slot) = light availability this tick.
        // Photosynthesis reactions reference catalyst=15 to be light-dependent.
        self.internal[LIGHT_SPECIES] = light;

        // === PHASE 3: INTRACELLULAR REACTIONS ===
        for rxn in &self.ruleset.reactions {
            if rxn.v_max.abs() < sim.active_reaction_threshold {
                continue;
            } // inactive slot
            if rxn.substrate == rxn.product {
                continue;
            } // no-op (e.g. energy→energy)

            let sub_idx = rxn.substrate as usize;
            let prod_idx = rxn.product as usize;
            let cat_idx = rxn.catalyst as usize;
            if sub_idx >= CHEMICAL_INT_SPECIES
                || prod_idx >= CHEMICAL_INT_SPECIES
                || cat_idx >= M_INT
            {
                continue;
            }

            let s = self.internal[sub_idx];
            let c = self.internal[cat_idx];

            // Michaelis-Menten with epsilon background rate
            let substrate_term = s / (rxn.k_m + s + f32::EPSILON);
            let catalyst_term = sim.epsilon + c / (rxn.k_cat + c + f32::EPSILON);
            let mut rate = rxn.v_max * substrate_term * catalyst_term;
            if !sim.stoich_enforcement.is_strict() {
                rate *= sim.reaction_descriptor_factor(sub_idx, prod_idx, rxn.cofactor);
            }

            // Optional cofactor
            if rxn.cofactor != 0xFF {
                let cof_idx = rxn.cofactor as usize;
                if cof_idx < CHEMICAL_INT_SPECIES {
                    let cof = self.internal[cof_idx];
                    rate *= cof / (1.0 + cof);
                } else {
                    continue;
                }
            }

            // Clamp flux to available substrate and product capacity so reactions
            // cannot create overflow or destroy substrate at a saturated product.
            let cof_idx = (rxn.cofactor != 0xFF).then_some(rxn.cofactor as usize);
            let mut max_flux = rate * dt;
            if let Some(cof_idx) = cof_idx {
                if cof_idx == sub_idx {
                    max_flux = max_flux.min(self.internal[sub_idx] / 1.5);
                } else {
                    max_flux = max_flux.min(self.internal[sub_idx]);
                    max_flux = max_flux.min(self.internal[cof_idx] / 0.5);
                }
            } else {
                max_flux = max_flux.min(self.internal[sub_idx]);
            }

            let product_headroom = (sim.c_max - self.internal[prod_idx]).max(0.0);
            let product_gain_per_flux = if cof_idx == Some(prod_idx) { 0.5 } else { 1.0 };
            let flux = max_flux.min(product_headroom / product_gain_per_flux);
            if flux <= 0.0 {
                continue;
            }
            let template = balanced_template_for_reaction(
                rxn.substrate,
                rxn.product,
                rxn.catalyst,
                rxn.cofactor,
            );
            if sim.stoich_enforcement.is_strict() && template.is_none() {
                if let Some(ledger) = stoich.as_deref_mut() {
                    ledger.record_strict_rejection(
                        StoichStage::Reactions,
                        flux,
                        self.lineage_id,
                        marl_config::stoich::STRICT_TEMPLATE_NONE,
                        keep_stoich_events,
                    );
                }
                continue;
            }

            self.internal[sub_idx] -= flux;
            if let Some(cof_idx) = cof_idx {
                self.internal[cof_idx] -= 0.5 * flux;
            }
            self.internal[prod_idx] += flux;
            if let Some(ledger) = stoich.as_deref_mut() {
                if sim.stoich_enforcement.is_strict() {
                    ledger.record_balanced_reaction(
                        template.expect("strict mode template checked above"),
                        flux,
                        self.lineage_id,
                        keep_stoich_events,
                    );
                } else {
                    ledger.record_legacy_reaction(
                        LegacyReactionAudit::new(
                            rxn.substrate,
                            rxn.product,
                            rxn.catalyst,
                            rxn.cofactor,
                            flux,
                        )
                        .actor(self.lineage_id),
                        keep_stoich_events,
                    );
                }
            }
        }

        // Maintenance energy drain: each tick, the cell loses lambda_maintenance
        // fraction of its current energy. During division prep, this cost is
        // multiplied — the cell is duplicating its genome, ribosomes, membranes.
        // Cells that evolve shorter prep times pay even more (rush penalty).
        let prep_multiplier = if self.prep_remaining > 0 {
            // Rush penalty: faster prep = more expensive per tick
            let evolved_prep = self.ruleset.fate.division_prep_ticks.max(1.0);
            let rush = (sim.base_division_prep - evolved_prep).max(0.0);
            sim.prep_maintenance_multiplier + rush * sim.rush_penalty_rate
        } else {
            1.0
        };
        let maintenance_fraction = (sim.lambda_maintenance * prep_multiplier * dt).clamp(0.0, 1.0);
        let energy_before_maintenance = self.internal[0];
        self.internal[0] *= 1.0 - maintenance_fraction;
        let maintenance_loss = energy_before_maintenance - self.internal[0];
        if record_full_stoich
            && maintenance_loss > 0.0
            && let Some(ledger) = stoich.as_deref_mut()
        {
            ledger.record(
                StoichRecord::new(
                    StoichStage::Maintenance,
                    StoichEventKind::Maintenance,
                    maintenance_loss,
                )
                .model_delta(internal_delta(0, -maintenance_loss))
                .balancing_reservoir(StoichReservoir::HeatSink)
                .actor(self.lineage_id)
                .species(0),
                keep_stoich_events,
            );
        }

        // Protein expression cost: each active enzyme requires transcription,
        // translation, and folding resources. No-op reactions (substrate == product)
        // are skipped — they don't encode a real enzyme.
        let active_rxn_count = self
            .ruleset
            .reactions
            .iter()
            .filter(|r| r.v_max.abs() > sim.active_reaction_threshold && r.substrate != r.product)
            .count();
        let expression_loss =
            (active_rxn_count as f32 * sim.reaction_maintenance * dt).min(self.internal[0]);
        self.internal[0] = (self.internal[0] - expression_loss).max(0.0);
        if record_full_stoich
            && expression_loss > 0.0
            && let Some(ledger) = stoich.as_deref_mut()
        {
            ledger.record(
                StoichRecord::new(
                    StoichStage::Maintenance,
                    StoichEventKind::ExpressionMaintenance,
                    expression_loss,
                )
                .model_delta(internal_delta(0, -expression_loss))
                .balancing_reservoir(StoichReservoir::HeatSink)
                .actor(self.lineage_id)
                .species(0),
                keep_stoich_events,
            );
        }

        // === PHASE 4: EFFECTOR PASS ===
        if !self.quiescent {
            for eff in &self.ruleset.effectors {
                let int_idx = eff.int_species as usize;
                let ext_idx = eff.ext_species as usize;
                if int_idx >= CHEMICAL_INT_SPECIES || ext_idx >= S_EXT {
                    continue;
                }
                if self.internal[int_idx] > eff.threshold {
                    let requested = eff.rate.max(0.0) * self.internal[int_idx]
                        / (1.0 + self.internal[int_idx])
                        * dt;
                    let amount = requested.min(self.internal[int_idx]);
                    if sim.stoich_enforcement.is_strict()
                        && amount > sim.active_reaction_threshold
                        && !transfer_delta(ext_idx, int_idx, -1.0).is_near_zero(1e-4)
                    {
                        if record_full_stoich && let Some(ledger) = stoich.as_deref_mut() {
                            ledger.record_strict_rejection(
                                StoichStage::Effectors,
                                amount,
                                self.lineage_id,
                                marl_config::stoich::STRICT_TEMPLATE_NONE,
                                keep_stoich_events,
                            );
                        }
                        continue;
                    }
                    self.internal[int_idx] -= amount;
                    field_deltas[ext_idx] += amount;
                    if record_full_stoich
                        && amount > 0.0
                        && let Some(ledger) = stoich.as_deref_mut()
                    {
                        ledger.record(
                            StoichRecord::new(
                                StoichStage::Effectors,
                                StoichEventKind::Effector,
                                amount,
                            )
                            .model_delta(transfer_delta(ext_idx, int_idx, -amount))
                            .actor(self.lineage_id)
                            .species(ext_idx),
                            keep_stoich_events,
                        );
                    }
                }
            }
        }

        // === PHASE 5: FATE DECISION ===
        let energy = self.internal[0];
        self.age += 1;

        // Death check first — hard floor is non-evolvable physics.
        let effective_death = self.ruleset.fate.death_energy.max(sim.hard_death_floor);

        let event = if energy < effective_death {
            // Dead — including cells that ran out of energy during prep.
            // Failed division attempts are a real biological phenomenon.
            self.prep_remaining = 0;
            CellEvent::Death
        } else if self.prep_remaining > 0 {
            // Currently in division prep — counting down.
            // The extra maintenance cost was already applied above.
            self.prep_remaining -= 1;
            if self.prep_remaining == 0 {
                // Prep complete — ready to divide!
                CellEvent::Division
            } else {
                CellEvent::None
            }
        } else if energy > self.ruleset.fate.division_energy {
            // Energy threshold reached — enter division prep phase.
            // The cell doesn't divide yet; it starts the costly prep countdown.
            let prep = self.ruleset.fate.division_prep_ticks.max(1.0) as u16;
            self.prep_remaining = prep;
            self.quiescent = false;
            CellEvent::None // division happens when prep_remaining hits 0
        } else if energy < self.ruleset.fate.quiescence_energy {
            self.quiescent = true;
            CellEvent::Quiescence
        } else {
            self.quiescent = false;
            CellEvent::None
        };

        (field_deltas, event)
    }
}

// ============================================================================
// MUTATION — the engine of evolution
// ============================================================================
//
// When a cell divides, its daughter's ruleset is mutated. Two kinds of
// mutation are implemented:
//
// 1. **Parametric mutation** (common): small Gaussian perturbations to
//    continuous parameters (rates, thresholds, kinetic constants). This
//    is like fine-tuning an enzyme's binding affinity.
//
// 2. **Structural mutation** (rare, 10x less likely): rewiring which
//    species a reaction acts on. This is like evolving a new enzyme
//    that catalyzes a completely different reaction.
//
// The mutation rate itself is evolvable (meta-evolution). Populations
// under strong selection pressure may evolve higher mutation rates to
// explore more of the fitness landscape. This is observed in real
// microbial populations (mutator phenotypes).

use rand::Rng;
use rand_distr::{Distribution, Normal};

impl Ruleset {
    fn sanitize_indices(&mut self) {
        for t in &mut self.transport {
            if t.ext_species as usize >= S_EXT || t.int_species as usize >= CHEMICAL_INT_SPECIES {
                t.uptake_rate = 0.0;
                t.secrete_rate = 0.0;
                t.ext_species = 0;
                t.int_species = 0;
            }
            if t.gate_receptor as usize >= S_RECEPTORS {
                t.gate_receptor = 0;
            }
            if t.gate_weight.is_finite() {
                t.gate_weight = t
                    .gate_weight
                    .clamp(-TRANSPORT_GATE_WEIGHT_LIMIT, TRANSPORT_GATE_WEIGHT_LIMIT);
            } else {
                t.gate_weight = 0.0;
            }
        }
        for r in &mut self.reactions {
            if r.substrate as usize >= CHEMICAL_INT_SPECIES
                || r.product as usize >= CHEMICAL_INT_SPECIES
                || r.catalyst as usize >= M_INT
            {
                r.v_max = 0.0;
                r.substrate = 0;
                r.product = 0;
                r.catalyst = 0;
            }
            if r.cofactor != 0xFF && r.cofactor as usize >= CHEMICAL_INT_SPECIES {
                r.cofactor = 0xFF;
            }
        }
        for e in &mut self.effectors {
            if e.int_species as usize >= CHEMICAL_INT_SPECIES || e.ext_species as usize >= S_EXT {
                e.rate = 0.0;
                e.int_species = 0;
                e.ext_species = 0;
            }
        }
    }

    /// Apply random mutations to all evolvable parameters.
    ///
    /// Called on the daughter cell's ruleset after division. Each parameter
    /// independently has a `mutation_rate` probability of being perturbed
    /// by a small Gaussian (mean=0, std=0.1). Parameters are clamped to
    /// non-negative values since rates and concentrations can't be negative.
    ///
    /// Structural mutations (rewiring reaction substrate/product/catalyst)
    /// happen 10x less frequently — they're more disruptive and usually lethal,
    /// but occasionally create entirely new metabolic capabilities.
    pub fn mutate(&mut self, rng: &mut impl Rng, sim: &SimulationConfig) {
        let rate = self.mutation_rate;
        let normal = Normal::new(0.0f32, sim.mutation_stddev).unwrap();

        // Helper: with probability `rate`, add a small Gaussian perturbation.
        fn maybe_mutate(val: &mut f32, rate: f32, normal: &Normal<f32>, rng: &mut impl Rng) {
            if rng.random::<f32>() < rate {
                *val += normal.sample(rng);
                *val = val.max(0.0);
            }
        }

        fn maybe_mutate_signed(
            val: &mut f32,
            rate: f32,
            normal: &Normal<f32>,
            rng: &mut impl Rng,
            limit: f32,
        ) {
            if rng.random::<f32>() < rate {
                *val += normal.sample(rng);
            }
            if val.is_finite() {
                *val = val.clamp(-limit, limit);
            } else {
                *val = 0.0;
            }
        }

        // Mutate receptor sensitivities
        for r in &mut self.receptors {
            maybe_mutate(&mut r.k_half, rate, &normal, rng);
            maybe_mutate(&mut r.n_hill, rate, &normal, rng);
            r.n_hill = r
                .n_hill
                .clamp(sim.hill_exponent_clamp_low, sim.hill_exponent_clamp_high);
            maybe_mutate(&mut r.gain, rate, &normal, rng);
        }

        // Mutate transport rates (+ rare structural: change which species)
        for t in &mut self.transport {
            maybe_mutate(&mut t.uptake_rate, rate, &normal, rng);
            maybe_mutate(&mut t.secrete_rate, rate, &normal, rng);
            maybe_mutate_signed(
                &mut t.gate_weight,
                rate,
                &normal,
                rng,
                TRANSPORT_GATE_WEIGHT_LIMIT,
            );
            if rng.random::<f32>() < rate * sim.structural_mutation_rate_mult {
                t.ext_species = rng.random_range(0..S_EXT as u8);
            }
            if rng.random::<f32>() < rate * sim.structural_mutation_rate_mult {
                t.int_species = rng.random_range(0..CHEMICAL_INT_SPECIES as u8);
            }
            if rng.random::<f32>() < rate * sim.structural_mutation_rate_mult {
                t.gate_receptor = rng.random_range(0..S_RECEPTORS as u8);
            }
        }

        // Mutate reaction kinetics (+ rare structural: gene duplication model)
        //
        // Gene duplication (Ohno, 1970) is the dominant mode of metabolic
        // innovation in real microbes: copy an existing working enzyme, then
        // let the copy diverge. We model this by copying substrate/product/
        // catalyst from a randomly-chosen ACTIVE reaction, rather than
        // assembling random species indices from scratch.
        //
        // First, collect indices of active reactions for duplication source.
        let active_indices: Vec<usize> = self
            .reactions
            .iter()
            .enumerate()
            .filter(|(_, r)| r.v_max.abs() > sim.active_reaction_threshold)
            .map(|(i, _)| i)
            .collect();

        for i in 0..R_MAX {
            maybe_mutate(&mut self.reactions[i].k_m, rate, &normal, rng);
            maybe_mutate(&mut self.reactions[i].v_max, rate, &normal, rng);
            maybe_mutate(&mut self.reactions[i].k_cat, rate, &normal, rng);

            if sim.stoich_enforcement.is_strict() {
                if rng.random::<f32>() < rate * sim.structural_mutation_rate_mult {
                    let template = BALANCED_REACTION_TEMPLATES
                        [rng.random_range(0..BALANCED_REACTION_TEMPLATES.len())];
                    self.reactions[i].substrate = template.substrate;
                    self.reactions[i].product = template.product;
                    self.reactions[i].catalyst = template.catalyst;
                    self.reactions[i].cofactor = template.cofactor;
                }
                continue;
            }

            // Structural mutation: copy topology from an existing active reaction
            // (gene duplication + divergence). Falls back to random if no active
            // reactions exist (shouldn't happen in practice).
            if rng.random::<f32>() < rate * sim.structural_mutation_rate_mult {
                if let Some(&donor) =
                    active_indices.get(rng.random_range(0..active_indices.len().max(1)))
                {
                    if donor != i {
                        // don't copy from self
                        self.reactions[i].substrate = self.reactions[donor].substrate;
                    }
                } else {
                    self.reactions[i].substrate = rng.random_range(0..CHEMICAL_INT_SPECIES as u8);
                }
            }
            if rng.random::<f32>() < rate * sim.structural_mutation_rate_mult {
                if let Some(&donor) =
                    active_indices.get(rng.random_range(0..active_indices.len().max(1)))
                {
                    if donor != i {
                        self.reactions[i].product = self.reactions[donor].product;
                    }
                } else {
                    self.reactions[i].product = rng.random_range(0..CHEMICAL_INT_SPECIES as u8);
                }
            }
            if rng.random::<f32>() < rate * sim.structural_mutation_rate_mult {
                if let Some(&donor) =
                    active_indices.get(rng.random_range(0..active_indices.len().max(1)))
                {
                    if donor != i {
                        self.reactions[i].catalyst = self.reactions[donor].catalyst;
                    }
                } else {
                    self.reactions[i].catalyst = rng.random_range(0..M_INT as u8);
                }
            }
        }

        // Mutate effector thresholds and rates
        for e in &mut self.effectors {
            maybe_mutate(&mut e.threshold, rate, &normal, rng);
            maybe_mutate(&mut e.rate, rate, &normal, rng);
        }

        // Mutate fate decision thresholds and cell cycle timing
        maybe_mutate(&mut self.fate.division_energy, rate, &normal, rng);
        maybe_mutate(&mut self.fate.death_energy, rate, &normal, rng);
        maybe_mutate(&mut self.fate.quiescence_energy, rate, &normal, rng);
        maybe_mutate(&mut self.fate.division_prep_ticks, rate, &normal, rng);
        self.fate.division_prep_ticks = self.fate.division_prep_ticks.max(1.0); // minimum 1 tick
        maybe_mutate(&mut self.hgt_propensity, rate, &normal, rng);

        // Meta-evolution: the mutation rate itself can mutate.
        // This happens at a fixed rate (not gated by the current mutation_rate)
        // to prevent runaway suppression of evolvability.
        if rng.random::<f32>() < sim.meta_mutation_rate {
            self.mutation_rate += normal.sample(rng) * sim.meta_mutation_rate;
            self.mutation_rate = self
                .mutation_rate
                .clamp(sim.meta_mutation_clamp_low, sim.meta_mutation_clamp_high);
        }

        self.sanitize_indices();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn inactive_receptor() -> ReceptorParams {
        ReceptorParams {
            k_half: 1.0,
            n_hill: 1.0,
            gain: 0.0,
        }
    }

    fn inactive_transport() -> TransportParams {
        TransportParams {
            uptake_rate: 0.0,
            secrete_rate: 0.0,
            ext_species: 0,
            int_species: 0,
            gate_receptor: 0,
            gate_weight: 0.0,
        }
    }

    fn inactive_reaction() -> Reaction {
        Reaction {
            substrate: 0,
            product: 0,
            catalyst: 0,
            cofactor: 0xFF,
            k_m: 1.0,
            v_max: 0.0,
            k_cat: 1.0,
        }
    }

    fn inactive_effector() -> EffectorParams {
        EffectorParams {
            threshold: f32::MAX,
            rate: 0.0,
            int_species: 0,
            ext_species: 0,
        }
    }

    fn test_ruleset() -> Ruleset {
        Ruleset {
            receptors: std::array::from_fn(|_| inactive_receptor()),
            transport: std::array::from_fn(|_| inactive_transport()),
            reactions: std::array::from_fn(|_| inactive_reaction()),
            effectors: std::array::from_fn(|_| inactive_effector()),
            fate: FateParams {
                division_energy: 100.0,
                death_energy: 0.0,
                quiescence_energy: 0.0,
                division_prep_ticks: 20.0,
            },
            hgt_propensity: 0.0,
            mutation_rate: 0.0,
        }
    }

    fn test_cell(ruleset: Ruleset) -> CellState {
        let mut internal = [0.0f32; M_INT];
        internal[0] = 10.0;
        CellState {
            pos: [1, 1, 1],
            lineage_id: 1,
            age: 0,
            internal,
            ruleset,
            quiescent: false,
            starter_type: 0,
            prep_remaining: 0,
        }
    }

    #[test]
    fn transport_uptake_is_limited_by_available_external_concentration() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let mut cell = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[2] = 0.1;

        let (deltas, _) = cell.tick(&ext, 0.0, &SimulationConfig::default());

        assert!((cell.internal[2] - 0.1).abs() < 1e-6);
        assert!((deltas[2] + 0.1).abs() < 1e-6);
    }

    #[test]
    fn accepted_transport_flux_pays_descriptor_energy_cost() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: EXT_REDUCTANT as u8,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let sim = SimulationConfig {
            lambda_maintenance: 0.0,
            ..SimulationConfig::default()
        };
        let mut cell = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[EXT_REDUCTANT] = 1.0;

        let (deltas, _) = cell.tick(&ext, 0.0, &sim);

        let expected_cost = sim.transport_energy_cost_per_unit(EXT_REDUCTANT);
        assert!((cell.internal[2] - 1.0).abs() < 1e-6);
        assert!((deltas[EXT_REDUCTANT] + 1.0).abs() < 1e-6);
        assert!((cell.internal[0] - (10.0 - expected_cost)).abs() < 1e-6);
    }

    #[test]
    fn zero_energy_cell_cannot_move_costed_nonenergy_flux() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: EXT_REDUCTANT as u8,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let sim = SimulationConfig {
            lambda_maintenance: 0.0,
            ..SimulationConfig::default()
        };
        let mut cell = test_cell(ruleset);
        cell.internal[0] = 0.0;
        let mut ext = [0.0f32; S_EXT];
        ext[EXT_REDUCTANT] = 1.0;

        let (deltas, _) = cell.tick(&ext, 0.0, &sim);

        assert_eq!(cell.internal[2], 0.0);
        assert_eq!(deltas[EXT_REDUCTANT], 0.0);
        assert_eq!(cell.internal[0], 0.0);
    }

    #[test]
    fn scarce_energy_scales_transport_flux_and_audit_amounts() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: EXT_REDUCTANT as u8,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let sim = SimulationConfig {
            lambda_maintenance: 0.0,
            ..SimulationConfig::default()
        };
        let full_cost = sim.transport_energy_cost_per_unit(EXT_REDUCTANT);
        let mut cell = test_cell(ruleset);
        cell.internal[0] = full_cost * 0.5;
        let mut ext = [0.0f32; S_EXT];
        ext[EXT_REDUCTANT] = 1.0;
        let mut ledger = StoichTickLedger::default();

        let (deltas, _) = cell.tick_with_stoich(&ext, 0.0, &sim, Some(&mut ledger), true, true);

        assert!((cell.internal[2] - 0.5).abs() < 1e-6);
        assert!((deltas[EXT_REDUCTANT] + 0.5).abs() < 1e-6);
        assert!(cell.internal[0].abs() < 1e-6);

        let transport_event = ledger
            .events
            .iter()
            .find(|event| event.kind == StoichEventKind::Transport)
            .expect("scaled transport event should be recorded");
        let cost_event = ledger
            .events
            .iter()
            .find(|event| event.kind == StoichEventKind::TransportEnergyCost)
            .expect("scaled transport cost event should be recorded");
        assert!((transport_event.amount - 0.5).abs() < 1e-6);
        assert!((cost_event.amount - full_cost * 0.5).abs() < 1e-6);
    }

    #[test]
    fn transport_energy_cost_is_audited_as_heat_loss() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: EXT_STRUCTURAL as u8,
            int_species: 7,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let sim = SimulationConfig {
            lambda_maintenance: 0.0,
            stoich_enforcement: marl_config::stoich::StoichEnforcement::Audit,
            ..SimulationConfig::default()
        };
        let mut cell = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[EXT_STRUCTURAL] = 1.0;
        let mut ledger = StoichTickLedger::default();

        cell.tick_with_stoich(&ext, 0.0, &sim, Some(&mut ledger), true, true);

        let cost_event = ledger
            .events
            .iter()
            .find(|event| event.kind == StoichEventKind::TransportEnergyCost)
            .expect("transport cost event should be recorded");
        let expected_cost = sim.transport_energy_cost_per_unit(EXT_STRUCTURAL);
        assert!((cost_event.amount - expected_cost).abs() < 1e-6);
        assert!(cost_event.balanced);
        assert!((cost_event.model_delta.energy + expected_cost).abs() < 1e-6);
        assert!((cost_event.reservoir_delta.energy - expected_cost).abs() < 1e-6);
    }

    #[test]
    fn transport_uptake_stops_when_internal_pool_is_full() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let sim = SimulationConfig::default();
        let mut cell = test_cell(ruleset);
        cell.internal[2] = sim.c_max;
        let mut ext = [0.0f32; S_EXT];
        ext[2] = 1.0;

        let (deltas, _) = cell.tick(&ext, 0.0, &sim);

        assert_eq!(cell.internal[2], sim.c_max);
        assert_eq!(deltas[2], 0.0);
    }

    #[test]
    fn duplicated_transporters_share_external_uptake_budget() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        ruleset.transport[1] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 3,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let mut cell = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[2] = 1.0;

        let (deltas, _) = cell.tick(&ext, 0.0, &SimulationConfig::default());

        assert!((cell.internal[2] + cell.internal[3] - 1.0).abs() < 1e-6);
        assert!((deltas[2] + 1.0).abs() < 1e-6);
    }

    #[test]
    fn duplicated_secretors_share_internal_pool_budget() {
        let mut ruleset = test_ruleset();
        ruleset.transport[0] = TransportParams {
            uptake_rate: 0.0,
            secrete_rate: 100.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        ruleset.transport[1] = TransportParams {
            uptake_rate: 0.0,
            secrete_rate: 100.0,
            ext_species: 3,
            int_species: 2,
            gate_receptor: 0,
            gate_weight: 0.0,
        };
        let mut cell = test_cell(ruleset);
        cell.internal[2] = 1.0;

        let (deltas, _) = cell.tick(&[0.0; S_EXT], 0.0, &SimulationConfig::default());

        assert_eq!(cell.internal[2], 0.0);
        assert!((deltas[2] + deltas[3] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn receptor_gate_can_increase_transport_rate() {
        let mut ruleset = test_ruleset();
        ruleset.receptors[2] = ReceptorParams {
            k_half: 1.0,
            n_hill: 1.0,
            gain: 1.0,
        };
        ruleset.transport[0] = TransportParams {
            uptake_rate: 1.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 2,
            gate_weight: 1.0,
        };
        let mut cell = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[2] = 1.0;

        let (deltas, _) = cell.tick(&ext, 0.0, &SimulationConfig::default());

        let expected = 1.5 * ext[2] / (1.0 + ext[2]);
        assert!((cell.internal[2] - expected).abs() < 1e-6);
        assert!((deltas[2] + expected).abs() < 1e-6);
    }

    #[test]
    fn receptor_gate_can_decrease_transport_rate() {
        let mut ruleset = test_ruleset();
        ruleset.receptors[2] = ReceptorParams {
            k_half: 1.0,
            n_hill: 1.0,
            gain: 1.0,
        };
        ruleset.transport[0] = TransportParams {
            uptake_rate: 1.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 2,
            gate_weight: -1.0,
        };
        let mut cell = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[2] = 1.0;

        let (deltas, _) = cell.tick(&ext, 0.0, &SimulationConfig::default());

        let expected = 0.5 * ext[2] / (1.0 + ext[2]);
        assert!((cell.internal[2] - expected).abs() < 1e-6);
        assert!((deltas[2] + expected).abs() < 1e-6);
    }

    #[test]
    fn receptor_gate_uses_selected_receptor_not_transport_species() {
        let mut ruleset = test_ruleset();
        ruleset.receptors[1] = ReceptorParams {
            k_half: 0.5,
            n_hill: 2.0,
            gain: 1.0,
        };
        ruleset.transport[0] = TransportParams {
            uptake_rate: 1.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 1,
            gate_weight: 1.0,
        };
        let mut without_signal = test_cell(ruleset.clone());
        let mut with_signal = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[2] = 1.0;

        without_signal.tick(&ext, 0.0, &SimulationConfig::default());
        ext[1] = 0.5;
        with_signal.tick(&ext, 0.0, &SimulationConfig::default());

        assert!(with_signal.internal[2] > without_signal.internal[2]);
    }

    #[test]
    fn degenerate_gate_is_bounded_and_can_shut_off_transport() {
        let mut ruleset = test_ruleset();
        ruleset.receptors[2] = ReceptorParams {
            k_half: 1e-12,
            n_hill: 100.0,
            gain: 1e9,
        };
        ruleset.transport[0] = TransportParams {
            uptake_rate: 100.0,
            secrete_rate: 0.0,
            ext_species: 2,
            int_species: 2,
            gate_receptor: 2,
            gate_weight: -100.0,
        };
        let mut cell = test_cell(ruleset);
        let mut ext = [0.0f32; S_EXT];
        ext[2] = 1.0;

        let (deltas, _) = cell.tick(&ext, 0.0, &SimulationConfig::default());

        assert_eq!(cell.internal[2], 0.0);
        assert_eq!(deltas[2], 0.0);
        assert!(cell.internal.iter().all(|value| value.is_finite()));
        assert!(deltas.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn effector_secretion_is_limited_by_internal_pool() {
        let mut ruleset = test_ruleset();
        ruleset.effectors[0] = EffectorParams {
            threshold: 0.0,
            rate: 100.0,
            int_species: 2,
            ext_species: 4,
        };
        let mut cell = test_cell(ruleset);
        cell.internal[2] = 0.05;

        let (deltas, _) = cell.tick(&[0.0; S_EXT], 0.0, &SimulationConfig::default());

        assert_eq!(cell.internal[2], 0.0);
        assert!((deltas[4] - 0.05).abs() < 1e-6);
    }

    #[test]
    fn saturated_reaction_product_does_not_destroy_substrate() {
        let mut ruleset = test_ruleset();
        ruleset.reactions[0] = Reaction {
            substrate: 2,
            product: 3,
            catalyst: 0,
            cofactor: 0xFF,
            k_m: 0.01,
            v_max: 100.0,
            k_cat: 0.01,
        };
        let sim = SimulationConfig::default();
        let mut cell = test_cell(ruleset);
        cell.internal[2] = 1.0;
        cell.internal[3] = sim.c_max;

        cell.tick(&[0.0; S_EXT], 0.0, &sim);

        assert!((cell.internal[2] - 1.0).abs() < 1e-6);
        assert_eq!(cell.internal[3], sim.c_max);
    }

    #[test]
    fn cofactor_consumption_tracks_actual_flux() {
        let mut ruleset = test_ruleset();
        ruleset.reactions[0] = Reaction {
            substrate: 2,
            product: 3,
            catalyst: 0,
            cofactor: 1,
            k_m: 0.01,
            v_max: 100.0,
            k_cat: 0.01,
        };
        let mut cell = test_cell(ruleset);
        cell.internal[1] = 10.0;
        cell.internal[2] = 0.1;

        cell.tick(&[0.0; S_EXT], 0.0, &SimulationConfig::default());

        assert!((cell.internal[2] - 0.0).abs() < 1e-6);
        assert!((cell.internal[3] - 0.1).abs() < 1e-6);
        assert!((cell.internal[1] - 9.95).abs() < 1e-6);
    }

    #[test]
    fn cofactor_same_as_substrate_cannot_overdraw_pool() {
        let mut ruleset = test_ruleset();
        ruleset.reactions[0] = Reaction {
            substrate: 2,
            product: 3,
            catalyst: 0,
            cofactor: 2,
            k_m: 0.01,
            v_max: 100.0,
            k_cat: 0.01,
        };
        let mut cell = test_cell(ruleset);
        cell.internal[2] = 0.15;

        cell.tick(&[0.0; S_EXT], 0.0, &SimulationConfig::default());

        assert!(cell.internal[2].abs() < 1e-6);
        assert!((cell.internal[3] - 0.1).abs() < 1e-6);
    }

    #[test]
    fn descriptor_coupling_scales_reaction_flux_but_strict_mode_skips_it() {
        let mut ruleset = test_ruleset();
        ruleset.reactions[0] = Reaction {
            substrate: EXT_CARBON as u8,
            product: EXT_ENERGY as u8,
            catalyst: LIGHT_SPECIES as u8,
            cofactor: 0xFF,
            k_m: 1.0,
            v_max: 0.05,
            k_cat: 1.0,
        };

        let legacy_sim = SimulationConfig {
            lambda_maintenance: 0.0,
            reaction_maintenance: 0.0,
            reaction_descriptor_coupling_strength: 0.0,
            ..SimulationConfig::default()
        };
        let coupled_sim = SimulationConfig {
            reaction_descriptor_coupling_strength: 1.0,
            ..legacy_sim.clone()
        };
        let strict_sim = SimulationConfig {
            stoich_enforcement: marl_config::stoich::StoichEnforcement::Strict,
            reaction_descriptor_coupling_strength: 1.0,
            ..legacy_sim.clone()
        };
        let mut legacy = test_cell(ruleset.clone());
        let mut coupled = test_cell(ruleset.clone());
        let mut strict = test_cell(ruleset);
        legacy.internal[EXT_CARBON] = 10.0;
        coupled.internal[EXT_CARBON] = 10.0;
        strict.internal[EXT_CARBON] = 10.0;

        legacy.tick(&[0.0; S_EXT], 1.0, &legacy_sim);
        coupled.tick(&[0.0; S_EXT], 1.0, &coupled_sim);
        strict.tick(&[0.0; S_EXT], 1.0, &strict_sim);

        let legacy_flux = legacy.internal[EXT_ENERGY] - 10.0;
        let coupled_flux = coupled.internal[EXT_ENERGY] - 10.0;
        let expected_factor = coupled_sim.reaction_descriptor_factor(EXT_CARBON, EXT_ENERGY, 0xFF);
        assert!((coupled_flux - legacy_flux * expected_factor).abs() < 1e-6);
        assert!((strict.internal[EXT_ENERGY] - legacy.internal[EXT_ENERGY]).abs() < 1e-6);
    }

    #[test]
    fn stoich_audit_mode_preserves_tick_behavior() {
        let mut ruleset = test_ruleset();
        ruleset.reactions[0] = Reaction {
            substrate: 2,
            product: 3,
            catalyst: 0,
            cofactor: 1,
            k_m: 0.01,
            v_max: 100.0,
            k_cat: 0.01,
        };
        let sim = SimulationConfig::default();
        let mut plain = test_cell(ruleset.clone());
        let mut audited = test_cell(ruleset);
        plain.internal[1] = 10.0;
        plain.internal[2] = 0.1;
        audited.internal[1] = 10.0;
        audited.internal[2] = 0.1;

        let (plain_deltas, plain_event) = plain.tick(&[0.0; S_EXT], 0.0, &sim);
        let mut ledger = StoichTickLedger::default();
        let (audit_deltas, audit_event) =
            audited.tick_with_stoich(&[0.0; S_EXT], 0.0, &sim, Some(&mut ledger), false, false);

        assert_eq!(plain.pos, audited.pos);
        assert_eq!(plain.lineage_id, audited.lineage_id);
        assert_eq!(plain.age, audited.age);
        assert_eq!(plain.quiescent, audited.quiescent);
        assert_eq!(plain.starter_type, audited.starter_type);
        assert_eq!(plain.prep_remaining, audited.prep_remaining);
        for (lhs, rhs) in plain.internal.iter().zip(audited.internal.iter()) {
            assert!((lhs - rhs).abs() < 1e-6);
        }
        for (lhs, rhs) in plain_deltas.iter().zip(audit_deltas.iter()) {
            assert!((lhs - rhs).abs() < 1e-6);
        }
        assert_eq!(
            std::mem::discriminant(&plain_event),
            std::mem::discriminant(&audit_event)
        );
        assert_eq!(ledger.reaction_count, 1);
        assert!((ledger.active_flux - 0.1).abs() < 1e-6);
    }

    #[test]
    fn mutation_keeps_light_as_catalyst_only_pseudo_species() {
        let mut ruleset = test_ruleset();
        ruleset.mutation_rate = 1.0;
        ruleset.transport[0].int_species = LIGHT_SPECIES as u8;
        ruleset.reactions[0] = Reaction {
            substrate: LIGHT_SPECIES as u8,
            product: LIGHT_SPECIES as u8,
            catalyst: M_INT as u8,
            cofactor: LIGHT_SPECIES as u8,
            k_m: 1.0,
            v_max: 1.0,
            k_cat: 1.0,
        };
        ruleset.effectors[0].int_species = LIGHT_SPECIES as u8;
        ruleset.effectors[0].ext_species = S_EXT as u8;
        let sim = SimulationConfig {
            structural_mutation_rate_mult: 1.0,
            ..SimulationConfig::default()
        };
        let mut rng = StdRng::seed_from_u64(7);

        ruleset.mutate(&mut rng, &sim);

        for t in &ruleset.transport {
            assert!((t.int_species as usize) < LIGHT_SPECIES);
        }
        for r in &ruleset.reactions {
            assert!((r.substrate as usize) < LIGHT_SPECIES);
            assert!((r.product as usize) < LIGHT_SPECIES);
            assert!((r.catalyst as usize) < M_INT);
            assert!(r.cofactor == 0xFF || (r.cofactor as usize) < LIGHT_SPECIES);
        }
        for e in &ruleset.effectors {
            assert!((e.int_species as usize) < LIGHT_SPECIES);
            assert!((e.ext_species as usize) < S_EXT);
        }
    }

    #[test]
    fn strict_mutation_uses_balanced_reaction_templates() {
        let mut ruleset = test_ruleset();
        ruleset.mutation_rate = 1.0;
        let sim = SimulationConfig {
            stoich_enforcement: marl_config::stoich::StoichEnforcement::Strict,
            structural_mutation_rate_mult: 1.0,
            ..SimulationConfig::default()
        };
        let mut rng = StdRng::seed_from_u64(17);

        ruleset.mutate(&mut rng, &sim);

        for reaction in &ruleset.reactions {
            assert!(
                balanced_template_for_reaction(
                    reaction.substrate,
                    reaction.product,
                    reaction.catalyst,
                    reaction.cofactor
                )
                .is_some()
            );
        }
    }

    #[test]
    fn strict_partial_structural_mutation_preserves_balanced_templates() {
        let mut ruleset = test_ruleset();
        ruleset.mutation_rate = 1.0;
        for (reaction, template) in ruleset
            .reactions
            .iter_mut()
            .zip(BALANCED_REACTION_TEMPLATES.iter().cycle())
        {
            reaction.substrate = template.substrate;
            reaction.product = template.product;
            reaction.catalyst = template.catalyst;
            reaction.cofactor = template.cofactor;
            reaction.v_max = 1.0;
        }
        let sim = SimulationConfig {
            stoich_enforcement: marl_config::stoich::StoichEnforcement::Strict,
            structural_mutation_rate_mult: 0.35,
            ..SimulationConfig::default()
        };
        let mut rng = StdRng::seed_from_u64(301);

        ruleset.mutate(&mut rng, &sim);

        for reaction in &ruleset.reactions {
            assert!(
                balanced_template_for_reaction(
                    reaction.substrate,
                    reaction.product,
                    reaction.catalyst,
                    reaction.cofactor
                )
                .is_some()
            );
        }
    }
}
