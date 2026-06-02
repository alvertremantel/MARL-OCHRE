use marl_cell::cell::CellState;
use marl_config::stoich::{
    StoichEventKind, StoichRecord, StoichReservoir, StoichStage, StoichTickLedger, external_delta,
};
use marl_config::{GRID_X, GRID_Y, GRID_Z, SimulationConfig};
use marl_field::field::Field;
use rand::Rng;
use std::collections::HashMap;

/// Initialize field with thin boundary layers only.
/// The bulk of the field starts empty — gradients build from boundary sources + diffusion.
pub fn init_field_boundaries(field: &mut Field, sim: &SimulationConfig) {
    init_field_boundaries_with_stoich(field, sim, None, false);
}

pub fn init_field_boundaries_with_stoich(
    field: &mut Field,
    sim: &SimulationConfig,
    mut stoich: Option<&mut StoichTickLedger>,
    keep_events: bool,
) {
    // Prime only the boundary faces (configurable layers deep) so initial cells
    // have a local substrate source but the bulk field is empty.
    let layers = sim.boundary_prime_layers.min(GRID_Z);
    for y in 0..GRID_Y {
        for x in 0..GRID_X {
            // Top layers: some oxidant and carbon (atmosphere analog)
            for z in 0..layers {
                field.set(x, y, z, 1, sim.boundary_prime_oxidant); // oxidant
                if let Some(ledger) = stoich.as_deref_mut() {
                    record_prime(ledger, 1, sim.boundary_prime_oxidant, keep_events);
                }
                field.set(x, y, z, 3, sim.boundary_prime_carbon); // carbon
                if let Some(ledger) = stoich.as_deref_mut() {
                    record_prime(ledger, 3, sim.boundary_prime_carbon, keep_events);
                }
            }
            // Bottom layers: some reductant (geological source analog)
            for z in (GRID_Z - layers)..GRID_Z {
                field.set(x, y, z, 2, sim.boundary_prime_reductant); // reductant
                if let Some(ledger) = stoich.as_deref_mut() {
                    record_prime(ledger, 2, sim.boundary_prime_reductant, keep_events);
                }
            }
        }
    }
}

fn record_prime(ledger: &mut StoichTickLedger, species: usize, amount: f32, keep_events: bool) {
    if amount <= 0.0 {
        return;
    }
    ledger.record(
        StoichRecord::new(
            StoichStage::BoundaryPriming,
            StoichEventKind::BoundaryPrime,
            amount,
        )
        .model_delta(external_delta(species, amount))
        .balancing_reservoir(StoichReservoir::BoundaryInput)
        .species(species),
        keep_events,
    );
}

#[allow(clippy::too_many_arguments)]
pub fn seed_cells(
    cells: &mut Vec<CellState>,
    cell_map: &mut HashMap<[u16; 3], usize>,
    rng: &mut impl Rng,
    count: usize,
    z_lo: u16,
    z_hi: u16,
    factory: fn([u16; 3], u64) -> CellState,
    sim: &SimulationConfig,
) {
    let margin = sim.seed_margin;
    let x_hi = (GRID_X as u16).saturating_sub(margin);
    let y_hi = (GRID_Y as u16).saturating_sub(margin);
    let z_hi = z_hi.min(GRID_Z as u16);
    let z_lo = z_lo.min(GRID_Z as u16);
    if margin >= x_hi || margin >= y_hi || z_lo >= z_hi {
        return;
    }

    let mut seeded = 0;
    for _ in 0..count.saturating_mul(8) {
        if seeded >= count {
            break;
        }
        let x = rng.random_range(margin..x_hi);
        let y = rng.random_range(margin..y_hi);
        let z = rng.random_range(z_lo..z_hi);
        let pos = [x, y, z];
        if cell_map.contains_key(&pos) {
            continue;
        }

        let cell = factory(pos, rng.random::<u64>());
        let idx = cells.len();
        cell_map.insert(pos, idx);
        cells.push(cell);
        seeded += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use marl_cell::cell::{
        CellState, EffectorParams, FateParams, Reaction, ReceptorParams, Ruleset, TransportParams,
    };
    use marl_config::M_INT;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn test_cell(pos: [u16; 3], lineage_id: u64) -> CellState {
        CellState {
            pos,
            lineage_id,
            age: 0,
            internal: [0.0; M_INT],
            ruleset: Ruleset {
                receptors: std::array::from_fn(|_| ReceptorParams {
                    k_half: 1.0,
                    n_hill: 1.0,
                    gain: 0.0,
                }),
                transport: std::array::from_fn(|_| TransportParams {
                    uptake_rate: 0.0,
                    secrete_rate: 0.0,
                    ext_species: 0,
                    int_species: 0,
                }),
                reactions: std::array::from_fn(|_| Reaction {
                    substrate: 0,
                    product: 0,
                    catalyst: 0,
                    cofactor: 0xFF,
                    k_m: 1.0,
                    v_max: 0.0,
                    k_cat: 1.0,
                }),
                effectors: std::array::from_fn(|_| EffectorParams {
                    threshold: 1.0,
                    rate: 0.0,
                    int_species: 0,
                    ext_species: 0,
                }),
                fate: FateParams {
                    division_energy: 1.0,
                    death_energy: 0.0,
                    quiescence_energy: 0.0,
                    division_prep_ticks: 1.0,
                },
                hgt_propensity: 0.0,
                mutation_rate: 0.0,
            },
            quiescent: false,
            starter_type: 0,
            prep_remaining: 0,
        }
    }

    #[test]
    fn boundary_priming_clamps_layers_to_grid_depth() {
        let sim = SimulationConfig {
            boundary_prime_layers: GRID_Z + 10,
            ..SimulationConfig::default()
        };
        let mut field = Field::new();

        init_field_boundaries(&mut field, &sim);

        assert_eq!(field.get(0, 0, 0, 1), sim.boundary_prime_oxidant);
        assert_eq!(field.get(0, 0, GRID_Z - 1, 2), sim.boundary_prime_reductant);
    }

    #[test]
    fn seeding_returns_without_panicking_for_empty_ranges() {
        let mut cells = Vec::new();
        let mut cell_map = HashMap::new();
        let mut rng = StdRng::seed_from_u64(11);
        let sim = SimulationConfig {
            seed_margin: GRID_X as u16,
            ..SimulationConfig::default()
        };

        seed_cells(
            &mut cells,
            &mut cell_map,
            &mut rng,
            10,
            GRID_Z as u16,
            GRID_Z as u16,
            test_cell,
            &sim,
        );

        assert!(cells.is_empty());
        assert!(cell_map.is_empty());
    }

    #[test]
    fn seeding_places_cells_inside_clamped_valid_ranges() {
        let mut cells = Vec::new();
        let mut cell_map = HashMap::new();
        let mut rng = StdRng::seed_from_u64(13);
        let sim = SimulationConfig {
            seed_margin: 1,
            ..SimulationConfig::default()
        };

        seed_cells(
            &mut cells,
            &mut cell_map,
            &mut rng,
            4,
            0,
            2,
            test_cell,
            &sim,
        );

        assert_eq!(cells.len(), 4);
        for cell in &cells {
            assert!(cell.pos[0] >= 1 && cell.pos[0] < GRID_X as u16 - 1);
            assert!(cell.pos[1] >= 1 && cell.pos[1] < GRID_Y as u16 - 1);
            assert!(cell.pos[2] < 2);
        }
    }
}
