use crate::cell::Ruleset;
use rand::Rng;

/// Transfer a random reaction rule from donor to recipient.
/// This transfers complete reaction rules (not individual parameters),
/// paralleling real bacterial HGT of metabolic operons.
pub fn transfer_reaction(donor: &Ruleset, recipient: &mut Ruleset, rng: &mut impl Rng) -> bool {
    // Reservoir-sample a non-trivial donor reaction without allocating.
    let mut donor_idx = None;
    let mut seen = 0usize;
    for (i, reaction) in donor.reactions.iter().enumerate() {
        if reaction.v_max.abs() <= 1e-9 {
            continue;
        }
        seen += 1;
        if rng.random_range(0..seen) == 0 {
            donor_idx = Some(i);
        }
    }

    let Some(donor_idx) = donor_idx else {
        return false;
    };
    let donated = donor.reactions[donor_idx].clone();

    // Find an inactive slot in recipient (or overwrite random slot)
    let target_idx = (0..recipient.reactions.len())
        .find(|&i| recipient.reactions[i].v_max.abs() < 1e-9)
        .unwrap_or_else(|| rng.random_range(0..recipient.reactions.len()));

    recipient.reactions[target_idx] = donated;
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::{
        EffectorParams, FateParams, Reaction, ReceptorParams, Ruleset, TransportParams,
    };
    use marl_config::{R_MAX, S_EFFECTORS, S_RECEPTORS, S_TRANSPORTERS};
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

    #[test]
    fn transfer_reaction_copies_complete_active_reaction() {
        let mut donor = test_ruleset();
        donor.reactions[3] = Reaction {
            substrate: 2,
            product: 5,
            catalyst: 7,
            cofactor: 1,
            k_m: 0.25,
            v_max: 1.5,
            k_cat: 0.75,
        };
        let mut recipient = test_ruleset();
        let mut rng = StdRng::seed_from_u64(11);

        assert!(transfer_reaction(&donor, &mut recipient, &mut rng));

        let transferred = recipient
            .reactions
            .iter()
            .find(|reaction| reaction.v_max == 1.5)
            .expect("active donor reaction should be copied");
        assert_eq!(transferred.substrate, 2);
        assert_eq!(transferred.product, 5);
        assert_eq!(transferred.catalyst, 7);
        assert_eq!(transferred.cofactor, 1);
        assert_eq!(transferred.k_m, 0.25);
        assert_eq!(transferred.k_cat, 0.75);
    }

    #[test]
    fn transfer_reaction_fails_without_active_donor_rule() {
        let donor = test_ruleset();
        let mut recipient = test_ruleset();
        let before = recipient
            .reactions
            .iter()
            .filter(|r| r.v_max != 0.0)
            .count();
        let mut rng = StdRng::seed_from_u64(12);

        assert!(!transfer_reaction(&donor, &mut recipient, &mut rng));
        assert_eq!(
            recipient
                .reactions
                .iter()
                .filter(|r| r.v_max != 0.0)
                .count(),
            before
        );
    }

    #[test]
    fn transfer_reaction_overwrites_when_recipient_has_no_empty_slots() {
        let mut donor = test_ruleset();
        donor.reactions[0].v_max = 9.0;
        donor.reactions[0].substrate = 4;
        let mut recipient = test_ruleset();
        for i in 0..R_MAX {
            recipient.reactions[i].v_max = i as f32 + 1.0;
            recipient.reactions[i].substrate = 1;
        }
        let mut rng = StdRng::seed_from_u64(13);

        assert!(transfer_reaction(&donor, &mut recipient, &mut rng));
        assert!(
            recipient
                .reactions
                .iter()
                .any(|reaction| reaction.v_max == 9.0 && reaction.substrate == 4)
        );
    }

    #[test]
    fn test_constants_match_ruleset_array_lengths() {
        let ruleset = test_ruleset();
        assert_eq!(ruleset.receptors.len(), S_RECEPTORS);
        assert_eq!(ruleset.transport.len(), S_TRANSPORTERS);
        assert_eq!(ruleset.reactions.len(), R_MAX);
        assert_eq!(ruleset.effectors.len(), S_EFFECTORS);
    }
}
