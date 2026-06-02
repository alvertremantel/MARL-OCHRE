use crate::cell::Ruleset;
use rand::Rng;

/// Transfer a random reaction rule from donor to recipient.
/// This transfers complete reaction rules (not individual parameters),
/// paralleling real bacterial HGT of metabolic operons.
pub fn transfer_reaction(donor: &Ruleset, recipient: &mut Ruleset, rng: &mut impl Rng) {
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
        return;
    };
    let donated = donor.reactions[donor_idx].clone();

    // Find an inactive slot in recipient (or overwrite random slot)
    let target_idx = (0..recipient.reactions.len())
        .find(|&i| recipient.reactions[i].v_max.abs() < 1e-9)
        .unwrap_or_else(|| rng.random_range(0..recipient.reactions.len()));

    recipient.reactions[target_idx] = donated;
}
