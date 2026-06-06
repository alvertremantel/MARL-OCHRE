# Toward Evolvable Chemistry

The current MARL chemistry is still mostly index chemistry: external species
and internal species are integer slots, and evolution mutates which slots a
cell transports, senses, secretes, and reacts over. That has already produced
interesting behavior, including an evolved extracellular species-0 energy pool,
but it also exposes the limit of purely slot-based chemistry. If we want richer
outcomes without hand-writing every future ecology, the slots need chemical
meaning that can influence physics and metabolism.

The adjusted plan is to model chemicals as abstract molecules, not literal
real-world molecules. A chemical species should have enough structure to answer
questions like:

- how much usable bond energy does this molecule tend to carry?
- what abstract building blocks is it made from?
- how stable is it outside a cell?
- how expensive or easy is it to move across membranes?
- does it behave like a storage molecule, signal, toxin, structural material,
  redox carrier, or work currency?

This should stay deliberately coarse. We do not need molecular geometry, atom
positions, stereochemistry, or physically exact reaction energetics. We need a
useful ecological abstraction where composition affects interaction with the
environment and where enzymes can evolve to make some transformations more
useful than others.

## Core Abstraction

Each chemical species gets a descriptor:

- `composition`: coarse amounts of abstract components such as carbon backbone,
  oxidizing group, reducing group, phosphate-like activation, lipid-like tail,
  signal group, structural group, and toxin group.
- `bond_energy`: scalar stored energy potential available to metabolism if an
  evolved reaction can couple it productively.
- `extracellular_stability`: how quickly the molecule tends to decay or
  hydrolyze outside cells.
- `membrane_permeability`: how naturally it crosses membranes before transporter
  specialization.
- `storage_density`: how suitable it is as compact reserve material.
- `work_coupling`: how suitable it is as an intracellular energy/work currency.

This differs from a simple moiety ledger. Composition matters, but bond energy
is not just another conserved component. It is a state of the molecule. A
protein/enzyme can make a transformation more or less useful by changing
coupling efficiency, leakage, specificity, and activation barriers.

## Why This Matters

An ATP-equivalent should become useful because it is a good intracellular work
currency: high work coupling, meaningful bond energy, controlled synthesis, and
poor long-term extracellular persistence. A fatty-acid-like molecule should
become useful for different reasons: reduced carbon richness, high storage
density, low diffusion, membrane affinity, and potentially structural value.

Those roles should not be assigned directly as fitness bonuses. They should
fall out of descriptor-driven physics plus evolved enzymes and transporters.

## Implementation Direction

1. Add a built-in chemical descriptor catalog for current external species.
2. Use descriptors in analysis and diagnostics first, so runs can be interpreted
   in terms of energy currencies, storage molecules, signals, toxins, and
   structural materials.
3. Gradually derive more behavior from descriptors:
   - decay and hydrolysis from extracellular stability
   - diffusion/permeability from molecule class
   - light absorption from structural/toxin/organic character
   - transport cost from permeability and size-like composition
4. Let mutation explore descriptor-compatible reaction space:
   - enzymes mutate substrate/product/cofactor as they do now
   - reaction viability is filtered or weighted by descriptor compatibility
   - coupling efficiency determines how much bond energy becomes internal work
   - leakage and byproducts allow cross-feeding and public-good chemistry
5. Eventually allow lineages to repurpose open species slots into newly
   discovered molecule classes, rather than requiring all future chemistry to be
   predeclared by hand.

## First Slice

The first implementation slice is intentionally small: add a descriptor catalog
for the current external species and pin the current decay profile to descriptor
defaults. This does not yet change simulation dynamics. It gives the codebase a
shared vocabulary for future behavior and analysis while keeping the current
model stable.

## Current Status

- Descriptor catalog: implemented in `marl-config`.
- Analysis/diagnostics: implemented in `marl-analysis` and `marl-analyze`; run
  reports now include role labels and descriptor summaries for tracked external
  species.
- Passive physics defaults: diffusion and decay defaults now come from the
  descriptor catalog.
- Membrane crossing cost: accepted transporter flux now pays an internal-energy
  cost derived from each external molecule's permeability and composition load.
  This is configurable through `transport_energy_cost_scale`,
  `transport_permeability_cost_weight`, and
  `transport_composition_cost_weight`; setting the scale to `0.0` restores the
  previous free-transport behavior for control runs.
- Reaction compatibility: non-strict reaction rates are now weighted by
  descriptor-derived substrate/product/cofactor compatibility. Bond energy,
  work-coupling, and composition overlap make some transformations more or less
  efficient without forbidding useless evolved reactions outright. Strict
  stoichiometry deliberately keeps exact template-balanced reaction dynamics and
  skips this descriptor multiplier until the strict template catalog itself is
  descriptor-aware.
- Reaction heat leakage: poorly coupled non-strict reactions now leak internal
  energy to heat according to descriptor mismatch. This makes bad chemistry
  costly without banning it, alongside the material byproduct routing below.
- Reaction byproducts: poorly coupled non-strict reactions now divert a bounded
  fraction of accepted product flux into extracellular byproducts. Routing is
  descriptor-based at a coarse level: signal-like chemistry goes to signal
  pools, structural/storage-heavy chemistry goes to structural material, and
  most carbon/lipid/toxin-bearing failures become organic waste.
- Byproduct analysis: `marl-analyze` now reads v2 stoichiometry events when
  available and reports reaction-byproduct and heat-leakage totals in terminal,
  Markdown, and JSON reports.
- Calibration support: `marl-analyze compare` now carries reaction-byproduct
  and heat-leakage totals across runs, and TOML parsing/validation rejects
  unknown keys and invalid values so typoed calibration settings fail before a
  run starts.
- Byproduct calibration proxies: `marl-analyze` now tracks all external species
  in field snapshots, estimates retained byproduct pools, connects produced
  byproduct species to evolved transporter uptake pressure, and flags
  cross-feeding/public-pool candidates for sweep interpretation.
- Pilot byproduct sweep: a 24x24x12, 160-tick ladder at strengths 0.0, 0.02,
  0.08, and 0.25 showed similar population growth across settings and no
  automated cross-feeding/public-pool candidates. Gross final organic pools were
  misleading because background organic dominated the field; compare reports now
  include zero-baseline-adjusted excess pools. In that pilot, excess retained
  fractions were about 0.325, 0.534, and 0.385 for low/default/high byproduct
  strengths, respectively, so the current default remains plausible but not
  settled.
- Reproducible experiments: `rng_seed` can now be set in `[simulation]`, or via
  `--rng-seed`, to make seeding, mutation, HGT, and division-placement
  randomness deterministic under the recorded run RNG algorithm. Unset seeds are
  resolved from entropy and recorded as concrete replay seeds in `summary.md`.
  Smoke runs confirmed identical tick logs, stoichiometry event logs, and final
  cell binaries for repeated same-seed runs, with different hashes for a
  neighboring seed.
- Replicated comparison support: `marl-analyze compare` now carries recorded
  `rng_seed` metadata into reports and, when multiple zero-byproduct controls
  are present, applies byproduct excess-pool adjustments against the matching
  same-seed baseline instead of the first control run.
- Seed-paired byproduct pilot: a 24x24x12, 241-tick pilot over seeds 41001 and
  41002 compared zero/default/high byproduct strengths in one multi-run compare.
  The new seed-paired adjustment produced low excess retained fractions:
  default 0.040 and 0.026, high 0.126 and 0.072. No automated cross-feeding or
  public-pool candidates were reported, and final populations stayed close to
  their same-seed zero controls. Gross retained fractions remained misleadingly
  high because background organic dominates the field pool.
- Calibration scaling issue: full `stoich_v2_events.csv` logging is now a
  practical bottleneck for replicated sweeps. One 24x24x12, 241-tick
  default-strength pilot run produced a roughly 515 MiB event log. Before the
  larger 4-seed 48x48x24 sweep, add a summary/counter mode or event filtering
  for calibration-relevant stoichiometry so replicated runs are not mostly I/O.
- Compact stoichiometry calibration counters: v2 summaries now include
  reaction-byproduct totals, per-external-species byproduct totals, and
  reaction-leakage totals. `marl-analyze` can use these compact counters when
  `stoich_v2_events.csv` is absent, preserving byproduct calibration reports
  without full row-level event logs. A summary-only 16x16x8, 81-tick smoke run
  confirmed byproduct/leakage analysis and calibration species reporting with no
  event CSV.

The next implementation target is the larger replicated byproduct sweep using
compact v2 summaries by default: run four paired seeds across zero, low,
default, and high byproduct strengths; inspect excess retention, public-pool and
cross-feeding candidates, paired population effects, and whether organic
byproduct accumulation remains low after longer evolutionary time.
