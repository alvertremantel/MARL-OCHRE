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

The next implementation target is reaction leakage and byproducts: descriptor
compatibility should not only change reaction speed, but also determine how
much potential becomes useful internal work versus heat, waste, or secreted
cross-feeding products.
