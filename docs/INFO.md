# INFO

## Overview

MARL is a 3D reaction-diffusion cellular automaton written in Rust. It is aimed at open-ended microbial evolution in a vertically structured environment rather than at a single fixed game or benchmark. The codebase is organized as a Cargo workspace with three crates. The current implementation is a CPU prototype (with an optional GPU diffusion path) of a Winogradsky-column-like system with sparse cells, a dense extracellular field, top-down light attenuation, and lineage-producing cell division.

Today the code is small enough to read end to end, and it is already coherent. It is also clearly mid-iteration: a few systems are fully implemented, a few are intentionally skeletal, and a few were started and then left disconnected when work moved elsewhere.

## High-Level Architecture

The project has a clean split between environment, cells, orchestration, and outputs.

- `crates/marl-engine/src/config.rs` defines compile-time physics constants and runtime run parameters.
- `crates/marl-engine/src/field.rs` owns the extracellular chemical field and the diffusion solver.
- `crates/marl-engine/src/cell.rs` owns the evolvable cell ruleset, internal state, per-tick update, and mutation logic.
- `crates/marl-engine/src/light.rs` computes a separate light availability field from the current chemistry and occupancy.
- `crates/marl-engine/src/main.rs` ties everything together: seeding, tick order, births, deaths, and output cadence.
- `crates/marl-output/src/data.rs` and `crates/marl-output/src/snapshot.rs` convert state into files for later analysis.
- `crates/marl-engine/src/hgt.rs` contains a horizontal gene transfer primitive that is currently not invoked.
- `crates/marl-output/src/binary_dump.rs` writes raw binary field/cell snapshots and `run_meta.json` for the viewer.
- `crates/marl-format/src/lib.rs` owns the shared binary schema (`RunMeta`, `ViewerCellRecord`, field layout constants).

Conceptually, the simulation loop is:

1. Add boundary source terms.
2. Diffuse the extracellular field with occupied voxels excluded.
3. Recompute the light field.
4. Tick each cell against the chemistry visible from neighboring empty voxels.
5. Apply births and deaths.
6. Log and snapshot.

One important caveat: cells are updated sequentially inside the tick, and each cell's field deltas are applied immediately. Later cells in the iteration therefore see a slightly newer extracellular state than earlier cells.

## State Representation

### Extracellular Field

The extracellular environment is a dense 3D field stored as a flat `Vec<f32>` in `crates/marl-engine/src/field.rs`. Each voxel stores `S_EXT = 12` external species. The current default grid is `128 x 128 x 64`, so the field is calibration-scale rather than the much larger target implied by earlier project notes.

Important external species in current use:

- `0`: unused placeholder in the field layout
- `1`: oxidant
- `2`: reductant
- `3`: carbon
- `4`: organic waste
- `5`, `6`: signal channels reserved but unused by starter metabolisms
- `7`: structural deposit that slows local diffusion
- `8..11`: spare capacity

Boundary sourcing is asymmetric by depth:

- top surface: oxidant and carbon
- bottom surface: reductant

That asymmetry is the main environmental driver of the current ecological setup.

### Cells

Cells live in `Vec<CellState>` storage plus a `HashMap<[u16; 3], usize>` for occupancy lookup. This gives contiguous iteration and constant-time spatial queries.

Each `CellState` contains:

- voxel position
- lineage ID
- age
- `16` internal chemical pools
- an evolvable `Ruleset`
- quiescence flag
- `starter_type` ancestry marker
- division prep countdown

Only one cell may occupy a voxel.

### Rulesets

`Ruleset` is the unit of evolution. It contains:

- receptors
- transporters
- reactions
- effectors
- fate thresholds
- HGT propensity
- mutation rate

This gives the code a good separation between cell identity and cell state: the ruleset expresses what a cell can do, while the internal pools express what it currently has.

## Tick Semantics

The main biological logic is in `CellState::tick`.

### 1. Receptor Pass

Receptors compute Hill-function activations from external concentrations. This machinery is present and documented, but the resulting activation vector is currently unused. In practice, this means the code has sensing primitives but not yet response gating.

This is one of the clearest signs that the project was in the middle of another iteration: the subsystem is not missing, but it is not connected to downstream behavior.

### 2. Transport Pass

Transporters move chemicals between extracellular species and internal pools. Uptake and secretion are saturating functions. Transport is unconditional right now because receptor activation is not yet wired in.

Cells do not read the chemistry in their own occupied voxel. Instead, they average the chemistry of empty face-neighbor voxels. This is an important and deliberate choice: chemicals live in extracellular space, not inside the body-occupied voxel.

### 3. Intracellular Reactions

Cells then run their catalytic network. Reactions are Michaelis-Menten-like, optionally use a cofactor, and include a small epsilon background rate to avoid evolutionary dead ends.

This is a pragmatic research choice rather than a strictly physical one. It makes the search space more navigable for mutation-driven discovery.

Light enters this phase by being written into the last internal slot each tick and then used as a catalyst by light-dependent reactions.

### 4. Maintenance And Effector Pass

After reactions, cells pay maintenance. They also pay a per-active-reaction cost, which is a useful pressure against gratuitously large catalytic networks.

Effectors then secrete selected internal species back into neighboring extracellular voxels, unless the cell is quiescent.

### 5. Fate

Fate decisions are energy-driven.

- too little energy: death
- enough energy: division prep
- sustained prep completion: division event
- low-but-not-dead energy: quiescence

There is also a non-evolvable hard death floor. Cells below `HARD_DEATH_FLOOR = 0.01` die even if their evolved `death_energy` would otherwise permit survival.

Division is intentionally not instantaneous. A cell enters a prep period and pays extra maintenance while preparing to divide. Shorter prep can evolve, but rushing is penalized.

## Spatial Model

The spatial model is one of the strongest parts of the current implementation.

### Occupied Voxels Exclude Diffusion

The field diffusion step treats occupied voxels as excluded space. Chemistry does not diffuse through cells. Occupied neighbors are treated like no-flux boundaries.

This has two major consequences:

- dense biomass creates nutrient shadows and interior starvation
- carrying capacity emerges from geometry and transport limits rather than from an imposed rule

This design is central to the project's behavior and is much more important than several of the more visible but still unwired features.

### Neighbor Exchange

Cells exchange chemistry only with empty face neighbors. If a cell is fully enclosed, it cannot access fresh resources and cannot release waste to the field. In that case it tends to starve.

That means the simulation's notion of crowding is not abstract. It is implemented directly in the geometry of exchange.

## Light Model

`crates/marl-engine/src/light.rs` computes a separate scalar field using a Beer-Lambert style top-down sweep.

Attenuation sources are currently:

- occupied voxels
- organic waste concentration

Light is spatially depth-structured and influences photosynthetic reactions by acting as a catalyst-like input. There is no direct free-energy injection from light in the current implementation. The code documents that design explicitly.

## Seeded Ecologies

The current run seeds three metabolisms in different depth bands.

### Phototrophs

- live near the surface
- use carbon plus light-linked reactions
- produce oxidant and waste

### Chemolithotrophs

- live near the redox interface
- oxidize reductant using oxidant as cofactor support
- start with a carbon reserve to survive while gradients develop

### Anaerobes

- live deeper in the column
- use reductant as their main energy source
- include an oxidant-toxicity mechanism that makes oxygenated environments hostile

These are encoded directly in `crates/marl-engine/src/main.rs` as starter factory functions, not as external data files.

## Evolution And Lineage

Division creates a daughter with:

- a fresh lineage ID
- a mutated copy of the parent ruleset
- half of every internal species pool

The full-pool split matters. It avoids a common exploit where only energy is duplicated and the rest of state is ignored.

Mutation has two levels:

- common parametric perturbations
- rarer structural rewiring

The structural rewiring logic is gene-duplication-inspired, but not a literal full-topology copy. Substrate, product, and catalyst are each sampled independently from active reactions during rare structural mutation, so new reactions are often chimeric combinations assembled from previously active parts rather than exact duplicates of a single donor reaction.

## Data And Outputs

The output side of the project is already quite useful.

### Binary Viewer Outputs

`crates/marl-output/src/binary_dump.rs` writes the high-fidelity binary outputs consumed by the standalone viewer:

- `run_meta.json` — grid dimensions, species counts, field byte length, and cell record stride
- `tick_<T>.field.bin.zst` — losslessly compressed little-endian `f32` field data in `[z][y][x][species]` order by default (`.bin` if compression is disabled)
- `tick_<T>.cells.bin.zst` — losslessly compressed packed 25-byte `ViewerCellRecord` array by default (`.bin` if compression is disabled)
- `tick_<T>.ruleset_layers.bin.zst` — optional per-z-layer averages of continuous ruleset parameters, written on an independent cadence
- `tick_<T>.rulesets.bin.zst` — optional per-cell deduplicated full ruleset dump with dictionary (dict section + per-cell references), written on an independent cadence

The shared schema for these files lives in `crates/marl-format/` so both engine and viewer can reference the same constants without code duplication.

### CSV And Markdown Outputs

`crates/marl-output/src/data.rs` writes:

- `ticks.csv`
- `chem_<tick>.csv`
- `cells_<tick>.csv`
- `reactions_<tick>.csv`
- `reaction_registry.csv`
- `summary.md`

The reaction registry is especially useful because it gives stable IDs to reaction topologies across the run, which makes later lineage and convergence analysis much more tractable.

### Image Outputs

`crates/marl-output/src/snapshot.rs` writes raw PPM images for:

- XZ chemical cross-sections
- XY carbon slices
- cell density slices
- ancestry-colored XZ views

This is intentionally simple and dependency-light.

## Current Technical Characterization

The project is in a good prototype state. It is not a toy, but it is also not yet a fully generalized platform.

### Clearly Implemented And Working

- 3D field storage and diffusion
- cell-body exclusion from diffusion
- light attenuation field
- cell tick loop with transport, reaction, secretion, death, and division prep
- mutation and lineage generation
- binary viewer snapshots (field, cells, metadata) and an interactive `wgpu`/`egui` 3D viewer
- run summaries and useful raw outputs
- a coherent seeded ecological scenario

### Present But Not Fully Integrated

- receptor activations are computed but not used
- HGT transfer logic exists but is not called
- signaling species exist in the chemistry space but are not meaningfully used by starters
- structural deposit species affects diffusion, but current starter metabolisms do not actively build a structural niche

### Stale Or Aspirational Elements

- older descriptions of much larger grid sizes no longer match the actual default configuration
- older GPU-facing intent from early project notes predates the current CPU-first architecture (the optional GPU diffusion path now exists as a prototype behind the `gpu` feature)

## Practical Caveats

Several current simplifications matter if this code is used for serious experimental interpretation.

- Bookkeeping conservation is enforced in the current cell update path: accepted uptake, secretion, reaction flux, cofactor consumption, and division splits are bounded by available pools.
- Physical stoichiometry has two output layers. The v1 `stoich_summary.json` and `stoich_ticks.csv` files remain a legacy-reaction audit. The v2 `stoich_v2_summary.json` and `stoich_v2_events.csv` files audit accepted full-system fluxes across boundary inputs, diffusion/decay, spatial exchange, transport, reactions, maintenance, effectors, division, death, and light availability.
- Strict stoichiometry is opt-in via `stoich_enforcement = "strict"`. Strict mode rejects unbalanced transport/effectors and active reactions that do not match a balanced template catalog, while explicit reservoirs close modeled sources and sinks without adding new species slots.
- Dead cells are removed; their internals are not lysed back into the field.
- Quiescence is partial rather than a deep dormancy mode.
- Cells are updated sequentially with immediate field writes inside each tick.
- Most physics, chemistry, and output parameters are runtime-configurable via TOML + CLI. Grid dimensions and species counts remain compile-time constants (they determine array sizes).
- Unit and integration tests exist for engine field diffusion, binary dump layout, GPU diffusion equivalence, viewer CLI/IO/camera/renderer/GUI, and the shared format crate. Run with `cargo test --workspace`.

These are not necessarily flaws for the present phase, but they define the boundary between prototype behavior and stronger scientific claims.

## Suggested Reading Order In `crates/marl-engine/src/`

If you want to reacquire context quickly, this is the best reading sequence:

1. `crates/marl-engine/src/config.rs`
2. `crates/marl-engine/src/field.rs`
3. `crates/marl-engine/src/cell.rs`
4. `crates/marl-engine/src/main.rs`
5. `crates/marl-engine/src/light.rs`
6. `crates/marl-output/src/data.rs`
7. `crates/marl-output/src/snapshot.rs`
8. `crates/marl-engine/src/hgt.rs`
9. `crates/marl-output/src/binary_dump.rs`
10. `crates/marl-format/src/lib.rs`

That order follows the dependency chain from assumptions, to field physics, to cell logic, to orchestration, then to outputs and unfinished extension points.

## Bottom Line

The codebase today is best described as a coherent CPU research prototype for spatial microbial evolution. Its strongest ideas are already in place: spatial exclusion, chemically mediated interaction, depth-structured ecology, lineage-producing division, and decent analysis outputs. Its most obvious unfinished step is moving from passive chemistry-following cells to cells whose sensing machinery actually modulates behavior, with HGT as a secondary unfinished branch.
