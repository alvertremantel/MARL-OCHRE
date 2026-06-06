# INFO

## Overview

MARL is a 3D reaction-diffusion cellular automaton written in Rust. It is aimed at open-ended microbial evolution in a vertically structured environment rather than at a single fixed game or benchmark. The codebase is organized as a multi-crate Cargo workspace. The current implementation is a CPU-first prototype, with optional GPU diffusion support behind a feature, of a Winogradsky-column-like system with sparse cells, a dense extracellular field, top-down light attenuation, lineage-producing cell division, receptor-gated transport, and optional local HGT.

Today the code is small enough to read end to end, and it is already coherent. It is also clearly mid-iteration: a few systems are fully implemented, a few are intentionally skeletal, and a few were started and then left disconnected when work moved elsewhere.

## High-Level Architecture

The project now has a deliberate split between configuration, field storage, cell biology, orchestration, outputs, analysis, and viewing.

- `crates/marl-engine/src/main.rs` is a thin binary entry point: it loads config and calls `marl_sim::run`.
- `crates/marl-config/src/lib.rs` defines runtime grid, physics, HGT, and output configuration plus fixed species/ruleset array sizes.
- `crates/marl-field/src/field.rs` owns the runtime-sized extracellular chemical field and CPU diffusion solver.
- `crates/marl-field/src/light.rs` computes a separate light availability field from current chemistry and occupancy.
- `crates/marl-cell/src/cell.rs` owns cell state, evolvable rulesets, receptor-gated transport, intracellular reactions, fate decisions, and mutation logic.
- `crates/marl-cell/src/hgt.rs` contains the reaction-rule HGT primitive used by the optional local HGT phase.
- `crates/marl-sim/src/lib.rs` orchestrates the tick loop: boundary sources, diffusion, light, cell ticks, births/deaths, HGT, logging, and snapshots.
- `crates/marl-sim/src/starter_metabolisms.rs`, `seeding.rs`, `spatial.rs`, and `stats.rs` own starter cell construction, initial placement, local exchange/neighborhood searches, and stdout summaries.
- `crates/marl-output/src/data.rs`, `binary_dump.rs`, and `snapshot.rs` convert state into CSV/Markdown, compressed binary viewer/analysis records, and raw PPM images.
- `crates/marl-format/src/lib.rs` owns the shared binary schema (`RunMeta`, `ViewerCellRecord`, field layout constants, ruleset dump constants).
- `crates/marl-analysis` and `crates/marl-analyze` provide the canonical headless analysis workflow for completed runs.
- `crates/marl-viewer-core`, `crates/marl-viewer-render`, and `crates/marl-viewer-rs` split viewer argument/IO types, renderer/GUI, and the standalone viewer binary.
- `crates/marl-gpu` contains the optional GPU diffusion prototype used by `marl-engine --features gpu -- --gpu-diffusion`.

Conceptually, the simulation loop is:

1. Add boundary source terms.
2. Diffuse the extracellular field with occupied voxels excluded.
3. Recompute the light field.
4. Tick each cell against the chemistry visible from neighboring empty voxels.
5. Apply births and deaths.
6. Run optional local HGT if enabled and due.
7. Log and snapshot.

One important caveat: cells are updated sequentially inside the tick, and each cell's field deltas are applied immediately. Later cells in the iteration therefore see a slightly newer extracellular state than earlier cells.

## State Representation

### Extracellular Field

The extracellular environment is a dense 3D field stored as a flat `Vec<f32>` in `crates/marl-field/src/field.rs`. Each voxel stores `S_EXT = 12` external species. Grid dimensions are runtime-configurable through `[grid]` in TOML; the default is `128 x 128 x 64`.

Important external species in current use:

- `0`: free-energy-like extracellular pool; not boundary-sourced by default and now decays faster than carbon
- `1`: oxidant
- `2`: reductant
- `3`: carbon
- `4`: organic waste
- `5`, `6`: signal channels reserved but unused by starter metabolisms
- `7`: structural deposit that slows local diffusion
- `8..11`: spare capacity

Boundary sourcing is asymmetric by depth by default:

- top surface: oxidant and carbon
- bottom surface: reductant

That asymmetry is the main environmental driver of the current ecological setup.
The TOML config can also provide explicit `boundary_sources` and
`boundary_primes` lists, allowing any external species to be sourced or primed
on the top or bottom face without adding new hardcoded Rust fields.

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

Receptors compute bounded Hill-function activations from external concentrations. Each transporter has an evolvable `gate_receptor` selector and signed `gate_weight`. The transport pass multiplies that transporter's uptake and secretion rates by:

```text
clamp(1 + gate_weight * activation[gate_receptor], 0, 4)
```

A zero `gate_weight` keeps transport unconditional. Positive weights amplify transport when the selected receptor is active; negative weights suppress it and can shut a transporter off. This is not a fitness bonus or survival guard. It only changes requested membrane flux before the existing external availability, internal-pool, and `c_max` caps are applied.

### 2. Transport Pass

Transporters move chemicals between extracellular species and internal pools. Uptake and secretion are saturating functions, then receptor-gated through the per-transporter factor above. Accepted membrane flux pays an internal-energy cost derived from the external molecule descriptor: poorly permeable and compositionally rich molecules are harder to move. Starter metabolisms use neutral gates, so they begin with the same unconditional transport behavior and must evolve useful gating through mutation.

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

Cells exchange chemistry only with empty face neighbors. If a cell is fully enclosed, it cannot access fresh resources and cannot release waste to the field. In that case it tends to starve. Optional HGT is separate from this chemistry exchange path: it uses a local cubic neighborhood, can include diagonal neighbors, and copies reaction rules rather than moving chemical mass.

That means the simulation's notion of crowding is not abstract. It is implemented directly in the geometry of exchange.

## Light Model

`crates/marl-field/src/light.rs` computes a separate scalar field using a Beer-Lambert style top-down sweep.

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

These are encoded directly in `crates/marl-sim/src/starter_metabolisms.rs` as starter factory functions, not as external data files.

## Evolution And Lineage

Division creates a daughter with:

- a fresh lineage ID
- a mutated copy of the parent ruleset
- half of every internal species pool

The full-pool split matters. It avoids a common exploit where only energy is duplicated and the rest of state is ignored.

Mutation has two levels:

- common parametric perturbations
- rarer structural rewiring

Transport rates, gate weights, receptors, reactions, effectors, fate thresholds, HGT propensity, and mutation rate can all evolve. Transporter structural mutation can rewire external species, internal species, and the gate receptor selector. Reaction structural mutation is gene-duplication-inspired, but not a literal full-topology copy in legacy mode: substrate, product, and catalyst are each sampled independently from active reactions during rare structural mutation, so new reactions are often chimeric combinations assembled from previously active parts rather than exact duplicates of a single donor reaction. Strict stoichiometry mode instead draws structural reaction mutations from balanced templates.

Optional HGT is a separate horizontal path. When enabled, the sim runs a local neighborhood phase on its own cadence. A recipient's finite, positive `hgt_propensity` is multiplied by `hgt_base_rate`; successful transfers copy one complete active reaction rule from a local donor ruleset snapshot into the recipient. HGT is capped per tick and per recipient, disabled by default, and logged as an event count in `ticks.csv`.

## Data And Outputs

The output side of the project is already quite useful.

### Binary Viewer Outputs

`crates/marl-output/src/binary_dump.rs` writes binary outputs consumed by the standalone viewer and analysis scripts. These files are data products, not restartable simulation checkpoints:

- `run_meta.json` — grid dimensions, species counts, field byte length, and cell record stride
- `tick_<T>.field.bin.zst` — compressed little-endian `f32` extracellular field data in `[z][y][x][species]` order by default (`.bin` if compression is disabled)
- `tick_<T>.cells.bin.zst` — compressed packed 25-byte `ViewerCellRecord` array by default (`.bin` if compression is disabled); records contain position, lineage ID, starter type, and energy only
- `tick_<T>.ruleset_layers.bin.zst` — optional per-z-layer slot-wise averages of continuous ruleset parameters, including transporter gate weights, written on an independent cadence; topology/species IDs are omitted, so averaged slots can mix different reaction or transport semantics
- `tick_<T>.rulesets.bin.zst` — optional per-cell deduplicated ruleset genotype dump with dictionary (dict section + per-cell position references), written on an independent cadence; current v2 transporter records include `uptake_rate`, `secrete_rate`, `ext_species`, `int_species`, `gate_receptor`, and `gate_weight`; it omits lineage, energy, and internal pools unless interpreted alongside the cell dump. The Rust analysis path and utility scripts also understand legacy v1/536-byte dumps, treating transporters as ungated.

The shared schema for these files lives in `crates/marl-format/` so both engine and viewer can reference the same constants without code duplication.

### CSV And Markdown Outputs

`crates/marl-output/src/data.rs` writes:

- `ticks.csv`
- `chem_<tick>.csv`
- `cells_<tick>.csv`
- `reactions_<tick>.csv`
- `reaction_registry.csv`
- `summary.md`

`ticks.csv` includes population, energy, active reaction averages, per-tick divisions/deaths, HGT event counts, and per-z-layer population counts. The reaction registry is especially useful because it gives stable IDs to reaction topologies across the run, which makes later lineage and convergence analysis much more tractable.

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
- cell tick loop with receptor-gated transport, reaction, secretion, death, and division prep
- mutation and lineage generation
- optional local HGT with runtime guardrails and tick logging
- binary viewer records (field, compact cells, metadata) and an interactive `wgpu`/`egui` 3D viewer
- headless analysis for trajectories, chemistry, zonation, genotype diversity, and transporter ecology
- run summaries and useful raw outputs
- a coherent seeded ecological scenario

### Present But Not Fully Integrated

- receptor activations gate transporter uptake/secretion but do not yet gate reactions
- optional HGT is wired into the tick loop, but disabled by default and still experimental
- signaling species exist in the chemistry space but are not meaningfully used by starters
- structural deposit species affects diffusion, but current starter metabolisms do not actively build a structural niche
- `marl-analyze` reads genotype/ruleset and transporter summaries, but does not yet surface HGT event counts from `ticks.csv`

### Stale Or Aspirational Elements

- older descriptions of much larger grid sizes no longer match the current runtime-configurable default
- older GPU-facing intent from early project notes predates the current CPU-first architecture (the optional GPU diffusion path now exists as a prototype behind the `gpu` feature)

## Practical Caveats

Several current simplifications matter if this code is used for serious experimental interpretation.

- Bookkeeping conservation is enforced in the current cell update path: accepted uptake, secretion, reaction flux, cofactor consumption, and division splits are bounded by available pools.
- Physical stoichiometry has two output layers. The v1 `stoich_summary.json` and `stoich_ticks.csv` files remain a legacy-reaction audit. The v2 `stoich_v2_summary.json` and `stoich_v2_events.csv` files audit accepted full-system fluxes across boundary inputs, diffusion/decay, spatial exchange, transport, reactions, maintenance, effectors, division, death, and light availability.
- Strict stoichiometry is opt-in via `stoich_enforcement = "strict"`. Strict mode rejects unbalanced transport/effectors and active reactions that do not match a balanced template catalog, while explicit reservoirs close modeled sources and sinks without adding new species slots.
- External species 0 is not a boundary input, but organisms can evolve to export and re-import internal energy through it. It decays by default, and headless analysis now reports when it accumulates. Treat it as a real modeled pool, not an unused spare.
- Dead cells are removed; their internals are not lysed back into the field.
- Quiescence is partial rather than a deep dormancy mode.
- Cells are updated sequentially with immediate field writes inside each tick.
- Grid dimensions, physics, chemistry, and output parameters are runtime-configurable via TOML + CLI. Species counts and ruleset slot counts remain compile-time constants because they determine fixed-size cell/ruleset arrays.
- Unit and integration tests exist for config parsing, field/light behavior, cell transport/reaction/HGT logic, sim orchestration, binary dump layout, analysis parsing, GPU diffusion equivalence, viewer CLI/IO/camera/renderer/GUI, and the shared format crate. Run with `cargo test --workspace`.

These are not necessarily flaws for the present phase, but they define the boundary between prototype behavior and stronger scientific claims.

## Suggested Reading Order

If you want to reacquire context quickly, this is the best reading sequence:

1. `crates/marl-config/src/lib.rs`
2. `crates/marl-field/src/field.rs`
3. `crates/marl-field/src/light.rs`
4. `crates/marl-cell/src/cell.rs`
5. `crates/marl-cell/src/hgt.rs`
6. `crates/marl-sim/src/lib.rs`
7. `crates/marl-sim/src/spatial.rs`
8. `crates/marl-sim/src/starter_metabolisms.rs`
9. `crates/marl-output/src/data.rs`
10. `crates/marl-output/src/binary_dump.rs`
11. `crates/marl-format/src/lib.rs`
12. `crates/marl-analysis/src/lib.rs`
13. `crates/marl-viewer-core/src/io.rs`

That order follows the dependency chain from assumptions, to field physics, to cell logic, orchestration, binary outputs, analysis, and viewer consumption.

## Bottom Line

The codebase today is best described as a coherent CPU research prototype for spatial microbial evolution. Its strongest ideas are already in place: spatial exclusion, chemically mediated interaction, depth-structured ecology, lineage-producing division, receptor-gated transport, optional local HGT, and decent analysis outputs. Its most obvious unfinished steps are extending regulatory wiring beyond transport, calibrating HGT, and making broader chemistry/signaling stories easier to analyze.
