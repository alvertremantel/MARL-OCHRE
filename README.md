# MARL

### Microbial Automata with Reaction-diffusion and Lineage

> marl, old term used to refer to an earthy mixture of fine-grained minerals. The term was applied to a great variety of sediments and rocks with a considerable range of composition. Calcareous marls grade into clays, by diminution in the amount of lime, and into clayey limestones. Greensand marls contain the green, potash-rich mica mineral glauconite; widely distributed along the Atlantic coast in the United States and Europe, they are used as water softeners.

-- *the encyclopedia brittanica, for some reason*

MARL is a CPU-based Rust research prototype for 3D reaction-diffusion cellular automata with simulated lineage. It models sparse microbial cells embedded in a continuous extracellular chemical field. Cells do not interact directly. They interact through transport, secretion, diffusion-limited access to neighboring empty voxels, and a depth-dependent light field.

The current codebase is a Winogradsky-column style simulation: oxidant and carbon are sourced at the top boundary, reductant is sourced at the bottom boundary, and three starter metabolisms are seeded at different depths. Birth, death, quiescence, and mutation emerge from local chemistry and per-cell rulesets rather than from any explicit fitness function.

## What The Simulation Does

- Maintains a dense 3D extracellular field with 12 external chemical species
- Maintains a sparse cell population with 16 internal species and up to 16 reactions per cell
- Diffuses chemistry with cell-body exclusion and local diffusion slowdown from structural deposits
- Computes a per-voxel light field from top-down Beer-Lambert attenuation
- Updates each cell through receptor-gated transport, reaction, effector, and fate phases
- Supports mutation of kinetic parameters and rare structural rewiring of reactions
- Writes compressed binary field arrays and compact viewer cell records, plus optional CSV/PPM diagnostics
- Can optionally write per-layer ruleset-parameter averages or per-cell genotype records on a slower cadence

## Key Crates

- [`marl-engine`](crates/marl-engine/) — simulation library, engine binary, and optional GPU diffusion prototype
- [`marl-analyze`](crates/marl-analyze/) — headless run/comparison analysis CLI
- [`marl-analysis`](crates/marl-analysis/) — reusable analysis/reporting library
- [`marl-viewer-rs`](crates/marl-viewer-rs/) — standalone `wgpu` 3D viewer with `egui` GUI
- [`marl-format`](crates/marl-format/) — shared binary metadata and cell-record schema for engine/viewer interop

A detailed architecture walkthrough lives in [`docs/INFO.md`](docs/INFO.md).

## Important Model Choices

- There is no explicit fitness function.
- Occupied voxels are excluded from diffusion, so dense clusters starve internally.
- Cells sample only empty face-neighbor voxels, not their own occupied voxel.
- Division splits all internal species 50/50 to avoid division as a free-energy exploit.
- Optional horizontal gene transfer copies complete active reaction rules between local neighbors,
  gated by runtime config and evolved `hgt_propensity`.
- Light is used as a catalyst signal inside cells, not as a free energy source.
- Selection is thermodynamic and spatial, driven by local chemistry and access constraints.

## Known In-Progress Or Partial Areas

- Receptors can gate transporter rates through evolvable per-transporter gate weights; reaction gating is still not implemented.
- Horizontal gene transfer is behaviorally wired but disabled by default; it remains intentionally
  stochastic and limited to reaction-rule chunks.
- The code includes spare external and internal species capacity for future chemistry expansion.
- GPU diffusion exists as an optional prototype behind the engine crate's `gpu` feature.

## State

- Language: Rust
- Execution model: CPU default, optional GPU field diffusion prototype
- Default grid: `128 × 128 × 64`
- Transport medium: diffusion only
- Cell occupancy: one cell per voxel
- Light model: Beer-Lambert attenuation from the surface
- Evolution path: vertical mutation during division
- Configuration: runtime TOML + CLI for grid dimensions, physics, and output parameters

This repository is a functional prototype, not a polished platform. The core simulation loop, field physics, lineage tracking, viewer-oriented dumps, run summaries, and interactive viewer are implemented. Several intended extensions are present only in partial form.

## Status Summary

The project is already a real simulation rather than a scaffold. Its current strengths are the field/cell split, the spatial exclusion model, the seeded ecological gradient, receptor-gated transport, optional local HGT, the 3D viewer, and the data products. Its main unfinished areas are reaction-level regulatory wiring, HGT calibration, and broader chemistry expansion.

## Documentation

- **[`docs/USAGE.md`](docs/USAGE.md)** — build, configure, run the engine and viewer, troubleshoot
- **[`docs/INFO.md`](docs/INFO.md)** — deep technical characterization, architecture, tick semantics, spatial model, evolution
- **[`docs/SCRIPTS.md`](docs/SCRIPTS.md)** — project utility scripts

## Quick Start

```bash
# Run a short simulation
cargo run -p marl-engine --release -- --ticks 5000 --stats 100 --snapshot 500

# View the results
cargo run -p marl-viewer-rs --release -- output/run_128x128x64 --tick 0

# Validate a binary viewer dump
python scripts/check_binary_dump.py output/run_128x128x64 0

# Generate headless biological interpretation reports
cargo run -p marl-analyze -- run output/run_128x128x64
```

See [`docs/USAGE.md`](docs/USAGE.md) for comprehensive build and run instructions.

## Citation

If you use this software in your scholarly work, please attribute it with this citation:

> Geosmin Jones. (2026). *Microbial Automata with Reaction-Diffusion and Lineage* (Version 0.2.0.1) [Desktop Software]. Github. https://github.com/alvertremantel/marl

If you're publishing in a big kid venue, you're free to shoot me an email at geosminjones@gmail.com for my real name and Google Scholar profile.
