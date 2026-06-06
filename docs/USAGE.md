# MARL Usage Guide

This document covers everything you need to build, configure, run, and inspect
the MARL simulation and its standalone viewer.

---

## Prerequisites

- **Rust toolchain** (stable, 1.85+). Install via [rustup](https://rustup.rs):
  ```bash
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  ```
- **Python 3.8+** (optional) — for the binary output validation script.
- **A GPU with Vulkan, Metal, or DX12 support** (optional) — for the `wgpu`
  viewer and the optional GPU diffusion path. Both work on CPU-only machines;
  the viewer will use a software adapter if available (`Vulkan` → `llvmpipe`).

## Building

### Full workspace (release)

```bash
cargo build --release --workspace
```

This builds the engine, viewer, analysis CLI, and all library crates under
`target/release/`.

### Engine only

```bash
cargo build -p marl-engine --release
```

To include the **experimental GPU diffusion** path:

```bash
cargo build -p marl-engine --release --features gpu
```

This compiles optional GPU compute shaders for field diffusion. The GPU path is
activated at runtime with `--gpu-diffusion`, not at compile time. Without the
flag, the engine always uses the CPU solver even when compiled with `--features
gpu`.

### Viewer only

```bash
cargo build -p marl-viewer-rs --release
```

### Running without building first

Cargo will build automatically if you use `cargo run`:

```bash
cargo run -p marl-engine --release -- --ticks 100 --stats 10
cargo run -p marl-viewer-rs --release -- output/run_128x128x64 --tick 0
```

---

## Running the Engine

### Quick start

```bash
cargo run -p marl-engine --release -- --ticks 5000 --stats 100 --snapshot 500 --images 500
```

This runs 5000 ticks with the default `128×128×64` grid and three seeded
microbial metabolisms. Stats print to stdout every 100 ticks, binary viewer
dumps write every 500 ticks, and PPM images (if enabled in TOML) write every
500 ticks.

### CLI flags

All run-control parameters can be set on the command line. These override both
built-in defaults and any TOML config file:

| Flag | Description | Default |
|------|-------------|---------|
| `--config <path>` | Path to TOML config file | `marl.toml` (CWD) |
| `--ticks <n>` | Total simulation ticks | 5000 |
| `--stats <n>` | Stdout stats interval (ticks) | 100 |
| `--snapshot <n>` | Binary (and optional CSV) snapshot interval | 500 |
| `--ruleset-interval <n>` | Ruleset sidecar binary interval | 1000 |
| `--images <n>` | PPM image snapshot interval | 500 |
| `--seed <n>` | Cells to seed per starter metabolism | 30 |
| `--rng-seed <n>` | Deterministic RNG seed for reproducible runs | entropy-derived, recorded with the RNG algorithm in `summary.md` |
| `--output <dir>` | Output directory | `output/run_128x128x64` |
| `--gpu-diffusion` | Use GPU diffusion (requires `--features gpu` at build time) | off |

### Grid dimensions

Grid dimensions are runtime-configurable in TOML:

```toml
[grid]
x = 64
y = 64
z = 32
```

Species counts (`S_EXT`, `M_INT`) and ruleset slot counts remain compile-time
because they determine fixed-size cell/ruleset arrays.

Suggested sizes:
- `64×64×32` — quick debug runs (~7 ticks/s)
- `128×128×64` — calibration runs (default)
- `256×256×128` — production runs (needs more than one thread; ~0.1 ticks/s estimated)

GPU diffusion currently uses a shader compiled for the default `128×128×64`
grid. Non-default grids automatically fall back to CPU diffusion.

---

## Runtime Configuration (TOML)

All physics, chemistry, biology, and output parameters are configurable via an
optional TOML file. A sample `marl.toml` with all defaults is included in the
repository root.

### Loading a custom config

```bash
cargo run -p marl-engine --release -- --config myrun.toml
```

Copy `marl.toml` to a new name, edit values, and run with `--config`. Missing
keys fall back to built-in defaults, so partial TOML files work — specify only
what you want to change.

### Configuration hierarchy

1. **Built-in defaults** (hardcoded in `crates/marl-config/src/lib.rs`)
2. **TOML file** (overrides defaults for any keys present)
3. **CLI flags** (override run-control fields: `--ticks`, `--stats`, etc.)

### `[simulation]` section

Core physics and biology parameters:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `dx` | f32 | 0.0001 | Voxel size (100 µm) |
| `dt` | f32 | 1.0 | Ticks per day |
| `diffusion_substeps` | usize | 10 | Diffusion sub-steps per tick |
| `rng_seed` | u64? | unset | Optional deterministic RNG seed; controls seeding, mutation, HGT, and division-placement randomness. When unset, the run derives and records a concrete entropy seed plus RNG algorithm. |
| `d_voxel` | [f32; 12] | [0, 1.5, 1.0, …] | Diffusion coefficients per external species |
| `lambda_decay` | [f32; 12] | [0.2, 0.01, …] | Decay rates per external species (fraction/tick); species 0 is free-energy-like and decays faster by default |
| `source_rate_oxidant` | f32 | 0.4 | Legacy oxidant top-boundary source rate, used only when `boundary_sources` is empty |
| `source_rate_carbon` | f32 | 0.15 | Legacy carbon top-boundary source rate, used only when `boundary_sources` is empty |
| `source_rate_reductant` | f32 | 0.5 | Legacy reductant bottom-boundary source rate, used only when `boundary_sources` is empty |
| `boundary_sources` | list | [] | Optional explicit source list: `{ species, face = "top"/"bottom", rate }`; replaces the legacy source-rate fields when non-empty |
| `epsilon` | f32 | 0.001 | Small background reaction rate |
| `c_max` | f32 | 10.0 | Saturation ceiling for transporter kinetics |
| `lambda_maintenance` | f32 | 0.12 | Base maintenance cost per tick |
| `hard_death_floor` | f32 | 0.01 | Energy below which cells die regardless of evolved threshold |
| `reaction_maintenance` | f32 | 0.003 | Per-active-reaction cost per tick |
| `transport_energy_cost_scale` | f32 | 0.02 | Descriptor-derived internal-energy cost per accepted membrane flux unit; set to `0.0` for legacy free transport |
| `transport_permeability_cost_weight` | f32 | 1.0 | Weight for the cost of moving poorly permeable external molecules |
| `transport_composition_cost_weight` | f32 | 0.5 | Weight for the size-like cost of moving compositionally richer molecules |
| `reaction_descriptor_coupling_strength` | f32 | 0.25 | Strength of descriptor-derived reaction-rate weighting; set to `0.0` for legacy slot-only reaction rates |
| `reaction_descriptor_min_factor` | f32 | 0.25 | Minimum descriptor multiplier for reaction rates |
| `reaction_descriptor_max_factor` | f32 | 1.5 | Maximum descriptor multiplier for reaction rates |
| `reaction_leakage_strength` | f32 | 0.15 | Heat leakage from poorly coupled non-strict reactions |
| `reaction_leakage_max_fraction` | f32 | 0.5 | Maximum leakage fraction of an accepted non-strict reaction flux |
| `reaction_byproduct_strength` | f32 | 0.08 | Product diversion into extracellular byproducts for poorly coupled non-strict reactions |
| `reaction_byproduct_max_fraction` | f32 | 0.25 | Maximum fraction of accepted non-strict reaction flux diverted into byproducts |
| `base_division_prep` | f32 | 20.0 | Tick count for full division prep |
| `prep_maintenance_multiplier` | f32 | 2.0 | Maintenance multiplier during division prep |
| `rush_penalty_rate` | f32 | 0.05 | Penalty for evolving shorter division prep |
| `alpha_eps` | f32 | 0.8 | Niche construction deposit efficiency |
| `k_eps` | f32 | 2.0 | Niche construction saturation constant |
| `light_efficiency` | f32 | 0.0 | Light-to-catalysis efficiency |
| `surface_intensity` | f32 | 1.0 | Light intensity at the top surface |
| `cell_absorption` | f32 | 0.3 | Light attenuation per occupied voxel |
| `chemical_absorption` | f32 | 0.05 | Light attenuation per organic waste unit |
| `light_floor` | f32 | 1e-7 | Minimum light value (prevents zero) |

**Mutation parameters:**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `mutation_stddev` | f32 | 0.1 | Standard deviation of parametric perturbation |
| `structural_mutation_rate_mult` | f32 | 0.1 | Multiplier on rare structural rewiring rate |
| `meta_mutation_rate` | f32 | 0.01 | Rate at which mutation rate itself mutates |
| `meta_mutation_clamp_low` | f32 | 0.001 | Minimum evolvable mutation rate |
| `meta_mutation_clamp_high` | f32 | 0.5 | Maximum evolvable mutation rate |
| `hill_exponent_clamp_low` | f32 | 0.5 | Minimum evolvable Hill coefficient |
| `hill_exponent_clamp_high` | f32 | 8.0 | Maximum evolvable Hill coefficient |
| `active_reaction_threshold` | f32 | 1e-9 | Flux threshold below which a reaction slot counts as inactive |

**Horizontal gene transfer:**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `hgt_enabled` | bool | false | Enables local horizontal gene transfer phase |
| `hgt_interval` | u32 | 10 | Tick cadence for HGT attempts; tick 0 is eligible when enabled; must be >0 when `hgt_enabled = true` |
| `hgt_radius` | u8 | 1 | Local cubic voxel radius for possible donors; 1 includes adjacent and diagonal neighbors; must be >0 when enabled |
| `hgt_base_rate` | f32 | 0.02 | Base per-recipient attempt probability multiplied by evolved `hgt_propensity` |
| `hgt_max_events_per_tick` | usize | 100 | Global cap on successful HGT transfers per tick; must be >0 when enabled |

**Seeding geometry:**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `seed_margin` | u16 | 5 | Minimum distance (voxels) from domain edges |
| `phototroph_z_lo` | f32 | 0.0 | Phototroph seeding band lower bound |
| `phototroph_z_hi` | f32 | 3.0 | Phototroph seeding band upper bound |
| `chemolithotroph_z_lo` | f32 | 80.0 | Chemolithotroph seeding band lower bound |
| `chemolithotroph_z_hi` | f32 | 130.0 | Chemolithotroph seeding band upper bound |
| `anaerobe_z_lo` | f32 | 120.0 | Anaerobe seeding band lower bound |
| `anaerobe_z_hi` | f32 | 180.0 | Anaerobe seeding band upper bound |
| `division_neighbor_distance` | u8 | 2 | Search radius for empty voxels during division |

**Boundary priming:**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `boundary_prime_layers` | usize | 2 | Number of z-layers to prime at boundaries |
| `boundary_prime_oxidant` | f32 | 0.5 | Legacy initial oxidant concentration in top primed layers, used only when `boundary_primes` is empty |
| `boundary_prime_carbon` | f32 | 0.3 | Legacy initial carbon concentration in top primed layers, used only when `boundary_primes` is empty |
| `boundary_prime_reductant` | f32 | 0.5 | Legacy initial reductant concentration in bottom primed layers, used only when `boundary_primes` is empty |
| `boundary_primes` | list | [] | Optional explicit boundary priming list: `{ species, face = "top"/"bottom", concentration }`; replaces the legacy prime fields when non-empty |
| `stoich_enforcement` | enum | `"off"` | `"off"`, `"audit"`, or `"strict"` full-system stoichiometry policy |

External species indices are currently interpreted as:

| Index | Name | Notes |
|---:|---|---|
| 0 | `free_energy` | Extracellular energy-like pool; can evolve through transport/effectors and decays at `0.2` by default |
| 1 | `oxidant` | Default top source |
| 2 | `reductant` | Default bottom source |
| 3 | `carbon` | Default top source |
| 4 | `organic` | Organic waste / absorber |
| 5-6 | `signal_a`, `signal_b` | Reserved signal channels |
| 7 | `structural` | Structural deposit that slows local diffusion |
| 8-11 | `spare_*` | Open chemistry capacity |

Example explicit boundary chemistry:

```toml
[simulation]
boundary_sources = [
  { species = 1, face = "top", rate = 0.4 },
  { species = 3, face = "top", rate = 0.15 },
  { species = 2, face = "bottom", rate = 0.5 },
]

boundary_primes = [
  { species = 1, face = "top", concentration = 0.5 },
  { species = 3, face = "top", concentration = 0.3 },
  { species = 2, face = "bottom", concentration = 0.5 },
]
```

### `[grid]` section

Grid dimensions for the 3D simulation domain:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `x` | usize | 128 | X dimension, in voxels |
| `y` | usize | 128 | Y dimension, in voxels |
| `z` | usize | 64 | Vertical depth, in voxels |

All dimensions must be nonzero and fit the simulation's coordinate math.

### `[output]` section

Output cadence, directories, and format toggles:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `max_ticks` | u32 | 5000 | Total ticks to simulate |
| `stats_interval` | u32 | 100 | Ticks between stdout stats lines |
| `snapshot_interval` | u32 | 500 | Ticks between binary field/cell dumps |
| `image_interval` | u32 | 500 | Ticks between PPM image dumps |
| `seed_count` | usize | 30 | Initial cells per starter metabolism |
| `output_dir` | string | `"output/run_128x128x64"` | Output root directory |
| `write_binary_field` | bool | true | Write `tick_<T>.field.bin.zst` by default |
| `write_binary_cells` | bool | true | Write `tick_<T>.cells.bin.zst` by default |
| `binary_compression` | string | `"zstd"` | Binary payload compression: `"none"` or `"zstd"` |
| `binary_compression_level` | i32 | 3 | Zstd compression level for binary payloads |
| `ruleset_interval` | u32 | 1000 | Ticks between ruleset binary sidecars |
| `ruleset_output_mode` | string | `"off"` | `"off"` / `"layer_averages"` / `"full"` / `"both"` |
| `write_tick_log` | bool | false | Write `ticks.csv` |
| `write_csv_snapshots` | bool | false | Write per-tick CSV snapshots |
| `write_stoich_summary` | bool | false | Write end-of-run physical stoichiometry audit summary |
| `write_stoich_tick_log` | bool | false | Write per-tick physical stoichiometry audit CSV |
| `write_stoich_v2_summary` | bool | false | Write versioned full-system stoichiometry summary |
| `write_stoich_v2_events` | bool | false | Write versioned full-system stoichiometry event CSV |
| `xz_snapshot_species` | [usize] | [] | Species indices for XZ cross-section PPMs |
| `xy_slice_depths_frac` | [f32] | [] | Fractional depths for XY slice PPMs |
| `write_ancestry_map` | bool | false | Write ancestry-colored XZ PPMs |
| `write_density_map` | bool | false | Write cell density PPMs |

---

## Output Files

Runs write into `output_dir` (default: `output/run_128x128x64`).

### Binary viewer outputs

| File | Format | Description |
|------|--------|-------------|
| `run_meta.json` | JSON | Grid dimensions, species counts, binary byte layouts, snapshot intervals, compression mode, and optional ruleset schema metadata |
| `tick_<T>.field.bin.zst` | zstd-compressed raw f32 LE | Full extracellular field in `[z][y][x][species]` order when `write_binary_field = true` |
| `tick_<T>.cells.bin.zst` | zstd-compressed packed binary | Sparse viewer cell records (pos, lineage_id, starter_type, energy) when `write_binary_cells = true`; not full cell state |
| `summary.md` | Markdown | End-of-run population, chemistry, and configuration summary |

### Opt-in via `marl.toml`

| File | Requires | Description |
|------|----------|-------------|
| `ticks.csv` | `write_tick_log = true` | Per-tick population, average energy/enzyme/activity summaries, divisions, deaths, HGT event counts, and z-layer counts |
| `stoich_summary.json` | `write_stoich_summary = true` | End-of-run audit of legacy intracellular reactions, including gross imbalance totals and net signed deltas |
| `stoich_ticks.csv` | `write_stoich_tick_log = true` | Per-tick stoichiometry audit totals with gross imbalance columns plus net signed deltas |
| `stoich_v2_summary.json` | `write_stoich_v2_summary = true` or stoich enforcement enabled | Versioned full-system ledger with stage and reservoir totals |
| `stoich_v2_events.csv` | `write_stoich_v2_events = true` | Per-event full-system stoichiometry rows |
| `tick_<T>.ruleset_layers.bin.zst` | `ruleset_output_mode = "layer_averages"` or `"both"` | Per-z-layer slot-wise averages of continuous ruleset parameters, including transporter gate weights; topology/species IDs are omitted, so slots can mix semantics |
| `tick_<T>.rulesets.bin.zst` | `ruleset_output_mode = "full"` or `"both"` | Deduplicated per-cell ruleset genotype dump with dictionary and positions; current v2 transporter entries include gate receptor/weight; omits lineage, energy, and internal pools unless combined with cell records |
| `chem_<tick>.csv` | `write_csv_snapshots = true` | Full field dump as CSV |
| `cells_<tick>.csv` | `write_csv_snapshots = true` | All cell states as CSV |
| `reactions_<tick>.csv` | `write_csv_snapshots = true` | All active reactions as CSV |
| `reaction_registry.csv` | `write_csv_snapshots = true` | Stable reaction topology IDs across the run |
| `*.ppm` | Various toggles | XZ cross-sections, XY slices, density maps, ancestry maps |

### Validating binary output

```bash
python scripts/check_binary_dump.py output/run_128x128x64 0
```

See [`docs/SCRIPTS.md`](SCRIPTS.md) for details.

### Headless biological analysis

`marl-analyze` is the canonical CLI for non-GUI interpretation of completed
runs. It reuses the Rust binary readers used by the viewer and writes durable
report artifacts instead of living under `scripts/`.

Analyze one run:

```bash
cargo run -p marl-analyze -- run output/run_128x128x64
```

This prints a terminal summary and writes `analysis/analysis.json` plus
`analysis/analysis.md` inside the run directory.

Compare runs:

```bash
cargo run -p marl-analyze -- compare output/run_a output/run_b --out-dir output/compare_a_b
```

When compared runs include v2 event rows, the comparison report includes
reaction-byproduct and reaction-leakage totals so chemistry calibration sweeps
can be inspected headlessly. It also reports byproduct calibration proxies:
final byproduct pool, retained fraction, cross-feeding candidates, and
public-pool candidates. A cross-feeding candidate means a produced byproduct
species has evolved uptake pressure and low final retention; a public-pool
candidate means it accumulated without detected uptake pressure.
The byproduct calibration section reports the field and ruleset ticks used for
those proxies. Aggregate retained fraction is omitted when any produced
byproduct species is missing from the field snapshot coverage. Single-run
retained fractions are gross field-pool ratios; comparison reports also compute
control-adjusted excess pools when a comparable zero-byproduct baseline run is
included. The baseline must have v2 event data, a chemistry snapshot, matching
grid dimensions, and the same calibration field tick.
For replicated calibration sweeps, keep the zero-byproduct baseline and each
byproduct setting paired by `rng_seed` so stochastic seeding and mutation are
controlled across the comparison.
If `rng_seed` is omitted for exploratory work, copy the resolved seed from
`summary.md` into follow-up configs before comparing variants, and keep the
recorded RNG algorithm with the run notes.

By default, analysis reads the full `ticks.csv` trajectory and samples the
first, middle, and latest binary snapshots. Use `--all-snapshots`,
`--latest-only`, or `--ticks 0,500,5000` to change that policy. Use
`--no-rulesets` when full ruleset sidecars are unavailable or not needed.
When `stoich_v2_summary.json` is present, the same reports include v2
stoichiometry metadata. Descriptor-driven `reaction_byproduct` and
`reaction_leakage` amounts require `stoich_v2_events.csv`.

When full ruleset sidecars are available, `marl-analyze` also reports
transporter ecology:

- population-weighted active transporter slots and active slots per cell
- uptake-dominant, secretion-dominant, and bidirectional active slots
- receptor-gated active slots and average absolute `gate_weight`
- common `(ext_species, int_species)` transporter pairs with average uptake and secretion rates
- dominant genotype transporter counts for genotype-level story hooks

Current full ruleset dumps are v2/576-byte payloads. `marl-analyze`,
`check_binary_dump.py`, and `inspect_rulesets.py` also support legacy
v1/536-byte payloads and report those transporters as ungated.

HGT event counts are currently written to `ticks.csv` as
`hgt_events_this_tick`. `marl-analyze` parses the trajectory but does not yet
surface HGT counts in its terminal or Markdown summaries.

### Stoichiometry audit outputs

- `stoich_ticks.csv` columns: `reaction_count`, `active_flux`, `imbalanced_reaction_count`, `unknown_species_flux`, `carbon_to_energy_flux`, `reductant_to_energy_flux`, `gross_material_abs`, `gross_total_abs`, then net signed `delta_*` budget columns.
- `stoich_summary.json` reports the same run-level ledger, with `material_abs_sum` as the gross material imbalance magnitude, `gross_total_abs_sum` as the gross material plus redox plus energy magnitude, and `net_material_abs_sum` as the final signed material imbalance after cancellation.
- v1 files cover legacy intracellular reactions. `stoich_v2_summary.json` and `stoich_v2_events.csv` cover boundary inputs, diffusion/decay, light availability, transport, spatial exchange, intracellular reactions, maintenance, effectors, division, and death.
- `marl-analyze` summarizes `reaction_byproduct` totals by external species, final retained byproduct pools from binary field snapshots, transporter uptake pressure from full ruleset sidecars, and `reaction_leakage` energy sent to heat when v2 event rows are available.
- `stoich_enforcement = "audit"` records the v2 ledger without changing dynamics. `stoich_enforcement = "strict"` keeps the feature opt-in, rejects unbalanced transport/effectors and untemplated active reactions, and closes modeled sources/sinks through explicit reservoirs.
- Strict mode also constrains structural reaction mutations to the balanced template catalog, so strict evolutionary proposal distributions are intentionally different from legacy runs.

### Compression and viewer compatibility

- Binary field/cell payloads are **zstd-compressed by default**.
- Set `binary_compression = "none"` to restore legacy raw `.bin` files.
- The viewer auto-detects both raw and `.zst` snapshots from `run_meta.json`.
- `field_byte_len` in `run_meta.json` always describes the **decompressed** field size.
- The viewer requires `write_binary_field = true`; cell overlays are skipped when `write_binary_cells = false`.
- Binary field, cell, and ruleset files are viewer/analysis records. They do not contain enough state to resume a simulation.

---

## Running the Viewer

The standalone viewer renders 3D isometric volumes with direct cell voxel
overlay, or legacy top-down field-only projections.

### Quick start

```bash
cargo run -p marl-viewer-rs --release -- output/run_128x128x64 --tick 0
```

### CLI flags

| Flag | Values | Default | Description |
|------|--------|---------|-------------|
| `--dir <path>` | path | (positional) | Output directory containing `run_meta.json` |
| `--tick <n>` | integer | 0 | Snapshot tick to load |
| `--species <n>` | integer | 1 | External chemical species to render |
| `--view <mode>` | `iso`, `top` | `iso` | Isometric volume or top-down projection |
| `--cells <mode>` | `off`, `starter`, `energy` | `starter` | Cell coloring mode |
| `--cell-alpha <f>` | (0, 1] | 0.95 | Opacity of cell voxel markers |
| `--scale <f>` | float | 2.0 | Concentration-to-density transfer function scale |
| `--exposure <f>` | float | 18.0 | Raymarch opacity multiplier |
| `--steps <n>` | integer | 160 | Raymarch sample count |

### Microbe coloring

Cell rendering uses the `starter_type` field from cell records — an ancestry
category from the initial seeding, **not** an inferred genotype-level species:

- **Red** — phototroph
- **Green** — chemolithotroph
- **Blue** — anaerobe
- **Magenta** — other/unknown

Cell rendering requires `write_binary_cells = true` in the output config
(enabled by default). The viewer still opens field-only snapshots when cell
output is disabled.

### Legacy top-down field-only rendering

```bash
cargo run -p marl-viewer-rs --release -- output/run_128x128x64 --tick 0 --view top --cells off --species 1
```

This renders a single chemical species as a top-down z-projection without cell
markers, matching the pre-isometric viewer behavior.

---

## Viewer GUI

The viewer includes an `egui` GUI shell overlaid on the 3D render. It opens
automatically — no extra flags needed.

### Controls

**Directory toolbar (top):**
- Text field — enter or edit the output directory path.
- `Open…` — native folder picker dialog (platform-dependent; falls back
  gracefully if unavailable).
- `Load Dir` — loads the entered directory, discovers available snapshot ticks,
  and opens the first available tick.
- `Reload` — rescans the current directory for new ticks (useful when a
  simulation is still running and writing new snapshots).

**Tick navigation:**
- Numeric tick entry + `Go` button.
- `First` / `Prev` / `Next` / `Last` buttons for stepping through available
  snapshot ticks.

**View Settings (collapsible side panel):**
- **Species** — which extracellular chemical species to render.
- **View mode** — isometric or top-down.
- **Cell mode** — off, starter-colored, or energy-colored.
- **Cell alpha** — opacity of cell voxel markers.
- **Density scale** — concentration-to-density mapping.
- **Exposure** — opacity multiplier for the raymarch.
- **Raymarch steps** — number of samples along each ray.
- `Apply` — reloads the snapshot with current settings.
- `Reset` — returns all settings to their defaults.

Changes are not live; click `Apply` to reload the snapshot with the new
settings.

### Startup behavior

If the viewer is launched with a missing or invalid output directory, it opens
a window with a 1×1×1 placeholder render and shows the GUI so you can navigate
to a valid directory. The window title will show "no snapshot loaded" in this
state.

---

## Common Workflows

### Parameter sweep

1. Copy `marl.toml`:
   ```bash
   cp marl.toml sweep_mutation.toml
   ```
2. Edit one value (e.g., `mutation_stddev = 0.05`).
3. Run with a distinct output directory:
   ```bash
   cargo run -p marl-engine --release -- --config sweep_mutation.toml --output output/sweep_mut_0.05
   ```

Repeat for different values. The `summary.md` in each output directory
records the final state.

### Running and viewing simultaneously

1. Start the engine in one terminal:
   ```bash
   cargo run -p marl-engine --release -- --ticks 5000 --snapshot 500 --output output/my_run
   ```
2. Open the viewer in another terminal:
   ```bash
   cargo run -p marl-viewer-rs --release -- output/my_run --tick 0
   ```
3. Use `Reload` in the viewer GUI to pick up new ticks as the engine writes
   them.

### Validating output quickly

```bash
python scripts/check_binary_dump.py output/run_128x128x64 0
python scripts/check_binary_dump.py output/run_128x128x64 500
```

### Inspecting chemistry as images

Enable a few XZ cross-sections in `marl.toml`:

```toml
[output]
xz_snapshot_species = [1, 3, 4]   # oxidant, carbon, waste
image_interval = 100
```

This produces `oxidant_xz_<tick>.ppm`, `carbon_xz_<tick>.ppm`, and
`organic_waste_xz_<tick>.ppm` every 100 ticks. PPM files open in most image
viewers or can be converted with ImageMagick:

```bash
convert oxidant_xz_500.ppm oxidant_xz_500.png
```

### Changing grid size

Set the `[grid]` section in your config file and choose a matching
`output_dir`, for example:

```toml
[grid]
x = 64
y = 64
z = 32

[output]
output_dir = "output/run_64x64x32"
```

Binary metadata and viewer dimensions are written from the runtime grid in
`run_meta.json`.

---

## Troubleshooting

### "Adapter not found" / black window in viewer

The viewer requires a `wgpu`-compatible backend (Vulkan, Metal, DX12). If your
system doesn't have GPU drivers, try:
- On Linux, install `mesa-vulkan-drivers` or `vulkan-tools` for software
  rendering via `llvmpipe`.
- Verify with `vulkaninfo | grep deviceName`.

### Viewer says "no snapshot loaded"

The output directory doesn't contain `run_meta.json`, `write_binary_field` is
disabled, or the metadata's `field_file_pattern` does not resolve to an
existing snapshot for the selected tick. Verify the path and that the engine has
written at least one binary field snapshot.
Use the GUI `Open…` button or text field to navigate to a valid directory.

### "field size mismatch"

The selected field file does not match the dimensions in `run_meta.json`. Verify
that the output directory and tick belong to the same run and were not mixed with
files from another grid size.

### `rfd` folder picker does nothing

The native file dialog (`rfd` crate) may not be available on all platforms or
desktop environments. Use the text field to type or paste the output directory
path, then click `Load Dir`.

### Snapshot ticks don't update while engine is running

Click `Reload` in the viewer GUI to rescan the output directory for new tick
files. The viewer does not watch the filesystem automatically.

### TOML config is rejected

Check your TOML syntax and section names. Unknown keys are rejected, including
fields in the wrong section, so typoed calibration settings cannot silently fall
back to defaults. Run with only the config change to isolate:

```bash
cargo run -p marl-engine --release -- --config my.toml --ticks 10 --stats 1
```

On TOML parse or validation failure the engine prints an error and exits before
starting the run.

### Mutation standard deviation is rejected

Setting `mutation_stddev = 0.0`, a negative value, or a non-finite value is
rejected before the run starts. Keep `mutation_stddev > 0.0`.
