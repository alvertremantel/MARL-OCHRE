# Seed-Paired Byproduct Calibration, 48x48x24

This experiment is the larger replicated calibration sweep called out in
`CHEMISTRY_PLAN.md`.

## Matrix

- Seeds: `41001`, `41002`, `41003`, `41004`
- Byproduct strengths:
  - `byp000`: `reaction_byproduct_strength = 0.0`
  - `byp002`: `reaction_byproduct_strength = 0.02`
  - `byp008`: `reaction_byproduct_strength = 0.08`
  - `byp025`: `reaction_byproduct_strength = 0.25`
- Grid: `48x48x24`
- Ticks: `801`, with binary field/cell/ruleset snapshots at `0`, `400`, and
  final tick `800`

Each seed has a zero-byproduct control. `marl-analyze compare` pairs each
byproduct run with the zero control that has the same recorded `rng_seed`.

## Output Policy

These configs intentionally set:

```toml
write_stoich_v2_summary = true
write_stoich_v2_events = false
```

The compact v2 summary contains byproduct/leakage totals and per-species
byproduct counters needed for calibration. Full `stoich_v2_events.csv` logs are
not written because the 24x24x12 pilot produced about 515 MiB for one 241-tick
run.

## Run

From the repository root:

```bash
for cfg in experiments/byproduct_calibration_48x24/configs/*.toml; do
  cargo run -p marl-engine --release -- --config "$cfg"
done
```

## Analyze

One all-seed comparison:

```bash
cargo run -p marl-analyze --release -- compare \
  output/byprod_calib_48x24/s41001_byp000 \
  output/byprod_calib_48x24/s41001_byp002 \
  output/byprod_calib_48x24/s41001_byp008 \
  output/byprod_calib_48x24/s41001_byp025 \
  output/byprod_calib_48x24/s41002_byp000 \
  output/byprod_calib_48x24/s41002_byp002 \
  output/byprod_calib_48x24/s41002_byp008 \
  output/byprod_calib_48x24/s41002_byp025 \
  output/byprod_calib_48x24/s41003_byp000 \
  output/byprod_calib_48x24/s41003_byp002 \
  output/byprod_calib_48x24/s41003_byp008 \
  output/byprod_calib_48x24/s41003_byp025 \
  output/byprod_calib_48x24/s41004_byp000 \
  output/byprod_calib_48x24/s41004_byp002 \
  output/byprod_calib_48x24/s41004_byp008 \
  output/byprod_calib_48x24/s41004_byp025 \
  --ticks 0,400,800 \
  --out-dir output/byprod_calib_48x24/reports/compare_all
```

Per-seed comparisons can be produced by passing each seed's four directories to
`marl-analyze compare` with the same `--ticks 0,400,800` flag.

## Primary Readouts

- `reaction_byproduct_amount`
- `reaction_leakage_energy_to_heat`
- `byproduct_excess_retained_fraction`
- `adjusted_uptake_pressure_candidates`
- `adjusted_public_pool_candidates`
- paired final population and growth versus `byp000`

The uptake-pressure candidate count is an indirect proxy: it means a produced
byproduct species has low same-seed baseline-adjusted retention and evolved
uptake pressure in the final rulesets. It is not direct proof of cross-feeding;
that will require flux and lineage attribution in a later analysis slice.
Trace byproduct species remain visible in the adjusted species table, but are
not counted as uptake-pressure or public-pool candidates.

Decision heuristics:

- `byp025` is too high if it repeatedly produces adjusted public-pool
  candidates or clipped excess retained fraction `>= 0.75`.
- `byp008` remains plausible if it produces measurable byproducts without
  repeated adjusted public-pool flags or strong paired population divergence.
- `byp002` is the low-effect guardrail if `byp008` is too noisy or too
  accumulative.
