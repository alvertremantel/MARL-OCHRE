# MARL Utility Scripts

This document describes utility scripts in the `scripts/` directory.

---

## `check_binary_dump.py`

**Purpose:** Sanity-checks MARL binary output files (field dump, cell dump,
optional ruleset-layer dump, optional full ruleset dump, and metadata) for a given tick.

**Requirements:** Python 3.6+; for compressed `.zst` snapshots, the `zstd` CLI
must also be available on `PATH`.

### Usage

```bash
python scripts/check_binary_dump.py <run_dir> <tick>
python scripts/check_binary_dump.py <run_dir> <tick> --require-rulesets
python scripts/check_binary_dump.py <run_dir> <tick> --require-full-rulesets
```

### Arguments

| Argument | Description |
|----------|-------------|
| `run_dir` | Path to the engine output directory containing `run_meta.json`, field snapshots, and cell snapshots when `write_binary_cells = true` |
| `tick` | Tick number to inspect (e.g., `0`, `500`, `1000`) |
| `--require-rulesets` | Also require and validate `tick_<N>.ruleset_layers.bin(.zst)` using metadata stride/count fields |
| `--require-full-rulesets` | Also require and validate `tick_<N>.rulesets.bin(.zst)` using metadata header/cell-ref/dict fields |

### Checks Performed

1. **Metadata:** Reads `run_meta.json` and validates the shared binary schema
   fields (`endianness`, `field_dtype`, `field_layout`, snapshot patterns, and
   `cell_record_stride`).
2. **Field file size:** Verifies that the decompressed `tick_<T>.field.bin(.zst)` payload has exactly `field_byte_len` bytes.
3. **First field value:** Reads the first `f32` (little-endian) from the field
   file and confirms it is a finite number (not `NaN` or `inf`).
4. **Cell file integrity:** When `write_binary_cells = true`, verifies that the decompressed `tick_<T>.cells.bin(.zst)` payload size is evenly divisible by `cell_record_stride`; otherwise reports zero cells without requiring a cell file.
5. **Ruleset layer integrity (optional):** When enabled in metadata or via `--require-rulesets`, verifies that `tick_<T>.ruleset_layers.bin(.zst)` contains exactly `grid_z` fixed-stride records.
6. **Full ruleset integrity (optional):** When enabled in metadata or via `--require-full-rulesets`, verifies:
   - Magic bytes (`MRSF`) and format version match
   - File size matches `header + dict_count × ruleset_byte_size + cell_count × cell_ref_stride`
   - All per-cell dict_id references are within dictionary bounds
   - All cell positions are within grid bounds

### Output

On success, prints a summary line and exits with code 0:

```
ok: first_f32=<value>, field_bytes=<n>, cell_count=<n>
```

On failure, exits non-zero with an error message describing the mismatch.

### Example

```bash
$ python scripts/check_binary_dump.py output/run_128x128x64 0
ok: first_f32=0, field_bytes=50331648, cell_count=90

$ python scripts/check_binary_dump.py output/run_128x128x64 500
ok: first_f32=0.0034521, field_bytes=50331648, cell_count=127

$ python scripts/check_binary_dump.py output/run_128x128x64 500 --require-rulesets
ok: first_f32=0.0034521, field_bytes=50331648, cell_count=127

# Error example — tick without snapshot
$ python scripts/check_binary_dump.py output/run_128x128x64 100
FileNotFoundError: [Errno 2] No such file or directory: 'output/run_128x128x64/tick_100.field.bin'
```

### Typical Use

Validate output after an engine run, or while a run is in progress to confirm
data integrity:

```bash
# Check the latest snapshot
ls -t output/run_128x128x64/tick_*.field.bin* | head -1
python scripts/check_binary_dump.py output/run_128x128x64 5000
```

---

## `inspect_rulesets.py`

**Purpose:** Inspect a MARL full deduplicated ruleset dump file
(`tick_<T>.rulesets.bin` or `tick_<T>.rulesets.bin.zst`). Prints a
human-readable summary of the on-disk structure, dictionary usage statistics,
and a decoded glimpse of one ruleset entry.

**Requirements:** Python 3.6+; for compressed `.zst` files, the `zstd` CLI
must also be available on `PATH`.

### Usage

```bash
python scripts/inspect_rulesets.py <path_to_ruleset_file>
python scripts/inspect_rulesets.py <path> --no-glimpse
python scripts/inspect_rulesets.py <path> --no-cell-sample
```

### Arguments

| Argument | Description |
|----------|-------------|
| `file` | Path to a `tick_<N>.rulesets.bin` or `.bin.zst` file |
| `--no-glimpse` | Skip the decoded ruleset entry glimpse (header + stats only) |
| `--no-cell-sample` | Skip the per-cell ref position/dict_id sampling |
| `--full-histogram` | Show complete dict usage histogram instead of top 10 |

### Checks Performed

1. **Magic bytes:** Verifies the file begins with `MRSF`.
2. **Format version:** Confirms version matches the expected value (currently 1).
3. **Ruleset byte size:** Confirms each canonical ruleset payload is exactly 536 bytes.
4. **File size:** Validates total size matches `header(24) + dict_count × ruleset_size + cell_count × cell_ref_stride(10)`.
5. **Dict ID range:** Every per-cell `dict_id` reference must be within dictionary bounds.

### Output

Prints a structured report containing:

- **Header:** magic, version, dict/cell count, ruleset byte size, uniqueness ratio
- **Layout:** byte breakdown by section (header, dictionary, per-cell refs)
- **Dict Usage:** frequency histogram showing which dict entries are most common
- **Cell refs:** sample of per-cell positions and their dict_id references
- **Ruleset Glimpse:** decoded first receptor, transport, reaction, effector, fate,
  hgt_propensity, and mutation_rate from the first dictionary entry

On failure, exits non-zero with an error message describing the mismatch.

### Example

```bash
$ python scripts/inspect_rulesets.py output/run_128x128x64/tick_1000.rulesets.bin.zst
Inspecting: output/run_128x128x64/tick_1000.rulesets.bin.zst
  Raw size: 45,818 bytes

  ── Header ──
  Magic:    MRSF
  Version:  1
  Flags:    0
  Dict entries: 82
  Cells:        527
  Ruleset byte size: 536
  Compression: 82 unique / 527 cells = 15.6% unique

  ── Layout ──
  Header:            24 B
  Dictionary:        43,952 B (82 × 536)
  Per-cell refs:     5,270 B (527 × 10)
  Total:             49,246 B

  ── Dict Usage ──
  Dict usage: 82 unique entries across 527 cell refs
  Most-used entry: dict[0] — 48 cells
  Singletons (used once): 8/82
  Distribution (top 10):
    dict[   0]:     48  ########################################
    dict[   1]:     40  ###################################
    ...

  ── Ruleset Glimpse ──
  Decoded glimpse (dict entry 0, 536 B):
    receptor[0]: {'k_half': 1.5, 'n_hill': 2.3, 'gain': 0.8}
    transport[0]: {'uptake_rate': 3.0, 'secrete_rate': 4.0, 'ext_species': 0, 'int_species': 0}
    reaction[0]: {'substrate': 0, 'product': 1, 'catalyst': 2, 'cofactor': 255, 'k_m': 5.0, 'v_max': 6.0, 'k_cat': 7.0}
    effector[0]: {'threshold': 8.0, 'rate': 9.0, 'int_species': 0, 'ext_species': 0}
    fate: {'division_energy': 10.0, 'death_energy': 11.0, 'quiescence_energy': 12.0, 'division_prep_ticks': 13.0}
    hgt_propensity: 14.0
    mutation_rate: 15.0

  ✓ Structure valid

# Error example — bad file
$ python scripts/inspect_rulesets.py /tmp/garbage.bin
Error: magic mismatch: got b'\x00\x01\x02\x03', expected b'MRSF'
```

### Typical Use

Quickly inspect a ruleset dump to understand diversity and verify format
integrity during development or after an engine run:

```bash
# Inspect a specific tick's ruleset dump
python scripts/inspect_rulesets.py output/run_128x128x64/tick_5000.rulesets.bin.zst

# Headers-only summary (no ruleset glimpse, no cell refs)
python scripts/inspect_rulesets.py --no-glimpse --no-cell-sample output/run_128x128x64/tick_5000.rulesets.bin.zst
```
