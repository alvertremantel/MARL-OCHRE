#!/usr/bin/env python3
"""Inspect MARL full deduplicated ruleset dump files (tick_<T>.rulesets.bin[.zst]).

Prints a human-readable summary of the on-disk format:
  - header metadata (magic, version, dict/cell count, ruleset byte size)
  - dictionary usage frequency histogram
  - a decoded glimpse of at least one dictionary entry
  - basic structural integrity validation

Works on both uncompressed .bin and zstd-compressed .bin.zst files.
If the file is compressed, requires the `zstd` CLI on PATH.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import struct
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Format constants (mirroring marl-format/src/lib.rs)
# ---------------------------------------------------------------------------

MAGIC = b"MRSF"
HEADER_SIZE = 24
CELL_REF_STRIDE = 10
CURRENT_RULESET_SIZE = 576
CURRENT_FORMAT_VERSION = 2
LEGACY_RULESET_SIZE_V1 = 536
LEGACY_FORMAT_VERSION_V1 = 1


@dataclass(frozen=True)
class RulesetLayout:
    version: int
    ruleset_size: int
    receptor_off: int
    receptor_count: int
    receptor_size: int
    transport_off: int
    transport_count: int
    transport_size: int
    transport_has_gate: bool
    reaction_off: int
    reaction_count: int
    reaction_size: int
    effector_off: int
    effector_count: int
    effector_size: int
    fate_off: int
    fate_size: int
    hgt_off: int
    mut_off: int


LAYOUT_V1 = RulesetLayout(
    version=LEGACY_FORMAT_VERSION_V1,
    ruleset_size=LEGACY_RULESET_SIZE_V1,
    receptor_off=0,
    receptor_count=8,
    receptor_size=12,
    transport_off=96,
    transport_count=8,
    transport_size=10,
    transport_has_gate=False,
    reaction_off=176,
    reaction_count=16,
    reaction_size=16,
    effector_off=432,
    effector_count=8,
    effector_size=10,
    fate_off=512,
    fate_size=16,
    hgt_off=528,
    mut_off=532,
)

LAYOUT_V2 = RulesetLayout(
    version=CURRENT_FORMAT_VERSION,
    ruleset_size=CURRENT_RULESET_SIZE,
    receptor_off=0,
    receptor_count=8,
    receptor_size=12,
    transport_off=96,
    transport_count=8,
    transport_size=15,
    transport_has_gate=True,
    reaction_off=216,
    reaction_count=16,
    reaction_size=16,
    effector_off=472,
    effector_count=8,
    effector_size=10,
    fate_off=552,
    fate_size=16,
    hgt_off=568,
    mut_off=572,
)

SUPPORTED_LAYOUTS = {
    (LAYOUT_V1.version, LAYOUT_V1.ruleset_size): LAYOUT_V1,
    (LAYOUT_V2.version, LAYOUT_V2.ruleset_size): LAYOUT_V2,
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def read_maybe_compressed(path: Path) -> bytes:
    """Read a file, decompressing with zstd if the extension is .zst."""
    if path.suffix == ".zst":
        try:
            return subprocess.check_output(["zstd", "-q", "-d", "-c", str(path)])
        except FileNotFoundError:
            sys.exit("Error: compressed file requires the `zstd` CLI on PATH")
        except subprocess.CalledProcessError as exc:
            sys.exit(f"Error: failed to decompress {path}: {exc}")
    return path.read_bytes()


def parse_tick(raw_tick: str) -> int:
    """Parse a CLI tick argument with a clearer error than argparse's int."""
    try:
        tick = int(raw_tick, 10)
    except ValueError:
        sys.exit(f"Error: malformed tick {raw_tick!r}; expected a non-negative integer")
    if tick < 0:
        sys.exit(f"Error: malformed tick {raw_tick!r}; expected a non-negative integer")
    return tick


def load_run_meta(run_dir: Path) -> dict:
    """Load run_meta.json from an engine output directory."""
    meta_path = run_dir / "run_meta.json"
    if not meta_path.exists():
        sys.exit(f"Error: run metadata not found: {meta_path}")
    try:
        return json.loads(meta_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        sys.exit(f"Error: invalid run metadata {meta_path}: {exc}")


def resolve_pattern(meta: dict, key: str, tick: int) -> Path:
    """Resolve a run_meta.json snapshot pattern by substituting the tick."""
    pattern = meta.get(key)
    if not isinstance(pattern, str) or not pattern:
        sys.exit(
            "Error: full ruleset output is off or missing in run_meta.json "
            f"({key} is absent)"
        )
    if "<T>" not in pattern:
        sys.exit(f"Error: {key} must contain <T> tick placeholder: {pattern!r}")
    return Path(pattern.replace("<T>", str(tick)))


def resolve_input_path(path_or_run_dir: Path, raw_tick: Optional[str]) -> Path:
    """Resolve either direct-file usage or RUN_DIR TICK usage."""
    if raw_tick is None:
        if path_or_run_dir.is_dir():
            sys.exit(
                "Error: direct-file usage requires a ruleset file path; "
                "for a run directory use: scripts/inspect_rulesets.py RUN_DIR TICK"
            )
        return path_or_run_dir

    tick = parse_tick(raw_tick)
    run_dir = path_or_run_dir
    if not run_dir.exists():
        sys.exit(f"Error: run directory not found: {run_dir}")
    if not run_dir.is_dir():
        sys.exit(f"Error: RUN_DIR is not a directory: {run_dir}")

    meta = load_run_meta(run_dir)
    return run_dir / resolve_pattern(meta, "ruleset_full_file_pattern", tick)


def validate_header(
    payload: bytes, allow_unsupported_version: bool
) -> tuple[int, int, int, RulesetLayout]:
    """Parse and validate header. Returns (flags, dict_count, cell_count, layout)."""
    if len(payload) < HEADER_SIZE:
        sys.exit(
            f"Error: file too small for header ({len(payload)} < {HEADER_SIZE} bytes)"
        )

    magic = payload[0:4]
    version = struct.unpack_from("<I", payload, 4)[0]
    flags = struct.unpack_from("<I", payload, 8)[0]
    dict_count = struct.unpack_from("<I", payload, 12)[0]
    cell_count = struct.unpack_from("<I", payload, 16)[0]
    ruleset_byte_size = struct.unpack_from("<I", payload, 20)[0]

    # Validate magic
    if magic != MAGIC:
        sys.exit(f"Error: magic mismatch: got {magic!r}, expected {MAGIC!r}")

    layout = SUPPORTED_LAYOUTS.get((version, ruleset_byte_size))
    if layout is None:
        message = (
            f"unsupported ruleset layout: version {version}, "
            f"ruleset_byte_size {ruleset_byte_size}; supported layouts are "
            f"v1/{LEGACY_RULESET_SIZE_V1}B and v2/{CURRENT_RULESET_SIZE}B"
        )
        if not allow_unsupported_version:
            sys.exit(f"Error: {message}")
        print(f"WARNING: {message}")
        layout = LAYOUT_V2

    # Validate total file size
    cell_refs_off = HEADER_SIZE + dict_count * ruleset_byte_size
    expected_total = cell_refs_off + cell_count * CELL_REF_STRIDE
    if len(payload) != expected_total:
        sys.exit(
            f"Error: file size mismatch: got {len(payload)} bytes, "
            f"expected {expected_total} "
            f"(header={HEADER_SIZE} + dict={dict_count}×{ruleset_byte_size} "
            f"+ refs={cell_count}×{CELL_REF_STRIDE})"
        )

    return flags, dict_count, cell_count, layout


def decode_receptor(payload: bytes, base: int, i: int, layout: RulesetLayout) -> dict:
    """Decode one receptor record from a dictionary entry."""
    off = base + layout.receptor_off + i * layout.receptor_size
    return {
        "k_half": struct.unpack_from("<f", payload, off)[0],
        "n_hill": struct.unpack_from("<f", payload, off + 4)[0],
        "gain": struct.unpack_from("<f", payload, off + 8)[0],
    }


def decode_transport(payload: bytes, base: int, i: int, layout: RulesetLayout) -> dict:
    """Decode one transport record from a dictionary entry."""
    off = base + layout.transport_off + i * layout.transport_size
    decoded = {
        "uptake_rate": struct.unpack_from("<f", payload, off)[0],
        "secrete_rate": struct.unpack_from("<f", payload, off + 4)[0],
        "ext_species": payload[off + 8],
        "int_species": payload[off + 9],
    }
    if layout.transport_has_gate:
        decoded["gate_receptor"] = payload[off + 10]
        decoded["gate_weight"] = struct.unpack_from("<f", payload, off + 11)[0]
    return decoded


def decode_reaction(payload: bytes, base: int, i: int, layout: RulesetLayout) -> dict:
    """Decode one reaction record from a dictionary entry."""
    off = base + layout.reaction_off + i * layout.reaction_size
    return {
        "substrate": payload[off],
        "product": payload[off + 1],
        "catalyst": payload[off + 2],
        "cofactor": payload[off + 3],
        "k_m": struct.unpack_from("<f", payload, off + 4)[0],
        "v_max": struct.unpack_from("<f", payload, off + 8)[0],
        "k_cat": struct.unpack_from("<f", payload, off + 12)[0],
    }


def decode_effector(payload: bytes, base: int, i: int, layout: RulesetLayout) -> dict:
    """Decode one effector record from a dictionary entry."""
    off = base + layout.effector_off + i * layout.effector_size
    return {
        "threshold": struct.unpack_from("<f", payload, off)[0],
        "rate": struct.unpack_from("<f", payload, off + 4)[0],
        "int_species": payload[off + 8],
        "ext_species": payload[off + 9],
    }


def decode_fate(payload: bytes, base: int, layout: RulesetLayout) -> dict:
    """Decode fate parameters from a dictionary entry."""
    off = base + layout.fate_off
    return {
        "division_energy": struct.unpack_from("<f", payload, off)[0],
        "death_energy": struct.unpack_from("<f", payload, off + 4)[0],
        "quiescence_energy": struct.unpack_from("<f", payload, off + 8)[0],
        "division_prep_ticks": struct.unpack_from("<f", payload, off + 12)[0],
    }


def decode_ruleset_glimpse(
    payload: bytes, dict_start: int, layout: RulesetLayout, entry_name: str = "entry 0"
) -> dict:
    """Decode a readable glimpse of one ruleset dictionary entry."""
    r = {}
    # First receptor
    r["receptor[0]"] = decode_receptor(payload, dict_start, 0, layout)
    # First transport
    r["transport[0]"] = decode_transport(payload, dict_start, 0, layout)
    # First reaction
    r["reaction[0]"] = decode_reaction(payload, dict_start, 0, layout)
    # First effector
    r["effector[0]"] = decode_effector(payload, dict_start, 0, layout)
    # Fate
    r["fate"] = decode_fate(payload, dict_start, layout)
    # hgt, mutation
    r["hgt_propensity"] = struct.unpack_from("<f", payload, dict_start + layout.hgt_off)[0]
    r["mutation_rate"] = struct.unpack_from("<f", payload, dict_start + layout.mut_off)[0]
    return r


def print_dict_usage_histogram(dict_refs: list[int]) -> None:
    """Print a small usage-frequency summary of dictionary entries."""
    freq = Counter(dict_refs)
    n = len(dict_refs)

    unique = len(freq)
    # Top 10 entries
    most_common = freq.most_common(10)

    print(f"  Dict usage: {unique} unique entries across {n} cell refs")
    if unique == 0:
        return

    max_refs = most_common[0][1]
    single_use = sum(1 for c in freq.values() if c == 1)

    print(f"  Most-used entry: dict[{most_common[0][0]}] — {max_refs} cells")
    print(f"  Singletons (used once): {single_use}/{unique}")
    if unique > 1:
        print(f"  Distribution (top {len(most_common)}):")
        for did, count in most_common:
            bar = "#" * min(count, 40)
            print(f"    dict[{did:>4}]: {count:>6}  {bar}")

    # Also show dict_ids sorted to check for contiguous 0..N-1
    ids_sorted = sorted(freq.keys())
    if ids_sorted == list(range(unique)):
        print(f"  dict_ids range: contiguous 0..{unique - 1}")
    else:
        print(
            f"  dict_ids range: {ids_sorted[0]}..{ids_sorted[-1]} "
            f"(non-contiguous, {len(ids_sorted)} entries)"
        )


def print_cell_ref_sampling(
    payload: bytes, cell_refs_off: int, cell_count: int, max_show: int = 10
) -> None:
    """Print a sample of per-cell position + dict_id references."""
    show = min(cell_count, max_show)
    print(f"\n  Cell refs (first {show} of {cell_count}):")
    for i in range(show):
        off = cell_refs_off + i * CELL_REF_STRIDE
        x = struct.unpack_from("<H", payload, off)[0]
        y = struct.unpack_from("<H", payload, off + 2)[0]
        z = struct.unpack_from("<H", payload, off + 4)[0]
        dict_id = struct.unpack_from("<I", payload, off + 6)[0]
        print(f"    [{i:>5}] pos=({x:>3},{y:>3},{z:>2})  dict_id={dict_id}")


def print_ruleset_glimpse(
    payload: bytes, dict_start: int, ruleset_size: int, layout: RulesetLayout
) -> None:
    """Print a decoded glimpse of the first dictionary entry."""
    glimpse = decode_ruleset_glimpse(payload, dict_start, layout)

    print(f"\n  Decoded glimpse (dict entry 0, {ruleset_size} B):")
    for key, val in glimpse.items():
        print(f"    {key}: {val}")


def collect_dict_refs(
    payload: bytes, cell_refs_off: int, cell_count: int, dict_count: int
) -> list[int]:
    """Collect all dict_id references and validate range. Returns the list."""
    refs = []
    for i in range(cell_count):
        off = cell_refs_off + i * CELL_REF_STRIDE
        dict_id = struct.unpack_from("<I", payload, off + 6)[0]
        if dict_id >= dict_count:
            sys.exit(
                f"Error: cell ref {i}: dict_id {dict_id} >= dict_count {dict_count}"
            )
        refs.append(dict_id)
    return refs


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect a MARL full ruleset dump file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python scripts/inspect_rulesets.py output/run_128x128x64 1000
  python scripts/inspect_rulesets.py output/run_128x128x64/tick_1000.rulesets.bin.zst
  python scripts/inspect_rulesets.py /tmp/opencode/marl_output_ruleset_full_test/tick_10.rulesets.bin.zst
  python scripts/inspect_rulesets.py output/run_128x128x64/tick_500.rulesets.bin
        """,
    )
    parser.add_argument(
        "path",
        type=Path,
        help="Path to a ruleset file, or a run directory when TICK is provided",
    )
    parser.add_argument(
        "tick",
        nargs="?",
        help="Tick to inspect when PATH is a run directory",
    )
    parser.add_argument(
        "--no-glimpse",
        action="store_true",
        help="Skip the decoded ruleset glimpse (just header + usage stats)",
    )
    parser.add_argument(
        "--no-cell-sample",
        action="store_true",
        help="Skip the per-cell ref sample",
    )
    parser.add_argument(
        "--full-histogram",
        action="store_true",
        help="Show complete dict usage histogram (not just top 10)",
    )
    parser.add_argument(
        "--allow-unsupported-version",
        action="store_true",
        help="Warn and continue when the file format version is unsupported",
    )
    args = parser.parse_args()

    path = resolve_input_path(args.path, args.tick)
    if not path.exists():
        sys.exit(f"Error: file not found: {path}")

    print(f"Inspecting: {path}")
    payload = read_maybe_compressed(path)
    print(f"  Raw size: {len(payload):,} bytes")

    # --- Header ---
    flags, dict_count, cell_count, layout = validate_header(
        payload, args.allow_unsupported_version
    )
    ruleset_size = struct.unpack_from("<I", payload, 20)[0]

    print()
    print("  ── Header ──")
    print(f"  Magic:    {MAGIC.decode('ascii')}")
    print(f"  Version:  {struct.unpack_from('<I', payload, 4)[0]}")
    print(f"  Flags:    {flags}")
    print(f"  Dict entries: {dict_count:,}")
    print(f"  Cells:        {cell_count:,}")
    print(f"  Ruleset byte size: {ruleset_size}")
    print(f"  Decoded layout: v{layout.version}")
    if cell_count > 0 and dict_count > 0:
        compress_ratio = dict_count / cell_count * 100
        print(
            f"  Compression: {dict_count} unique / {cell_count} cells = "
            f"{compress_ratio:.1f}% unique"
        )

    # --- Offsets ---
    cell_refs_off = HEADER_SIZE + dict_count * ruleset_size
    total_dict_bytes = dict_count * ruleset_size
    total_refs_bytes = cell_count * CELL_REF_STRIDE

    print()
    print("  ── Layout ──")
    print(f"  Header:            {HEADER_SIZE} B")
    print(
        f"  Dictionary:        {total_dict_bytes:,} B ({dict_count} × {ruleset_size})"
    )
    print(
        f"  Per-cell refs:     {total_refs_bytes:,} B ({cell_count} × {CELL_REF_STRIDE})"
    )
    print(
        f"  Total:             {HEADER_SIZE + total_dict_bytes + total_refs_bytes:,} B"
    )

    # --- Cell refs ---
    if cell_count > 0:
        dict_refs = collect_dict_refs(payload, cell_refs_off, cell_count, dict_count)
        print()
        print("  ── Dict Usage ──")
        print_dict_usage_histogram(dict_refs)

        if not args.no_cell_sample:
            print_cell_ref_sampling(payload, cell_refs_off, cell_count)

    # --- Ruleset glimpse ---
    if dict_count > 0 and not args.no_glimpse:
        print()
        print("  ── Ruleset Glimpse ──")
        print_ruleset_glimpse(payload, HEADER_SIZE, ruleset_size, layout)

    print()
    print("  ✓ Structure valid")


if __name__ == "__main__":
    main()
