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
import struct
import subprocess
import sys
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------------
# Format constants (mirroring marl-format/src/lib.rs)
# ---------------------------------------------------------------------------

MAGIC = b"MRSF"
HEADER_SIZE = 24
CELL_REF_STRIDE = 10
CANONICAL_RULESET_SIZE = 536
FORMAT_VERSION = 1

# Ruleset payload layout (byte offsets and sizes of each section)
RECEPTORS_OFF = 0
RECEPTORS_COUNT = 8
RECEPTOR_SIZE = 12  # k_half:f32, n_hill:f32, gain:f32

TRANSPORT_OFF = RECEPTORS_OFF + RECEPTORS_COUNT * RECEPTOR_SIZE  # 96
TRANSPORT_COUNT = 8
TRANSPORT_SIZE = 10  # uptake_rate:f32, secrete_rate:f32, ext_species:u8, int_species:u8

REACTIONS_OFF = TRANSPORT_OFF + TRANSPORT_COUNT * TRANSPORT_SIZE  # 176
REACTIONS_COUNT = 16
REACTION_SIZE = 16  # substrate:u8, product:u8, catalyst:u8, cofactor:u8, k_m:f32, v_max:f32, k_cat:f32

EFFECTORS_OFF = REACTIONS_OFF + REACTIONS_COUNT * REACTION_SIZE  # 432
EFFECTORS_COUNT = 8
EFFECTOR_SIZE = 10  # threshold:f32, rate:f32, int_species:u8, ext_species:u8

FATE_OFF = EFFECTORS_OFF + EFFECTORS_COUNT * EFFECTOR_SIZE  # 512
FATE_SIZE = 16  # division_energy:f32, death_energy:f32, quiescence_energy:f32, division_prep_ticks:f32

HGT_OFF = FATE_OFF + FATE_SIZE  # 528
HGT_SIZE = 4  # hgt_propensity:f32

MUT_OFF = HGT_OFF + HGT_SIZE  # 532
MUT_SIZE = 4  # mutation_rate:f32

assert MUT_OFF + MUT_SIZE == CANONICAL_RULESET_SIZE, (
    f"Payload layout mismatch: expected {CANONICAL_RULESET_SIZE}, "
    f"computed {MUT_OFF + MUT_SIZE}"
)

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


def validate_header(payload: bytes) -> tuple[int, int, int]:
    """Parse and validate header. Returns (flags, dict_count, cell_count)."""
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

    # Validate version
    if version != FORMAT_VERSION:
        print(f"WARNING: version mismatch: got {version}, expected {FORMAT_VERSION}")

    # Validate ruleset byte size
    if ruleset_byte_size != CANONICAL_RULESET_SIZE:
        sys.exit(
            f"Error: ruleset_byte_size mismatch: got {ruleset_byte_size}, "
            f"expected {CANONICAL_RULESET_SIZE}"
        )

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

    return flags, dict_count, cell_count


def decode_receptor(payload: bytes, base: int, i: int) -> dict:
    """Decode one receptor record from a dictionary entry."""
    off = base + RECEPTORS_OFF + i * RECEPTOR_SIZE
    return {
        "k_half": struct.unpack_from("<f", payload, off)[0],
        "n_hill": struct.unpack_from("<f", payload, off + 4)[0],
        "gain": struct.unpack_from("<f", payload, off + 8)[0],
    }


def decode_transport(payload: bytes, base: int, i: int) -> dict:
    """Decode one transport record from a dictionary entry."""
    off = base + TRANSPORT_OFF + i * TRANSPORT_SIZE
    return {
        "uptake_rate": struct.unpack_from("<f", payload, off)[0],
        "secrete_rate": struct.unpack_from("<f", payload, off + 4)[0],
        "ext_species": payload[off + 8],
        "int_species": payload[off + 9],
    }


def decode_reaction(payload: bytes, base: int, i: int) -> dict:
    """Decode one reaction record from a dictionary entry."""
    off = base + REACTIONS_OFF + i * REACTION_SIZE
    return {
        "substrate": payload[off],
        "product": payload[off + 1],
        "catalyst": payload[off + 2],
        "cofactor": payload[off + 3],
        "k_m": struct.unpack_from("<f", payload, off + 4)[0],
        "v_max": struct.unpack_from("<f", payload, off + 8)[0],
        "k_cat": struct.unpack_from("<f", payload, off + 12)[0],
    }


def decode_effector(payload: bytes, base: int, i: int) -> dict:
    """Decode one effector record from a dictionary entry."""
    off = base + EFFECTORS_OFF + i * EFFECTOR_SIZE
    return {
        "threshold": struct.unpack_from("<f", payload, off)[0],
        "rate": struct.unpack_from("<f", payload, off + 4)[0],
        "int_species": payload[off + 8],
        "ext_species": payload[off + 9],
    }


def decode_fate(payload: bytes, base: int) -> dict:
    """Decode fate parameters from a dictionary entry."""
    off = base + FATE_OFF
    return {
        "division_energy": struct.unpack_from("<f", payload, off)[0],
        "death_energy": struct.unpack_from("<f", payload, off + 4)[0],
        "quiescence_energy": struct.unpack_from("<f", payload, off + 8)[0],
        "division_prep_ticks": struct.unpack_from("<f", payload, off + 12)[0],
    }


def decode_ruleset_glimpse(
    payload: bytes, dict_start: int, entry_name: str = "entry 0"
) -> dict:
    """Decode a readable glimpse of one ruleset dictionary entry."""
    r = {}
    # First receptor
    r["receptor[0]"] = decode_receptor(payload, dict_start, 0)
    # First transport
    r["transport[0]"] = decode_transport(payload, dict_start, 0)
    # First reaction
    r["reaction[0]"] = decode_reaction(payload, dict_start, 0)
    # First effector
    r["effector[0]"] = decode_effector(payload, dict_start, 0)
    # Fate
    r["fate"] = decode_fate(payload, dict_start)
    # hgt, mutation
    r["hgt_propensity"] = struct.unpack_from("<f", payload, dict_start + HGT_OFF)[0]
    r["mutation_rate"] = struct.unpack_from("<f", payload, dict_start + MUT_OFF)[0]
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


def print_ruleset_glimpse(payload: bytes, dict_start: int, ruleset_size: int) -> None:
    """Print a decoded glimpse of the first dictionary entry."""
    glimpse = decode_ruleset_glimpse(payload, dict_start)

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
  python scripts/inspect_rulesets.py output/run_128x128x64/tick_1000.rulesets.bin.zst
  python scripts/inspect_rulesets.py /tmp/opencode/marl_output_ruleset_full_test/tick_10.rulesets.bin.zst
  python scripts/inspect_rulesets.py output/run_128x128x64/tick_500.rulesets.bin
        """,
    )
    parser.add_argument(
        "file",
        type=Path,
        help="Path to a tick_<T>.rulesets.bin or .bin.zst file",
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
    args = parser.parse_args()

    path = args.file
    if not path.exists():
        sys.exit(f"Error: file not found: {path}")

    print(f"Inspecting: {path}")
    payload = read_maybe_compressed(path)
    print(f"  Raw size: {len(payload):,} bytes")

    # --- Header ---
    flags, dict_count, cell_count = validate_header(payload)
    ruleset_size = struct.unpack_from("<I", payload, 20)[0]

    print()
    print("  ── Header ──")
    print(f"  Magic:    {MAGIC.decode('ascii')}")
    print(f"  Version:  {struct.unpack_from('<I', payload, 4)[0]}")
    print(f"  Flags:    {flags}")
    print(f"  Dict entries: {dict_count:,}")
    print(f"  Cells:        {cell_count:,}")
    print(f"  Ruleset byte size: {ruleset_size}")
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
        print_ruleset_glimpse(payload, HEADER_SIZE, ruleset_size)

    print()
    print("  ✓ Structure valid")


if __name__ == "__main__":
    main()
