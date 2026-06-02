#!/usr/bin/env python3
"""Sanity-check MARL binary viewer output files."""

from __future__ import annotations

import argparse
import json
import math
import struct
import subprocess
from pathlib import Path


def resolve_pattern(meta: dict, key: str, tick: int, default: str) -> Path:
    pattern = str(meta.get(key, default))
    if "<T>" not in pattern:
        raise SystemExit(f"{key} must contain <T> tick placeholder: {pattern!r}")
    return Path(pattern.replace("<T>", str(tick)))


def read_maybe_compressed(path: Path) -> bytes:
    if path.suffix == ".zst":
        try:
            return subprocess.check_output(["zstd", "-q", "-d", "-c", str(path)])
        except FileNotFoundError as exc:
            raise SystemExit(
                "compressed snapshots require the `zstd` CLI on PATH"
            ) from exc
        except subprocess.CalledProcessError as exc:
            raise SystemExit(f"failed to decompress {path}: {exc}") from exc
    return path.read_bytes()


def check_full_ruleset(run_dir: Path, meta: dict, tick: int) -> None:
    """Validate the full deduplicated ruleset sidecar file."""
    ruleset_path = run_dir / resolve_pattern(
        meta,
        "ruleset_full_file_pattern",
        tick,
        "tick_<T>.rulesets.bin",
    )
    if not ruleset_path.exists():
        raise SystemExit(f"full ruleset file not found: {ruleset_path}")

    payload = read_maybe_compressed(ruleset_path)
    if len(payload) < 24:
        raise SystemExit(
            f"full ruleset file too small for header: {len(payload)} bytes"
        )

    # Parse header
    magic = payload[0:4]
    version = struct.unpack_from("<I", payload, 4)[0]
    flags = struct.unpack_from("<I", payload, 8)[0]
    dict_count = struct.unpack_from("<I", payload, 12)[0]
    cell_count = struct.unpack_from("<I", payload, 16)[0]
    ruleset_byte_size = struct.unpack_from("<I", payload, 20)[0]

    expected_magic = meta.get("ruleset_full_magic_ascii", "MRSF").encode("ascii")
    if magic != expected_magic:
        raise SystemExit(
            f"full ruleset magic mismatch: got {magic!r}, expected {expected_magic!r}"
        )

    expected_version = int(meta.get("ruleset_full_format_version", 1))
    if version != expected_version:
        raise SystemExit(
            f"full ruleset version mismatch: got {version}, expected {expected_version}"
        )
    if flags != 0:
        raise SystemExit(f"full ruleset flags mismatch: got {flags}, expected 0")

    expected_ruleset_size = int(meta.get("ruleset_full_ruleset_byte_size", 536))
    if ruleset_byte_size != expected_ruleset_size:
        raise SystemExit(
            f"full ruleset byte size mismatch: got {ruleset_byte_size}, expected {expected_ruleset_size}"
        )

    expected_cell_ref_stride = int(meta.get("ruleset_full_cell_ref_stride", 10))
    header_size = int(meta.get("ruleset_full_header_size", 24))
    expected_total = (
        header_size
        + dict_count * ruleset_byte_size
        + cell_count * expected_cell_ref_stride
    )
    if len(payload) != expected_total:
        raise SystemExit(
            f"full ruleset file size mismatch: got {len(payload)}, expected {expected_total} "
            f"(dict_count={dict_count}, cell_count={cell_count}, ruleset_byte_size={ruleset_byte_size})"
        )

    # Validate dict IDs in cell refs are in range
    cell_refs_off = header_size + dict_count * ruleset_byte_size
    for i in range(cell_count):
        off = cell_refs_off + i * expected_cell_ref_stride
        x = struct.unpack_from("<H", payload, off)[0]
        y = struct.unpack_from("<H", payload, off + 2)[0]
        z = struct.unpack_from("<H", payload, off + 4)[0]
        dict_id = struct.unpack_from("<I", payload, off + 6)[0]

        grid_x = int(meta["grid_x"])
        grid_y = int(meta["grid_y"])
        grid_z = int(meta["grid_z"])
        if x >= grid_x or y >= grid_y or z >= grid_z:
            raise SystemExit(
                f"cell ref {i}: position ({x},{y},{z}) out of bounds ({grid_x},{grid_y},{grid_z})"
            )
        if dict_id >= dict_count:
            raise SystemExit(
                f"cell ref {i}: dict_id {dict_id} >= dict_count {dict_count}"
            )

    print(
        f"  rulesets: OK  dict({dict_count}) cells({cell_count}) each({ruleset_byte_size}B)"
    )


def require_meta_value(meta: dict, key: str, expected: object) -> None:
    actual = meta.get(key)
    if actual != expected:
        raise SystemExit(f"metadata {key} mismatch: got {actual!r}, expected {expected!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="Output run directory")
    parser.add_argument("tick", type=int, help="Tick number to inspect")
    parser.add_argument(
        "--require-rulesets",
        action="store_true",
        help="Require and validate per-layer ruleset average sidecars",
    )
    parser.add_argument(
        "--require-full-rulesets",
        action="store_true",
        help="Require and validate full deduplicated ruleset sidecars",
    )
    args = parser.parse_args()

    meta_path = args.run_dir / "run_meta.json"

    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    require_meta_value(meta, "endianness", "little")
    require_meta_value(meta, "field_dtype", "f32")
    require_meta_value(meta, "field_layout", "z_y_x_species")
    field_path = args.run_dir / resolve_pattern(
        meta, "field_file_pattern", args.tick, "tick_<T>.field.bin"
    )

    if not bool(meta.get("write_binary_field", True)):
        raise SystemExit("write_binary_field is false; no field snapshot is expected")
    expected_field_bytes = int(meta["field_byte_len"])
    field_payload = read_maybe_compressed(field_path)
    field_bytes = len(field_payload)
    if field_bytes != expected_field_bytes:
        raise SystemExit(
            f"field size mismatch: got {field_bytes}, expected {expected_field_bytes}"
        )
    if field_bytes < 4:
        raise SystemExit(f"field file too small to contain one f32: {field_bytes} bytes")

    first = struct.unpack("<f", field_payload[:4])[0]
    if not math.isfinite(first):
        raise SystemExit(f"first field value is not finite: {first!r}")

    stride = int(meta["cell_record_stride"])
    if stride != 25:
        raise SystemExit(f"cell_record_stride mismatch: got {stride}, expected 25")
    if bool(meta.get("write_binary_cells", True)):
        cells_path = args.run_dir / resolve_pattern(
            meta, "cell_file_pattern", args.tick, "tick_<T>.cells.bin"
        )
        cell_payload = read_maybe_compressed(cells_path)
        cell_bytes = len(cell_payload)
        if cell_bytes % stride != 0:
            raise SystemExit(
                f"cell file size {cell_bytes} is not divisible by stride {stride}"
            )
    else:
        cell_bytes = 0

    ruleset_mode = str(meta.get("ruleset_output_mode", "off"))
    if args.require_rulesets or ruleset_mode in ("layer_averages", "both"):
        ruleset_path = args.run_dir / resolve_pattern(
            meta,
            "ruleset_layer_file_pattern",
            args.tick,
            "tick_<T>.ruleset_layers.bin",
        )
        ruleset_stride = int(meta.get("ruleset_layer_record_stride", 0))
        if ruleset_stride <= 0:
            raise SystemExit("ruleset_layer_record_stride missing or invalid")
        ruleset_payload = read_maybe_compressed(ruleset_path)
        if len(ruleset_payload) % ruleset_stride != 0:
            raise SystemExit(
                f"ruleset file size {len(ruleset_payload)} is not divisible by stride {ruleset_stride}"
            )
        grid_z = int(meta["grid_z"])
        layer_records = len(ruleset_payload) // ruleset_stride
        if layer_records != grid_z:
            raise SystemExit(
                f"ruleset layer record count mismatch: got {layer_records}, expected {grid_z}"
            )

    if args.require_full_rulesets or ruleset_mode in ("full", "both"):
        check_full_ruleset(args.run_dir, meta, args.tick)

    print(
        f"ok: first_f32={first:.6g}, field_bytes={field_bytes}, "
        f"cell_count={cell_bytes // stride}"
    )


if __name__ == "__main__":
    main()
