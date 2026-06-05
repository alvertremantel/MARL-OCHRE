#!/usr/bin/env python3
"""Sanity-check MARL binary viewer output files."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import struct
import subprocess
from pathlib import Path


POSITION_INTEGER_TOLERANCE = 0.001
CURRENT_RULESET_SIZE = 576
CURRENT_FORMAT_VERSION = 2
LEGACY_RULESET_SIZE_V1 = 536
LEGACY_FORMAT_VERSION_V1 = 1


@dataclass(frozen=True)
class RulesetLayout:
    version: int
    ruleset_size: int
    float_sections: tuple[tuple[int, int, int, tuple[int, ...], str], ...]


LAYOUT_V1 = RulesetLayout(
    version=LEGACY_FORMAT_VERSION_V1,
    ruleset_size=LEGACY_RULESET_SIZE_V1,
    float_sections=(
        (0, 8, 12, (0, 4, 8), "receptor"),
        (96, 8, 10, (0, 4), "transport"),
        (176, 16, 16, (4, 8, 12), "reaction"),
        (432, 8, 10, (0, 4), "effector"),
        (512, 1, 16, (0, 4, 8, 12), "fate"),
        (528, 1, 4, (0,), "hgt_propensity"),
        (532, 1, 4, (0,), "mutation_rate"),
    ),
)

LAYOUT_V2 = RulesetLayout(
    version=CURRENT_FORMAT_VERSION,
    ruleset_size=CURRENT_RULESET_SIZE,
    float_sections=(
        (0, 8, 12, (0, 4, 8), "receptor"),
        (96, 8, 15, (0, 4, 11), "transport"),
        (216, 16, 16, (4, 8, 12), "reaction"),
        (472, 8, 10, (0, 4), "effector"),
        (552, 1, 16, (0, 4, 8, 12), "fate"),
        (568, 1, 4, (0,), "hgt_propensity"),
        (572, 1, 4, (0,), "mutation_rate"),
    ),
)

SUPPORTED_RULESET_LAYOUTS = {
    (LAYOUT_V1.version, LAYOUT_V1.ruleset_size): LAYOUT_V1,
    (LAYOUT_V2.version, LAYOUT_V2.ruleset_size): LAYOUT_V2,
}


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


def cadence_due(tick: int, max_ticks: int, interval: int) -> bool:
    return tick + 1 == max_ticks or (interval > 0 and tick % interval == 0)


def ruleset_sidecar_due(meta: dict, tick: int) -> bool:
    return cadence_due(tick, int(meta["max_ticks"]), int(meta["ruleset_interval"]))


def validate_field_floats(payload: bytes) -> None:
    if len(payload) % 4 != 0:
        raise SystemExit(
            f"field file size {len(payload)} is not divisible by f32 width 4"
        )
    for i, (value,) in enumerate(struct.iter_unpack("<f", payload)):
        if not math.isfinite(value):
            raise SystemExit(f"field value {i} is not finite: {value!r}")


def validate_cell_payload(payload: bytes, stride: int, meta: dict) -> None:
    if len(payload) % stride != 0:
        raise SystemExit(
            f"cell file size {len(payload)} is not divisible by stride {stride}"
        )

    bounds = (
        ("x", int(meta["grid_x"])),
        ("y", int(meta["grid_y"])),
        ("z", int(meta["grid_z"])),
    )
    for index in range(len(payload) // stride):
        off = index * stride
        positions = struct.unpack_from("<fff", payload, off)
        energy = struct.unpack_from("<f", payload, off + 21)[0]

        for axis, value in enumerate(positions):
            name, limit = bounds[axis]
            if not math.isfinite(value) or value < 0.0:
                raise SystemExit(
                    f"cell {index}: position {name}={value!r} is not a finite non-negative float"
                )
            rounded = round(value)
            if abs(value - rounded) > POSITION_INTEGER_TOLERANCE:
                raise SystemExit(
                    f"cell {index}: position {name}={value!r} is not close to an integer voxel index"
                )
            if rounded >= limit:
                raise SystemExit(
                    f"cell {index}: position {name}={value!r} out of bounds (grid_{name}={limit})"
                )

        if not math.isfinite(energy):
            raise SystemExit(f"cell {index}: energy={energy!r} is not finite")


def validate_ruleset_layer_payload(payload: bytes, stride: int, meta: dict) -> None:
    if len(payload) % stride != 0:
        raise SystemExit(
            f"ruleset file size {len(payload)} is not divisible by stride {stride}"
        )
    if stride < 8 or (stride - 8) % 4 != 0:
        raise SystemExit(f"ruleset_layer_record_stride has invalid layout: {stride}")

    grid_x = int(meta["grid_x"])
    grid_y = int(meta["grid_y"])
    grid_z = int(meta["grid_z"])
    max_layer_cells = grid_x * grid_y
    layer_records = len(payload) // stride
    if layer_records != grid_z:
        raise SystemExit(
            f"ruleset layer record count mismatch: got {layer_records}, expected {grid_z}"
        )

    avg_count = (stride - 8) // 4
    for index in range(layer_records):
        off = index * stride
        z = struct.unpack_from("<H", payload, off)[0]
        reserved = struct.unpack_from("<H", payload, off + 2)[0]
        cell_count = struct.unpack_from("<I", payload, off + 4)[0]

        if z != index:
            raise SystemExit(
                f"ruleset layer {index}: z index mismatch: got {z}, expected {index}"
            )
        if z >= grid_z:
            raise SystemExit(f"ruleset layer {index}: z index {z} >= grid_z {grid_z}")
        if reserved != 0:
            raise SystemExit(
                f"ruleset layer {index}: reserved field is {reserved}, expected 0"
            )
        if cell_count > max_layer_cells:
            raise SystemExit(
                f"ruleset layer {index}: cell_count {cell_count} exceeds layer capacity {max_layer_cells}"
            )

        for avg_index in range(avg_count):
            value = struct.unpack_from("<f", payload, off + 8 + avg_index * 4)[0]
            if not math.isfinite(value):
                raise SystemExit(
                    f"ruleset layer {index}: average {avg_index} is not finite: {value!r}"
                )


def validate_ruleset_dictionary_floats(
    payload: bytes, dict_off: int, dict_count: int, layout: RulesetLayout
) -> None:
    for dict_id in range(dict_count):
        base = dict_off + dict_id * layout.ruleset_size
        for section_off, count, stride, float_offsets, section_name in layout.float_sections:
            for entry_index in range(count):
                entry_base = base + section_off + entry_index * stride
                for rel_off in float_offsets:
                    value = struct.unpack_from("<f", payload, entry_base + rel_off)[0]
                    if not math.isfinite(value):
                        raise SystemExit(
                            f"full ruleset dict[{dict_id}] {section_name}[{entry_index}] "
                            f"float@+{rel_off} is not finite: {value!r}"
                        )


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

    expected_version = int(meta.get("ruleset_full_format_version", version))
    if version != expected_version:
        raise SystemExit(
            f"full ruleset version mismatch: got {version}, expected {expected_version}"
        )
    if flags != 0:
        raise SystemExit(f"full ruleset flags mismatch: got {flags}, expected 0")

    expected_ruleset_size = int(meta.get("ruleset_full_ruleset_byte_size", ruleset_byte_size))
    if ruleset_byte_size != expected_ruleset_size:
        raise SystemExit(
            f"full ruleset byte size mismatch: got {ruleset_byte_size}, expected {expected_ruleset_size}"
        )
    layout = SUPPORTED_RULESET_LAYOUTS.get((version, ruleset_byte_size))
    if layout is None:
        raise SystemExit(
            f"unsupported full ruleset layout: version {version}, size {ruleset_byte_size}; "
            f"supported layouts are v1/{LEGACY_RULESET_SIZE_V1}B and v2/{CURRENT_RULESET_SIZE}B"
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

    dict_off = header_size
    validate_ruleset_dictionary_floats(payload, dict_off, dict_count, layout)

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
    validate_field_floats(field_payload)

    stride = int(meta["cell_record_stride"])
    if stride != 25:
        raise SystemExit(f"cell_record_stride mismatch: got {stride}, expected 25")
    if bool(meta.get("write_binary_cells", True)):
        cells_path = args.run_dir / resolve_pattern(
            meta, "cell_file_pattern", args.tick, "tick_<T>.cells.bin"
        )
        cell_payload = read_maybe_compressed(cells_path)
        cell_bytes = len(cell_payload)
        validate_cell_payload(cell_payload, stride, meta)
    else:
        cell_bytes = 0

    ruleset_mode = str(meta.get("ruleset_output_mode", "off"))
    auto_ruleset_due = (
        ruleset_mode in ("layer_averages", "both")
        and ruleset_sidecar_due(meta, args.tick)
    )
    if args.require_rulesets or auto_ruleset_due:
        ruleset_path = args.run_dir / resolve_pattern(
            meta,
            "ruleset_layer_file_pattern",
            args.tick,
            "tick_<T>.ruleset_layers.bin",
        )
        if not ruleset_path.exists():
            raise SystemExit(f"ruleset layer file not found: {ruleset_path}")
        ruleset_stride = int(meta.get("ruleset_layer_record_stride", 0))
        if ruleset_stride <= 0:
            raise SystemExit("ruleset_layer_record_stride missing or invalid")
        ruleset_payload = read_maybe_compressed(ruleset_path)
        validate_ruleset_layer_payload(ruleset_payload, ruleset_stride, meta)

    auto_full_ruleset_due = ruleset_mode in ("full", "both") and ruleset_sidecar_due(
        meta, args.tick
    )
    if args.require_full_rulesets or auto_full_ruleset_due:
        check_full_ruleset(args.run_dir, meta, args.tick)

    print(
        f"ok: first_f32={first:.6g}, field_bytes={field_bytes}, "
        f"cell_count={cell_bytes // stride}"
    )


if __name__ == "__main__":
    main()
