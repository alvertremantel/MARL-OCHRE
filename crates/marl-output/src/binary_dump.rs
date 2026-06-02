//! Binary outputs for high-throughput viewer and analysis ingestion.

use marl_cell::cell::{CellState, Ruleset};
use marl_config::*;
use marl_field::field::Field;
use marl_format::{
    CELL_FILE_PATTERN_RAW, CELL_RECORD_STRIDE, FIELD_FILE_PATTERN_RAW, RULESET_FULL_CANONICAL_SIZE,
    RULESET_FULL_CELL_REF_STRIDE, RULESET_FULL_FILE_PATTERN_RAW, RULESET_FULL_FORMAT_VERSION,
    RULESET_FULL_HEADER_SIZE, RULESET_FULL_MAGIC, RULESET_FULL_PAYLOAD_LAYOUT,
    RULESET_LAYER_FILE_PATTERN_RAW, RULESET_LAYER_RECORD_LAYOUT, RunMeta, ViewerCellRecord,
    field_byte_len, ruleset_layer_record_stride,
};
use serde::Serialize;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Error, ErrorKind, Write};
use std::mem;
use std::path::{Path, PathBuf};

#[cfg(not(target_endian = "little"))]
compile_error!("binary dumps declare little-endian layout and require a little-endian target");

const RULESET_LAYER_AVG_F32_COUNT: usize =
    3 * S_RECEPTORS + 2 * S_TRANSPORTERS + 3 * R_MAX + 2 * S_EFFECTORS + 6;
const RULESET_LAYER_RECORD_STRIDE: usize = 8 + 4 * RULESET_LAYER_AVG_F32_COUNT;

fn as_bytes<T>(slice: &[T]) -> &[u8] {
    let len = mem::size_of_val(slice);
    let ptr = slice.as_ptr().cast::<u8>();
    // SAFETY: `slice` is valid for `len` bytes, and u8 has alignment 1.
    unsafe { std::slice::from_raw_parts(ptr, len) }
}

fn viewer_cell_from(cell: &CellState) -> ViewerCellRecord {
    ViewerCellRecord {
        pos: [cell.pos[0] as f32, cell.pos[1] as f32, cell.pos[2] as f32],
        lineage_id: cell.lineage_id,
        starter_type: cell.starter_type,
        energy: cell.internal[0],
    }
}

fn compressed_pattern(base: &str, compression: BinaryCompression) -> String {
    match compression {
        BinaryCompression::None => base.to_string(),
        BinaryCompression::Zstd => format!("{base}.zst"),
    }
}

fn field_file_pattern(out: &OutputConfig) -> String {
    compressed_pattern(FIELD_FILE_PATTERN_RAW, out.binary_compression)
}

fn cell_file_pattern(out: &OutputConfig) -> String {
    compressed_pattern(CELL_FILE_PATTERN_RAW, out.binary_compression)
}

fn ruleset_layer_file_pattern(out: &OutputConfig) -> String {
    compressed_pattern(RULESET_LAYER_FILE_PATTERN_RAW, out.binary_compression)
}

fn ruleset_full_file_pattern(out: &OutputConfig) -> String {
    compressed_pattern(RULESET_FULL_FILE_PATTERN_RAW, out.binary_compression)
}

fn binary_file_path(out: &OutputConfig, tick: u64, stem: &str) -> PathBuf {
    let suffix = match out.binary_compression {
        BinaryCompression::None => format!("tick_{tick}.{stem}.bin"),
        BinaryCompression::Zstd => format!("tick_{tick}.{stem}.bin.zst"),
    };
    Path::new(&out.output_dir).join(suffix)
}

fn temp_path_for(path: &Path) -> PathBuf {
    let mut temp = path.as_os_str().to_os_string();
    temp.push(".tmp");
    PathBuf::from(temp)
}

fn write_binary_payload(path: &Path, payload: &[u8], out: &OutputConfig) -> std::io::Result<()> {
    let temp_path = temp_path_for(path);
    let result = (|| -> std::io::Result<()> {
        let file = File::create(&temp_path)?;
        match out.binary_compression {
            BinaryCompression::None => {
                let mut writer = BufWriter::new(file);
                writer.write_all(payload)?;
                writer.flush()?;
            }
            BinaryCompression::Zstd => {
                let writer = BufWriter::new(file);
                let mut encoder =
                    zstd::stream::write::Encoder::new(writer, out.binary_compression_level)?;
                encoder.write_all(payload)?;
                let mut writer = encoder.finish()?;
                writer.flush()?;
            }
        }
        commit_temp_path(&temp_path, path)?;
        Ok(())
    })();

    if result.is_err() {
        let _ = fs::remove_file(&temp_path);
    }
    result
}

fn commit_temp_path(temp_path: &Path, path: &Path) -> std::io::Result<()> {
    match fs::rename(temp_path, path) {
        Ok(()) => Ok(()),
        Err(err) if err.kind() == ErrorKind::AlreadyExists => {
            fs::remove_file(path)?;
            fs::rename(temp_path, path)
        }
        Err(err) => Err(err),
    }
}

pub fn write_field_dump(field: &Field, tick: u64, out: &OutputConfig) -> std::io::Result<()> {
    fs::create_dir_all(&out.output_dir)?;
    let path = binary_file_path(out, tick, "field");
    write_binary_payload(&path, as_bytes(&field.data), out)
}

pub fn write_cell_dump(cells: &[CellState], tick: u64, out: &OutputConfig) -> std::io::Result<()> {
    fs::create_dir_all(&out.output_dir)?;
    let viewer_cells: Vec<ViewerCellRecord> = cells.iter().map(viewer_cell_from).collect();
    let path = binary_file_path(out, tick, "cells");
    write_binary_payload(&path, as_bytes(&viewer_cells), out)
}

fn append_ruleset_layer_values(buf: &mut [f64; RULESET_LAYER_AVG_F32_COUNT], ruleset: &Ruleset) {
    let mut i = 0;
    for receptor in &ruleset.receptors {
        buf[i] += receptor.k_half as f64;
        buf[i + 1] += receptor.n_hill as f64;
        buf[i + 2] += receptor.gain as f64;
        i += 3;
    }
    for transport in &ruleset.transport {
        buf[i] += transport.uptake_rate as f64;
        buf[i + 1] += transport.secrete_rate as f64;
        i += 2;
    }
    for reaction in &ruleset.reactions {
        buf[i] += reaction.k_m as f64;
        buf[i + 1] += reaction.v_max as f64;
        buf[i + 2] += reaction.k_cat as f64;
        i += 3;
    }
    for effector in &ruleset.effectors {
        buf[i] += effector.threshold as f64;
        buf[i + 1] += effector.rate as f64;
        i += 2;
    }
    buf[i] += ruleset.fate.division_energy as f64;
    buf[i + 1] += ruleset.fate.death_energy as f64;
    buf[i + 2] += ruleset.fate.quiescence_energy as f64;
    buf[i + 3] += ruleset.fate.division_prep_ticks as f64;
    buf[i + 4] += ruleset.hgt_propensity as f64;
    buf[i + 5] += ruleset.mutation_rate as f64;
}

pub fn write_ruleset_layer_dump(
    cells: &[CellState],
    tick: u64,
    out: &OutputConfig,
) -> std::io::Result<()> {
    fs::create_dir_all(&out.output_dir)?;

    let mut bytes = Vec::with_capacity(GRID_Z * RULESET_LAYER_RECORD_STRIDE);
    for z in 0..GRID_Z {
        let mut sums = [0.0f64; RULESET_LAYER_AVG_F32_COUNT];
        let mut count: u32 = 0;

        for cell in cells.iter().filter(|cell| cell.pos[2] as usize == z) {
            append_ruleset_layer_values(&mut sums, &cell.ruleset);
            count += 1;
        }

        bytes.extend_from_slice(&(z as u16).to_le_bytes());
        bytes.extend_from_slice(&0u16.to_le_bytes());
        bytes.extend_from_slice(&count.to_le_bytes());

        let denom = if count > 0 { count as f64 } else { 1.0 };
        for sum in sums {
            let avg = if count > 0 { (sum / denom) as f32 } else { 0.0 };
            bytes.extend_from_slice(&avg.to_le_bytes());
        }
    }

    let path = binary_file_path(out, tick, "ruleset_layers");
    write_binary_payload(&path, &bytes, out)
}

// ---------------------------------------------------------------------------
// Full per-cell deduplicated ruleset dump
// ---------------------------------------------------------------------------

/// Serialize a `Ruleset` into its canonical fixed-size byte representation.
///
/// All fields are written in a deterministic order so that byte-for-byte
/// comparison produces exact deduplication. The layout matches
/// [`RULESET_FULL_PAYLOAD_LAYOUT`] and the total output is exactly
/// [`RULESET_FULL_CANONICAL_SIZE`] bytes.
fn ruleset_to_canonical_bytes(ruleset: &Ruleset) -> [u8; RULESET_FULL_CANONICAL_SIZE as usize] {
    let mut buf = [0u8; RULESET_FULL_CANONICAL_SIZE as usize];
    let mut off: usize = 0;

    // receptors: 8 × {k_half:f32, n_hill:f32, gain:f32}
    for r in &ruleset.receptors {
        buf[off..off + 4].copy_from_slice(&r.k_half.to_le_bytes());
        buf[off + 4..off + 8].copy_from_slice(&r.n_hill.to_le_bytes());
        buf[off + 8..off + 12].copy_from_slice(&r.gain.to_le_bytes());
        off += 12;
    }

    // transport: 8 × {uptake_rate:f32, secrete_rate:f32, ext_species:u8, int_species:u8}
    for t in &ruleset.transport {
        buf[off..off + 4].copy_from_slice(&t.uptake_rate.to_le_bytes());
        buf[off + 4..off + 8].copy_from_slice(&t.secrete_rate.to_le_bytes());
        buf[off + 8] = t.ext_species;
        buf[off + 9] = t.int_species;
        off += 10;
    }

    // reactions: 16 × {substrate:u8, product:u8, catalyst:u8, cofactor:u8, k_m:f32, v_max:f32, k_cat:f32}
    for rxn in &ruleset.reactions {
        buf[off] = rxn.substrate;
        buf[off + 1] = rxn.product;
        buf[off + 2] = rxn.catalyst;
        buf[off + 3] = rxn.cofactor;
        buf[off + 4..off + 8].copy_from_slice(&rxn.k_m.to_le_bytes());
        buf[off + 8..off + 12].copy_from_slice(&rxn.v_max.to_le_bytes());
        buf[off + 12..off + 16].copy_from_slice(&rxn.k_cat.to_le_bytes());
        off += 16;
    }

    // effectors: 8 × {threshold:f32, rate:f32, int_species:u8, ext_species:u8}
    for e in &ruleset.effectors {
        buf[off..off + 4].copy_from_slice(&e.threshold.to_le_bytes());
        buf[off + 4..off + 8].copy_from_slice(&e.rate.to_le_bytes());
        buf[off + 8] = e.int_species;
        buf[off + 9] = e.ext_species;
        off += 10;
    }

    // fate: {division_energy:f32, death_energy:f32, quiescence_energy:f32, division_prep_ticks:f32}
    buf[off..off + 4].copy_from_slice(&ruleset.fate.division_energy.to_le_bytes());
    buf[off + 4..off + 8].copy_from_slice(&ruleset.fate.death_energy.to_le_bytes());
    buf[off + 8..off + 12].copy_from_slice(&ruleset.fate.quiescence_energy.to_le_bytes());
    buf[off + 12..off + 16].copy_from_slice(&ruleset.fate.division_prep_ticks.to_le_bytes());
    off += 16;

    // hgt_propensity, mutation_rate
    buf[off..off + 4].copy_from_slice(&ruleset.hgt_propensity.to_le_bytes());
    buf[off + 4..off + 8].copy_from_slice(&ruleset.mutation_rate.to_le_bytes());

    buf
}

/// Write a deduplicated full ruleset dump for all cells at a tick.
///
/// Identical rulesets are stored once in a dictionary; each cell references
/// its dictionary entry by ID. The file layout:
///
/// ```text
/// HEADER (24 bytes):
///   magic:        [u8; 4] = b"MRSF"
///   version:      u32 le  (currently 1)
///   flags:        u32 le  (reserved, write 0)
///   dict_count:   u32 le
///   cell_count:   u32 le
///   ruleset_byte_size: u32 le
///
/// DICTIONARY:  dict_count × ruleset_byte_size bytes
///   Sequential canonical ruleset payloads
///
/// PER-CELL REFS: cell_count × 10 bytes each
///   x:       u16 le
///   y:       u16 le
///   z:       u16 le
///   dict_id: u32 le
/// ```
pub fn write_ruleset_full_dump(
    cells: &[CellState],
    tick: u64,
    out: &OutputConfig,
) -> std::io::Result<()> {
    fs::create_dir_all(&out.output_dir)?;

    let ruleset_byte_size = RULESET_FULL_CANONICAL_SIZE;
    let cell_count = u32::try_from(cells.len()).map_err(|_| {
        Error::new(
            ErrorKind::InvalidInput,
            "too many cells to encode in full ruleset dump",
        )
    })?;

    // Build deduplicated dictionary: canonical bytes → dict_id
    let mut dict: HashMap<[u8; RULESET_FULL_CANONICAL_SIZE as usize], u32> = HashMap::new();
    let cell_refs_capacity = (cell_count as usize)
        .checked_mul(RULESET_FULL_CELL_REF_STRIDE as usize)
        .ok_or_else(|| Error::new(ErrorKind::InvalidInput, "cell refs size overflow"))?;
    let mut cell_refs: Vec<u8> = Vec::with_capacity(cell_refs_capacity);

    for cell in cells {
        let canonical = ruleset_to_canonical_bytes(&cell.ruleset);
        let next_id = dict.len() as u32;
        let dict_id = *dict.entry(canonical).or_insert(next_id);

        cell_refs.extend_from_slice(&cell.pos[0].to_le_bytes());
        cell_refs.extend_from_slice(&cell.pos[1].to_le_bytes());
        cell_refs.extend_from_slice(&cell.pos[2].to_le_bytes());
        cell_refs.extend_from_slice(&dict_id.to_le_bytes());
    }

    let dict_count = dict.len() as u32;

    // Assemble payload: header + dictionary + cell refs
    let dict_storage_size = dict_count
        .checked_mul(ruleset_byte_size)
        .ok_or_else(|| Error::new(ErrorKind::InvalidInput, "ruleset dictionary size overflow"))?
        as usize;
    let total = (RULESET_FULL_HEADER_SIZE as usize)
        .checked_add(dict_storage_size)
        .and_then(|v| v.checked_add(cell_refs.len()))
        .ok_or_else(|| {
            Error::new(
                ErrorKind::InvalidInput,
                "full ruleset payload size overflow",
            )
        })?;
    let mut payload = Vec::with_capacity(total);

    // Header
    payload.extend_from_slice(&RULESET_FULL_MAGIC);
    payload.extend_from_slice(&RULESET_FULL_FORMAT_VERSION.to_le_bytes());
    payload.extend_from_slice(&0u32.to_le_bytes()); // flags (reserved)
    payload.extend_from_slice(&dict_count.to_le_bytes());
    payload.extend_from_slice(&cell_count.to_le_bytes());
    payload.extend_from_slice(&ruleset_byte_size.to_le_bytes());

    // Dictionary: write entries in dict_id order (0..dict_count-1)
    let mut entries: Vec<([u8; RULESET_FULL_CANONICAL_SIZE as usize], u32)> =
        dict.into_iter().collect();
    entries.sort_by_key(|(_, id)| *id);
    for (bytes, _id) in &entries {
        payload.extend_from_slice(bytes);
    }

    // Per-cell refs (already in cell iteration order)
    payload.extend_from_slice(&cell_refs);

    let path = binary_file_path(out, tick, "rulesets");
    write_binary_payload(&path, &payload, out)
}

#[derive(Serialize)]
struct RunMetaFile {
    #[serde(flatten)]
    core: RunMeta,
    snapshot_interval: u32,
    max_ticks: u32,
    field_count: usize,
    cell_header_byte_len: u32,
    cell_count_source: &'static str,
    cell_pos_units: &'static str,
    cell_record_layout: &'static str,
    ruleset_interval: u32,
    ruleset_output_mode: &'static str,
    ruleset_layer_file_pattern: String,
    ruleset_layer_record_stride: u32,
    ruleset_layer_record_layout: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    ruleset_full_file_pattern: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ruleset_full_header_size: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ruleset_full_ruleset_byte_size: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ruleset_full_cell_ref_stride: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ruleset_full_magic_ascii: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ruleset_full_format_version: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ruleset_full_payload_layout: Option<&'static str>,
}

pub fn write_run_meta(out: &OutputConfig) -> std::io::Result<()> {
    fs::create_dir_all(&out.output_dir)?;
    let path = Path::new(&out.output_dir).join("run_meta.json");
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    let field_count = GRID_X * GRID_Y * GRID_Z * S_EXT;
    let field_byte_len =
        field_byte_len(GRID_X as u32, GRID_Y as u32, GRID_Z as u32, S_EXT as u32).unwrap_or(0);
    let ruleset_stride = ruleset_layer_record_stride(
        S_RECEPTORS as u32,
        S_TRANSPORTERS as u32,
        R_MAX as u32,
        S_EFFECTORS as u32,
    )
    .unwrap_or(0);

    let mode = out.ruleset_output_mode;
    let writes_full = mode.writes_full_dump();

    let mut core = RunMeta::new(
        GRID_X as u32,
        GRID_Y as u32,
        GRID_Z as u32,
        S_EXT as u32,
        M_INT as u32,
        out.write_binary_field,
        out.write_binary_cells,
    );
    core.binary_compression = out.binary_compression.as_str().to_string();
    core.binary_compression_level = match out.binary_compression {
        BinaryCompression::None => 0,
        BinaryCompression::Zstd => out.binary_compression_level,
    };
    core.field_file_pattern = field_file_pattern(out);
    core.cell_file_pattern = cell_file_pattern(out);
    core.field_byte_len = field_byte_len;
    core.cell_record_stride = CELL_RECORD_STRIDE;

    let meta = RunMetaFile {
        core,
        snapshot_interval: out.snapshot_interval,
        max_ticks: out.max_ticks,
        field_count,
        cell_header_byte_len: 0,
        cell_count_source: "file_size_divided_by_cell_record_stride",
        cell_pos_units: "grid_voxel_indices",
        cell_record_layout: "pos:f32[3],lineage_id:u64,starter_type:u8,energy:f32",
        ruleset_interval: out.ruleset_interval,
        ruleset_output_mode: mode.as_str(),
        ruleset_layer_file_pattern: ruleset_layer_file_pattern(out),
        ruleset_layer_record_stride: ruleset_stride,
        ruleset_layer_record_layout: RULESET_LAYER_RECORD_LAYOUT,
        ruleset_full_file_pattern: if writes_full {
            Some(ruleset_full_file_pattern(out))
        } else {
            None
        },
        ruleset_full_header_size: writes_full.then_some(RULESET_FULL_HEADER_SIZE),
        ruleset_full_ruleset_byte_size: writes_full.then_some(RULESET_FULL_CANONICAL_SIZE),
        ruleset_full_cell_ref_stride: writes_full.then_some(RULESET_FULL_CELL_REF_STRIDE),
        ruleset_full_magic_ascii: writes_full.then_some("MRSF"),
        ruleset_full_format_version: writes_full.then_some(RULESET_FULL_FORMAT_VERSION),
        ruleset_full_payload_layout: writes_full.then_some(RULESET_FULL_PAYLOAD_LAYOUT),
    };

    serde_json::to_writer_pretty(writer, &meta)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::field_reassign_with_default)]

    use super::*;
    use marl_cell::cell::{
        EffectorParams, FateParams, Reaction, ReceptorParams, Ruleset, TransportParams,
    };
    use std::io::Read;

    fn test_output_dir(name: &str) -> String {
        let dir = format!("/tmp/opencode/{name}");
        let _ = fs::remove_dir_all(&dir);
        dir
    }

    fn read_maybe_zstd(path: &Path) -> Vec<u8> {
        let mut bytes = Vec::new();
        File::open(path).unwrap().read_to_end(&mut bytes).unwrap();
        if path.extension().and_then(|s| s.to_str()) == Some("zst") {
            zstd::stream::decode_all(&bytes[..]).unwrap()
        } else {
            bytes
        }
    }

    fn test_ruleset(seed: f32) -> Ruleset {
        Ruleset {
            receptors: std::array::from_fn(|i| ReceptorParams {
                k_half: seed + i as f32,
                n_hill: seed + 0.5 + i as f32,
                gain: seed + 1.0 + i as f32,
            }),
            transport: std::array::from_fn(|i| TransportParams {
                uptake_rate: seed + 2.0 + i as f32,
                secrete_rate: seed + 3.0 + i as f32,
                ext_species: i as u8,
                int_species: i as u8,
            }),
            reactions: std::array::from_fn(|i| Reaction {
                substrate: i as u8,
                product: (i + 1) as u8,
                catalyst: (i + 2) as u8,
                cofactor: 0xFF,
                k_m: seed + 4.0 + i as f32,
                v_max: seed + 5.0 + i as f32,
                k_cat: seed + 6.0 + i as f32,
            }),
            effectors: std::array::from_fn(|i| EffectorParams {
                threshold: seed + 7.0 + i as f32,
                rate: seed + 8.0 + i as f32,
                int_species: i as u8,
                ext_species: i as u8,
            }),
            fate: FateParams {
                division_energy: seed + 9.0,
                death_energy: seed + 10.0,
                quiescence_energy: seed + 11.0,
                division_prep_ticks: seed + 12.0,
            },
            hgt_propensity: seed + 13.0,
            mutation_rate: seed + 14.0,
        }
    }

    fn test_cell(pos: [u16; 3], lineage_id: u64, seed: f32) -> CellState {
        let mut internal = [0.0f32; M_INT];
        internal[0] = seed;
        CellState {
            pos,
            lineage_id,
            age: 0,
            internal,
            ruleset: test_ruleset(seed),
            quiescent: false,
            starter_type: 1,
            prep_remaining: 0,
        }
    }

    #[test]
    fn viewer_cell_record_layout_is_stable() {
        assert_eq!(mem::size_of::<ViewerCellRecord>(), 25);
        assert_eq!(mem::align_of::<ViewerCellRecord>(), 1);
    }

    #[test]
    fn viewer_cell_record_bytes_match_metadata_layout() {
        let cell = ViewerCellRecord {
            pos: [1.0, 2.0, 3.0],
            lineage_id: 0x0102_0304_0506_0708,
            starter_type: 2,
            energy: 4.5,
        };
        let bytes = as_bytes(std::slice::from_ref(&cell));

        assert_eq!(bytes.len(), 25);
        assert_eq!(&bytes[0..4], &1.0f32.to_le_bytes());
        assert_eq!(&bytes[4..8], &2.0f32.to_le_bytes());
        assert_eq!(&bytes[8..12], &3.0f32.to_le_bytes());
        assert_eq!(&bytes[12..20], &0x0102_0304_0506_0708u64.to_le_bytes());
        assert_eq!(bytes[20], 2);
        assert_eq!(&bytes[21..25], &4.5f32.to_le_bytes());
    }

    #[test]
    fn commit_temp_path_replaces_existing_file() {
        let out_dir = test_output_dir("marl_output_commit_replace_test");
        fs::create_dir_all(&out_dir).unwrap();
        let path = Path::new(&out_dir).join("payload.bin");
        let temp_path = Path::new(&out_dir).join("payload.bin.tmp");
        fs::write(&path, b"old").unwrap();
        fs::write(&temp_path, b"new").unwrap();

        commit_temp_path(&temp_path, &path).unwrap();

        assert_eq!(fs::read(&path).unwrap(), b"new");
        assert!(!temp_path.exists());
    }

    #[test]
    fn write_run_meta_includes_compression_and_ruleset_fields() {
        let out_dir = test_output_dir("marl_output_meta_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();
        out.binary_compression = BinaryCompression::Zstd;
        out.binary_compression_level = 3;
        out.ruleset_output_mode = RulesetOutputMode::LayerAverages;
        out.ruleset_interval = 777;

        write_run_meta(&out).unwrap();

        let meta_bytes = fs::read(Path::new(&out_dir).join("run_meta.json")).unwrap();
        let meta_json: serde_json::Value = serde_json::from_slice(&meta_bytes).unwrap();
        let meta: RunMeta = serde_json::from_slice(&meta_bytes).unwrap();
        assert!(meta.validate().is_ok());
        assert_eq!(
            meta.binary_compression,
            marl_format::BINARY_COMPRESSION_ZSTD
        );
        assert_eq!(meta.binary_compression_level, 3);
        assert_eq!(meta.field_file_pattern, "tick_<T>.field.bin.zst");
        assert_eq!(meta.cell_file_pattern, "tick_<T>.cells.bin.zst");
        assert_eq!(meta_json["ruleset_output_mode"], "layer_averages");
        assert_eq!(meta_json["ruleset_interval"], 777);
        assert_eq!(
            meta_json["ruleset_layer_file_pattern"],
            "tick_<T>.ruleset_layers.bin.zst"
        );
        assert_eq!(meta_json["ruleset_layer_record_stride"], 448);
    }

    #[test]
    fn write_cell_dump_compresses_when_requested() {
        let out_dir = test_output_dir("marl_output_cells_zstd_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();
        out.binary_compression = BinaryCompression::Zstd;
        out.binary_compression_level = 3;

        let cells = vec![test_cell([1, 2, 3], 42, 4.5)];
        write_cell_dump(&cells, 9, &out).unwrap();

        let path = Path::new(&out_dir).join("tick_9.cells.bin.zst");
        assert!(path.exists());
        let decoded = read_maybe_zstd(&path);
        assert_eq!(decoded.len(), CELL_RECORD_STRIDE as usize);
    }

    #[test]
    fn write_ruleset_layer_dump_writes_expected_stride_and_average() {
        let out_dir = test_output_dir("marl_output_ruleset_layer_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();

        let cells = vec![
            test_cell([0, 0, 0], 1, 1.0),
            test_cell([1, 1, 0], 2, 3.0),
            test_cell([2, 2, 1], 3, 5.0),
        ];
        write_ruleset_layer_dump(&cells, 7, &out).unwrap();

        let path = Path::new(&out_dir).join("tick_7.ruleset_layers.bin.zst");
        let bytes = read_maybe_zstd(&path);
        assert_eq!(bytes.len(), GRID_Z * RULESET_LAYER_RECORD_STRIDE);

        let first = &bytes[..RULESET_LAYER_RECORD_STRIDE];
        assert_eq!(u16::from_le_bytes(first[0..2].try_into().unwrap()), 0);
        assert_eq!(u32::from_le_bytes(first[4..8].try_into().unwrap()), 2);
        let avg_first_k_half = f32::from_le_bytes(first[8..12].try_into().unwrap());
        assert!((avg_first_k_half - 2.0).abs() < 1e-6);

        let second = &bytes[RULESET_LAYER_RECORD_STRIDE..2 * RULESET_LAYER_RECORD_STRIDE];
        assert_eq!(u16::from_le_bytes(second[0..2].try_into().unwrap()), 1);
        assert_eq!(u32::from_le_bytes(second[4..8].try_into().unwrap()), 1);
        let second_first_k_half = f32::from_le_bytes(second[8..12].try_into().unwrap());
        assert!((second_first_k_half - 5.0).abs() < 1e-6);

        let third = &bytes[2 * RULESET_LAYER_RECORD_STRIDE..3 * RULESET_LAYER_RECORD_STRIDE];
        assert_eq!(u32::from_le_bytes(third[4..8].try_into().unwrap()), 0);
        let zero_k_half = f32::from_le_bytes(third[8..12].try_into().unwrap());
        assert_eq!(zero_k_half, 0.0);
    }

    #[test]
    fn canonical_bytes_identical_rulesets_produce_identical_bytes() {
        let a = test_ruleset(1.0);
        let b = test_ruleset(1.0);
        let ba = ruleset_to_canonical_bytes(&a);
        let bb = ruleset_to_canonical_bytes(&b);
        assert_eq!(ba.len(), RULESET_FULL_CANONICAL_SIZE as usize);
        assert_eq!(ba, bb);
    }

    #[test]
    fn canonical_bytes_different_rulesets_produce_different_bytes() {
        let a = test_ruleset(1.0);
        let b = test_ruleset(2.0);
        let ba = ruleset_to_canonical_bytes(&a);
        let bb = ruleset_to_canonical_bytes(&b);
        assert_ne!(ba, bb);
    }

    #[test]
    fn canonical_bytes_roundtrip_first_receptor() {
        let ruleset = test_ruleset(1.0);
        let bytes = ruleset_to_canonical_bytes(&ruleset);
        let k_half = f32::from_le_bytes(bytes[0..4].try_into().unwrap());
        assert!((k_half - 1.0).abs() < 1e-6);
        let n_hill = f32::from_le_bytes(bytes[4..8].try_into().unwrap());
        assert!((n_hill - 1.5).abs() < 1e-6);
        let gain = f32::from_le_bytes(bytes[8..12].try_into().unwrap());
        assert!((gain - 2.0).abs() < 1e-6);
    }

    #[test]
    fn write_ruleset_full_dump_deduplicates() {
        let out_dir = test_output_dir("marl_output_ruleset_full_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();

        // 5 cells: 2 share ruleset A, 2 share ruleset B, 1 unique C
        let cells = vec![
            test_cell([0, 0, 0], 1, 1.0), // A
            test_cell([1, 0, 0], 2, 1.0), // A (same ruleset)
            test_cell([2, 0, 0], 3, 3.0), // B
            test_cell([3, 0, 0], 4, 3.0), // B (same ruleset)
            test_cell([4, 0, 0], 5, 5.0), // C (unique)
        ];
        write_ruleset_full_dump(&cells, 10, &out).unwrap();

        let path = Path::new(&out_dir).join("tick_10.rulesets.bin.zst");
        let bytes = read_maybe_zstd(&path);

        // Parse header
        assert_eq!(&bytes[0..4], b"MRSF");
        let version = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
        assert_eq!(version, 1);
        let _flags = u32::from_le_bytes(bytes[8..12].try_into().unwrap());
        let dict_count = u32::from_le_bytes(bytes[12..16].try_into().unwrap());
        let cell_count = u32::from_le_bytes(bytes[16..20].try_into().unwrap());
        let ruleset_size = u32::from_le_bytes(bytes[20..24].try_into().unwrap());

        assert_eq!(dict_count, 3); // 3 unique rulesets
        assert_eq!(cell_count, 5);
        assert_eq!(ruleset_size, RULESET_FULL_CANONICAL_SIZE);

        let header_size = 24usize;
        let dict_end = header_size + dict_count as usize * ruleset_size as usize;
        let cell_refs_end = dict_end + cell_count as usize * RULESET_FULL_CELL_REF_STRIDE as usize;
        assert_eq!(bytes.len(), cell_refs_end);

        // Check cell refs (after dictionary)
        // Cell 0: (0,0,0) with seed=1.0 → dict A (which should be dict_id 0)
        let ref0 = &bytes[dict_end..dict_end + 10];
        assert_eq!(u16::from_le_bytes(ref0[0..2].try_into().unwrap()), 0); // x
        assert_eq!(u16::from_le_bytes(ref0[2..4].try_into().unwrap()), 0); // y
        assert_eq!(u16::from_le_bytes(ref0[4..6].try_into().unwrap()), 0); // z
        let id0 = u32::from_le_bytes(ref0[6..10].try_into().unwrap());
        assert!(id0 < dict_count);

        // Cell 1: same ruleset as cell 0 → same dict_id
        let ref1 = &bytes[dict_end + 10..dict_end + 20];
        let id1 = u32::from_le_bytes(ref1[6..10].try_into().unwrap());
        assert_eq!(id0, id1, "identical rulesets must share dict_id");

        // Cell 2: seed=3.0 → different dict_id
        let ref2 = &bytes[dict_end + 20..dict_end + 30];
        let id2 = u32::from_le_bytes(ref2[6..10].try_into().unwrap());
        assert_ne!(id0, id2, "different rulesets must have different dict_id");

        // Cell 3: same as cell 2
        let ref3 = &bytes[dict_end + 30..dict_end + 40];
        let id3 = u32::from_le_bytes(ref3[6..10].try_into().unwrap());
        assert_eq!(id2, id3);

        // Cell 4: unique
        let ref4 = &bytes[dict_end + 40..dict_end + 50];
        let id4 = u32::from_le_bytes(ref4[6..10].try_into().unwrap());
        assert_ne!(id4, id0);
        assert_ne!(id4, id2);
    }

    #[test]
    fn write_ruleset_full_dump_empty_cells() {
        let out_dir = test_output_dir("marl_output_ruleset_full_empty_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();

        let cells: Vec<CellState> = Vec::new();
        write_ruleset_full_dump(&cells, 11, &out).unwrap();

        let path = Path::new(&out_dir).join("tick_11.rulesets.bin.zst");
        let bytes = read_maybe_zstd(&path);

        assert_eq!(&bytes[0..4], b"MRSF");
        let dict_count = u32::from_le_bytes(bytes[12..16].try_into().unwrap());
        let cell_count = u32::from_le_bytes(bytes[16..20].try_into().unwrap());
        assert_eq!(dict_count, 0);
        assert_eq!(cell_count, 0);
        assert_eq!(bytes.len(), 24); // just header
    }

    #[test]
    fn write_run_meta_includes_full_ruleset_fields_when_full_mode() {
        let out_dir = test_output_dir("marl_output_meta_full_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();
        out.ruleset_output_mode = RulesetOutputMode::Full;
        out.ruleset_interval = 500;

        write_run_meta(&out).unwrap();

        let meta_bytes = fs::read(Path::new(&out_dir).join("run_meta.json")).unwrap();
        let meta_json: serde_json::Value = serde_json::from_slice(&meta_bytes).unwrap();

        assert_eq!(meta_json["ruleset_output_mode"], "full");
        assert_eq!(
            meta_json["ruleset_full_file_pattern"],
            "tick_<T>.rulesets.bin.zst"
        );
        assert_eq!(meta_json["ruleset_full_header_size"], 24);
        assert_eq!(
            meta_json["ruleset_full_ruleset_byte_size"],
            RULESET_FULL_CANONICAL_SIZE
        );
        assert_eq!(meta_json["ruleset_full_cell_ref_stride"], 10);
        assert_eq!(meta_json["ruleset_full_magic_ascii"], "MRSF");
        assert_eq!(meta_json["ruleset_full_format_version"], 1);
        assert!(meta_json["ruleset_full_payload_layout"].is_string());
    }

    #[test]
    fn write_run_meta_omits_full_ruleset_fields_when_off() {
        let out_dir = test_output_dir("marl_output_meta_off_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();
        out.ruleset_output_mode = RulesetOutputMode::Off;

        write_run_meta(&out).unwrap();

        let meta_bytes = fs::read(Path::new(&out_dir).join("run_meta.json")).unwrap();
        let meta_json: serde_json::Value = serde_json::from_slice(&meta_bytes).unwrap();

        assert_eq!(meta_json["ruleset_output_mode"], "off");
        assert!(meta_json.get("ruleset_full_file_pattern").is_none());
        assert!(meta_json.get("ruleset_full_header_size").is_none());
    }

    #[test]
    fn write_ruleset_full_dump_with_zstd_compression() {
        let out_dir = test_output_dir("marl_output_ruleset_full_zstd_test");
        let mut out = OutputConfig::default();
        out.output_dir = out_dir.clone();
        out.binary_compression = BinaryCompression::Zstd;
        out.binary_compression_level = 3;

        let cells = vec![test_cell([0, 0, 0], 1, 1.0)];
        write_ruleset_full_dump(&cells, 12, &out).unwrap();

        let path = Path::new(&out_dir).join("tick_12.rulesets.bin.zst");
        assert!(path.exists());

        let path_no_zst = Path::new(&out_dir).join("tick_12.rulesets.bin");
        assert!(!path_no_zst.exists());

        let decoded = read_maybe_zstd(&path);
        assert_eq!(&decoded[0..4], b"MRSF");
        assert_eq!(u32::from_le_bytes(decoded[12..16].try_into().unwrap()), 1); // dict_count
        assert_eq!(u32::from_le_bytes(decoded[16..20].try_into().unwrap()), 1); // cell_count
    }
}
