use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use marl_format::RunMeta;

use crate::args::{CellMode, ViewerArgs};

// ---------------------------------------------------------------------------
// Loaded cell representation
// ---------------------------------------------------------------------------

/// A validated cell record unpacked from the binary dump.
#[derive(Debug, Clone)]
pub struct LoadedCell {
    pub pos: [u32; 3],
    pub lineage_id: u64,
    pub starter_type: u8,
    pub energy: f32,
}

// ---------------------------------------------------------------------------
// Snapshot payload
// ---------------------------------------------------------------------------

pub struct SnapshotPayload {
    pub meta: RunMeta,
    pub field_bytes: Vec<u8>,
    pub cells: Vec<LoadedCell>,
    pub tick: u64,
    pub species: u32,
    pub exposure: f32,
    pub density_scale: f32,
    pub steps: u32,
    pub cell_mode: CellMode,
    pub cell_alpha: f32,
}

// ---------------------------------------------------------------------------
// Snapshot loading
// ---------------------------------------------------------------------------

/// Load and validate `run_meta.json` from an output directory.
pub fn load_run_meta(output_dir: &std::path::Path) -> Result<RunMeta, Box<dyn Error>> {
    let meta_path = output_dir.join("run_meta.json");
    let meta_bytes =
        fs::read(&meta_path).map_err(|e| format!("failed to read {}: {e}", meta_path.display()))?;
    let meta: RunMeta = serde_json::from_slice(&meta_bytes)
        .map_err(|e| format!("failed to parse {}: {e}", meta_path.display()))?;
    meta.validate()
        .map_err(|e| format!("invalid run_meta.json: {e}"))?;
    Ok(meta)
}

pub fn load_snapshot(args: &ViewerArgs) -> Result<SnapshotPayload, Box<dyn Error>> {
    let meta = load_run_meta(&args.output_dir)?;

    if args.species >= meta.s_ext {
        return Err(format!(
            "species {} is out of range for {} external species",
            args.species, meta.s_ext
        )
        .into());
    }

    // --- tick_<T>.field.bin ---
    let field_path = snapshot_path(&args.output_dir, &meta.field_file_pattern, args.tick)?;
    let field_bytes = read_binary_payload(&field_path)?;
    if field_bytes.len() as u64 != meta.field_byte_len {
        return Err(format!(
            "{} has {} bytes, expected {} from run_meta.json",
            field_path.display(),
            field_bytes.len(),
            meta.field_byte_len
        )
        .into());
    }

    // --- tick_<T>.cells.bin ---
    let cells = if args.cell_mode == CellMode::Off {
        Vec::new()
    } else {
        load_cell_records(args, &meta)?
    };

    Ok(SnapshotPayload {
        meta,
        field_bytes,
        cells,
        tick: args.tick,
        species: args.species,
        exposure: args.exposure,
        density_scale: args.density_scale,
        steps: args.steps,
        cell_mode: args.cell_mode,
        cell_alpha: args.cell_alpha,
    })
}

// ---------------------------------------------------------------------------
// Tick discovery
// ---------------------------------------------------------------------------

/// Discover available snapshot ticks in an output directory.
///
/// Scans for files matching `tick_<digits>.field.bin`, parses the tick
/// number, sorts ascending, and deduplicates. Returns an empty `Vec` if
/// the directory is readable but contains no matching files.
pub fn discover_field_ticks(output_dir: &std::path::Path) -> Result<Vec<u64>, Box<dyn Error>> {
    let entries = fs::read_dir(output_dir)
        .map_err(|e| format!("failed to read directory {}: {e}", output_dir.display()))?;

    let mut ticks: Vec<u64> = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|e| {
            Box::<dyn Error>::from(format!(
                "failed to read entry in {}: {e}",
                output_dir.display()
            ))
        })?;
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if let Some(tick) = parse_field_tick_file_name(&name_str) {
            ticks.push(tick);
        }
    }

    ticks.sort_unstable();
    ticks.dedup();
    Ok(ticks)
}

/// Parse a tick number from a filename like `tick_42.field.bin`.
///
/// Returns `Some(tick)` if the name starts with `tick_`, ends with
/// `.field.bin`, and the middle portion is a valid `u64`.
fn parse_field_tick_file_name(name: &str) -> Option<u64> {
    let without_prefix = name.strip_prefix("tick_")?;
    let digits = without_prefix
        .strip_suffix(".field.bin.zst")
        .or_else(|| without_prefix.strip_suffix(".field.bin"))?;
    // Reject empty digit string or strings with non-digit characters
    if digits.is_empty() {
        return None;
    }
    if !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    digits.parse::<u64>().ok()
}

// ---------------------------------------------------------------------------
// Cell record parsing
// ---------------------------------------------------------------------------

fn load_cell_records(args: &ViewerArgs, meta: &RunMeta) -> Result<Vec<LoadedCell>, Box<dyn Error>> {
    let cells_path = snapshot_path(&args.output_dir, &meta.cell_file_pattern, args.tick)?;

    if !meta.write_binary_cells {
        eprintln!(
            "[viewer] warning: cell output was disabled for this run; rendering without cells"
        );
        return Ok(Vec::new());
    }
    if !cells_path.exists() {
        eprintln!(
            "[viewer] warning: cell file {} not found; rendering without cells",
            cells_path.display()
        );
        return Ok(Vec::new());
    }

    let raw = read_binary_payload(&cells_path)?;

    let stride = marl_format::CELL_RECORD_STRIDE as usize;
    if raw.len() % stride != 0 {
        return Err(format!(
            "{} length {} is not a multiple of cell record stride {}",
            cells_path.display(),
            raw.len(),
            stride
        )
        .into());
    }

    let count = raw.len() / stride;
    let mut cells = Vec::with_capacity(count);

    for i in 0..count {
        let offset = i * stride;
        let record = &raw[offset..offset + stride];
        let cell = parse_one_cell_record(record, i, &cells_path, meta)?;
        cells.push(cell);
    }

    Ok(cells)
}

fn snapshot_path(output_dir: &Path, pattern: &str, tick: u64) -> Result<PathBuf, Box<dyn Error>> {
    if !pattern.contains("<T>") {
        return Err(format!("snapshot pattern {pattern:?} does not contain <T>").into());
    }
    Ok(output_dir.join(pattern.replace("<T>", &tick.to_string())))
}

fn read_binary_payload(path: &Path) -> Result<Vec<u8>, Box<dyn Error>> {
    let raw = fs::read(path).map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    if path.extension().and_then(|ext| ext.to_str()) == Some("zst") {
        zstd::stream::decode_all(&raw[..])
            .map_err(|e| format!("failed to decompress {}: {e}", path.display()).into())
    } else {
        Ok(raw)
    }
}

/// Parse a single 25-byte cell record manually with from_le_bytes.
///
/// Does not transmute or take references to packed fields.
fn parse_one_cell_record(
    record: &[u8],
    index: usize,
    path: &PathBuf,
    meta: &RunMeta,
) -> Result<LoadedCell, Box<dyn Error>> {
    // Layout: pos_x:f32, pos_y:f32, pos_z:f32, lineage_id:u64, starter_type:u8, energy:f32

    let pos_x = f32::from_le_bytes(record[0..4].try_into().unwrap());
    let pos_y = f32::from_le_bytes(record[4..8].try_into().unwrap());
    let pos_z = f32::from_le_bytes(record[8..12].try_into().unwrap());
    let lineage_id = u64::from_le_bytes(record[12..20].try_into().unwrap());
    let starter_type = record[20];
    let energy = f32::from_le_bytes(record[21..25].try_into().unwrap());

    // Validate position: finite, non-negative, close to integer, within bounds
    for (val, dim, name) in [
        (pos_x, meta.grid_x, "x"),
        (pos_y, meta.grid_y, "y"),
        (pos_z, meta.grid_z, "z"),
    ] {
        if !val.is_finite() || val < 0.0 {
            return Err(format!(
                "{} record {}: position {name}={val} is not a finite non-negative float",
                path.display(),
                index
            )
            .into());
        }
        let rounded = val.round();
        if (val - rounded).abs() > 0.001 {
            return Err(format!(
                "{} record {}: position {name}={val} is not close to an integer voxel index",
                path.display(),
                index
            )
            .into());
        }
        if rounded as u32 >= dim {
            return Err(format!(
                "{} record {}: position {name}={val} out of bounds (grid_{name}={dim})",
                path.display(),
                index,
                dim = name
            )
            .into());
        }
    }

    if !energy.is_finite() {
        return Err(format!(
            "{} record {}: energy={energy} is not finite",
            path.display(),
            index
        )
        .into());
    }

    Ok(LoadedCell {
        pos: [pos_x as u32, pos_y as u32, pos_z as u32],
        lineage_id,
        starter_type,
        energy,
    })
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic 25-byte cell record for testing.
    fn make_record(x: f32, y: f32, z: f32, lineage: u64, stype: u8, energy: f32) -> Vec<u8> {
        let mut buf = Vec::with_capacity(25);
        buf.extend_from_slice(&x.to_le_bytes());
        buf.extend_from_slice(&y.to_le_bytes());
        buf.extend_from_slice(&z.to_le_bytes());
        buf.extend_from_slice(&lineage.to_le_bytes());
        buf.push(stype);
        buf.extend_from_slice(&energy.to_le_bytes());
        buf
    }

    fn test_meta() -> RunMeta {
        RunMeta::new(128, 128, 64, 12, 16, true, true)
    }

    #[test]
    fn parse_one_valid_record() {
        let record = make_record(10.0, 20.0, 5.0, 42, 1, 3.14);
        let meta = test_meta();
        let cell =
            parse_one_cell_record(&record, 0, &PathBuf::from("test.cells.bin"), &meta).unwrap();
        assert_eq!(cell.pos, [10, 20, 5]);
        assert_eq!(cell.lineage_id, 42);
        assert_eq!(cell.starter_type, 1);
        assert!((cell.energy - 3.14).abs() < 0.001);
    }

    #[test]
    fn parse_multiple_records() {
        let mut buf = Vec::new();
        buf.extend(&make_record(0.0, 0.0, 0.0, 1, 0, 0.5));
        buf.extend(&make_record(127.0, 127.0, 63.0, 2, 2, 1.0));
        buf.extend(&make_record(1.0, 2.0, 3.0, 3, 1, 2.0));

        let meta = test_meta();
        let cells =
            load_cell_records_from_bytes(&buf, &PathBuf::from("t.cells.bin"), &meta).unwrap();
        assert_eq!(cells.len(), 3);
        assert_eq!(cells[0].pos, [0, 0, 0]);
        assert_eq!(cells[1].pos, [127, 127, 63]);
        assert_eq!(cells[2].pos, [1, 2, 3]);
    }

    #[test]
    fn reject_bad_byte_length() {
        let buf = vec![0u8; 27]; // 27 % 25 != 0
        let meta = test_meta();
        let err = load_cell_records_from_bytes(&buf, &PathBuf::from("bad.bin"), &meta).unwrap_err();
        assert!(err.to_string().contains("not a multiple"));
    }

    #[test]
    fn reject_oob_position() {
        let record = make_record(128.0, 0.0, 0.0, 1, 0, 0.5);
        let meta = test_meta();
        let err = parse_one_cell_record(&record, 0, &PathBuf::from("oob.bin"), &meta).unwrap_err();
        assert!(err.to_string().contains("out of bounds"));
    }

    #[test]
    fn reject_non_integral_position() {
        let record = make_record(10.5, 0.0, 0.0, 1, 0, 0.5);
        let meta = test_meta();
        let err =
            parse_one_cell_record(&record, 0, &PathBuf::from("nonint.bin"), &meta).unwrap_err();
        assert!(err.to_string().contains("not close to an integer"));
    }

    #[test]
    fn reject_non_finite_energy() {
        let record = make_record(0.0, 0.0, 0.0, 1, 0, f32::NAN);
        let meta = test_meta();
        let err = parse_one_cell_record(&record, 0, &PathBuf::from("nan.bin"), &meta).unwrap_err();
        assert!(err.to_string().contains("energy"));
    }

    // --- Tick discovery tests ---

    #[test]
    fn parse_field_tick_file_name_valid() {
        assert_eq!(parse_field_tick_file_name("tick_0.field.bin"), Some(0));
        assert_eq!(parse_field_tick_file_name("tick_0.field.bin.zst"), Some(0));
        assert_eq!(parse_field_tick_file_name("tick_42.field.bin"), Some(42));
        assert_eq!(
            parse_field_tick_file_name("tick_18446744073709551615.field.bin"),
            Some(18446744073709551615)
        ); // u64::MAX
    }

    #[test]
    fn parse_field_tick_file_name_invalid() {
        assert_eq!(parse_field_tick_file_name("tick_.field.bin"), None);
        assert_eq!(parse_field_tick_file_name("field.bin"), None);
        assert_eq!(parse_field_tick_file_name("tick_0"), None);
        assert_eq!(parse_field_tick_file_name("tick_abc.field.bin"), None);
        assert_eq!(parse_field_tick_file_name("tick_1.cells.bin"), None);
        assert_eq!(parse_field_tick_file_name("tick_1.cells.bin.zst"), None);
        assert_eq!(parse_field_tick_file_name("tick_-1.field.bin"), None);
        assert_eq!(parse_field_tick_file_name("extra_tick_1.field.bin"), None);
        assert_eq!(parse_field_tick_file_name(""), None);
    }

    #[test]
    fn discover_field_ticks_sorts_dedups() {
        // Use a directory under /tmp/opencode for testing
        let dir = std::path::Path::new("/tmp/opencode/io_test_discover");
        let _ = fs::remove_dir_all(dir);
        fs::create_dir_all(dir).unwrap();
        // Create files in non-sorted order with duplicate tick
        for name in &[
            "tick_100.field.bin",
            "tick_3.field.bin",
            "tick_42.field.bin",
            "tick_42.field.bin.zst",
            "tick_3.field.bin", // duplicate — overwrites, same tick
        ] {
            fs::write(dir.join(name), b"data").unwrap();
        }
        // Also create some non-matching files
        fs::write(dir.join("tick_1.cells.bin"), b"data").unwrap();
        fs::write(dir.join("run_meta.json"), b"{}").unwrap();

        let ticks = discover_field_ticks(dir).unwrap();
        assert_eq!(ticks, vec![3, 42, 100]);

        // Cleanup
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn discover_field_ticks_empty_dir() {
        let dir = std::path::Path::new("/tmp/opencode/io_test_empty");
        let _ = fs::remove_dir_all(dir);
        fs::create_dir_all(dir).unwrap();
        let ticks = discover_field_ticks(dir).unwrap();
        assert!(ticks.is_empty());
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn discover_field_ticks_unreadable_dir() {
        let result = discover_field_ticks(std::path::Path::new("/nonexistent/path/12345"));
        assert!(result.is_err());
    }

    #[test]
    fn snapshot_path_uses_run_meta_pattern() {
        let path = snapshot_path(
            std::path::Path::new("/tmp/opencode/meta_pattern_test"),
            "tick_<T>.field.bin.zst",
            42,
        )
        .unwrap();
        assert_eq!(
            path,
            PathBuf::from("/tmp/opencode/meta_pattern_test/tick_42.field.bin.zst")
        );
    }

    #[test]
    fn read_binary_payload_decompresses_zstd() {
        let dir = std::path::Path::new("/tmp/opencode/io_test_zstd");
        let _ = fs::remove_dir_all(dir);
        fs::create_dir_all(dir).unwrap();
        let path = dir.join("sample.bin.zst");
        let payload = b"hello compressed world";
        fs::write(&path, zstd::stream::encode_all(&payload[..], 3).unwrap()).unwrap();

        let decoded = read_binary_payload(&path).unwrap();
        assert_eq!(decoded, payload);

        let _ = fs::remove_dir_all(dir);
    }

    /// Helper: parse cells from in-memory bytes for testing.
    fn load_cell_records_from_bytes(
        raw: &[u8],
        path: &PathBuf,
        meta: &RunMeta,
    ) -> Result<Vec<LoadedCell>, Box<dyn Error>> {
        let stride = marl_format::CELL_RECORD_STRIDE as usize;
        if raw.len() % stride != 0 {
            return Err(format!(
                "{} length {} is not a multiple of cell record stride {}",
                path.display(),
                raw.len(),
                stride
            )
            .into());
        }
        let count = raw.len() / stride;
        let mut cells = Vec::with_capacity(count);
        for i in 0..count {
            let offset = i * stride;
            let record = &raw[offset..offset + stride];
            cells.push(parse_one_cell_record(record, i, path, meta)?);
        }
        Ok(cells)
    }
}
