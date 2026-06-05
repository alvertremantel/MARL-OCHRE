//! `marl-format` - Binary schema shared by the MARL engine and viewer.
//!
//! This crate owns the durable on-disk format constants, metadata structs,
//! and validation helpers for the engine's binary field/cell dumps and
//! `run_meta.json`.
//!
//! # Constants
//!
//! - [`ENDIANNESS`]: data endianness (`"little"`)
//! - [`FIELD_DTYPE`]: field element dtype (`"f32"`)
//! - [`FIELD_LAYOUT`]: field memory layout (`"z_y_x_species"`)
//! - [`CELL_RECORD_STRIDE`]: size of one packed cell record in bytes (`25`)
//!
//! # Types
//!
//! - [`RunMeta`]: serializable metadata written to `run_meta.json`
//! - [`ViewerCellRecord`]: packed 25-byte viewer record for binary cell dumps
//! - [`FormatError`]: validation error type

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Data endianness for binary field and cell dumps.
pub const ENDIANNESS: &str = "little";

/// Field element dtype as written to disk.
pub const FIELD_DTYPE: &str = "f32";

/// Field memory layout: outer dimension is z, then y, then x, then species.
pub const FIELD_LAYOUT: &str = "z_y_x_species";

/// Size of one packed cell record in bytes.
pub const CELL_RECORD_STRIDE: u32 = 25;

/// Current binary format version written by the engine.
pub const BINARY_FORMAT_VERSION: u32 = 2;

/// Compression label for uncompressed binary payloads.
pub const BINARY_COMPRESSION_NONE: &str = "none";

/// Compression label for zstd-compressed binary payloads.
pub const BINARY_COMPRESSION_ZSTD: &str = "zstd";

/// Legacy/raw field snapshot filename pattern.
pub const FIELD_FILE_PATTERN_RAW: &str = "tick_<T>.field.bin";

/// Legacy/raw cell snapshot filename pattern.
pub const CELL_FILE_PATTERN_RAW: &str = "tick_<T>.cells.bin";

/// Per-layer ruleset average sidecar filename pattern.
pub const RULESET_LAYER_FILE_PATTERN_RAW: &str = "tick_<T>.ruleset_layers.bin";

/// Binary layout string for per-layer ruleset averages.
pub const RULESET_LAYER_RECORD_LAYOUT: &str = "z:u16,reserved:u16,cell_count:u32,receptors:{k_half:f32,n_hill:f32,gain:f32}[s_receptors],transport:{uptake_rate:f32,secrete_rate:f32,gate_weight:f32}[s_transporters],reactions:{k_m:f32,v_max:f32,k_cat:f32}[r_max],effectors:{threshold:f32,rate:f32}[s_effectors],fate:f32[4],hgt_propensity:f32,mutation_rate:f32";

/// Full per-cell deduplicated ruleset dump filename pattern.
pub const RULESET_FULL_FILE_PATTERN_RAW: &str = "tick_<T>.rulesets.bin";

/// Magic bytes for the full ruleset dump file header: "MRSF" in ASCII.
pub const RULESET_FULL_MAGIC: [u8; 4] = [b'M', b'R', b'S', b'F'];

/// Version of the full ruleset dump binary format.
pub const RULESET_FULL_FORMAT_VERSION: u32 = 2;

/// Byte size of the full ruleset file header (magic + version + flags + dict_count + cell_count + ruleset_byte_size).
pub const RULESET_FULL_HEADER_SIZE: u32 = 24;

/// Byte stride of one per-cell reference record in the full ruleset dump.
/// x:u16, y:u16, z:u16, dict_id:u32 = 10 bytes.
pub const RULESET_FULL_CELL_REF_STRIDE: u32 = 10;

/// Canonical byte size of a single deduplicated ruleset payload.
///
/// Layout (all little-endian; u8 fields are 1 byte, f32 fields are 4 bytes):
///   receptors: 8 × {k_half:f32, n_hill:f32, gain:f32}           = 96 B
///   transport: 8 × {uptake_rate:f32, secrete_rate:f32, ext_species:u8, int_species:u8, gate_receptor:u8, gate_weight:f32} = 120 B
///   reactions: 16 × {substrate:u8, product:u8, catalyst:u8, cofactor:u8, k_m:f32, v_max:f32, k_cat:f32} = 256 B
///   effectors: 8 × {threshold:f32, rate:f32, int_species:u8, ext_species:u8} = 80 B
///   fate: {division_energy:f32, death_energy:f32, quiescence_energy:f32, division_prep_ticks:f32} = 16 B
///   hgt_propensity: f32 = 4 B
///   mutation_rate: f32 = 4 B
///   total = 576
pub const RULESET_FULL_CANONICAL_SIZE: u32 = 576;

/// Binary layout string for one canonical ruleset payload.
pub const RULESET_FULL_PAYLOAD_LAYOUT: &str = "\
receptors:{k_half:f32,n_hill:f32,gain:f32}[s_receptors],\
transport:{uptake_rate:f32,secrete_rate:f32,ext_species:u8,int_species:u8,gate_receptor:u8,gate_weight:f32}[s_transporters],\
reactions:{substrate:u8,product:u8,catalyst:u8,cofactor:u8,k_m:f32,v_max:f32,k_cat:f32}[r_max],\
effectors:{threshold:f32,rate:f32,int_species:u8,ext_species:u8}[s_effectors],\
fate:{division_energy:f32,death_energy:f32,quiescence_energy:f32,division_prep_ticks:f32},\
hgt_propensity:f32,mutation_rate:f32";

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Validation error returned by [`RunMeta::validate`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormatError {
    pub message: String,
}

impl FormatError {
    fn new(msg: impl Into<String>) -> Self {
        Self {
            message: msg.into(),
        }
    }
}

impl std::fmt::Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FormatError: {}", self.message)
    }
}

impl std::error::Error for FormatError {}

// ---------------------------------------------------------------------------
// Run metadata
// ---------------------------------------------------------------------------

/// Metadata written to `run_meta.json` at startup when binary output is enabled.
///
/// All field names and value shapes are preserved for compatibility with
/// downstream binary-dump consumers.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RunMeta {
    #[serde(rename = "grid_x")]
    pub grid_x: u32,
    #[serde(rename = "grid_y")]
    pub grid_y: u32,
    #[serde(rename = "grid_z")]
    pub grid_z: u32,
    #[serde(rename = "s_ext")]
    pub s_ext: u32,
    #[serde(default, rename = "m_int")]
    pub m_int: u32,
    #[serde(rename = "field_dtype")]
    pub field_dtype: String,
    #[serde(rename = "field_layout")]
    pub field_layout: String,
    #[serde(rename = "field_byte_len")]
    pub field_byte_len: u64,
    #[serde(rename = "cell_record_stride")]
    pub cell_record_stride: u32,
    #[serde(rename = "endianness")]
    pub endianness: String,
    #[serde(
        default = "default_binary_format_version",
        rename = "binary_format_version"
    )]
    pub binary_format_version: u32,
    #[serde(default = "default_binary_compression", rename = "binary_compression")]
    pub binary_compression: String,
    #[serde(default, rename = "binary_compression_level")]
    pub binary_compression_level: i32,
    #[serde(rename = "write_binary_field")]
    pub write_binary_field: bool,
    #[serde(rename = "write_binary_cells")]
    pub write_binary_cells: bool,
    #[serde(default = "default_field_file_pattern", rename = "field_file_pattern")]
    pub field_file_pattern: String,
    #[serde(default = "default_cell_file_pattern", rename = "cell_file_pattern")]
    pub cell_file_pattern: String,
}

fn default_binary_format_version() -> u32 {
    1
}

fn default_binary_compression() -> String {
    BINARY_COMPRESSION_NONE.to_string()
}

fn default_field_file_pattern() -> String {
    FIELD_FILE_PATTERN_RAW.to_string()
}

fn default_cell_file_pattern() -> String {
    CELL_FILE_PATTERN_RAW.to_string()
}

impl RunMeta {
    /// Build a new `RunMeta` from grid dimensions, species counts, and output toggles.
    pub fn new(
        grid_x: u32,
        grid_y: u32,
        grid_z: u32,
        s_ext: u32,
        m_int: u32,
        write_binary_field: bool,
        write_binary_cells: bool,
    ) -> Self {
        let field_byte_len = field_byte_len(grid_x, grid_y, grid_z, s_ext).unwrap_or(0);
        Self {
            grid_x,
            grid_y,
            grid_z,
            s_ext,
            m_int,
            field_dtype: FIELD_DTYPE.to_string(),
            field_layout: FIELD_LAYOUT.to_string(),
            field_byte_len,
            cell_record_stride: CELL_RECORD_STRIDE,
            endianness: ENDIANNESS.to_string(),
            binary_format_version: BINARY_FORMAT_VERSION,
            binary_compression: BINARY_COMPRESSION_NONE.to_string(),
            binary_compression_level: 0,
            write_binary_field,
            write_binary_cells,
            field_file_pattern: FIELD_FILE_PATTERN_RAW.to_string(),
            cell_file_pattern: CELL_FILE_PATTERN_RAW.to_string(),
        }
    }

    /// Validate that stored schema constants match the loaded metadata.
    ///
    /// Returns `Ok(())` if all constants are consistent, or a [`FormatError`]
    /// describing the first mismatch.
    pub fn validate(&self) -> Result<(), FormatError> {
        if self.endianness != ENDIANNESS {
            return Err(FormatError::new(format!(
                "endianness mismatch: expected {}, got {}",
                ENDIANNESS, self.endianness
            )));
        }
        if self.field_dtype != FIELD_DTYPE {
            return Err(FormatError::new(format!(
                "field_dtype mismatch: expected {}, got {}",
                FIELD_DTYPE, self.field_dtype
            )));
        }
        if self.field_layout != FIELD_LAYOUT {
            return Err(FormatError::new(format!(
                "field_layout mismatch: expected {}, got {}",
                FIELD_LAYOUT, self.field_layout
            )));
        }
        if self.cell_record_stride != CELL_RECORD_STRIDE {
            return Err(FormatError::new(format!(
                "cell_record_stride mismatch: expected {}, got {}",
                CELL_RECORD_STRIDE, self.cell_record_stride
            )));
        }
        if self.binary_format_version == 0 || self.binary_format_version > BINARY_FORMAT_VERSION {
            return Err(FormatError::new(format!(
                "binary_format_version mismatch: expected 1..={}, got {}",
                BINARY_FORMAT_VERSION, self.binary_format_version
            )));
        }
        if self.binary_compression != BINARY_COMPRESSION_NONE
            && self.binary_compression != BINARY_COMPRESSION_ZSTD
        {
            return Err(FormatError::new(format!(
                "binary_compression mismatch: expected {} or {}, got {}",
                BINARY_COMPRESSION_NONE, BINARY_COMPRESSION_ZSTD, self.binary_compression
            )));
        }
        validate_snapshot_pattern("field_file_pattern", &self.field_file_pattern)?;
        validate_snapshot_pattern("cell_file_pattern", &self.cell_file_pattern)?;
        validate_compression_pattern(
            "field_file_pattern",
            &self.field_file_pattern,
            &self.binary_compression,
        )?;
        validate_compression_pattern(
            "cell_file_pattern",
            &self.cell_file_pattern,
            &self.binary_compression,
        )?;
        let expected_len = field_byte_len(self.grid_x, self.grid_y, self.grid_z, self.s_ext)
            .ok_or_else(|| {
                FormatError::new(
                    "field dimensions must be non-zero and field_byte_len must not overflow",
                )
            })?;
        if self.field_byte_len != expected_len {
            return Err(FormatError::new(format!(
                "field_byte_len mismatch: expected {}, got {}",
                expected_len, self.field_byte_len
            )));
        }
        Ok(())
    }
}

fn validate_snapshot_pattern(name: &str, pattern: &str) -> Result<(), FormatError> {
    if pattern.is_empty() {
        return Err(FormatError::new(format!("{name} must not be empty")));
    }
    if pattern.matches("<T>").count() != 1 {
        return Err(FormatError::new(format!(
            "{name} must contain exactly one <T> tick placeholder"
        )));
    }
    for component in std::path::Path::new(pattern).components() {
        match component {
            std::path::Component::Prefix(_) | std::path::Component::RootDir => {
                return Err(FormatError::new(format!(
                    "{name} must be relative to output_dir"
                )));
            }
            std::path::Component::ParentDir => {
                return Err(FormatError::new(format!(
                    "{name} must not contain parent directory components"
                )));
            }
            std::path::Component::CurDir | std::path::Component::Normal(_) => {}
        }
    }
    Ok(())
}

fn validate_compression_pattern(
    name: &str,
    pattern: &str,
    compression: &str,
) -> Result<(), FormatError> {
    let has_zst_suffix = pattern.ends_with(".zst");
    match compression {
        BINARY_COMPRESSION_ZSTD if !has_zst_suffix => Err(FormatError::new(format!(
            "{name} must end with .zst when binary_compression is {BINARY_COMPRESSION_ZSTD}"
        ))),
        BINARY_COMPRESSION_NONE if has_zst_suffix => Err(FormatError::new(format!(
            "{name} must not end with .zst when binary_compression is {BINARY_COMPRESSION_NONE}"
        ))),
        _ => Ok(()),
    }
}

// ---------------------------------------------------------------------------
// Field byte length helper
// ---------------------------------------------------------------------------

/// Compute the expected byte length of a raw field dump.
///
/// Returns `None` if any argument is zero (avoids panics from multiply-add).
#[inline]
pub fn field_byte_len(grid_x: u32, grid_y: u32, grid_z: u32, s_ext: u32) -> Option<u64> {
    if grid_x == 0 || grid_y == 0 || grid_z == 0 || s_ext == 0 {
        return None;
    }
    let count = u64::from(grid_x)
        .checked_mul(u64::from(grid_y))?
        .checked_mul(u64::from(grid_z))?
        .checked_mul(u64::from(s_ext))?;
    count.checked_mul(4)
}

/// Compute the number of `f32` averages stored in one per-layer ruleset record.
#[inline]
pub fn ruleset_layer_value_count(
    s_receptors: u32,
    s_transporters: u32,
    r_max: u32,
    s_effectors: u32,
) -> Option<u32> {
    s_receptors
        .checked_mul(3)?
        .checked_add(s_transporters.checked_mul(3)?)?
        .checked_add(r_max.checked_mul(3)?)?
        .checked_add(s_effectors.checked_mul(2)?)?
        .checked_add(6)
}

/// Compute the byte stride of one per-layer ruleset average record.
#[inline]
pub fn ruleset_layer_record_stride(
    s_receptors: u32,
    s_transporters: u32,
    r_max: u32,
    s_effectors: u32,
) -> Option<u32> {
    let values = ruleset_layer_value_count(s_receptors, s_transporters, r_max, s_effectors)?;
    8u32.checked_add(values.checked_mul(4)?)
}

// ---------------------------------------------------------------------------
// Packed viewer cell record
// ---------------------------------------------------------------------------

/// Packed 25-byte viewer cell record written to `tick_<T>.cells.bin`.
///
/// Each record contains position (3 × f32), lineage_id (u64), starter_type (u8),
/// and energy (f32) in little-endian byte order. It is not a full cell-state
/// checkpoint.
///
/// # Layout
///
/// | Offset | Size | Type       | Name          |
/// |--------|------|------------|---------------|
/// | 0      | 12   | f32[3]     | pos           |
/// | 12     | 8    | u64        | lineage_id    |
/// | 20     | 1    | u8         | starter_type  |
/// | 21     | 4    | f32        | energy        |
/// | 25     | —    | —          | total = 25    |
///
/// The struct is `#[repr(C, packed)]` to match the binary file layout.
/// Do not take references to multi-byte fields of a packed record;
/// copy values by value instead.
#[repr(C, packed)]
#[derive(Clone, Copy)]
pub struct ViewerCellRecord {
    /// World position of the cell (x, y, z).
    pub pos: [f32; 3],
    /// Unique lineage identifier.
    pub lineage_id: u64,
    /// Starter metabolism type encoded as a small integer.
    pub starter_type: u8,
    /// Current energy reserve.
    pub energy: f32,
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_byte_len_basic() {
        // 128 * 128 * 64 * 12 * 4 bytes
        let len = field_byte_len(128, 128, 64, 12);
        assert_eq!(len, Some(50_331_648));
    }

    #[test]
    fn test_field_byte_len_zero_arg() {
        assert_eq!(field_byte_len(0, 128, 64, 12), None);
        assert_eq!(field_byte_len(128, 0, 64, 12), None);
        assert_eq!(field_byte_len(128, 128, 0, 12), None);
        assert_eq!(field_byte_len(128, 128, 64, 0), None);
    }

    #[test]
    fn test_field_byte_len_overflow() {
        // very large values should not panic
        assert_eq!(field_byte_len(u32::MAX, u32::MAX, u32::MAX, u32::MAX), None);
    }

    #[test]
    fn test_run_meta_new_round_trip() {
        let meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        assert_eq!(meta.grid_x, 128);
        assert_eq!(meta.grid_y, 128);
        assert_eq!(meta.grid_z, 64);
        assert_eq!(meta.s_ext, 12);
        assert_eq!(meta.m_int, 8);
        assert_eq!(meta.field_dtype, "f32");
        assert_eq!(meta.field_layout, "z_y_x_species");
        assert_eq!(meta.endianness, "little");
        assert_eq!(meta.binary_format_version, BINARY_FORMAT_VERSION);
        assert_eq!(meta.binary_compression, "none");
        assert_eq!(meta.binary_compression_level, 0);
        assert_eq!(meta.cell_record_stride, 25);
        assert!(meta.write_binary_field);
        assert!(meta.write_binary_cells);
        assert_eq!(meta.field_byte_len, 50_331_648);
        assert_eq!(meta.field_file_pattern, FIELD_FILE_PATTERN_RAW);
        assert_eq!(meta.cell_file_pattern, CELL_FILE_PATTERN_RAW);
    }

    #[test]
    fn test_run_meta_validate_ok() {
        let meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        assert!(meta.validate().is_ok());
    }

    #[test]
    fn test_run_meta_validate_bad_endianness() {
        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.endianness = "big".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("endianness"));
    }

    #[test]
    fn test_run_meta_validate_bad_layout() {
        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.field_layout = "x_y_z_species".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("field_layout"));
    }

    #[test]
    fn test_run_meta_validate_bad_byte_len() {
        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.field_byte_len = 1;
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("field_byte_len"));
    }

    #[test]
    fn test_run_meta_validate_rejects_invalid_dimensions() {
        let meta = RunMeta::new(0, 128, 64, 12, 8, true, true);
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("dimensions"), "got: {err}");

        let mut meta = RunMeta::new(u32::MAX, u32::MAX, u32::MAX, u32::MAX, 8, true, true);
        meta.field_byte_len = u64::MAX;
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("overflow"), "got: {err}");
    }

    #[test]
    fn test_run_meta_validate_bad_compression() {
        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.binary_compression = "brotli".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("binary_compression"));
    }

    #[test]
    fn test_run_meta_validate_bad_snapshot_pattern() {
        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.field_file_pattern = "latest.field.bin".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("field_file_pattern"));

        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.cell_file_pattern = "latest.cells.bin".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("cell_file_pattern"));

        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.field_file_pattern = "tick_<T>_<T>.field.bin".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("exactly one <T>"));
    }

    #[test]
    fn test_run_meta_validate_rejects_escaping_snapshot_patterns() {
        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.field_file_pattern = "/tmp/tick_<T>.field.bin".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("relative to output_dir"), "got: {err}");

        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.cell_file_pattern = "../tick_<T>.cells.bin".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("parent directory"), "got: {err}");

        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.cell_file_pattern = "cells/../tick_<T>.cells.bin".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("parent directory"), "got: {err}");
    }

    #[test]
    fn test_run_meta_validate_compression_matches_patterns() {
        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.binary_compression = BINARY_COMPRESSION_ZSTD.to_string();
        meta.field_file_pattern = "tick_<T>.field.bin.zst".to_string();
        meta.cell_file_pattern = "tick_<T>.cells.bin.zst".to_string();
        assert!(meta.validate().is_ok());

        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.binary_compression = BINARY_COMPRESSION_ZSTD.to_string();
        meta.field_file_pattern = "tick_<T>.field.bin.zst".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("cell_file_pattern"), "got: {err}");
        assert!(err.message.contains(".zst"), "got: {err}");

        let mut meta = RunMeta::new(128, 128, 64, 12, 8, true, true);
        meta.field_file_pattern = "tick_<T>.field.bin.zst".to_string();
        let err = meta.validate().unwrap_err();
        assert!(err.message.contains("field_file_pattern"), "got: {err}");
        assert!(err.message.contains("must not end with .zst"), "got: {err}");
    }

    #[test]
    fn test_run_meta_legacy_defaults_still_validate() {
        let json = r#"{
            "grid_x": 128,
            "grid_y": 128,
            "grid_z": 64,
            "s_ext": 12,
            "m_int": 8,
            "field_dtype": "f32",
            "field_layout": "z_y_x_species",
            "field_byte_len": 50331648,
            "cell_record_stride": 25,
            "endianness": "little",
            "write_binary_field": true,
            "write_binary_cells": true
        }"#;
        let meta: RunMeta = serde_json::from_str(json).unwrap();
        assert_eq!(meta.binary_format_version, 1);
        assert_eq!(meta.binary_compression, BINARY_COMPRESSION_NONE);
        assert_eq!(meta.field_file_pattern, FIELD_FILE_PATTERN_RAW);
        assert_eq!(meta.cell_file_pattern, CELL_FILE_PATTERN_RAW);
        assert!(meta.validate().is_ok());
    }

    #[test]
    fn test_run_meta_serde_round_trip() {
        let meta = RunMeta::new(128, 128, 64, 12, 8, true, false);
        let json = serde_json::to_string_pretty(&meta).unwrap();
        let back: RunMeta = serde_json::from_str(&json).unwrap();
        assert_eq!(meta.grid_x, back.grid_x);
        assert_eq!(meta.grid_y, back.grid_y);
        assert_eq!(meta.grid_z, back.grid_z);
        assert_eq!(meta.s_ext, back.s_ext);
        assert_eq!(meta.m_int, back.m_int);
        assert_eq!(meta.write_binary_field, back.write_binary_field);
        assert_eq!(meta.write_binary_cells, back.write_binary_cells);
        assert_eq!(meta.binary_format_version, back.binary_format_version);
        assert_eq!(meta.binary_compression, back.binary_compression);
        assert_eq!(meta.field_file_pattern, back.field_file_pattern);
        assert_eq!(meta.cell_file_pattern, back.cell_file_pattern);
        // field names match JSON
        assert!(json.contains("\"grid_x\""));
        assert!(json.contains("\"field_dtype\""));
        assert!(json.contains("\"cell_record_stride\""));
        assert!(json.contains("\"binary_compression\""));
    }

    #[test]
    fn test_ruleset_layer_stride_helpers() {
        assert_eq!(ruleset_layer_value_count(8, 8, 16, 8), Some(118));
        assert_eq!(ruleset_layer_record_stride(8, 8, 16, 8), Some(480));
    }

    #[test]
    fn test_viewer_cell_record_size() {
        use std::mem::size_of;
        assert_eq!(size_of::<ViewerCellRecord>(), 25);
    }

    #[test]
    fn test_constants() {
        assert_eq!(ENDIANNESS, "little");
        assert_eq!(FIELD_DTYPE, "f32");
        assert_eq!(FIELD_LAYOUT, "z_y_x_species");
        assert_eq!(CELL_RECORD_STRIDE, 25);
        assert_eq!(FIELD_FILE_PATTERN_RAW, "tick_<T>.field.bin");
        assert_eq!(CELL_FILE_PATTERN_RAW, "tick_<T>.cells.bin");
        assert_eq!(
            RULESET_LAYER_FILE_PATTERN_RAW,
            "tick_<T>.ruleset_layers.bin"
        );
    }
}
