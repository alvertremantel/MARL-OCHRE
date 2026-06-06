// snapshot.rs — Visual snapshot module for field cross-sections and cell density maps.
//
// Writes PPM (Portable PixMap) image files that can be opened by most image
// viewers and converted to PNG with ImageMagick (`convert file.ppm file.png`).
//
// PPM P6 binary format is dead simple:
//   Line 1: "P6"
//   Line 2: "{width} {height}"
//   Line 3: "255"              (max color value)
//   Then:   width * height * 3 raw bytes (R, G, B per pixel, left-to-right,
//           top-to-bottom)
//
// No external image library needed — we write raw bytes directly.

use marl_config::*;
use marl_field::field::Field;
use marl_field::light::LightField;
use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

// ---------------------------------------------------------------------------
// Species naming
// ---------------------------------------------------------------------------

/// Map external species index (0..S_EXT) to a human-readable name suitable
/// for filenames. These match the comments in config.rs.
pub fn species_name(s: usize) -> &'static str {
    external_species_name(s)
}

// ---------------------------------------------------------------------------
// Color mapping
// ---------------------------------------------------------------------------

/// Map a scalar value to an RGB color using a scientific "jet-like" colormap.
///
/// The colormap goes through five anchor colors in order:
///   blue (0,0,255) -> cyan (0,255,255) -> green (0,255,0)
///     -> yellow (255,255,0) -> red (255,0,0)
///
/// This gives good perceptual separation across the range, making it easy to
/// distinguish low, medium, and high concentrations at a glance.
///
/// `val` is the raw concentration; `max_val` is the top of the scale.
/// Values are clamped to [0, max_val]. If max_val <= 0, returns black.
pub fn value_to_color(val: f32, max_val: f32) -> [u8; 3] {
    if max_val <= 0.0 {
        return [0, 0, 0];
    }

    // Normalize to 0..1
    let t = (val / max_val).clamp(0.0, 1.0);

    // Linearly interpolate between anchor colors in four segments.
    let (r, g, b) = if t < 0.25 {
        // Blue -> Cyan: red stays 0, green rises, blue stays 255
        let f = t / 0.25;
        (0.0, 255.0 * f, 255.0)
    } else if t < 0.50 {
        // Cyan -> Green: red stays 0, green stays 255, blue falls
        let f = (t - 0.25) / 0.25;
        (0.0, 255.0, 255.0 * (1.0 - f))
    } else if t < 0.75 {
        // Green -> Yellow: red rises, green stays 255, blue stays 0
        let f = (t - 0.50) / 0.25;
        (255.0 * f, 255.0, 0.0)
    } else {
        // Yellow -> Red: red stays 255, green falls, blue stays 0
        let f = (t - 0.75) / 0.25;
        (255.0, 255.0 * (1.0 - f), 0.0)
    };

    [r as u8, g as u8, b as u8]
}

// ---------------------------------------------------------------------------
// PPM writer helper
// ---------------------------------------------------------------------------

/// Write raw RGB pixel data as a P6 (binary) PPM file.
///
/// `pixels` must contain exactly `width * height` entries, stored row-major
/// (left-to-right, top-to-bottom).
fn write_ppm(path: &Path, width: usize, height: usize, pixels: &[[u8; 3]]) -> std::io::Result<()> {
    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);

    // PPM P6 header: magic, dimensions, max value, then raw RGB bytes.
    write!(w, "P6\n{} {}\n255\n", width, height)?;

    for &[r, g, b] in pixels {
        w.write_all(&[r, g, b])?;
    }

    w.flush()
}

// ---------------------------------------------------------------------------
// XZ cross-section (vertical slice through the column midpoint)
// ---------------------------------------------------------------------------

/// Write a PPM image of a chemical species concentration in the XZ plane.
///
/// Physical interpretation: this is a vertical cross-section through the
/// middle of the domain (at y = grid.y / 2). The X axis runs horizontally
/// (left to right), and the Z axis runs vertically with z=0 at the **top**
/// of the image — matching the physical layout where light and nutrients
/// enter from the surface (z=0) and depth increases downward.
///
/// The color scale auto-ranges to the maximum value found in this slice,
/// so even low-concentration species produce useful images.
pub fn write_xz_cross_section(
    field: &Field,
    species: usize,
    tick: u64,
    out_dir: &str,
) -> std::io::Result<()> {
    let grid = field.grid();
    let y_mid = grid.y / 2;
    let width = grid.x;
    let height = grid.z;

    // First pass: find the maximum concentration in this slice for auto-scaling.
    let mut max_val: f32 = 0.0;
    for z in 0..grid.z {
        for x in 0..grid.x {
            let v = field.get(x, y_mid, z, species);
            if v > max_val {
                max_val = v;
            }
        }
    }

    // Second pass: map concentrations to colors.
    let mut pixels = Vec::with_capacity(width * height);
    for z in 0..grid.z {
        for x in 0..grid.x {
            let v = field.get(x, y_mid, z, species);
            pixels.push(value_to_color(v, max_val));
        }
    }

    let name = species_name(species);
    let filename = format!("xz_{}_t{}.ppm", name, tick);
    let path = Path::new(out_dir).join(filename);
    write_ppm(&path, width, height, &pixels)
}

// ---------------------------------------------------------------------------
// XY cross-section (horizontal slice at a given depth z)
// ---------------------------------------------------------------------------

/// Write a PPM image of a chemical species concentration in the XY plane
/// at a given z-depth.
///
/// Physical interpretation: this is a top-down horizontal slice at depth `z`.
/// X runs left-to-right, Y runs top-to-bottom. Think of it as looking down
/// into the column from above and seeing the concentration at one specific
/// depth layer.
pub fn write_xy_cross_section(
    field: &Field,
    species: usize,
    z: usize,
    tick: u64,
    out_dir: &str,
) -> std::io::Result<()> {
    let grid = field.grid();
    let width = grid.x;
    let height = grid.y;

    // First pass: find max for auto-scaling.
    let mut max_val: f32 = 0.0;
    for y in 0..grid.y {
        for x in 0..grid.x {
            let v = field.get(x, y, z, species);
            if v > max_val {
                max_val = v;
            }
        }
    }

    // Second pass: render.
    let mut pixels = Vec::with_capacity(width * height);
    for y in 0..grid.y {
        for x in 0..grid.x {
            let v = field.get(x, y, z, species);
            pixels.push(value_to_color(v, max_val));
        }
    }

    let name = species_name(species);
    let filename = format!("xy_z{}_{}_t{}.ppm", z, name, tick);
    let path = Path::new(out_dir).join(filename);
    write_ppm(&path, width, height, &pixels)
}

// ---------------------------------------------------------------------------
// Cell density heatmap (XZ plane)
// ---------------------------------------------------------------------------

/// Write a grayscale PPM image showing cell density in the XZ plane at
/// y = grid.y / 2.
///
/// Currently each voxel holds at most one cell, so this is effectively a
/// binary image (black = empty, white = occupied). However, the code counts
/// all cells with matching (x, z) at the target y-slice, so it will work
/// correctly if we later allow multiple cells per voxel or do coarse-grained
/// averaging across y-neighbors.
///
/// Grayscale is achieved by setting R = G = B to the same intensity value.
pub fn write_cell_density_xz(
    grid: GridDims,
    cells: &HashMap<[u16; 3], usize>,
    tick: u64,
    out_dir: &str,
) -> std::io::Result<()> {
    let y_mid = grid.y / 2;
    let width = grid.x;
    let height = grid.z;

    // Count cells at each (x, z) position in this y-slice.
    let mut counts = vec![0u32; width * height];
    for pos in cells.keys() {
        if pos[1] as usize == y_mid {
            let x = pos[0] as usize;
            let z = pos[2] as usize;
            if x < grid.x && z < grid.z {
                counts[z * width + x] += 1;
            }
        }
    }

    // Find max count for scaling (usually 1, but future-proofed).
    let max_count = counts.iter().copied().max().unwrap_or(0);

    let mut pixels = Vec::with_capacity(width * height);
    for &c in &counts {
        let intensity = if max_count > 0 {
            ((c as f32 / max_count as f32) * 255.0) as u8
        } else {
            0
        };
        pixels.push([intensity, intensity, intensity]);
    }

    let filename = format!("density_xz_t{}.ppm", tick);
    let path = Path::new(out_dir).join(filename);
    write_ppm(&path, width, height, &pixels)
}

// ---------------------------------------------------------------------------
// Convenience: write a complete snapshot set for one tick
// ---------------------------------------------------------------------------

/// Write a full set of diagnostic images for the current simulation state.
///
/// This produces:
///   - XZ cross-sections (vertical slices at y=grid.y/2) for the configured
///     species (default: oxidant (1), reductant (2), carbon (3), organic (4)).
///     These show the vertical gradient structure that emerges from
///     top-sourced oxidant/carbon vs. bottom-sourced reductant.
///
///   - XY cross-sections (horizontal slices) for carbon (species 3) at the
///     configured depth fractions.
///     These reveal lateral heterogeneity and niche partitioning at each
///     depth layer.
///
///   - Cell density XZ heatmap showing where cells are physically located
///     in the vertical column (if enabled).
///
///   - Ancestry map showing starter-type origin of cells (if enabled).
///
/// All files are written to `out.output_dir`, which is created if it does not exist.
/// Filenames include the tick number for easy time-series assembly.
pub fn write_all_snapshots(
    field: &Field,
    _light: &LightField,
    cells: &HashMap<[u16; 3], usize>,
    cell_vec: &[marl_cell::cell::CellState],
    tick: u64,
    out: &OutputConfig,
    _sim: &SimulationConfig,
) -> std::io::Result<()> {
    let out_dir = &out.output_dir;
    let grid = field.grid();
    // Ensure the output directory exists.
    fs::create_dir_all(out_dir)?;

    // --- XZ cross-sections for key chemical species ---
    for &species in &out.xz_snapshot_species {
        if species < S_EXT {
            write_xz_cross_section(field, species, tick, out_dir)?;
        }
    }

    // --- XY cross-sections at configured depths ---
    for &frac in &out.xy_slice_depths_frac {
        let z = (frac.clamp(0.0, 1.0) * (grid.z - 1) as f32).round() as usize;
        write_xy_cross_section(field, 3, z, tick, out_dir)?;
    }

    // --- Cell density heatmap ---
    if out.write_density_map {
        write_cell_density_xz(grid, cells, tick, out_dir)?;
    }

    // --- Ancestry map (red=photo, green=chemo, blue=anaerobe) ---
    if out.write_ancestry_map {
        write_ancestry_xz(grid, cell_vec, cells, tick, out_dir)?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Ancestry map — color cells by their original starter metabolism
// ---------------------------------------------------------------------------

/// Write an XZ cross-section where cell color indicates ancestry:
///   Red   = phototroph descendants (starter_type 0)
///   Green = chemolithotroph descendants (starter_type 1)
///   Blue  = anaerobe descendants (starter_type 2)
///
/// Brightness scales with energy — bright = thriving, dim = starving.
/// Empty voxels are black. This immediately shows which original
/// metabolism colonized each depth zone and whether the middle was
/// invaded from above or below.
pub fn write_ancestry_xz(
    grid: GridDims,
    cells: &[marl_cell::cell::CellState],
    cell_map: &HashMap<[u16; 3], usize>,
    tick: u64,
    out_dir: &str,
) -> std::io::Result<()> {
    let y_mid = grid.y / 2;
    let width = grid.x;
    let height = grid.z;

    let mut pixels = vec![[0u8; 3]; width * height];

    for z in 0..grid.z {
        for x in 0..grid.x {
            let pos = [x as u16, y_mid as u16, z as u16];
            if let Some(&idx) = cell_map.get(&pos) {
                let cell = &cells[idx];
                // Brightness from energy: 0.0→dim(40), 1.0+→full(255)
                let brightness = ((cell.internal[0] / 0.5).clamp(0.0, 1.0) * 215.0 + 40.0) as u8;
                let color = match cell.starter_type {
                    0 => [brightness, brightness / 5, brightness / 5], // red = phototroph
                    1 => [brightness / 5, brightness, brightness / 5], // green = chemolithotroph
                    2 => [brightness / 5, brightness / 5, brightness], // blue = anaerobe
                    _ => [brightness; 3],                              // white = unknown
                };
                pixels[z * width + x] = color;
            }
            // else: black (empty voxel), already [0,0,0]
        }
    }

    let filename = format!("ancestry_xz_t{}.ppm", tick);
    let path = std::path::Path::new(out_dir).join(filename);
    write_ppm(&path, width, height, &pixels)
}
