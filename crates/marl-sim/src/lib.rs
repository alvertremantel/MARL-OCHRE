pub mod seeding;
pub mod spatial;
pub mod starter_metabolisms;
pub mod stats;

use marl_cell::cell::*;
use marl_config::stoich::{StoichRunLedger, StoichTickLedger};
use marl_config::*;
use marl_field::field::{Field, validate_diffusion_config};
use marl_field::light::LightField;
#[cfg(feature = "gpu")]
use marl_gpu::GpuFieldDiffuser;
use marl_output::binary_dump;
use marl_output::data::DataLogger;
use marl_output::snapshot;

use crate::seeding::{init_field_boundaries, seed_cells};
use crate::spatial::{
    apply_deltas_to_neighbors, find_empty_neighbor_avoiding, read_neighbor_environment,
};
use crate::starter_metabolisms::{make_anaerobe, make_chemolithotroph, make_phototroph};
use crate::stats::{print_stats, print_z_profile};

use rand::Rng;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

fn cadence_due(tick: u32, max_ticks: u32, interval: u32) -> bool {
    tick + 1 == max_ticks || (interval > 0 && tick.is_multiple_of(interval))
}

fn validate_run_config(cfg: &Config) -> Result<(), String> {
    validate_diffusion_config(&cfg.simulation)
}

/// Run the full MARL simulation tick loop.
///
/// This is the main orchestration function — all phases of the simulation
/// (boundary sources, diffusion, light, cell updates, fate processing,
/// logging, snapshots) happen here.
///
/// # Arguments
/// - `cfg` — fully-parsed simulation + output configuration
/// - `use_gpu_diffusion` — whether to attempt GPU-accelerated diffusion
pub fn run(cfg: Config, _use_gpu_diffusion: bool) {
    if let Err(e) = validate_run_config(&cfg) {
        eprintln!("Invalid simulation configuration: {e}");
        return;
    }

    let mut rng = rand::rng();

    let mut field = Field::new();
    let mut light = LightField::new();

    let writes_stoich = cfg.output.write_stoich_summary || cfg.output.write_stoich_tick_log;

    // Create the data logger for optional CSV diagnostics and summaries.
    let mut logger = DataLogger::new(
        &cfg.output.output_dir,
        cfg.output.write_tick_log,
        cfg.output.write_stoich_tick_log,
    )
    .expect("Failed to create data logger / output directory");
    let writes_binary = cfg.output.write_binary_field
        || cfg.output.write_binary_cells
        || cfg.output.ruleset_output_mode.is_enabled();
    if writes_binary {
        binary_dump::write_run_meta(&cfg.output).expect("Failed to write run metadata");
    }

    // Start with empty field — let boundary sources build gradients organically.
    // Pre-load only a thin boundary layer so initial cells can bootstrap.
    init_field_boundaries(&mut field, &cfg.simulation);

    // Cell storage: Vec for contiguous iteration + HashMap for O(1) spatial lookup.
    // The map stores position -> index into the Vec.
    let mut cells: Vec<CellState> = Vec::new();
    let mut cell_map: HashMap<[u16; 3], usize> = HashMap::new();

    // Seed three metabolisms — small populations at appropriate depths.
    // z_scale maps the "canonical" 200-layer depth to our actual grid depth,
    // so metabolisms land at the right relative positions regardless of GRID_Z.
    let z_scale = GRID_Z as f32 / 200.0;
    let sim = &cfg.simulation;

    // Phototrophs: surface
    let photo_lo = (sim.phototroph_z_lo * z_scale) as u16;
    let photo_hi = (sim.phototroph_z_hi * z_scale).max(photo_lo as f32 + 1.0) as u16;
    seed_cells(
        &mut cells,
        &mut cell_map,
        &mut rng,
        cfg.output.seed_count,
        photo_lo,
        photo_hi,
        make_phototroph,
        sim,
    );

    // Chemolithotrophs: chemocline — oxidize reductant using oxidant at the interface
    let chemo_lo = (sim.chemolithotroph_z_lo * z_scale) as u16;
    let chemo_hi = (sim.chemolithotroph_z_hi * z_scale).max(chemo_lo as f32 + 3.0) as u16;
    seed_cells(
        &mut cells,
        &mut cell_map,
        &mut rng,
        cfg.output.seed_count,
        chemo_lo,
        chemo_hi,
        make_chemolithotroph,
        sim,
    );

    // Anaerobes: deep zone — use reductant, killed by oxidant
    let ana_lo = (sim.anaerobe_z_lo * z_scale) as u16;
    let ana_hi = (sim.anaerobe_z_hi * z_scale).max(ana_lo as f32 + 3.0) as u16;
    seed_cells(
        &mut cells,
        &mut cell_map,
        &mut rng,
        cfg.output.seed_count,
        ana_lo,
        ana_hi,
        make_anaerobe,
        sim,
    );

    println!("MARL v0.3 — CPU Prototype (Winogradsky)");
    println!(
        "Grid: {}x{}x{} ({:.1}M voxels), Species: {} ext / {} int",
        GRID_X,
        GRID_Y,
        GRID_Z,
        (GRID_X * GRID_Y * GRID_Z) as f64 / 1e6,
        S_EXT,
        M_INT
    );
    println!("Seeded {} cells (photo/chemo/anaerobe)", cells.len());
    println!("Output: {}", cfg.output.output_dir);
    println!(
        "Plan: {} ticks, stats every {}, snapshots every {}, images every {}",
        cfg.output.max_ticks,
        cfg.output.stats_interval,
        cfg.output.snapshot_interval,
        cfg.output.image_interval
    );
    #[cfg(feature = "gpu")]
    println!(
        "Diffusion: {}",
        if _use_gpu_diffusion { "GPU" } else { "CPU" }
    );
    println!("---");

    #[cfg(feature = "gpu")]
    let mut gpu_diffuser = if _use_gpu_diffusion {
        match GpuFieldDiffuser::new() {
            Ok(diffuser) => Some(diffuser),
            Err(e) => {
                eprintln!(
                    "Warning: GPU diffusion unavailable ({e}); falling back to CPU diffusion"
                );
                None
            }
        }
    } else {
        None
    };

    let mut total_divisions: u64 = 0;
    let mut total_deaths: u64 = 0;
    let start = Instant::now();
    let writes_images = !cfg.output.xz_snapshot_species.is_empty()
        || !cfg.output.xy_slice_depths_frac.is_empty()
        || cfg.output.write_density_map
        || cfg.output.write_ancestry_map;

    // Track per-tick division/death counts for the data logger
    let mut tick_divisions: u64;
    let mut tick_deaths: u64;
    let mut occupancy = vec![false; GRID_X * GRID_Y * GRID_Z];
    let mut stoich_run = StoichRunLedger::default();

    for tick in 0..cfg.output.max_ticks {
        tick_divisions = 0;
        tick_deaths = 0;

        // === STEP 1: Boundary sources ===
        // Inject oxidant + carbon at top, reductant at bottom — the only
        // external energy inputs. Everything else is recycled by cells.
        field.apply_boundary_sources(sim);

        // === STEP 2: Diffusion ===
        // Sub-stepped forward Euler on 3D Laplacian. CFL-stable because
        // D * dt_sub < 1/6 for all species (see config.rs).
        // Build occupancy grid so the diffusion solver knows where cells are.
        // Occupied voxels are fully excluded from diffusion.
        occupancy.fill(false);
        for pos in cell_map.keys() {
            let idx =
                pos[2] as usize * GRID_Y * GRID_X + pos[1] as usize * GRID_X + pos[0] as usize;
            occupancy[idx] = true;
        }
        #[cfg(feature = "gpu")]
        if let Some(diffuser) = gpu_diffuser.as_mut() {
            if let Err(e) = diffuser.diffuse_tick_with_cells(&mut field, &occupancy, sim) {
                eprintln!(
                    "Warning: GPU diffusion failed at tick {tick} ({e}); falling back to CPU diffusion"
                );
                gpu_diffuser = None;
                field.diffuse_tick_with_cells(&occupancy, sim);
            }
        } else {
            field.diffuse_tick_with_cells(&occupancy, sim);
        }
        #[cfg(not(feature = "gpu"))]
        field.diffuse_tick_with_cells(&occupancy, sim);

        // === STEP 3: Light attenuation ===
        // Beer-Lambert top-down sweep. Light enters at z=0, attenuated by
        // cells and chemical absorbers. Stored per-voxel so photosynthesis
        // reactions can reference it as a catalyst.
        light.update(&field, &cell_map, sim);

        // === STEP 4: Cell update pass ===
        // Each cell runs the 5-phase tick (receptor, transport, reactions,
        // effector, fate) and returns field deltas + a fate event.
        let mut events: Vec<(usize, CellEvent)> = Vec::with_capacity(cells.len());
        let mut tick_stoich = StoichTickLedger::default();

        for (i, cell) in cells.iter_mut().enumerate() {
            let p = cell.pos;
            // Cells sense the extracellular medium via empty neighbors,
            // not their own voxel (which is excluded from diffusion).
            let ext = read_neighbor_environment(p, &field, &cell_map);
            let l = light.get(p[0] as usize, p[1] as usize, p[2] as usize);

            let (deltas, event) = if writes_stoich {
                cell.tick_with_stoich(&ext, l, sim, Some(&mut tick_stoich))
            } else {
                cell.tick(&ext, l, sim)
            };
            // Secretion/consumption distributed to neighboring empty voxels
            apply_deltas_to_neighbors(p, &mut field, &cell_map, &deltas);
            events.push((i, event));
        }

        if writes_stoich {
            stoich_run.add_tick(&tick_stoich);
            if let Err(e) = logger.log_stoich_tick(tick as u64, &tick_stoich) {
                eprintln!("Warning: failed to log stoichiometry tick {}: {}", tick, e);
            }
        }

        // === STEP 5: Process fate events ===
        let mut births: Vec<CellState> = Vec::new();
        let mut deaths: Vec<usize> = Vec::new();
        let mut reserved_birth_positions: HashSet<[u16; 3]> = HashSet::new();

        for (i, event) in &events {
            match event {
                CellEvent::Division => {
                    let parent = &cells[*i];
                    if let Some(daughter_pos) = find_empty_neighbor_avoiding(
                        parent.pos,
                        &cell_map,
                        Some(&reserved_birth_positions),
                        &mut rng,
                        sim,
                    ) {
                        reserved_birth_positions.insert(daughter_pos);
                        let mut daughter = parent.clone();
                        daughter.pos = daughter_pos;
                        daughter.age = 0;
                        daughter.prep_remaining = 0; // daughter starts fresh, not in prep
                        daughter.lineage_id = rng.random::<u64>();
                        // CRITICAL: split ALL 16 internal species, not just energy.
                        // This prevents division from being a free-energy exploit.
                        for k in 0..M_INT {
                            daughter.internal[k] *= 0.5;
                            cells[*i].internal[k] *= 0.5;
                        }
                        daughter.ruleset.mutate(&mut rng, sim);
                        births.push(daughter);
                        tick_divisions += 1;
                    }
                }
                CellEvent::Death => {
                    deaths.push(*i);
                    tick_deaths += 1;
                }
                _ => {}
            }
        }

        total_divisions += tick_divisions;
        total_deaths += tick_deaths;

        // Remove dead cells (reverse order to preserve indices during swap_remove)
        deaths.sort_unstable();
        deaths.dedup();
        for &i in deaths.iter().rev() {
            let pos = cells[i].pos;
            cell_map.remove(&pos);
            cells.swap_remove(i);
            if i < cells.len() {
                cell_map.insert(cells[i].pos, i);
            }
        }

        // Add newborns
        for cell in births {
            let pos = cell.pos;
            if let std::collections::hash_map::Entry::Vacant(entry) = cell_map.entry(pos) {
                let idx = cells.len();
                entry.insert(idx);
                cells.push(cell);
            }
        }

        // === STEP 6: Data logging and periodic output ===

        // Optionally log every tick to ticks.csv (lightweight — just one CSV row)
        if let Err(e) = logger.log_tick(tick as u64, &cells, tick_divisions, tick_deaths) {
            eprintln!("Warning: failed to log tick {}: {}", tick, e);
        }

        // Print human-readable stats to stdout at configured interval
        if cadence_due(tick, cfg.output.max_ticks, cfg.output.stats_interval) {
            print_stats(
                tick,
                &cells,
                &field,
                &light,
                total_divisions,
                total_deaths,
                &start,
            );
        }

        // Write raw binary snapshots for viewer ingestion.
        if cadence_due(tick, cfg.output.max_ticks, cfg.output.snapshot_interval) {
            let t = tick as u64;
            if cfg.output.write_binary_field
                && let Err(e) = binary_dump::write_field_dump(&field, t, &cfg.output)
            {
                eprintln!(
                    "Warning: failed to write binary field snapshot at tick {}: {}",
                    tick, e
                );
            }
            if cfg.output.write_binary_cells
                && let Err(e) = binary_dump::write_cell_dump(&cells, t, &cfg.output)
            {
                eprintln!(
                    "Warning: failed to write binary cell snapshot at tick {}: {}",
                    tick, e
                );
            }

            // Optional legacy CSV snapshots (chemistry profiles + cell dumps + reactions).
            if cfg.output.write_csv_snapshots {
                if let Err(e) = logger.snapshot_chemistry(t, &field, &light) {
                    eprintln!(
                        "Warning: failed to write chemistry snapshot at tick {}: {}",
                        tick, e
                    );
                }
                if let Err(e) = logger.snapshot_cells(t, &cells) {
                    eprintln!(
                        "Warning: failed to write cell snapshot at tick {}: {}",
                        tick, e
                    );
                }
                if let Err(e) = logger.snapshot_reactions(t, &cells) {
                    eprintln!(
                        "Warning: failed to write reaction snapshot at tick {}: {}",
                        tick, e
                    );
                }
            }
        }

        let mode = cfg.output.ruleset_output_mode;
        if mode.is_enabled() && cadence_due(tick, cfg.output.max_ticks, cfg.output.ruleset_interval)
        {
            let t = tick as u64;
            if mode.writes_layer_averages()
                && let Err(e) = binary_dump::write_ruleset_layer_dump(&cells, t, &cfg.output)
            {
                eprintln!(
                    "Warning: failed to write binary ruleset layer snapshot at tick {}: {}",
                    tick, e
                );
            }
            if mode.writes_full_dump()
                && let Err(e) = binary_dump::write_ruleset_full_dump(&cells, t, &cfg.output)
            {
                eprintln!(
                    "Warning: failed to write binary ruleset full dump at tick {}: {}",
                    tick, e
                );
            }
        }

        // Write PPM image snapshots (cross-sections, density maps)
        if writes_images
            && cadence_due(tick, cfg.output.max_ticks, cfg.output.image_interval)
            && let Err(e) = snapshot::write_all_snapshots(
                &field,
                &light,
                &cell_map,
                &cells,
                tick as u64,
                &cfg.output,
                sim,
            )
        {
            eprintln!(
                "Warning: failed to write image snapshots at tick {}: {}",
                tick, e
            );
        }
    }

    println!("\n=== FINAL Z-LAYER PROFILE ===");
    print_z_profile(&cells, &field, &light);

    let runtime = start.elapsed().as_secs_f32();
    println!(
        "\nDone. {} ticks in {:.1}s, final pop={}, div={}, death={}",
        cfg.output.max_ticks,
        runtime,
        cells.len(),
        total_divisions,
        total_deaths
    );

    // Write the post-run summary (lab notebook entry for this run)
    if let Err(e) = logger.write_summary(
        cfg.output.max_ticks,
        runtime,
        &cells,
        &field,
        &light,
        total_divisions,
        total_deaths,
        sim,
    ) {
        eprintln!("Warning: failed to write summary: {}", e);
    } else {
        println!("Summary written to {}/summary.md", cfg.output.output_dir);
    }

    // Write the reaction registry — maps IDs back to topologies for the CLI tool
    if cfg.output.write_csv_snapshots {
        if let Err(e) = logger.write_registry() {
            eprintln!("Warning: failed to write reaction registry: {}", e);
        } else {
            println!(
                "Reaction registry: {} unique topologies observed",
                logger.registry.count()
            );
        }
    }

    if cfg.output.write_stoich_summary {
        if let Err(e) = logger.write_stoich_summary(cfg.output.max_ticks, &stoich_run) {
            eprintln!("Warning: failed to write stoichiometry summary: {}", e);
        } else {
            println!(
                "Stoichiometry summary written to {}/stoich_summary.json",
                cfg.output.output_dir
            );
        }
    }

    // Write ancestry-colored XZ cross-section (red=photo, green=chemo, blue=anaerobe)
    if cfg.output.write_ancestry_map
        && let Err(e) = snapshot::write_ancestry_xz(
            &cells,
            &cell_map,
            cfg.output.max_ticks as u64,
            &cfg.output.output_dir,
        )
    {
        eprintln!("Warning: failed to write ancestry map: {}", e);
    }
}

#[cfg(test)]
mod tests {
    use super::{cadence_due, validate_run_config};
    use marl_config::{Config, SimulationConfig};

    #[test]
    fn zero_interval_disables_periodic_cadence_but_keeps_final_tick() {
        assert!(!cadence_due(0, 10, 0));
        assert!(cadence_due(9, 10, 0));
    }

    #[test]
    fn nonzero_interval_matches_periodic_or_final_tick() {
        assert!(cadence_due(0, 10, 5));
        assert!(cadence_due(5, 10, 5));
        assert!(!cadence_due(6, 10, 5));
        assert!(cadence_due(9, 10, 5));
    }

    #[test]
    fn invalid_diffusion_config_is_rejected_before_run_loop() {
        let cfg = Config {
            simulation: SimulationConfig {
                k_eps: 0.0,
                ..SimulationConfig::default()
            },
            ..Config::default()
        };

        let err = validate_run_config(&cfg).unwrap_err();
        assert!(err.contains("k_eps"));
    }
}
