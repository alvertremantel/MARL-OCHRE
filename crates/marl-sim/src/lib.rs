pub mod seeding;
pub mod spatial;
pub mod starter_metabolisms;
pub mod stats;

use marl_cell::cell::*;
use marl_config::stoich::{
    StoichEnforcement, StoichEventKind, StoichRecord, StoichReservoir, StoichRunLedger,
    StoichStage, StoichTickLedger, external_delta, internal_pool_delta,
};
use marl_config::*;
use marl_field::field::{Field, validate_diffusion_config};
use marl_field::light::LightField;
#[cfg(feature = "gpu")]
use marl_gpu::GpuFieldDiffuser;
use marl_output::binary_dump;
use marl_output::data::DataLogger;
use marl_output::snapshot;

use crate::seeding::{init_field_boundaries_with_stoich, seed_cells};
use crate::spatial::{
    apply_deltas_to_neighbors_with_stoich, find_empty_neighbor_avoiding, nearby_cell_indices,
    read_neighbor_environment,
};
use crate::starter_metabolisms::{make_anaerobe, make_chemolithotroph, make_phototroph};
use crate::stats::{print_stats, print_z_profile};

use marl_cell::hgt::transfer_reaction;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha12Rng;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

const RUN_RNG_ALGORITHM: &str = "chacha12";

fn cadence_due(tick: u32, max_ticks: u32, interval: u32) -> bool {
    tick + 1 == max_ticks || (interval > 0 && tick.is_multiple_of(interval))
}

fn validate_run_config(cfg: &Config) -> Result<(), String> {
    cfg.grid.validate()?;
    validate_diffusion_config(&cfg.simulation)?;
    cfg.simulation.validate_chemistry()?;
    validate_mutation_config(&cfg.simulation)?;
    validate_hgt_config(&cfg.simulation)
}

fn validate_mutation_config(sim: &SimulationConfig) -> Result<(), String> {
    if !sim.mutation_stddev.is_finite() || sim.mutation_stddev <= 0.0 {
        return Err(format!(
            "mutation_stddev must be finite and > 0.0, got {}",
            sim.mutation_stddev
        ));
    }
    Ok(())
}

fn validate_hgt_config(sim: &SimulationConfig) -> Result<(), String> {
    if !sim.hgt_base_rate.is_finite() {
        return Err(format!(
            "hgt_base_rate must be finite, got {}",
            sim.hgt_base_rate
        ));
    }
    if sim.hgt_base_rate < 0.0 {
        return Err(format!(
            "hgt_base_rate must be nonnegative, got {}",
            sim.hgt_base_rate
        ));
    }
    if sim.hgt_enabled {
        if sim.hgt_interval == 0 {
            return Err("hgt_interval must be > 0 when HGT is enabled".to_string());
        }
        if sim.hgt_radius == 0 {
            return Err("hgt_radius must be > 0 when HGT is enabled".to_string());
        }
        if sim.hgt_max_events_per_tick == 0 {
            return Err("hgt_max_events_per_tick must be > 0 when HGT is enabled".to_string());
        }
    }
    Ok(())
}

fn hgt_due(tick: u32, sim: &SimulationConfig) -> bool {
    sim.hgt_enabled
        && sim.hgt_interval > 0
        && sim.hgt_base_rate > 0.0
        && sim.hgt_max_events_per_tick > 0
        && tick.is_multiple_of(sim.hgt_interval)
}

fn hgt_accept_probability(base_rate: f32, propensity: f32) -> f32 {
    if !base_rate.is_finite() || !propensity.is_finite() || base_rate <= 0.0 || propensity <= 0.0 {
        return 0.0;
    }
    (base_rate * propensity).min(1.0)
}

fn run_hgt_phase(
    tick: u32,
    grid: GridDims,
    cells: &mut [CellState],
    cell_map: &HashMap<[u16; 3], usize>,
    sim: &SimulationConfig,
    rng: &mut impl Rng,
) -> u64 {
    if !hgt_due(tick, sim) || cells.len() < 2 {
        return 0;
    }

    let mut received = vec![false; cells.len()];
    let donor_rulesets: Vec<_> = cells.iter().map(|cell| cell.ruleset.clone()).collect();
    let mut events = 0u64;
    for recipient_idx in 0..cells.len() {
        if events as usize >= sim.hgt_max_events_per_tick {
            break;
        }
        if received[recipient_idx] {
            continue;
        }

        let accept_probability = hgt_accept_probability(
            sim.hgt_base_rate,
            cells[recipient_idx].ruleset.hgt_propensity,
        );
        if accept_probability <= 0.0 || rng.random::<f32>() >= accept_probability {
            continue;
        }

        let recipient_pos = cells[recipient_idx].pos;
        let donors = nearby_cell_indices(grid, recipient_pos, cell_map, sim.hgt_radius);
        if donors.is_empty() {
            continue;
        }
        let donor_idx = donors[rng.random_range(0..donors.len())];
        if donor_idx == recipient_idx {
            continue;
        }

        if transfer_reaction(
            &donor_rulesets[donor_idx],
            &mut cells[recipient_idx].ruleset,
            rng,
        ) {
            received[recipient_idx] = true;
            events += 1;
        }
    }

    events
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
pub fn run(cfg: Config, _use_gpu_diffusion: bool) -> Result<(), String> {
    validate_run_config(&cfg)?;
    let grid = cfg.grid;

    let (mut rng, rng_seed) = make_run_rng(&cfg.simulation);

    let mut field = Field::new(grid);
    let mut light = LightField::new(grid);

    let writes_stoich_v1 = cfg.output.write_stoich_summary || cfg.output.write_stoich_tick_log;
    let writes_stoich_v2 = cfg.simulation.stoich_enforcement.is_enabled()
        || cfg.output.write_stoich_v2_summary
        || cfg.output.write_stoich_v2_events;
    let writes_stoich = writes_stoich_v1 || writes_stoich_v2;
    let effective_stoich_enforcement =
        if writes_stoich_v2 && cfg.simulation.stoich_enforcement == StoichEnforcement::Off {
            StoichEnforcement::Audit
        } else {
            cfg.simulation.stoich_enforcement
        };

    // Create the data logger for optional CSV diagnostics and summaries.
    let mut logger = DataLogger::new(
        grid,
        &cfg.output.output_dir,
        cfg.output.write_tick_log,
        cfg.output.write_stoich_tick_log,
        cfg.output.write_stoich_v2_events,
    )
    .expect("Failed to create data logger / output directory");
    let writes_binary = cfg.output.write_binary_field
        || cfg.output.write_binary_cells
        || cfg.output.ruleset_output_mode.is_enabled();
    if writes_binary {
        binary_dump::write_run_meta(grid, &cfg.output).expect("Failed to write run metadata");
    }

    // Start with empty field — let boundary sources build gradients organically.
    // Pre-load only a thin boundary layer so initial cells can bootstrap.
    let mut initial_stoich = writes_stoich_v2.then(StoichTickLedger::default);
    init_field_boundaries_with_stoich(
        &mut field,
        &cfg.simulation,
        initial_stoich.as_mut(),
        cfg.output.write_stoich_v2_events,
    );

    // Cell storage: Vec for contiguous iteration + HashMap for O(1) spatial lookup.
    // The map stores position -> index into the Vec.
    let mut cells: Vec<CellState> = Vec::new();
    let mut cell_map: HashMap<[u16; 3], usize> = HashMap::new();

    // Seed three metabolisms — small populations at appropriate depths.
    // z_scale maps the "canonical" 200-layer depth to our actual grid depth,
    // so metabolisms land at the right relative positions regardless of grid depth.
    let z_scale = grid.z as f32 / 200.0;
    let sim = &cfg.simulation;

    // Phototrophs: surface
    let photo_lo = (sim.phototroph_z_lo * z_scale) as u16;
    let photo_hi = (sim.phototroph_z_hi * z_scale).max(photo_lo as f32 + 1.0) as u16;
    seed_cells(
        grid,
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
        grid,
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
        grid,
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
        grid.x,
        grid.y,
        grid.z,
        grid.voxel_count().unwrap_or(0) as f64 / 1e6,
        S_EXT,
        M_INT
    );
    println!("Seeded {} cells (photo/chemo/anaerobe)", cells.len());
    println!("Output: {}", cfg.output.output_dir);
    println!("RNG: {RUN_RNG_ALGORITHM}, seed: {rng_seed}");
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
        if grid != GridDims::default() {
            eprintln!(
                "Warning: GPU diffusion currently requires the default grid; falling back to CPU diffusion"
            );
            None
        } else {
            match GpuFieldDiffuser::new() {
                Ok(diffuser) => Some(diffuser),
                Err(e) => {
                    eprintln!(
                        "Warning: GPU diffusion unavailable ({e}); falling back to CPU diffusion"
                    );
                    None
                }
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
    let mut occupancy = vec![false; grid.voxel_count().expect("valid grid voxel count")];
    let mut stoich_run = StoichRunLedger::default();

    for tick in 0..cfg.output.max_ticks {
        tick_divisions = 0;
        tick_deaths = 0;

        // === STEP 1: Boundary sources ===
        // Inject oxidant + carbon at top, reductant at bottom — the only
        // external energy inputs. Everything else is recycled by cells.
        let mut tick_stoich = if tick == 0 {
            initial_stoich.take().unwrap_or_default()
        } else {
            StoichTickLedger::default()
        };

        if writes_stoich_v2 {
            field.apply_boundary_sources_with_stoich(
                sim,
                Some(&mut tick_stoich),
                cfg.output.write_stoich_v2_events,
            );
        } else {
            field.apply_boundary_sources(sim);
        }

        // === STEP 2: Diffusion ===
        // Sub-stepped forward Euler on 3D Laplacian. CFL-stable because
        // D * dt_sub < 1/6 for all species (see config.rs).
        // Build occupancy grid so the diffusion solver knows where cells are.
        // Occupied voxels are fully excluded from diffusion.
        occupancy.fill(false);
        for pos in cell_map.keys() {
            let idx =
                pos[2] as usize * grid.y * grid.x + pos[1] as usize * grid.x + pos[0] as usize;
            occupancy[idx] = true;
        }
        let field_totals_before_diffusion = writes_stoich_v2.then(|| field.species_totals());
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
        if let Some(before) = field_totals_before_diffusion {
            record_diffusion_losses(
                before,
                field.species_totals(),
                &mut tick_stoich,
                cfg.output.write_stoich_v2_events,
            );
        }

        // === STEP 3: Light attenuation ===
        // Beer-Lambert top-down sweep. Light enters at z=0, attenuated by
        // cells and chemical absorbers. Stored per-voxel so photosynthesis
        // reactions can reference it as a catalyst.
        light.update(&field, &cell_map, sim);
        if writes_stoich_v2 {
            let total_light = light.data.iter().copied().sum::<f32>();
            tick_stoich.record(
                StoichRecord::new(
                    StoichStage::Light,
                    StoichEventKind::LightAvailability,
                    total_light,
                ),
                cfg.output.write_stoich_v2_events,
            );
        }

        // === STEP 4: Cell update pass ===
        // Each cell runs the 5-phase tick (receptor, transport, reactions,
        // effector, fate) and returns field deltas + a fate event.
        let mut events: Vec<(usize, CellEvent)> = Vec::with_capacity(cells.len());
        for (i, cell) in cells.iter_mut().enumerate() {
            let p = cell.pos;
            // Cells sense the extracellular medium via empty neighbors,
            // not their own voxel (which is excluded from diffusion).
            let ext = read_neighbor_environment(p, &field, &cell_map);
            let l = light.get(p[0] as usize, p[1] as usize, p[2] as usize);

            let (deltas, event) = if writes_stoich {
                cell.tick_with_stoich(
                    &ext,
                    l,
                    sim,
                    Some(&mut tick_stoich),
                    writes_stoich_v2,
                    cfg.output.write_stoich_v2_events,
                )
            } else {
                cell.tick(&ext, l, sim)
            };
            // Secretion/consumption distributed to neighboring empty voxels
            if writes_stoich_v2 {
                apply_deltas_to_neighbors_with_stoich(
                    p,
                    &mut field,
                    &cell_map,
                    &deltas,
                    Some(&mut tick_stoich),
                    cfg.output.write_stoich_v2_events,
                    cell.lineage_id,
                );
            } else {
                apply_deltas_to_neighbors_with_stoich(
                    p, &mut field, &cell_map, &deltas, None, false, 0,
                );
            }
            events.push((i, event));
        }

        // === STEP 5: Process fate events ===
        let mut births: Vec<CellState> = Vec::new();
        let mut deaths: Vec<usize> = Vec::new();
        let mut reserved_birth_positions: HashSet<[u16; 3]> = HashSet::new();

        for (i, event) in &events {
            match event {
                CellEvent::Division => {
                    let parent = &cells[*i];
                    let parent_lineage = parent.lineage_id;
                    if let Some(daughter_pos) = find_empty_neighbor_avoiding(
                        grid,
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
                        if writes_stoich_v2 {
                            tick_stoich.record(
                                StoichRecord::new(
                                    StoichStage::Division,
                                    StoichEventKind::DivisionSplit,
                                    0.0,
                                )
                                .actor(parent_lineage),
                                cfg.output.write_stoich_v2_events,
                            );
                        }
                        daughter.ruleset.mutate(&mut rng, sim);
                        births.push(daughter);
                        tick_divisions += 1;
                    }
                }
                CellEvent::Death => {
                    if writes_stoich_v2 {
                        let cell = &cells[*i];
                        let removed = internal_pool_delta(&cell.internal, -1.0);
                        tick_stoich.record(
                            StoichRecord::new(
                                StoichStage::Death,
                                StoichEventKind::DeathRemoval,
                                removed.total_abs_sum(),
                            )
                            .model_delta(removed)
                            .balancing_reservoir(StoichReservoir::RemovedBiomass)
                            .actor(cell.lineage_id),
                            cfg.output.write_stoich_v2_events,
                        );
                    }
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

        let tick_hgt_events = run_hgt_phase(tick, grid, &mut cells, &cell_map, sim, &mut rng);

        // === STEP 6: Data logging and periodic output ===
        if writes_stoich {
            stoich_run.add_tick(&tick_stoich);
            if let Err(e) = logger.log_stoich_tick(tick as u64, &tick_stoich) {
                eprintln!("Warning: failed to log stoichiometry tick {}: {}", tick, e);
            }
            if let Err(e) = logger.log_stoich_v2_events(tick as u64, &tick_stoich) {
                eprintln!(
                    "Warning: failed to log stoichiometry v2 events at tick {}: {}",
                    tick, e
                );
            }
        }

        // Optionally log every tick to ticks.csv (lightweight — just one CSV row)
        if let Err(e) = logger.log_tick(
            tick as u64,
            &cells,
            tick_divisions,
            tick_deaths,
            tick_hgt_events,
        ) {
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
                && let Err(e) = binary_dump::write_cell_dump(grid, &cells, t, &cfg.output)
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
                && let Err(e) = binary_dump::write_ruleset_layer_dump(grid, &cells, t, &cfg.output)
            {
                eprintln!(
                    "Warning: failed to write binary ruleset layer snapshot at tick {}: {}",
                    tick, e
                );
            }
            if mode.writes_full_dump()
                && let Err(e) = binary_dump::write_ruleset_full_dump(grid, &cells, t, &cfg.output)
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
        rng_seed,
        RUN_RNG_ALGORITHM,
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

    if cfg.output.write_stoich_v2_summary || cfg.simulation.stoich_enforcement.is_enabled() {
        if let Err(e) = logger.write_stoich_v2_summary(
            cfg.output.max_ticks,
            effective_stoich_enforcement,
            &stoich_run,
        ) {
            eprintln!("Warning: failed to write stoichiometry v2 summary: {}", e);
        } else {
            println!(
                "Stoichiometry v2 summary written to {}/stoich_v2_summary.json",
                cfg.output.output_dir
            );
        }
    }

    // Write ancestry-colored XZ cross-section (red=photo, green=chemo, blue=anaerobe)
    if cfg.output.write_ancestry_map
        && let Err(e) = snapshot::write_ancestry_xz(
            grid,
            &cells,
            &cell_map,
            cfg.output.max_ticks as u64,
            &cfg.output.output_dir,
        )
    {
        eprintln!("Warning: failed to write ancestry map: {}", e);
    }

    Ok(())
}

fn make_run_rng(sim: &SimulationConfig) -> (ChaCha12Rng, u64) {
    let seed = sim.rng_seed.unwrap_or_else(|| {
        let mut entropy = rand::rng();
        entropy.random::<u64>()
    });
    (ChaCha12Rng::seed_from_u64(seed), seed)
}

fn record_diffusion_losses(
    before: [f32; S_EXT],
    after: [f32; S_EXT],
    ledger: &mut StoichTickLedger,
    keep_events: bool,
) {
    for species in 0..S_EXT {
        let loss = (before[species] - after[species]).max(0.0);
        if loss <= f32::EPSILON {
            continue;
        }
        ledger.record(
            StoichRecord::new(
                StoichStage::DiffusionDecay,
                StoichEventKind::DiffusionDecay,
                loss,
            )
            .model_delta(external_delta(species, -loss))
            .balancing_reservoir(StoichReservoir::AbioticDecaySink)
            .species(species),
            keep_events,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::{
        cadence_due, hgt_accept_probability, hgt_due, make_run_rng, run, run_hgt_phase,
        validate_run_config,
    };
    use crate::starter_metabolisms::{make_chemolithotroph, make_phototroph};
    use marl_config::stoich::StoichEnforcement;
    use marl_config::{Config, GridDims, SimulationConfig};
    use rand::Rng;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rand_chacha::ChaCha12Rng;
    use std::collections::HashMap;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn test_output_dir(name: &str) -> String {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let unique = format!("{}_{}_{}", name, std::process::id(), nanos);
        let dir = std::env::temp_dir().join(unique);
        let _ = fs::remove_dir_all(&dir);
        dir.to_string_lossy().into_owned()
    }

    fn replay_smoke_config(output_dir: String, rng_seed: Option<u64>) -> Config {
        let mut cfg = Config::default();
        cfg.grid = GridDims { x: 8, y: 8, z: 4 };
        cfg.simulation.rng_seed = rng_seed;
        cfg.simulation.diffusion_substeps = 1;
        cfg.simulation.seed_margin = 1;
        cfg.simulation.phototroph_z_lo = 0.0;
        cfg.simulation.phototroph_z_hi = 1.0;
        cfg.simulation.chemolithotroph_z_lo = 1.0;
        cfg.simulation.chemolithotroph_z_hi = 3.0;
        cfg.simulation.anaerobe_z_lo = 2.0;
        cfg.simulation.anaerobe_z_hi = 3.0;
        cfg.simulation.boundary_prime_layers = 1;
        cfg.output.max_ticks = 3;
        cfg.output.stats_interval = 3;
        cfg.output.snapshot_interval = 3;
        cfg.output.ruleset_interval = 3;
        cfg.output.image_interval = 0;
        cfg.output.seed_count = 1;
        cfg.output.output_dir = output_dir;
        cfg.output.write_tick_log = true;
        cfg.output.write_binary_field = false;
        cfg.output.write_binary_cells = true;
        cfg
    }

    fn summary_rng_seed(summary: &str) -> u64 {
        summary
            .lines()
            .find_map(|line| line.strip_prefix("- RNG seed: "))
            .expect("summary includes RNG seed")
            .parse()
            .expect("summary RNG seed is u64")
    }

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
    fn configured_rng_seed_reproduces_random_sequence() {
        let seeded = SimulationConfig {
            rng_seed: Some(12345),
            ..SimulationConfig::default()
        };
        let other_seed = SimulationConfig {
            rng_seed: Some(54321),
            ..SimulationConfig::default()
        };

        let (mut rng_a, seed_a) = make_run_rng(&seeded);
        let (mut rng_b, seed_b) = make_run_rng(&seeded);
        let (mut rng_c, seed_c) = make_run_rng(&other_seed);
        let sequence_a = (0..4).map(|_| rng_a.random::<u64>()).collect::<Vec<_>>();
        let sequence_b = (0..4).map(|_| rng_b.random::<u64>()).collect::<Vec<_>>();
        let sequence_c = (0..4).map(|_| rng_c.random::<u64>()).collect::<Vec<_>>();

        assert_eq!(seed_a, 12345);
        assert_eq!(seed_b, 12345);
        assert_eq!(seed_c, 54321);
        assert_eq!(sequence_a, sequence_b);
        assert_ne!(sequence_a, sequence_c);
    }

    #[test]
    fn entropy_rng_seed_is_resolved_and_replayable() {
        let unseeded = SimulationConfig {
            rng_seed: None,
            ..SimulationConfig::default()
        };

        let (mut rng, resolved_seed) = make_run_rng(&unseeded);
        let mut replay = ChaCha12Rng::seed_from_u64(resolved_seed);

        let sequence = (0..4).map(|_| rng.random::<u64>()).collect::<Vec<_>>();
        let replay_sequence = (0..4).map(|_| replay.random::<u64>()).collect::<Vec<_>>();

        assert_eq!(sequence, replay_sequence);
    }

    #[test]
    fn entropy_run_can_be_replayed_from_summary_seed() {
        let entropy_dir = test_output_dir("marl_sim_entropy_replay_entropy");
        let replay_dir = test_output_dir("marl_sim_entropy_replay_seeded");

        run(replay_smoke_config(entropy_dir.clone(), None), false).unwrap();
        let summary = fs::read_to_string(PathBuf::from(&entropy_dir).join("summary.md")).unwrap();
        let seed = summary_rng_seed(&summary);
        assert!(summary.contains("- RNG algorithm: chacha12"));
        assert!(summary.contains(&format!("- Replay TOML: `rng_seed = {seed}`")));

        run(replay_smoke_config(replay_dir.clone(), Some(seed)), false).unwrap();

        let entropy_path = PathBuf::from(&entropy_dir);
        let replay_path = PathBuf::from(&replay_dir);
        assert_eq!(
            fs::read(entropy_path.join("ticks.csv")).unwrap(),
            fs::read(replay_path.join("ticks.csv")).unwrap()
        );
        assert_eq!(
            fs::read(entropy_path.join("tick_2.cells.bin.zst")).unwrap(),
            fs::read(replay_path.join("tick_2.cells.bin.zst")).unwrap()
        );

        let _ = fs::remove_dir_all(&entropy_dir);
        let _ = fs::remove_dir_all(&replay_dir);
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

    #[test]
    fn invalid_hgt_base_rate_is_rejected_before_run_loop() {
        for bad_rate in [f32::NAN, f32::INFINITY, -0.1] {
            let cfg = Config {
                simulation: SimulationConfig {
                    hgt_base_rate: bad_rate,
                    ..SimulationConfig::default()
                },
                ..Config::default()
            };

            let err = validate_run_config(&cfg).unwrap_err();
            assert!(err.contains("hgt_base_rate"));
        }
    }

    #[test]
    fn invalid_mutation_stddev_is_rejected_before_run_loop() {
        for mutation_stddev in [f32::NAN, f32::INFINITY, 0.0, -0.1] {
            let cfg = Config {
                simulation: SimulationConfig {
                    mutation_stddev,
                    ..SimulationConfig::default()
                },
                ..Config::default()
            };

            let err = validate_run_config(&cfg).unwrap_err();
            assert!(err.contains("mutation_stddev"));
        }
    }

    #[test]
    fn invalid_enabled_hgt_config_is_rejected_before_run_loop() {
        let invalid = [
            SimulationConfig {
                hgt_enabled: true,
                hgt_interval: 0,
                ..SimulationConfig::default()
            },
            SimulationConfig {
                hgt_enabled: true,
                hgt_radius: 0,
                ..SimulationConfig::default()
            },
            SimulationConfig {
                hgt_enabled: true,
                hgt_max_events_per_tick: 0,
                ..SimulationConfig::default()
            },
        ];

        for simulation in invalid {
            let cfg = Config {
                simulation,
                ..Config::default()
            };
            assert!(validate_run_config(&cfg).is_err());
        }
    }

    #[test]
    fn hgt_probability_treats_nan_and_negative_values_as_zero() {
        assert_eq!(hgt_accept_probability(f32::NAN, 1.0), 0.0);
        assert_eq!(hgt_accept_probability(1.0, f32::NAN), 0.0);
        assert_eq!(hgt_accept_probability(-0.1, 1.0), 0.0);
        assert_eq!(hgt_accept_probability(1.0, -0.1), 0.0);
        assert_eq!(hgt_accept_probability(0.75, 4.0), 1.0);
    }

    #[test]
    fn hgt_cadence_includes_tick_zero_when_enabled() {
        let sim = SimulationConfig {
            hgt_enabled: true,
            hgt_interval: 10,
            hgt_base_rate: 0.02,
            hgt_radius: 1,
            hgt_max_events_per_tick: 1,
            ..SimulationConfig::default()
        };

        assert!(hgt_due(0, &sim));
        assert!(!hgt_due(1, &sim));
        assert!(hgt_due(10, &sim));
    }

    #[test]
    fn disabled_hgt_phase_is_inactive() {
        let grid = GridDims { x: 4, y: 4, z: 4 };
        let mut recipient = make_phototroph([1, 1, 1], 1);
        for reaction in &mut recipient.ruleset.reactions {
            reaction.v_max = 0.0;
        }
        recipient.ruleset.hgt_propensity = 1.0;
        let donor = make_chemolithotroph([2, 1, 1], 2);
        let mut cells = vec![recipient, donor];
        let cell_map = HashMap::from([([1, 1, 1], 0), ([2, 1, 1], 1)]);
        let sim = SimulationConfig {
            hgt_enabled: false,
            hgt_interval: 1,
            hgt_base_rate: 1.0,
            hgt_radius: 1,
            hgt_max_events_per_tick: 1,
            ..SimulationConfig::default()
        };
        let mut rng = StdRng::seed_from_u64(21);

        let events = run_hgt_phase(0, grid, &mut cells, &cell_map, &sim, &mut rng);

        assert_eq!(events, 0);
        assert!(
            cells[0]
                .ruleset
                .reactions
                .iter()
                .all(|reaction| reaction.v_max == 0.0)
        );
    }

    #[test]
    fn nonfinite_or_negative_hgt_propensity_does_not_transfer() {
        for propensity in [f32::NAN, f32::INFINITY, -1.0] {
            let grid = GridDims { x: 4, y: 4, z: 4 };
            let mut recipient = make_phototroph([1, 1, 1], 1);
            for reaction in &mut recipient.ruleset.reactions {
                reaction.v_max = 0.0;
            }
            recipient.ruleset.hgt_propensity = propensity;
            let donor = make_chemolithotroph([2, 1, 1], 2);
            let mut cells = vec![recipient, donor];
            let cell_map = HashMap::from([([1, 1, 1], 0), ([2, 1, 1], 1)]);
            let sim = SimulationConfig {
                hgt_enabled: true,
                hgt_interval: 1,
                hgt_base_rate: 1.0,
                hgt_radius: 1,
                hgt_max_events_per_tick: 1,
                ..SimulationConfig::default()
            };
            let mut rng = StdRng::seed_from_u64(23);

            let events = run_hgt_phase(0, grid, &mut cells, &cell_map, &sim, &mut rng);

            assert_eq!(events, 0);
            assert!(
                cells[0]
                    .ruleset
                    .reactions
                    .iter()
                    .all(|reaction| reaction.v_max == 0.0)
            );
        }
    }

    #[test]
    fn enabled_hgt_phase_transfers_local_reaction() {
        let grid = GridDims { x: 4, y: 4, z: 4 };
        let mut recipient = make_phototroph([1, 1, 1], 1);
        for reaction in &mut recipient.ruleset.reactions {
            reaction.v_max = 0.0;
        }
        recipient.ruleset.hgt_propensity = 1.0;
        let donor = make_chemolithotroph([2, 1, 1], 2);
        let donor_reactions = donor.ruleset.reactions.clone();
        let mut cells = vec![recipient, donor];
        let cell_map = HashMap::from([([1, 1, 1], 0), ([2, 1, 1], 1)]);
        let sim = SimulationConfig {
            hgt_enabled: true,
            hgt_interval: 1,
            hgt_base_rate: 1.0,
            hgt_radius: 1,
            hgt_max_events_per_tick: 1,
            ..SimulationConfig::default()
        };
        let mut rng = StdRng::seed_from_u64(22);

        let events = run_hgt_phase(0, grid, &mut cells, &cell_map, &sim, &mut rng);

        assert_eq!(events, 1);
        assert!(cells[0].ruleset.reactions.iter().any(|recipient_reaction| {
            recipient_reaction.v_max != 0.0
                && donor_reactions.iter().any(|donor_reaction| {
                    donor_reaction.substrate == recipient_reaction.substrate
                        && donor_reaction.product == recipient_reaction.product
                        && donor_reaction.catalyst == recipient_reaction.catalyst
                        && donor_reaction.cofactor == recipient_reaction.cofactor
                        && donor_reaction.v_max == recipient_reaction.v_max
                })
        }));
    }

    #[test]
    fn stoich_output_flags_write_expected_files() {
        let out_dir = test_output_dir("marl_sim_stoich_output_test");
        let mut cfg = Config::default();
        cfg.output.output_dir = out_dir.clone();
        cfg.output.max_ticks = 1;
        cfg.output.stats_interval = 0;
        cfg.output.snapshot_interval = 0;
        cfg.output.image_interval = 0;
        cfg.output.seed_count = 1;
        cfg.output.write_binary_field = false;
        cfg.output.write_binary_cells = false;
        cfg.output.write_stoich_summary = true;
        cfg.output.write_stoich_tick_log = true;
        cfg.output.write_ancestry_map = false;
        cfg.output.write_density_map = false;

        run(cfg, false).unwrap();

        let dir = PathBuf::from(&out_dir);
        assert!(dir.join("stoich_summary.json").exists());
        assert!(dir.join("stoich_ticks.csv").exists());

        let summary = fs::read_to_string(dir.join("stoich_summary.json")).unwrap();
        assert!(summary.contains("\"gross_total_abs_sum\""));
        let summary_json: serde_json::Value = serde_json::from_str(&summary).unwrap();
        assert!(summary_json["ledger"]["reservoir_deltas"].is_null());
        assert!(summary_json["ledger"]["stage_summaries"].is_null());
        assert!(summary_json["ledger"]["strict_rejection_count"].is_null());
        let ticks = fs::read_to_string(dir.join("stoich_ticks.csv")).unwrap();
        assert!(ticks.contains("gross_material_abs,gross_total_abs"));

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn audit_stoich_v2_events_include_legacy_reactions() {
        let out_dir = test_output_dir("marl_sim_audit_stoich_v2_event_test");
        let mut cfg = Config::default();
        cfg.simulation.stoich_enforcement = StoichEnforcement::Audit;
        cfg.output.output_dir = out_dir.clone();
        cfg.output.max_ticks = 1;
        cfg.output.stats_interval = 0;
        cfg.output.snapshot_interval = 0;
        cfg.output.image_interval = 0;
        cfg.output.seed_count = 1;
        cfg.output.write_binary_field = false;
        cfg.output.write_binary_cells = false;
        cfg.output.write_stoich_v2_summary = true;
        cfg.output.write_stoich_v2_events = true;
        cfg.output.write_ancestry_map = false;
        cfg.output.write_density_map = false;

        run(cfg, false).unwrap();

        let dir = PathBuf::from(&out_dir);
        let events = fs::read_to_string(dir.join("stoich_v2_events.csv")).unwrap();
        assert!(events.contains(",reactions,legacy_reaction,"));

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn strict_stoich_run_writes_v2_outputs() {
        let out_dir = test_output_dir("marl_sim_strict_stoich_output_test");
        let mut cfg = Config::default();
        cfg.simulation.stoich_enforcement = StoichEnforcement::Strict;
        cfg.output.output_dir = out_dir.clone();
        cfg.output.max_ticks = 1;
        cfg.output.stats_interval = 0;
        cfg.output.snapshot_interval = 0;
        cfg.output.image_interval = 0;
        cfg.output.seed_count = 1;
        cfg.output.write_binary_field = false;
        cfg.output.write_binary_cells = false;
        cfg.output.write_stoich_v2_summary = false;
        cfg.output.write_stoich_v2_events = true;
        cfg.output.write_ancestry_map = false;
        cfg.output.write_density_map = false;

        run(cfg, false).unwrap();

        let dir = PathBuf::from(&out_dir);
        let summary = fs::read_to_string(dir.join("stoich_v2_summary.json")).unwrap();
        assert!(summary.contains("\"enforcement\": \"strict\""));
        assert!(summary.contains("\"schema_version\": 2"));
        let events = fs::read_to_string(dir.join("stoich_v2_events.csv")).unwrap();
        assert!(events.contains("tick,stage,kind"));

        let _ = fs::remove_dir_all(&dir);
    }
}
