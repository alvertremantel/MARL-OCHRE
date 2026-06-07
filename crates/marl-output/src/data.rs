//! Optional CSV diagnostics and end-of-run summaries.
//!
//! This module writes two kinds of output:
//!
//! 1. **ticks.csv** — one row per simulation tick with population-level summary
//!    statistics (population count, average energy, enzyme levels, spatial
//!    distribution by z-layer, etc.). Appended every tick when enabled.
//!
//! 2. **Periodic snapshots** — detailed per-cell and per-voxel chemistry dumps
//!    written at configurable intervals when enabled:
//!    - `chem_<tick>.csv`  — z-layer-averaged chemical concentrations
//!    - `cells_<tick>.csv` — full dump of every living cell's state
//!
//! All output uses plain CSV with headers. No external crate dependencies —
//! just `std::fs` and `std::io`. BufWriter is used for performance so that
//! thousands of small writes don't each trigger a syscall.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Result, Write};
use std::path::{Path, PathBuf};

use marl_cell::cell::*;
use marl_config::stoich::{
    STOICH_RESERVOIR_COUNT, STOICH_STAGE_COUNT, StoichEnforcement, StoichEvent, StoichReservoir,
    StoichRunLedger, StoichStage, StoichTickLedger,
};
use marl_config::*;
use marl_field::field::Field;
use marl_field::light::LightField;

// ============================================================================
// REACTION REGISTRY — stable IDs for every unique reaction topology
// ============================================================================

/// A reaction's identity is its topology: which species it connects.
/// Kinetic parameters (v_max, k_m, k_cat) vary continuously and don't
/// define the reaction — they're expression levels, not enzyme identity.
/// Two cells with the same ReactionTopology have the "same enzyme,"
/// possibly with different expression.
#[derive(Hash, Eq, PartialEq, Clone, Debug)]
pub struct ReactionTopology {
    pub substrate: u8,
    pub product: u8,
    pub catalyst: u8,
    pub cofactor: u8,
}

impl ReactionTopology {
    pub fn from_reaction(rxn: &Reaction) -> Self {
        Self {
            substrate: rxn.substrate,
            product: rxn.product,
            catalyst: rxn.catalyst,
            cofactor: rxn.cofactor,
        }
    }
}

/// Maps every unique reaction topology observed during the run to a
/// stable integer ID. The first time a topology is seen, it gets the
/// next available ID. IDs are permanent and never reused.
///
/// At simulation start, the ~15 reactions across the 3 starter
/// metabolisms are registered as IDs 0..N. As mutations create new
/// topologies, they get sequential IDs. Convergent evolution (two
/// independent lineages arriving at the same topology) produces the
/// same ID — that's the whole point.
pub struct ReactionRegistry {
    map: HashMap<ReactionTopology, u32>,
    next_id: u32,
}

impl ReactionRegistry {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
            next_id: 0,
        }
    }

    /// Register a topology and return its stable ID.
    /// If already seen, returns the existing ID.
    pub fn register(&mut self, topo: &ReactionTopology) -> u32 {
        if let Some(&id) = self.map.get(topo) {
            id
        } else {
            let id = self.next_id;
            self.map.insert(topo.clone(), id);
            self.next_id += 1;
            id
        }
    }

    /// How many unique topologies have been observed so far.
    pub fn count(&self) -> u32 {
        self.next_id
    }

    /// Write the full registry to `reaction_registry.csv` in the output dir.
    /// Called at end of run so the CLI tool can decode IDs back to topologies.
    pub fn write_registry(&self, output_dir: &Path) -> Result<()> {
        let path = output_dir.join("reaction_registry.csv");
        let file = File::create(&path)?;
        let mut w = BufWriter::new(file);
        writeln!(w, "reaction_id,substrate,product,catalyst,cofactor")?;
        // Sort by ID for readability
        let mut entries: Vec<_> = self.map.iter().collect();
        entries.sort_by_key(|(_, id)| **id);
        for (topo, id) in entries {
            writeln!(
                w,
                "{},{},{},{},{}",
                id, topo.substrate, topo.product, topo.catalyst, topo.cofactor
            )?;
        }
        w.flush()
    }
}

impl Default for ReactionRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Manages file handles and the output directory for all data logging.
///
/// Create one of these at simulation startup. Call `log_tick()` every tick,
/// and `snapshot_chemistry()` / `snapshot_cells()` at whatever interval you
/// want detailed dumps (e.g., every 100 ticks).
pub struct DataLogger {
    /// Path to the output directory (created on construction).
    output_dir: PathBuf,

    grid: GridDims,

    /// Buffered writer for ticks.csv — kept open for the entire run so we
    /// only pay the file-open cost once.
    ticks_writer: Option<BufWriter<File>>,

    /// Optional compact stoichiometry audit tick log.
    stoich_writer: Option<BufWriter<File>>,

    /// Optional full-system stoichiometry v2 event log.
    stoich_v2_event_writer: Option<BufWriter<File>>,

    /// Registry mapping reaction topologies to stable integer IDs.
    /// Populated incrementally as new topologies are observed.
    pub registry: ReactionRegistry,
}

impl DataLogger {
    /// Create a new DataLogger.
    ///
    /// - Creates `output_dir` if it doesn't already exist.
    /// - Optionally opens `ticks.csv` inside that directory and writes the header row.
    ///
    /// # Errors
    /// Returns `std::io::Error` if directory creation or file opening fails.
    pub fn new(
        grid: GridDims,
        output_dir: &str,
        write_tick_log: bool,
        write_stoich_tick_log: bool,
        write_stoich_v2_events: bool,
    ) -> Result<Self> {
        // Create the output directory (and any missing parents).
        let dir = PathBuf::from(output_dir);
        fs::create_dir_all(&dir)?;

        let ticks_writer = if write_tick_log {
            // Open ticks.csv and write the header line.
            let ticks_path = dir.join("ticks.csv");
            let file = File::create(&ticks_path)?;
            let mut writer = BufWriter::new(file);

            // Build header: fixed columns first, then one column per z-layer.
            // The z-layer columns let you see vertical population stratification
            // over time — important for checking whether cells self-organize into
            // depth-dependent niches (the Winogradsky hypothesis).
            write!(
                writer,
                "tick,population,avg_energy,avg_enzyme_a,avg_enzyme_b,avg_active_rxns,divisions_this_tick,deaths_this_tick,hgt_events_this_tick"
            )?;
            for z in 0..grid.z {
                write!(writer, ",z{}_cells", z)?;
            }
            writeln!(writer)?;
            writer.flush()?;
            Some(writer)
        } else {
            None
        };

        let stoich_writer = if write_stoich_tick_log {
            let stoich_path = dir.join("stoich_ticks.csv");
            let file = File::create(&stoich_path)?;
            let mut writer = BufWriter::new(file);
            writeln!(
                writer,
                "tick,reaction_count,active_flux,imbalanced_reaction_count,unknown_species_flux,carbon_to_energy_flux,reductant_to_energy_flux,gross_material_abs,gross_total_abs,delta_c,delta_h,delta_o,delta_s,delta_redox,delta_energy"
            )?;
            writer.flush()?;
            Some(writer)
        } else {
            None
        };

        let stoich_v2_event_writer = if write_stoich_v2_events {
            let events_path = dir.join("stoich_v2_events.csv");
            let file = File::create(&events_path)?;
            let mut writer = BufWriter::new(file);
            writeln!(
                writer,
                "tick,stage,kind,actor_id,species_index,template_id,amount,model_c,model_h,model_o,model_s,model_redox,model_energy,reservoir_c,reservoir_h,reservoir_o,reservoir_s,reservoir_redox,reservoir_energy,residual_abs,balanced"
            )?;
            writer.flush()?;
            Some(writer)
        } else {
            None
        };

        Ok(Self {
            output_dir: dir,
            grid,
            ticks_writer,
            stoich_writer,
            stoich_v2_event_writer,
            registry: ReactionRegistry::new(),
        })
    }

    /// Write one summary row to ticks.csv.
    ///
    /// Call this once per tick. It computes population-level averages from the
    /// cell list and writes them as a single CSV row.
    ///
    /// # Arguments
    /// - `tick` — current simulation tick number
    /// - `cells` — slice of all living cells this tick
    /// - `divisions` — number of division events that happened this tick
    /// - `deaths` — number of death events that happened this tick
    ///
    /// # Why these columns?
    /// - `avg_energy` tracks whether the population is thriving or starving.
    /// - `avg_enzyme_a/b` (internal[5]/[6]) show metabolic specialization.
    /// - `avg_active_rxns` measures catalytic network complexity.
    /// - Per-z-layer counts reveal vertical niche structure.
    pub fn log_tick(
        &mut self,
        tick: u64,
        cells: &[CellState],
        divisions: u64,
        deaths: u64,
        hgt_events: u64,
    ) -> Result<()> {
        let Some(writer) = self.ticks_writer.as_mut() else {
            return Ok(());
        };

        let pop = cells.len() as f64;

        // Accumulate sums for averaging. Use f64 to avoid precision loss
        // when summing thousands of f32 values.
        let mut sum_energy: f64 = 0.0;
        let mut sum_enzyme_a: f64 = 0.0;
        let mut sum_enzyme_b: f64 = 0.0;
        let mut sum_active_rxns: f64 = 0.0;

        // One counter per z-layer for spatial distribution.
        let mut z_counts = vec![0u64; self.grid.z];

        for cell in cells {
            sum_energy += cell.internal[0] as f64;
            sum_enzyme_a += cell.internal[5] as f64;
            sum_enzyme_b += cell.internal[6] as f64;

            // Count active reactions: those with |v_max| above a tiny threshold.
            // This tells us how complex each cell's metabolism actually is,
            // as opposed to how many reaction slots are allocated.
            let active = cell
                .ruleset
                .reactions
                .iter()
                .filter(|r| r.v_max.abs() > 1e-9)
                .count();
            sum_active_rxns += active as f64;

            // Tally which z-layer this cell lives in.
            let z = cell.pos[2] as usize;
            if z < self.grid.z {
                z_counts[z] += 1;
            }
        }

        // Compute averages (guard against division by zero when population = 0).
        let (avg_e, avg_ea, avg_eb, avg_rxn) = if pop > 0.0 {
            (
                sum_energy / pop,
                sum_enzyme_a / pop,
                sum_enzyme_b / pop,
                sum_active_rxns / pop,
            )
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };

        // Write the fixed columns.
        write!(
            writer,
            "{},{},{:.6},{:.6},{:.6},{:.4},{},{},{}",
            tick,
            cells.len(),
            avg_e,
            avg_ea,
            avg_eb,
            avg_rxn,
            divisions,
            deaths,
            hgt_events
        )?;

        // Write per-z-layer cell counts.
        for count in z_counts.iter().take(self.grid.z) {
            write!(writer, ",{count}")?;
        }
        writeln!(writer)?;

        // Flush periodically so that if the process crashes we don't lose
        // too many ticks of data. BufWriter batches the small writes above
        // into fewer syscalls, so this flush is the only one that hits disk.
        writer.flush()?;

        Ok(())
    }

    pub fn log_stoich_tick(&mut self, tick: u64, ledger: &StoichTickLedger) -> Result<()> {
        let Some(writer) = self.stoich_writer.as_mut() else {
            return Ok(());
        };

        writeln!(
            writer,
            "{},{},{:.6},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
            tick,
            ledger.reaction_count,
            ledger.active_flux,
            ledger.imbalanced_reaction_count,
            ledger.unknown_species_flux,
            ledger.carbon_to_energy_flux,
            ledger.reductant_to_energy_flux,
            ledger.gross_material_abs,
            ledger.gross_total_abs,
            ledger.delta.c,
            ledger.delta.h,
            ledger.delta.o,
            ledger.delta.s,
            ledger.delta.redox,
            ledger.delta.energy,
        )?;
        writer.flush()
    }

    pub fn log_stoich_v2_events(&mut self, tick: u64, ledger: &StoichTickLedger) -> Result<()> {
        let Some(writer) = self.stoich_v2_event_writer.as_mut() else {
            return Ok(());
        };

        for event in &ledger.events {
            write_stoich_event_row(writer, tick, event)?;
        }
        writer.flush()
    }

    /// Write a per-z-layer chemistry snapshot to `chem_<tick>.csv`.
    ///
    /// For each z-layer, averages each of the S_EXT external chemical species
    /// across the entire XY plane, plus the average light intensity. This gives
    /// a 1D depth profile of the chemical environment — the vertical structure
    /// that cells experience and (hopefully) create through niche construction.
    ///
    /// # Arguments
    /// - `tick` — current tick (used in the filename)
    /// - `field` — the chemical concentration field
    /// - `light` — the light attenuation field
    pub fn snapshot_chemistry(&self, tick: u64, field: &Field, light: &LightField) -> Result<()> {
        let grid = field.grid();
        let path = self.output_dir.join(format!("chem_{}.csv", tick));
        let file = File::create(&path)?;
        let mut w = BufWriter::new(file);

        // Header: z, then one column per external species, then light.
        write!(w, "z")?;
        for s in 0..S_EXT {
            write!(w, ",species_{}_avg", s)?;
        }
        writeln!(w, ",light_avg")?;

        // Number of voxels in one XY plane — the denominator for averaging.
        let xy_count = (grid.x * grid.y) as f64;

        for z in 0..grid.z {
            write!(w, "{}", z)?;

            // Sum each species across the XY plane for this z-layer.
            // Using f64 accumulators because large XY planes can otherwise lose precision.
            for s in 0..S_EXT {
                let mut sum: f64 = 0.0;
                for y in 0..grid.y {
                    for x in 0..grid.x {
                        sum += field.get(x, y, z, s) as f64;
                    }
                }
                write!(w, ",{:.6}", sum / xy_count)?;
            }

            // Average light intensity for this z-layer.
            let mut light_sum: f64 = 0.0;
            for y in 0..grid.y {
                for x in 0..grid.x {
                    light_sum += light.get(x, y, z) as f64;
                }
            }
            writeln!(w, ",{:.6}", light_sum / xy_count)?;
        }

        w.flush()?;
        Ok(())
    }

    /// Write a full cell dump to `cells_<tick>.csv`.
    ///
    /// One row per living cell with its position, key internal concentrations,
    /// reaction network complexity, lineage, age, and quiescence state. This is
    /// the raw data you need for phylogenetic analysis, spatial plots, and
    /// checking whether distinct metabolic strategies have emerged.
    ///
    /// # Arguments
    /// - `tick` — current tick (used in the filename)
    /// - `cells` — slice of all living cells
    pub fn snapshot_cells(&self, tick: u64, cells: &[CellState]) -> Result<()> {
        let path = self.output_dir.join(format!("cells_{}.csv", tick));
        let file = File::create(&path)?;
        let mut w = BufWriter::new(file);

        // Header row.
        writeln!(
            w,
            "x,y,z,energy,enzyme_a,enzyme_b,active_rxns,lineage_id,age,quiescent,starter_type"
        )?;

        for cell in cells {
            // Count active reactions (same threshold as log_tick).
            let active_rxns = cell
                .ruleset
                .reactions
                .iter()
                .filter(|r| r.v_max.abs() > 1e-9)
                .count();

            writeln!(
                w,
                "{},{},{},{:.6},{:.6},{:.6},{},{},{},{},{}",
                cell.pos[0],
                cell.pos[1],
                cell.pos[2],
                cell.internal[0], // energy
                cell.internal[5], // enzyme_a
                cell.internal[6], // enzyme_b
                active_rxns,
                cell.lineage_id,
                cell.age,
                cell.quiescent as u8,
                cell.starter_type, // 0=photo, 1=chemo, 2=anaerobe
            )?;
        }

        w.flush()?;
        Ok(())
    }

    /// Write a reaction snapshot to `reactions_<tick>.csv`.
    ///
    /// For each living cell, records its position, ancestry (starter_type),
    /// and the registry ID of every active reaction. This is the primary
    /// data for tracking metabolic evolution — which reaction topologies
    /// exist, where, and in which lineages.
    ///
    /// The registry maps IDs back to (substrate, product, catalyst, cofactor)
    /// tuples. See `reaction_registry.csv` written at end of run.
    ///
    /// Inactive reactions (v_max ≈ 0) are omitted to keep file size manageable.
    pub fn snapshot_reactions(&mut self, tick: u64, cells: &[CellState]) -> Result<()> {
        let path = self.output_dir.join(format!("reactions_{}.csv", tick));
        let file = File::create(&path)?;
        let mut w = BufWriter::new(file);

        writeln!(w, "x,y,z,starter_type,lineage_id,reaction_ids")?;

        for cell in cells {
            // Collect IDs of all active reactions for this cell
            let mut ids: Vec<u32> = Vec::new();
            for rxn in &cell.ruleset.reactions {
                if rxn.v_max.abs() > 1e-9 {
                    let topo = ReactionTopology::from_reaction(rxn);
                    let id = self.registry.register(&topo);
                    ids.push(id);
                }
            }

            // Write as semicolon-separated IDs (commas are the CSV delimiter)
            let id_str: Vec<String> = ids.iter().map(|id| id.to_string()).collect();
            writeln!(
                w,
                "{},{},{},{},{},{}",
                cell.pos[0],
                cell.pos[1],
                cell.pos[2],
                cell.starter_type,
                cell.lineage_id,
                id_str.join(";"),
            )?;
        }

        w.flush()?;
        Ok(())
    }

    /// Write the reaction registry to disk. Call at end of run.
    pub fn write_registry(&self) -> Result<()> {
        self.registry.write_registry(&self.output_dir)
    }

    pub fn write_stoich_summary(&self, total_ticks: u32, ledger: &StoichRunLedger) -> Result<()> {
        #[derive(serde::Serialize)]
        struct LegacyStoichLedger {
            tick_count: u64,
            reaction_count: u64,
            active_flux: f32,
            imbalanced_reaction_count: u64,
            unknown_species_flux: f32,
            carbon_to_energy_flux: f32,
            reductant_to_energy_flux: f32,
            gross_material_abs: f32,
            gross_total_abs: f32,
            delta: marl_config::stoich::StoichBudgetDelta,
        }

        #[derive(serde::Serialize)]
        struct StoichSummary {
            total_ticks: u32,
            ledger: LegacyStoichLedger,
            material_abs_sum: f32,
            net_material_abs_sum: f32,
            gross_total_abs_sum: f32,
            notes: [&'static str; 3],
        }

        let path = self.output_dir.join("stoich_summary.json");
        let file = File::create(&path)?;
        let writer = BufWriter::new(file);
        let summary = StoichSummary {
            total_ticks,
            ledger: LegacyStoichLedger {
                tick_count: ledger.tick_count,
                reaction_count: ledger.reaction_count,
                active_flux: ledger.active_flux,
                imbalanced_reaction_count: ledger.imbalanced_reaction_count,
                unknown_species_flux: ledger.unknown_species_flux,
                carbon_to_energy_flux: ledger.carbon_to_energy_flux,
                reductant_to_energy_flux: ledger.reductant_to_energy_flux,
                gross_material_abs: ledger.gross_material_abs,
                gross_total_abs: ledger.gross_total_abs,
                delta: ledger.delta,
            },
            material_abs_sum: ledger.material_abs_sum(),
            net_material_abs_sum: ledger.net_material_abs_sum(),
            gross_total_abs_sum: ledger.total_abs_sum(),
            notes: [
                "legacy reactions are audited, not enforced",
                "energy and light are non-material bookkeeping slots",
                "inactive species contribute to unknown_species_flux",
            ],
        };
        serde_json::to_writer_pretty(writer, &summary).map_err(std::io::Error::other)
    }

    pub fn write_stoich_v2_summary(
        &self,
        total_ticks: u32,
        enforcement: StoichEnforcement,
        ledger: &StoichRunLedger,
    ) -> Result<()> {
        #[derive(serde::Serialize)]
        struct NamedStage {
            stage: StoichStage,
            summary: marl_config::stoich::StoichStageSummary,
        }

        #[derive(serde::Serialize)]
        struct NamedReservoir {
            reservoir: StoichReservoir,
            delta: marl_config::stoich::StoichBudgetDelta,
        }

        #[derive(serde::Serialize)]
        struct StoichV2Summary<'a> {
            schema_version: u32,
            total_ticks: u32,
            enforcement: &'static str,
            tolerance: f32,
            ledger: &'a StoichRunLedger,
            stages: Vec<NamedStage>,
            reservoirs: Vec<NamedReservoir>,
            gross_residual_abs_sum: f32,
            notes: [&'static str; 4],
        }

        let stages = (0..STOICH_STAGE_COUNT)
            .map(|i| NamedStage {
                stage: StoichStage::ALL[i],
                summary: ledger.stage_summaries[i],
            })
            .collect();
        let reservoirs = (0..STOICH_RESERVOIR_COUNT)
            .map(|i| NamedReservoir {
                reservoir: StoichReservoir::ALL[i],
                delta: ledger.reservoir_deltas[i],
            })
            .collect();
        let summary = StoichV2Summary {
            schema_version: 2,
            total_ticks,
            enforcement: enforcement.as_str(),
            tolerance: marl_config::stoich::STOICH_TOLERANCE,
            ledger,
            stages,
            reservoirs,
            gross_residual_abs_sum: ledger.gross_residual_abs(),
            notes: [
                "v2 covers full-system accepted flux accounting",
                "strict mode is opt-in and uses balanced reaction templates",
                "reservoir deltas close modeled sources and sinks without adding species slots",
                "v1 stoich files remain legacy-reaction audit outputs",
            ],
        };

        let path = self.output_dir.join("stoich_v2_summary.json");
        let file = File::create(&path)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, &summary).map_err(std::io::Error::other)
    }

    // ========================================================================
    // POST-RUN SUMMARY
    // ========================================================================

    /// Write a human-readable Markdown summary of an entire simulation run.
    ///
    /// Call this once, after the main tick loop finishes. It produces
    /// `summary.md` in the output directory — a quick-glance report that
    /// scientists can open in any Markdown viewer (or plain text editor) to
    /// see whether the run produced interesting dynamics.
    ///
    /// # What goes in the summary
    /// - Grid and timing parameters (so you know *which* run this was)
    /// - Final population + division/death totals
    /// - Vertical zonation by z-layer thirds (surface / middle / deep)
    /// - Chemical depth profiles sampled at the center column
    /// - Metabolic diversity statistics (active reactions, enzyme levels)
    /// - A brief automated assessment paragraph
    ///
    /// # Arguments
    /// - `total_ticks`     — how many ticks the run actually executed
    /// - `runtime_secs`    — wall-clock time for the run (seconds)
    /// - `cells`           — the final living cell population
    /// - `field`           — the chemical concentration field at end of run
    /// - `light`           — the light attenuation field at end of run
    /// - `total_divisions` — cumulative division events over the whole run
    /// - `total_deaths`    — cumulative death events over the whole run
    #[allow(clippy::too_many_arguments)]
    pub fn write_summary(
        &self,
        total_ticks: u32,
        runtime_secs: f32,
        cells: &[CellState],
        field: &Field,
        light: &LightField,
        total_divisions: u64,
        total_deaths: u64,
        sim: &SimulationConfig,
        rng_seed: u64,
        rng_algorithm: &str,
    ) -> Result<()> {
        let path = self.output_dir.join("summary.md");
        let file = File::create(&path)?;
        let mut w = BufWriter::new(file);

        // -- Derived constants used in several sections below ----------------
        let grid = field.grid();
        let voxel_count = grid.voxel_count().unwrap_or(0);
        let pop = cells.len();
        let grid_capacity = voxel_count; // one cell per voxel maximum
        let pop_pct = if grid_capacity > 0 {
            (pop as f64 / grid_capacity as f64) * 100.0
        } else {
            0.0
        };
        let ticks_per_sec = if runtime_secs > 0.0 {
            total_ticks as f64 / runtime_secs as f64
        } else {
            0.0
        };

        // ====================================================================
        // HEADER
        // ====================================================================
        writeln!(w, "# Run Summary")?;
        writeln!(w)?;

        // ====================================================================
        // PARAMETERS — so you can reconstruct which run produced this file
        // ====================================================================
        writeln!(w, "## Parameters")?;
        writeln!(w)?;
        writeln!(
            w,
            "- Grid: {}x{}x{} ({} voxels)",
            grid.x, grid.y, grid.z, voxel_count
        )?;
        writeln!(w, "- Ticks: {}", total_ticks)?;
        writeln!(w, "- RNG algorithm: {rng_algorithm}")?;
        writeln!(w, "- RNG seed: {rng_seed}")?;
        writeln!(w, "- Replay TOML: `rng_seed = {rng_seed}`")?;
        writeln!(
            w,
            "- Runtime: {:.1}s ({:.2} ticks/sec)",
            runtime_secs, ticks_per_sec
        )?;
        writeln!(w, "- Maintenance rate: {}", sim.lambda_maintenance)?;
        let sources = sim
            .effective_boundary_sources()
            .into_iter()
            .filter(|source| source.rate > 0.0)
            .map(|source| {
                format!(
                    "{}@{}={}",
                    external_species_name(source.species),
                    source.face.as_str(),
                    source.rate
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(w, "- Boundary sources: {sources}")?;
        writeln!(w)?;

        // ====================================================================
        // POPULATION — final counts and lifetime turnover
        // ====================================================================
        writeln!(w, "## Population")?;
        writeln!(w)?;
        writeln!(w, "- Final population: {} ({:.1}% capacity)", pop, pop_pct)?;
        writeln!(w, "- Total divisions: {}", total_divisions)?;
        writeln!(w, "- Total deaths: {}", total_deaths)?;
        // Turnover ratio: deaths per division. A ratio near 1.0 means the
        // population is roughly in steady state; >> 1 means net decline.
        let turnover = if total_divisions > 0 {
            total_deaths as f64 / total_divisions as f64
        } else {
            0.0
        };
        writeln!(w, "- Turnover ratio: {:.2}", turnover)?;
        writeln!(w)?;

        // ====================================================================
        // VERTICAL ZONATION — do cells stratify by depth?
        // ====================================================================
        // Count cells per z-layer, then aggregate into thirds.
        writeln!(w, "## Vertical Zonation")?;
        writeln!(w)?;

        let mut z_counts = vec![0u64; grid.z];
        for cell in cells {
            let z = cell.pos[2] as usize;
            if z < grid.z {
                z_counts[z] += 1;
            }
        }

        // Report per-layer counts as a compact table.
        writeln!(w, "Per-layer cell counts:")?;
        writeln!(w)?;
        writeln!(w, "| z | cells |")?;
        writeln!(w, "|---|-------|")?;
        for (z, count) in z_counts.iter().enumerate().take(grid.z) {
            writeln!(w, "| {z} | {count} |")?;
        }
        writeln!(w)?;

        // Aggregate into three zones: surface, middle, deep.
        let third = grid.z / 3;
        let two_thirds = 2 * grid.z / 3;

        let surface_count: u64 = z_counts[..third].iter().sum();
        let middle_count: u64 = z_counts[third..two_thirds].iter().sum();
        let deep_count: u64 = z_counts[two_thirds..].iter().sum();

        writeln!(w, "- Surface zone (z < {}): {} cells", third, surface_count)?;
        writeln!(
            w,
            "- Middle zone ({} <= z < {}): {} cells",
            third, two_thirds, middle_count
        )?;
        writeln!(w, "- Deep zone (z >= {}): {} cells", two_thirds, deep_count)?;

        // Zonation ratio: max zone / min zone. A ratio >> 1 means cells
        // strongly prefer certain depths — evidence of niche partitioning.
        let zone_max = surface_count.max(middle_count).max(deep_count);
        let zone_min = surface_count.min(middle_count).min(deep_count);
        let zonation_ratio = if zone_min > 0 {
            zone_max as f64 / zone_min as f64
        } else if zone_max > 0 {
            f64::INFINITY
        } else {
            1.0 // no cells at all, ratio is meaningless
        };
        writeln!(
            w,
            "- Zonation ratio (max_zone / min_zone): {:.1}",
            zonation_ratio
        )?;
        writeln!(w)?;

        // ====================================================================
        // CHEMICAL GRADIENTS — depth profile at center column
        // ====================================================================
        // Sample the center column at three depths
        // to see whether the run produced the expected redox gradient.
        writeln!(w, "## Chemical Gradients (center column)")?;
        writeln!(w)?;

        let cx = grid.x / 2;
        let cy = grid.y / 2;
        let z_surface = 0;
        let z_mid = grid.z / 2;
        let z_deep = grid.z - 1;

        // Species indices (from config/field layout):
        //   1 = oxidant, 2 = reductant, 3 = carbon, 4 = organic waste
        let oxidant_s = field.get(cx, cy, z_surface, 1);
        let oxidant_m = field.get(cx, cy, z_mid, 1);
        let oxidant_d = field.get(cx, cy, z_deep, 1);

        let reductant_s = field.get(cx, cy, z_surface, 2);
        let reductant_m = field.get(cx, cy, z_mid, 2);
        let reductant_d = field.get(cx, cy, z_deep, 2);

        let carbon_s = field.get(cx, cy, z_surface, 3);
        let carbon_m = field.get(cx, cy, z_mid, 3);
        let carbon_d = field.get(cx, cy, z_deep, 3);

        let waste_s = field.get(cx, cy, z_surface, 4);
        let waste_m = field.get(cx, cy, z_mid, 4);
        let waste_d = field.get(cx, cy, z_deep, 4);

        writeln!(
            w,
            "- Oxidant: surface={:.3}, mid={:.3}, deep={:.3}",
            oxidant_s, oxidant_m, oxidant_d
        )?;
        writeln!(
            w,
            "- Reductant: surface={:.3}, mid={:.3}, deep={:.3}",
            reductant_s, reductant_m, reductant_d
        )?;
        writeln!(
            w,
            "- Carbon: surface={:.3}, mid={:.3}, deep={:.3}",
            carbon_s, carbon_m, carbon_d
        )?;
        writeln!(
            w,
            "- Organic waste: surface={:.3}, mid={:.3}, deep={:.3}",
            waste_s, waste_m, waste_d
        )?;

        // Oxidant penetration depth: scan from surface downward, find first z
        let oxidant_pen = (0..grid.z)
            .find(|&z| field.get(cx, cy, z, 1) < 0.01)
            .unwrap_or(grid.z);
        writeln!(
            w,
            "- Oxidant penetration depth: z={} (first z where oxidant < 0.01)",
            oxidant_pen
        )?;

        // Reductant penetration depth: scan from bottom upward, find first z
        // (from bottom) where reductant drops below 0.01. Reports depth
        let reductant_pen = (0..grid.z)
            .rev()
            .find(|&z| field.get(cx, cy, z, 2) < 0.01)
            .unwrap_or(0);
        writeln!(
            w,
            "- Reductant penetration depth: z={} (first z from bottom where reductant < 0.01)",
            reductant_pen
        )?;

        // Light intensity at the three sample depths — shows how quickly
        // light attenuates through the column (affected by cell density and EPS).
        let light_s = light.get(cx, cy, z_surface);
        let light_m = light.get(cx, cy, z_mid);
        let light_d = light.get(cx, cy, z_deep);
        writeln!(
            w,
            "- Light: surface={:.3}, mid={:.3}, deep={:.3}",
            light_s, light_m, light_d
        )?;
        writeln!(w)?;

        // ====================================================================
        // METABOLIC DIVERSITY — how complex and varied are the cells?
        // ====================================================================
        writeln!(w, "## Metabolic Diversity")?;
        writeln!(w)?;

        if pop > 0 {
            let mut sum_active: f64 = 0.0;
            let mut sum_energy: f64 = 0.0;
            let mut sum_enzyme_a: f64 = 0.0;
            let mut sum_enzyme_b: f64 = 0.0;

            for cell in cells {
                let active = cell
                    .ruleset
                    .reactions
                    .iter()
                    .filter(|r| r.v_max.abs() > 1e-9)
                    .count();
                sum_active += active as f64;
                sum_energy += cell.internal[0] as f64;
                sum_enzyme_a += cell.internal[5] as f64;
                sum_enzyme_b += cell.internal[6] as f64;
            }

            let n = pop as f64;
            writeln!(w, "- Average active reactions: {:.1}", sum_active / n)?;
            writeln!(w, "- Average energy: {:.3}", sum_energy / n)?;
            writeln!(
                w,
                "- Average enzyme_A (internal[5]): {:.4}",
                sum_enzyme_a / n
            )?;
            writeln!(
                w,
                "- Average enzyme_B (internal[6]): {:.4}",
                sum_enzyme_b / n
            )?;
        } else {
            writeln!(w, "- (no living cells — all metrics are zero)")?;
        }
        writeln!(w)?;

        // ====================================================================
        // ASSESSMENT — brief automated interpretation of the data above
        // ====================================================================
        // This is NOT a fitness function — it's just a convenience for the
        // scientist to quickly see whether the run is worth investigating.
        writeln!(w, "## Assessment")?;
        writeln!(w)?;

        if pop_pct > 80.0 {
            writeln!(
                w,
                "- **Grid saturated** — consider larger grid or higher maintenance"
            )?;
        } else if pop_pct < 5.0 {
            writeln!(
                w,
                "- **Population sparse** — dynamics may be resource-limited"
            )?;
        } else {
            writeln!(
                w,
                "- Population at {:.1}% capacity — within normal operating range",
                pop_pct
            )?;
        }

        if zonation_ratio > 2.0 {
            writeln!(
                w,
                "- **Clear vertical zonation detected** (ratio {:.1}x)",
                zonation_ratio
            )?;
        } else {
            writeln!(
                w,
                "- Weak or no vertical zonation (ratio {:.1}x)",
                zonation_ratio
            )?;
        }

        if (oxidant_pen as f64) < (grid.z as f64 / 2.0) {
            writeln!(
                w,
                "- **Anoxic zone established** — oxidant penetrates to z={}",
                oxidant_pen
            )?;
        } else {
            writeln!(
                w,
                "- No clear anoxic zone — oxidant penetrates to z={}",
                oxidant_pen
            )?;
        }

        if turnover > 0.0 && turnover < 0.5 {
            writeln!(
                w,
                "- Low turnover ({:.2}) — population is growing",
                turnover
            )?;
        } else if turnover > 1.5 {
            writeln!(
                w,
                "- High turnover ({:.2}) — population is declining",
                turnover
            )?;
        } else if total_divisions > 0 {
            writeln!(
                w,
                "- Turnover near unity ({:.2}) — approximate steady state",
                turnover
            )?;
        }

        writeln!(w)?;

        w.flush()?;
        Ok(())
    }
}

fn write_stoich_event_row(
    writer: &mut BufWriter<File>,
    tick: u64,
    event: &StoichEvent,
) -> Result<()> {
    writeln!(
        writer,
        "{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{}",
        tick,
        event.stage.as_str(),
        event.kind.as_str(),
        event.actor_id,
        event.species_index,
        event.template_id,
        event.amount,
        event.model_delta.c,
        event.model_delta.h,
        event.model_delta.o,
        event.model_delta.s,
        event.model_delta.redox,
        event.model_delta.energy,
        event.reservoir_delta.c,
        event.reservoir_delta.h,
        event.reservoir_delta.o,
        event.reservoir_delta.s,
        event.reservoir_delta.redox,
        event.reservoir_delta.energy,
        event.residual_abs,
        event.balanced as u8,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use marl_config::stoich::LegacyReactionAudit;

    fn test_output_dir(name: &str) -> String {
        let dir = std::env::temp_dir().join(name);
        let _ = fs::remove_dir_all(&dir);
        dir.to_string_lossy().into_owned()
    }

    #[test]
    fn stoich_tick_log_and_summary_are_written() {
        let out_dir = test_output_dir("marl_output_stoich_audit_test");
        let mut logger = DataLogger::new(GridDims::default(), &out_dir, false, true, true).unwrap();
        let mut tick = StoichTickLedger::default();
        tick.record_legacy_reaction(
            LegacyReactionAudit::new(3, 0, 15, 0xFF, 2.0).actor(123),
            true,
        );
        tick.record_transport_uptake(3, 0, 1.5);
        tick.record_transport_secretion(4, 1, 0.25);

        logger.log_stoich_tick(7, &tick).unwrap();
        logger.log_stoich_v2_events(7, &tick).unwrap();
        let mut run = StoichRunLedger::default();
        run.add_tick(&tick);
        logger.write_stoich_summary(8, &run).unwrap();
        logger
            .write_stoich_v2_summary(8, marl_config::stoich::StoichEnforcement::Audit, &run)
            .unwrap();

        let dir = PathBuf::from(&out_dir);
        let ticks = fs::read_to_string(dir.join("stoich_ticks.csv")).unwrap();
        assert!(ticks.contains("tick,reaction_count"));
        assert!(ticks.contains("gross_material_abs,gross_total_abs"));
        assert!(ticks.contains("7,1,2.000000,1,0.000000,2.000000,0.000000,6.000000,16.000000"));

        let summary = fs::read_to_string(dir.join("stoich_summary.json")).unwrap();
        assert!(summary.contains("\"total_ticks\": 8"));
        assert!(summary.contains("\"carbon_to_energy_flux\": 2.0"));
        assert!(summary.contains("\"material_abs_sum\": 6.0"));
        assert!(summary.contains("\"gross_total_abs_sum\": 16.0"));

        let v2_summary = fs::read_to_string(dir.join("stoich_v2_summary.json")).unwrap();
        assert!(v2_summary.contains("\"schema_version\": 2"));
        assert!(v2_summary.contains("\"enforcement\": \"audit\""));
        assert!(v2_summary.contains("\"transport_flux_by_species\""));
        assert!(v2_summary.contains("\"transport_flux_by_species_starter\""));
        assert!(v2_summary.contains("\"reaction_byproduct_by_species_starter\""));
        assert!(v2_summary.contains("\"uptake_amount\": 1.5"));
        assert!(v2_summary.contains("\"secretion_amount\": 0.25"));
        let v2_events = fs::read_to_string(dir.join("stoich_v2_events.csv")).unwrap();
        assert!(v2_events.contains("tick,stage,kind"));
        assert!(v2_events.contains("reactions,legacy_reaction,123"));

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn tick_log_includes_hgt_event_count() {
        let out_dir = test_output_dir("marl_output_tick_hgt_test");
        let mut logger =
            DataLogger::new(GridDims { x: 2, y: 2, z: 2 }, &out_dir, true, false, false).unwrap();

        logger.log_tick(3, &[], 1, 2, 4).unwrap();

        let dir = PathBuf::from(&out_dir);
        let ticks = fs::read_to_string(dir.join("ticks.csv")).unwrap();
        assert!(ticks.contains("deaths_this_tick,hgt_events_this_tick,z0_cells"));
        assert!(ticks.contains("3,0,0.000000,0.000000,0.000000,0.0000,1,2,4,0,0"));

        let _ = fs::remove_dir_all(&dir);
    }
}
