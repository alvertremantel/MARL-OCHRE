use std::collections::HashMap;
use std::error::Error;
use std::f64::consts::E;
use std::fs;
use std::path::{Path, PathBuf};

use marl_format::{
    RULESET_FULL_CANONICAL_SIZE, RULESET_FULL_CELL_REF_STRIDE, RULESET_FULL_FORMAT_VERSION,
    RULESET_FULL_HEADER_SIZE, RULESET_FULL_MAGIC, RunMeta,
};
use marl_viewer_core::io::{
    LoadedCell, discover_field_ticks, load_cell_records, load_field_bytes, load_run_meta,
    read_binary_payload_bounded, snapshot_path,
};
use serde::Serialize;

pub type AnalysisResult<T> = Result<T, Box<dyn Error>>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScanMode {
    Sampled,
    All,
    LatestOnly,
    Explicit(Vec<u64>),
}

#[derive(Debug, Clone)]
pub struct AnalysisConfig {
    pub scan_mode: ScanMode,
    pub include_rulesets: bool,
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            scan_mode: ScanMode::Sampled,
            include_rulesets: true,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct RunAnalysis {
    pub run_dir: PathBuf,
    pub grid: [u32; 3],
    pub available_ticks: Vec<u64>,
    pub sampled_ticks: Vec<u64>,
    pub trajectory: Option<TrajectorySummary>,
    pub zonation: Option<ZonationSummary>,
    pub cells: Option<CellSummary>,
    pub cell_timeline: Vec<CellSummary>,
    pub chemistry: Vec<SnapshotChemistry>,
    pub rulesets: Option<RulesetSummary>,
    pub ruleset_timeline: Vec<RulesetSummary>,
    pub transporter_pair_timeline: Vec<TransportPairTimeline>,
    pub findings: Vec<Finding>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ComparisonAnalysis {
    pub runs: Vec<RunComparisonEntry>,
    pub findings: Vec<Finding>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunComparisonEntry {
    pub name: String,
    pub run_dir: PathBuf,
    pub final_population: Option<u64>,
    pub max_population: Option<u64>,
    pub growth_factor: Option<f64>,
    pub final_avg_energy: Option<f64>,
    pub top_fraction: Option<f64>,
    pub middle_fraction: Option<f64>,
    pub deep_fraction: Option<f64>,
    pub genotype_unique_count: Option<u32>,
    pub genotype_dominant_fraction: Option<f64>,
    pub transporter_active_per_cell: Option<f64>,
    pub transporter_gated_slots: Option<u64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TrajectorySummary {
    pub rows: usize,
    pub initial_tick: u64,
    pub final_tick: u64,
    pub initial_population: u64,
    pub final_population: u64,
    pub min_population: u64,
    pub max_population: u64,
    pub growth_factor: f64,
    pub total_divisions: u64,
    pub total_deaths: u64,
    pub final_avg_energy: f64,
    pub max_avg_energy: f64,
    pub min_avg_energy: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ZonationSummary {
    pub z_layers: usize,
    pub top_count: u64,
    pub middle_count: u64,
    pub deep_count: u64,
    pub top_fraction: f64,
    pub middle_fraction: f64,
    pub deep_fraction: f64,
    pub occupied_z_min: Option<usize>,
    pub occupied_z_max: Option<usize>,
    pub entropy: f64,
    pub z_counts: Vec<u64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CellSummary {
    pub tick: u64,
    pub count: usize,
    pub avg_energy: f64,
    pub min_energy: f32,
    pub max_energy: f32,
    pub starter_counts: [u64; 4],
    pub starter_ancestry: Vec<StarterAncestrySummary>,
}

#[derive(Debug, Clone, Serialize)]
pub struct StarterAncestrySummary {
    pub starter_type: u8,
    pub label: String,
    pub count: u64,
    pub fraction: f64,
    pub avg_energy: f64,
    pub top_count: u64,
    pub middle_count: u64,
    pub deep_count: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct SnapshotChemistry {
    pub tick: u64,
    pub nonfinite_values: u64,
    pub negative_values: u64,
    pub oxidant_penetration_z: Option<usize>,
    pub reductant_penetration_z: Option<usize>,
    pub redox_overlap_layers: usize,
    pub species_profiles: Vec<SpeciesProfile>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SpeciesProfile {
    pub species: usize,
    pub name: String,
    pub total_concentration: f64,
    pub max_value: f64,
    pub surface_mean: f64,
    pub middle_mean: f64,
    pub deep_mean: f64,
    pub per_z_mean: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RulesetSummary {
    pub tick: u64,
    pub unique_count: u32,
    pub cell_count: u32,
    pub dominant_dict_id: Option<u32>,
    pub dominant_count: u32,
    pub dominant_fraction: f64,
    pub shannon_diversity: f64,
    pub transporters: TransporterSummary,
}

#[derive(Debug, Clone, Serialize)]
pub struct TransporterSummary {
    pub active_slots: u64,
    pub avg_active_per_cell: f64,
    pub uptake_dominant_slots: u64,
    pub secretion_dominant_slots: u64,
    pub bidirectional_slots: u64,
    pub gated_slots: u64,
    pub avg_abs_gate_weight: f64,
    pub uptake_rate_sum: f64,
    pub secretion_rate_sum: f64,
    pub common_pairs: Vec<TransportPairSummary>,
    pub dominant_genotype: Option<GenotypeTransportSummary>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TransportPairSummary {
    pub ext_species: u8,
    pub int_species: u8,
    pub active_slots: u64,
    pub avg_uptake_rate: f64,
    pub avg_secrete_rate: f64,
    pub gated_slots: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct GenotypeTransportSummary {
    pub dict_id: u32,
    pub active_slots: u32,
    pub gated_slots: u32,
    pub uptake_dominant_slots: u32,
    pub secretion_dominant_slots: u32,
}

#[derive(Debug, Clone, Serialize)]
pub struct TransportPairTimeline {
    pub ext_species: u8,
    pub int_species: u8,
    pub total_active_slots: u64,
    pub first_tick: Option<u64>,
    pub latest_tick: Option<u64>,
    pub points: Vec<TransportPairTimelinePoint>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TransportPairTimelinePoint {
    pub tick: u64,
    pub active_slots: u64,
    pub avg_uptake_rate: f64,
    pub avg_secrete_rate: f64,
    pub gated_slots: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct Finding {
    pub id: String,
    pub level: FindingLevel,
    pub title: String,
    pub narrative: String,
    pub evidence: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum FindingLevel {
    Info,
    Interesting,
    Warning,
}

#[derive(Debug, Clone)]
struct TickRow {
    tick: u64,
    population: u64,
    avg_energy: f64,
    divisions: u64,
    deaths: u64,
    z_counts: Vec<u64>,
}

const RULESET_FULL_FORMAT_VERSION_V1: u32 = 1;
const RULESET_FULL_CANONICAL_SIZE_V1: u32 = 536;
const TRANSPORTER_COUNT: usize = 8;
const RECEPTOR_PAYLOAD_BYTES: usize = 8 * 12;
const TRANSPORT_V1_STRIDE: usize = 10;
const TRANSPORT_V2_STRIDE: usize = 15;
const ACTIVE_TRANSPORT_THRESHOLD: f32 = 1e-9;
const MAX_COMMON_TRANSPORT_PAIRS: usize = 8;

#[derive(Debug, Clone, Copy)]
struct ParsedTransporter {
    uptake_rate: f32,
    secrete_rate: f32,
    ext_species: u8,
    int_species: u8,
    gate_weight: f32,
}

#[derive(Debug, Clone, Copy, Default)]
struct GenotypeTransportStats {
    active_slots: u32,
    gated_slots: u32,
    uptake_dominant_slots: u32,
    secretion_dominant_slots: u32,
}

#[derive(Debug, Clone, Default)]
struct TransportPairAccumulator {
    active_slots: u64,
    uptake_rate_sum: f64,
    secrete_rate_sum: f64,
    gated_slots: u64,
}

pub fn analyze_run(run_dir: impl AsRef<Path>, cfg: &AnalysisConfig) -> AnalysisResult<RunAnalysis> {
    let run_dir = run_dir.as_ref();
    let meta = load_run_meta(run_dir)?;
    let available_ticks = discover_field_ticks(run_dir)?;
    let sampled_ticks = choose_ticks(&available_ticks, &cfg.scan_mode);
    let mut warnings = Vec::new();

    let tick_rows = match read_ticks_csv(run_dir) {
        Ok(rows) => rows,
        Err(err) => {
            warnings.push(format!("ticks.csv unavailable or invalid: {err}"));
            Vec::new()
        }
    };

    let trajectory = summarize_trajectory(&tick_rows);
    let mut zonation = tick_rows
        .last()
        .map(|row| summarize_zonation(&row.z_counts));

    let mut chemistry = Vec::new();
    for &tick in &sampled_ticks {
        match load_field_bytes(run_dir, tick, &meta)
            .and_then(|bytes| summarize_field(tick, &meta, &bytes))
        {
            Ok(summary) => chemistry.push(summary),
            Err(err) => warnings.push(format!("tick {tick}: field analysis skipped: {err}")),
        }
    }

    let mut cell_timeline = Vec::new();
    for &tick in &sampled_ticks {
        match load_cell_records(run_dir, tick, &meta) {
            Ok(cells) => {
                if zonation.is_none() && sampled_ticks.last().copied() == Some(tick) {
                    zonation = Some(summarize_zonation_from_cells(meta.grid_z as usize, &cells));
                }
                cell_timeline.push(summarize_cells(tick, meta.grid_z as usize, &cells));
            }
            Err(err) => {
                warnings.push(format!("tick {tick}: cell analysis skipped: {err}"));
            }
        }
    }
    let latest_tick = sampled_ticks
        .last()
        .copied()
        .or_else(|| available_ticks.last().copied());
    let cells = cell_timeline.last().cloned();

    let mut ruleset_timeline = Vec::new();
    if cfg.include_rulesets && latest_tick.is_some() {
        if full_rulesets_enabled(run_dir)? {
            for &tick in &sampled_ticks {
                match summarize_rulesets(run_dir, tick, &meta) {
                    Ok(summary) => ruleset_timeline.push(summary),
                    Err(err) => {
                        warnings.push(format!("tick {tick}: ruleset analysis skipped: {err}"));
                    }
                }
            }
        } else {
            warnings.push(
                "ruleset analysis skipped: full ruleset dumps were not enabled for this run"
                    .to_string(),
            );
        }
    }
    let rulesets = ruleset_timeline.last().cloned();
    let transporter_pair_timeline = summarize_transporter_pair_timeline(&ruleset_timeline);

    let mut analysis = RunAnalysis {
        run_dir: run_dir.to_path_buf(),
        grid: [meta.grid_x, meta.grid_y, meta.grid_z],
        available_ticks,
        sampled_ticks,
        trajectory,
        zonation,
        cells,
        cell_timeline,
        chemistry,
        rulesets,
        ruleset_timeline,
        transporter_pair_timeline,
        findings: Vec::new(),
        warnings,
    };
    analysis.findings = classify_run_findings(&analysis);
    Ok(analysis)
}

pub fn compare_runs(
    run_dirs: &[PathBuf],
    cfg: &AnalysisConfig,
) -> AnalysisResult<ComparisonAnalysis> {
    let mut runs = Vec::new();
    let mut warnings = Vec::new();
    for run_dir in run_dirs {
        match analyze_run(run_dir, cfg) {
            Ok(analysis) => {
                warnings.extend(
                    analysis
                        .warnings
                        .iter()
                        .map(|warning| format!("{}: {warning}", run_name(&analysis.run_dir))),
                );
                runs.push(RunComparisonEntry::from_analysis(&analysis));
            }
            Err(err) => warnings.push(format!("{}: analysis failed: {err}", run_name(run_dir))),
        }
    }
    let findings = classify_comparison_findings(&runs);
    Ok(ComparisonAnalysis {
        runs,
        findings,
        warnings,
    })
}

pub fn write_run_reports(
    analysis: &RunAnalysis,
    out_dir: impl AsRef<Path>,
    write_json: bool,
    write_markdown: bool,
) -> AnalysisResult<()> {
    if !write_json && !write_markdown {
        return Ok(());
    }
    let out_dir = out_dir.as_ref();
    fs::create_dir_all(out_dir)?;
    if write_json {
        let path = out_dir.join("analysis.json");
        fs::write(path, serde_json::to_string_pretty(analysis)?)?;
    }
    if write_markdown {
        let path = out_dir.join("analysis.md");
        fs::write(path, render_run_markdown(analysis))?;
    }
    Ok(())
}

pub fn write_comparison_reports(
    analysis: &ComparisonAnalysis,
    out_dir: impl AsRef<Path>,
    write_json: bool,
    write_markdown: bool,
) -> AnalysisResult<()> {
    if !write_json && !write_markdown {
        return Ok(());
    }
    let out_dir = out_dir.as_ref();
    fs::create_dir_all(out_dir)?;
    if write_json {
        let path = out_dir.join("compare.json");
        fs::write(path, serde_json::to_string_pretty(analysis)?)?;
    }
    if write_markdown {
        let path = out_dir.join("compare.md");
        fs::write(path, render_comparison_markdown(analysis))?;
    }
    Ok(())
}

pub fn render_run_terminal(analysis: &RunAnalysis) -> String {
    let mut out = String::new();
    out.push_str(&format!("{}\n", run_name(&analysis.run_dir)));
    out.push_str(&format!(
        "  grid={}x{}x{}, snapshots={} sampled={:?}\n",
        analysis.grid[0],
        analysis.grid[1],
        analysis.grid[2],
        analysis.available_ticks.len(),
        analysis.sampled_ticks
    ));
    if let Some(traj) = &analysis.trajectory {
        out.push_str(&format!(
            "  population: {} -> {} (max {}, growth {:.2}x), div={}, death={}, final E={:.3}\n",
            traj.initial_population,
            traj.final_population,
            traj.max_population,
            traj.growth_factor,
            traj.total_divisions,
            traj.total_deaths,
            traj.final_avg_energy
        ));
    }
    if let Some(zonation) = &analysis.zonation {
        out.push_str(&format!(
            "  zones: top {:.1}%, mid {:.1}%, deep {:.1}%, occupied z={:?}..{:?}\n",
            zonation.top_fraction * 100.0,
            zonation.middle_fraction * 100.0,
            zonation.deep_fraction * 100.0,
            zonation.occupied_z_min,
            zonation.occupied_z_max
        ));
    }
    if let Some(cells) = &analysis.cells {
        let ancestry = cells
            .starter_ancestry
            .iter()
            .map(|starter| {
                format!(
                    "{}:{} ({:.1}%)",
                    starter.label,
                    starter.count,
                    starter.fraction * 100.0
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        out.push_str(&format!("  ancestry: {ancestry}\n"));
    }
    if let Some(rulesets) = &analysis.rulesets {
        out.push_str(&format!(
            "  rulesets: {} unique / {} cells, dominant {:.1}%, H={:.2}\n",
            rulesets.unique_count,
            rulesets.cell_count,
            rulesets.dominant_fraction * 100.0,
            rulesets.shannon_diversity
        ));
        out.push_str(&format!(
            "  transporters: active {:.2}/cell, gated {} slots, uptake_dom={}, secretion_dom={}\n",
            rulesets.transporters.avg_active_per_cell,
            rulesets.transporters.gated_slots,
            rulesets.transporters.uptake_dominant_slots,
            rulesets.transporters.secretion_dominant_slots
        ));
        if !rulesets.transporters.common_pairs.is_empty() {
            let pairs = rulesets
                .transporters
                .common_pairs
                .iter()
                .take(3)
                .map(|pair| {
                    format!(
                        "ext{}->int{}:{}",
                        pair.ext_species, pair.int_species, pair.active_slots
                    )
                })
                .collect::<Vec<_>>()
                .join(", ");
            out.push_str(&format!("  common transporter pairs: {pairs}\n"));
        }
    }
    if analysis.ruleset_timeline.len() > 1 {
        out.push_str(&format!(
            "  ruleset timeline: {} sampled ticks, active transport {:.2}->{:.2}/cell, gated {:.2}->{:.2}/cell\n",
            analysis.ruleset_timeline.len(),
            analysis
                .ruleset_timeline
                .first()
                .map(|rulesets| rulesets.transporters.avg_active_per_cell)
                .unwrap_or(0.0),
            analysis
                .ruleset_timeline
                .last()
                .map(|rulesets| rulesets.transporters.avg_active_per_cell)
                .unwrap_or(0.0),
            analysis
                .ruleset_timeline
                .first()
                .map(|rulesets| rulesets.transporters.gated_slots as f64 / rulesets.cell_count.max(1) as f64)
                .unwrap_or(0.0),
            analysis
                .ruleset_timeline
                .last()
                .map(|rulesets| rulesets.transporters.gated_slots as f64 / rulesets.cell_count.max(1) as f64)
                .unwrap_or(0.0)
        ));
    }
    if let (Some(first), Some(latest)) = (
        leading_transport_pair(analysis.ruleset_timeline.first()),
        leading_transport_pair(analysis.ruleset_timeline.last()),
    ) && (first.ext_species, first.int_species) != (latest.ext_species, latest.int_species)
    {
        out.push_str(&format!(
            "  leading transporter pair changed: ext{}->int{} (tick {}) -> ext{}->int{} (tick {})\n",
            first.ext_species,
            first.int_species,
            analysis.ruleset_timeline.first().map(|rulesets| rulesets.tick).unwrap_or(0),
            latest.ext_species,
            latest.int_species,
            analysis.ruleset_timeline.last().map(|rulesets| rulesets.tick).unwrap_or(0)
        ));
    }
    for finding in &analysis.findings {
        out.push_str(&format!(
            "  [{:?}] {}: {}\n",
            finding.level, finding.title, finding.narrative
        ));
    }
    if !analysis.warnings.is_empty() {
        out.push_str("  warnings:\n");
        for warning in &analysis.warnings {
            out.push_str(&format!("    - {warning}\n"));
        }
    }
    out
}

pub fn render_comparison_terminal(analysis: &ComparisonAnalysis) -> String {
    let mut out = String::new();
    out.push_str("MARL run comparison\n");
    for run in &analysis.runs {
        out.push_str(&format!(
            "  {}: final_pop={:?}, growth={:?}, top={:?}, dominant_genotype={:?}, active_transporters={:?}\n",
            run.name,
            run.final_population,
            run.growth_factor.map(|v| format!("{v:.2}x")),
            run.top_fraction.map(|v| format!("{:.1}%", v * 100.0)),
            run.genotype_dominant_fraction
                .map(|v| format!("{:.1}%", v * 100.0)),
            run.transporter_active_per_cell
                .map(|v| format!("{v:.2}/cell"))
        ));
    }
    for finding in &analysis.findings {
        out.push_str(&format!(
            "  [{:?}] {}: {}\n",
            finding.level, finding.title, finding.narrative
        ));
    }
    if !analysis.warnings.is_empty() {
        out.push_str("  warnings:\n");
        for warning in &analysis.warnings {
            out.push_str(&format!("    - {warning}\n"));
        }
    }
    out
}

fn leading_transport_pair(rulesets: Option<&RulesetSummary>) -> Option<&TransportPairSummary> {
    rulesets
        .and_then(|rulesets| rulesets.transporters.common_pairs.first())
        .filter(|pair| pair.active_slots > 0)
}

pub fn render_run_markdown(analysis: &RunAnalysis) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "# MARL Analysis: {}\n\n",
        run_name(&analysis.run_dir)
    ));
    out.push_str("## Overview\n\n");
    out.push_str(&format!(
        "- Grid: {}x{}x{}\n- Available snapshots: {}\n- Analyzed snapshots: {:?}\n",
        analysis.grid[0],
        analysis.grid[1],
        analysis.grid[2],
        analysis.available_ticks.len(),
        analysis.sampled_ticks
    ));
    if let Some(traj) = &analysis.trajectory {
        out.push_str(&format!(
            "- Population: {} -> {} (max {}, growth {:.2}x)\n- Total divisions/deaths: {}/{}\n- Final average energy: {:.3}\n",
            traj.initial_population,
            traj.final_population,
            traj.max_population,
            traj.growth_factor,
            traj.total_divisions,
            traj.total_deaths,
            traj.final_avg_energy
        ));
    }
    if let Some(zonation) = &analysis.zonation {
        out.push_str("\n## Zonation\n\n");
        out.push_str(&format!(
            "- Surface: {} cells ({:.1}%)\n- Middle: {} cells ({:.1}%)\n- Deep: {} cells ({:.1}%)\n- Occupied z span: {:?}..{:?}\n- Layer entropy: {:.3}\n",
            zonation.top_count,
            zonation.top_fraction * 100.0,
            zonation.middle_count,
            zonation.middle_fraction * 100.0,
            zonation.deep_count,
            zonation.deep_fraction * 100.0,
            zonation.occupied_z_min,
            zonation.occupied_z_max,
            zonation.entropy
        ));
    }
    if let Some(cells) = &analysis.cells
        && !cells.starter_ancestry.is_empty()
    {
        out.push_str("\n## Starter Ancestry\n\n");
        out.push_str(&format!(
            "- Latest cell snapshot: tick {} ({} cells)\n",
            cells.tick, cells.count
        ));
        out.push_str("| starter | cells | fraction | avg_energy | top | middle | deep |\n");
        out.push_str("|---|---:|---:|---:|---:|---:|---:|\n");
        for starter in &cells.starter_ancestry {
            out.push_str(&format!(
                "| {} | {} | {:.1}% | {:.3} | {} | {} | {} |\n",
                starter.label,
                starter.count,
                starter.fraction * 100.0,
                starter.avg_energy,
                starter.top_count,
                starter.middle_count,
                starter.deep_count
            ));
        }
    }
    if analysis.cell_timeline.len() > 1 {
        out.push_str("\n### Ancestry Timeline\n\n");
        out.push_str("| tick | cells | phototroph | chemolithotroph | anaerobe | other |\n");
        out.push_str("|---:|---:|---:|---:|---:|---:|\n");
        for cells in &analysis.cell_timeline {
            out.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} |\n",
                cells.tick,
                cells.count,
                cells.starter_counts[0],
                cells.starter_counts[1],
                cells.starter_counts[2],
                cells.starter_counts[3]
            ));
        }
    }
    if !analysis.chemistry.is_empty() {
        out.push_str("\n## Chemistry\n\n");
        out.push_str("| tick | oxidant_pen_z | reductant_pen_z | redox_overlap_layers | nonfinite | negative |\n");
        out.push_str("|---:|---:|---:|---:|---:|---:|\n");
        for chem in &analysis.chemistry {
            out.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} |\n",
                chem.tick,
                display_opt_usize(chem.oxidant_penetration_z),
                display_opt_usize(chem.reductant_penetration_z),
                chem.redox_overlap_layers,
                chem.nonfinite_values,
                chem.negative_values
            ));
        }
    }
    if let Some(rulesets) = &analysis.rulesets {
        out.push_str("\n## Rulesets\n\n");
        out.push_str(&format!(
            "- Tick: {}\n- Unique genotypes: {}\n- Cells with rulesets: {}\n- Dominant genotype: {:?} ({} cells, {:.1}%)\n- Shannon diversity: {:.3}\n",
            rulesets.tick,
            rulesets.unique_count,
            rulesets.cell_count,
            rulesets.dominant_dict_id,
            rulesets.dominant_count,
            rulesets.dominant_fraction * 100.0,
            rulesets.shannon_diversity
        ));
        out.push_str(&format!(
            "- Active transporter slots: {} ({:.2} per cell)\n- Gated active slots: {} (avg |gate_weight| {:.3})\n- Uptake/secretion dominant active slots: {}/{}\n- Bidirectional active slots: {}\n",
            rulesets.transporters.active_slots,
            rulesets.transporters.avg_active_per_cell,
            rulesets.transporters.gated_slots,
            rulesets.transporters.avg_abs_gate_weight,
            rulesets.transporters.uptake_dominant_slots,
            rulesets.transporters.secretion_dominant_slots,
            rulesets.transporters.bidirectional_slots
        ));
        if let Some(dominant) = &rulesets.transporters.dominant_genotype {
            out.push_str(&format!(
                "- Dominant genotype transport: {} active slots, {} gated, uptake/secretion dominant {}/{}\n",
                dominant.active_slots,
                dominant.gated_slots,
                dominant.uptake_dominant_slots,
                dominant.secretion_dominant_slots
            ));
        }
        if !rulesets.transporters.common_pairs.is_empty() {
            out.push_str("\n### Common Transporter Pairs\n\n");
            out.push_str(
                "| ext | int | active_slots | avg_uptake | avg_secretion | gated_slots |\n",
            );
            out.push_str("|---:|---:|---:|---:|---:|---:|\n");
            for pair in &rulesets.transporters.common_pairs {
                out.push_str(&format!(
                    "| {} | {} | {} | {:.3} | {:.3} | {} |\n",
                    pair.ext_species,
                    pair.int_species,
                    pair.active_slots,
                    pair.avg_uptake_rate,
                    pair.avg_secrete_rate,
                    pair.gated_slots
                ));
            }
        }
    }
    if analysis.ruleset_timeline.len() > 1 {
        out.push_str("\n### Ruleset Timeline\n\n");
        out.push_str("| tick | cells | unique | dominant% | shannon | active_transport/cell | gated_transport/cell | avg_abs_gate |\n");
        out.push_str("|---:|---:|---:|---:|---:|---:|---:|---:|\n");
        for rulesets in &analysis.ruleset_timeline {
            out.push_str(&format!(
                "| {} | {} | {} | {:.1}% | {:.3} | {:.2} | {:.2} | {:.3} |\n",
                rulesets.tick,
                rulesets.cell_count,
                rulesets.unique_count,
                rulesets.dominant_fraction * 100.0,
                rulesets.shannon_diversity,
                rulesets.transporters.avg_active_per_cell,
                rulesets.transporters.gated_slots as f64 / rulesets.cell_count.max(1) as f64,
                rulesets.transporters.avg_abs_gate_weight
            ));
        }
    }
    if !analysis.transporter_pair_timeline.is_empty() {
        out.push_str("\n### Transporter Pair Timeline\n\n");
        out.push_str("| tick |");
        for pair in &analysis.transporter_pair_timeline {
            out.push_str(&format!(
                " ext{}->int{} |",
                pair.ext_species, pair.int_species
            ));
        }
        out.push('\n');
        out.push_str("|---:|");
        for _ in &analysis.transporter_pair_timeline {
            out.push_str("---:|");
        }
        out.push('\n');
        for point_index in 0..analysis.transporter_pair_timeline[0].points.len() {
            let tick = analysis.transporter_pair_timeline[0].points[point_index].tick;
            out.push_str(&format!("| {tick} |"));
            for pair in &analysis.transporter_pair_timeline {
                let point = &pair.points[point_index];
                if point.active_slots == 0 {
                    out.push_str(" - |");
                } else if point.gated_slots > 0 {
                    out.push_str(&format!(
                        " {} (g{}) |",
                        point.active_slots, point.gated_slots
                    ));
                } else {
                    out.push_str(&format!(" {} |", point.active_slots));
                }
            }
            out.push('\n');
        }
    }
    out.push_str("\n## Findings\n\n");
    if analysis.findings.is_empty() {
        out.push_str("- No major automated findings.\n");
    } else {
        for finding in &analysis.findings {
            out.push_str(&format!(
                "- **{:?}: {}** — {}\n",
                finding.level, finding.title, finding.narrative
            ));
            if !finding.evidence.is_empty() {
                out.push_str(&format!("  Evidence: {}.\n", finding.evidence.join("; ")));
            }
        }
    }
    if !analysis.warnings.is_empty() {
        out.push_str("\n## Warnings\n\n");
        for warning in &analysis.warnings {
            out.push_str(&format!("- {warning}\n"));
        }
    }
    out
}

pub fn render_comparison_markdown(analysis: &ComparisonAnalysis) -> String {
    let mut out = String::new();
    out.push_str("# MARL Run Comparison\n\n");
    out.push_str("| run | final_pop | max_pop | growth | final_energy | top% | mid% | deep% | unique_genotypes | dominant_genotype% | active_transporters/cell | gated_slots |\n");
    out.push_str("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for run in &analysis.runs {
        out.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            run.name,
            display_opt_u64(run.final_population),
            display_opt_u64(run.max_population),
            display_opt_f64(run.growth_factor),
            display_opt_f64(run.final_avg_energy),
            display_opt_pct(run.top_fraction),
            display_opt_pct(run.middle_fraction),
            display_opt_pct(run.deep_fraction),
            run.genotype_unique_count
                .map(|v| v.to_string())
                .unwrap_or_else(|| "n/a".to_string()),
            display_opt_pct(run.genotype_dominant_fraction),
            display_opt_f64(run.transporter_active_per_cell),
            run.transporter_gated_slots
                .map(|v| v.to_string())
                .unwrap_or_else(|| "n/a".to_string())
        ));
    }
    out.push_str("\n## Findings\n\n");
    if analysis.findings.is_empty() {
        out.push_str("- No major automated comparison findings.\n");
    } else {
        for finding in &analysis.findings {
            out.push_str(&format!(
                "- **{:?}: {}** — {}\n",
                finding.level, finding.title, finding.narrative
            ));
        }
    }
    if !analysis.warnings.is_empty() {
        out.push_str("\n## Warnings\n\n");
        for warning in &analysis.warnings {
            out.push_str(&format!("- {warning}\n"));
        }
    }
    out
}

fn choose_ticks(available: &[u64], mode: &ScanMode) -> Vec<u64> {
    match mode {
        ScanMode::Sampled => sample_ticks(available),
        ScanMode::All => available.to_vec(),
        ScanMode::LatestOnly => available.last().copied().into_iter().collect(),
        ScanMode::Explicit(ticks) => {
            let mut ticks = ticks.clone();
            ticks.sort_unstable();
            ticks.dedup();
            ticks
        }
    }
}

fn sample_ticks(available: &[u64]) -> Vec<u64> {
    if available.is_empty() {
        return Vec::new();
    }
    let mut ticks = vec![
        available[0],
        available[available.len() / 2],
        *available.last().unwrap(),
    ];
    ticks.sort_unstable();
    ticks.dedup();
    ticks
}

fn read_ticks_csv(run_dir: &Path) -> AnalysisResult<Vec<TickRow>> {
    let path = run_dir.join("ticks.csv");
    let text = fs::read_to_string(&path)?;
    parse_ticks_csv(&text)
}

fn parse_ticks_csv(text: &str) -> AnalysisResult<Vec<TickRow>> {
    let mut lines = text.lines();
    let header = lines.next().ok_or("ticks.csv is empty")?;
    let columns: Vec<&str> = header.split(',').collect();
    let index = |name: &str| -> AnalysisResult<usize> {
        columns
            .iter()
            .position(|column| *column == name)
            .ok_or_else(|| format!("ticks.csv missing column {name}").into())
    };
    let tick_i = index("tick")?;
    let pop_i = index("population")?;
    let energy_i = index("avg_energy")?;
    let div_i = index("divisions_this_tick")?;
    let death_i = index("deaths_this_tick")?;
    let z_indices: Vec<usize> = columns
        .iter()
        .enumerate()
        .filter_map(|(i, column)| {
            (column.starts_with('z') && column.ends_with("_cells")).then_some(i)
        })
        .collect();

    let mut rows = Vec::new();
    for (line_no, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let values: Vec<&str> = line.split(',').collect();
        let get = |i: usize| -> AnalysisResult<&str> {
            values.get(i).copied().ok_or_else(|| {
                format!("ticks.csv line {} missing column {}", line_no + 2, i).into()
            })
        };
        rows.push(TickRow {
            tick: get(tick_i)?.parse()?,
            population: get(pop_i)?.parse()?,
            avg_energy: get(energy_i)?.parse()?,
            divisions: get(div_i)?.parse()?,
            deaths: get(death_i)?.parse()?,
            z_counts: z_indices
                .iter()
                .map(|&i| get(i).and_then(|value| value.parse::<u64>().map_err(|e| e.into())))
                .collect::<AnalysisResult<Vec<_>>>()?,
        });
    }
    Ok(rows)
}

fn summarize_trajectory(rows: &[TickRow]) -> Option<TrajectorySummary> {
    let first = rows.first()?;
    let last = rows.last()?;
    let min_population = rows.iter().map(|row| row.population).min().unwrap_or(0);
    let max_population = rows.iter().map(|row| row.population).max().unwrap_or(0);
    let min_avg_energy = rows
        .iter()
        .map(|row| row.avg_energy)
        .fold(f64::INFINITY, f64::min);
    let max_avg_energy = rows
        .iter()
        .map(|row| row.avg_energy)
        .fold(f64::NEG_INFINITY, f64::max);
    Some(TrajectorySummary {
        rows: rows.len(),
        initial_tick: first.tick,
        final_tick: last.tick,
        initial_population: first.population,
        final_population: last.population,
        min_population,
        max_population,
        growth_factor: if first.population > 0 {
            last.population as f64 / first.population as f64
        } else {
            0.0
        },
        total_divisions: rows.iter().map(|row| row.divisions).sum(),
        total_deaths: rows.iter().map(|row| row.deaths).sum(),
        final_avg_energy: last.avg_energy,
        max_avg_energy,
        min_avg_energy,
    })
}

fn summarize_zonation(z_counts: &[u64]) -> ZonationSummary {
    let z_layers = z_counts.len();
    let total: u64 = z_counts.iter().sum();
    let third = z_layers / 3;
    let two_thirds = 2 * z_layers / 3;
    let top_count: u64 = z_counts[..third].iter().sum();
    let middle_count: u64 = z_counts[third..two_thirds].iter().sum();
    let deep_count: u64 = z_counts[two_thirds..].iter().sum();
    let occupied_z_min = z_counts.iter().position(|count| *count > 0);
    let occupied_z_max = z_counts.iter().rposition(|count| *count > 0);
    let entropy = if total > 0 {
        z_counts
            .iter()
            .filter(|count| **count > 0)
            .map(|count| {
                let p = *count as f64 / total as f64;
                -p * p.log(E)
            })
            .sum()
    } else {
        0.0
    };
    ZonationSummary {
        z_layers,
        top_count,
        middle_count,
        deep_count,
        top_fraction: fraction(top_count, total),
        middle_fraction: fraction(middle_count, total),
        deep_fraction: fraction(deep_count, total),
        occupied_z_min,
        occupied_z_max,
        entropy,
        z_counts: z_counts.to_vec(),
    }
}

fn summarize_zonation_from_cells(z_layers: usize, cells: &[LoadedCell]) -> ZonationSummary {
    let mut z_counts = vec![0; z_layers];
    for cell in cells {
        if let Some(count) = z_counts.get_mut(cell.pos[2] as usize) {
            *count += 1;
        }
    }
    summarize_zonation(&z_counts)
}

fn summarize_cells(tick: u64, z_layers: usize, cells: &[LoadedCell]) -> CellSummary {
    let mut starter_counts = [0; 4];
    let mut starter_energy_sums = [0.0; 4];
    let mut starter_zone_counts = [[0u64; 3]; 4];
    let mut sum = 0.0;
    let mut min_energy = f32::INFINITY;
    let mut max_energy = f32::NEG_INFINITY;
    let third = z_layers / 3;
    let two_thirds = 2 * z_layers / 3;
    for cell in cells {
        let starter = match cell.starter_type {
            0..=2 => cell.starter_type as usize,
            _ => 3,
        };
        starter_counts[starter] += 1;
        starter_energy_sums[starter] += cell.energy as f64;
        let z = cell.pos[2] as usize;
        let zone = if z < third {
            0
        } else if z < two_thirds {
            1
        } else {
            2
        };
        starter_zone_counts[starter][zone] += 1;
        sum += cell.energy as f64;
        min_energy = min_energy.min(cell.energy);
        max_energy = max_energy.max(cell.energy);
    }
    let total = cells.len() as u64;
    let starter_ancestry = starter_counts
        .iter()
        .enumerate()
        .filter(|(_, count)| **count > 0)
        .map(|(starter, count)| StarterAncestrySummary {
            starter_type: starter as u8,
            label: starter_label(starter as u8).to_string(),
            count: *count,
            fraction: fraction(*count, total),
            avg_energy: starter_energy_sums[starter] / *count as f64,
            top_count: starter_zone_counts[starter][0],
            middle_count: starter_zone_counts[starter][1],
            deep_count: starter_zone_counts[starter][2],
        })
        .collect();
    CellSummary {
        tick,
        count: cells.len(),
        avg_energy: if cells.is_empty() {
            0.0
        } else {
            sum / cells.len() as f64
        },
        min_energy: if cells.is_empty() { 0.0 } else { min_energy },
        max_energy: if cells.is_empty() { 0.0 } else { max_energy },
        starter_counts,
        starter_ancestry,
    }
}

fn summarize_field(tick: u64, meta: &RunMeta, bytes: &[u8]) -> AnalysisResult<SnapshotChemistry> {
    if !bytes.len().is_multiple_of(4) {
        return Err("field byte length is not divisible by four".into());
    }
    let x = meta.grid_x as usize;
    let y = meta.grid_y as usize;
    let z = meta.grid_z as usize;
    let s_ext = meta.s_ext as usize;
    let expected_floats = x
        .checked_mul(y)
        .and_then(|v| v.checked_mul(z))
        .and_then(|v| v.checked_mul(s_ext))
        .ok_or("field dimensions overflow")?;
    if bytes.len() / 4 != expected_floats {
        return Err(format!(
            "field has {} floats, expected {expected_floats}",
            bytes.len() / 4
        )
        .into());
    }

    let tracked_species = [0usize, 1, 2, 3, 4];
    let mut sums = vec![vec![0.0; z]; tracked_species.len()];
    let mut totals = vec![0.0; tracked_species.len()];
    let mut max_values = vec![f64::NEG_INFINITY; tracked_species.len()];
    let mut nonfinite_values = 0;
    let mut negative_values = 0;
    for (i, chunk) in bytes.chunks_exact(4).enumerate() {
        let value = f32::from_le_bytes(chunk.try_into().expect("chunk size checked"));
        if !value.is_finite() {
            nonfinite_values += 1;
            continue;
        }
        if value < 0.0 {
            negative_values += 1;
        }
        let species = i % s_ext;
        if let Some(species_index) = tracked_species.iter().position(|s| *s == species) {
            let voxel = i / s_ext;
            let z_index = voxel / (x * y);
            let value = value as f64;
            sums[species_index][z_index] += value;
            totals[species_index] += value;
            max_values[species_index] = max_values[species_index].max(value);
        }
    }

    let layer_voxels = (x * y) as f64;
    let profiles: Vec<SpeciesProfile> = tracked_species
        .iter()
        .enumerate()
        .filter(|(_, species)| **species < s_ext)
        .map(|(i, species)| {
            let per_z_mean: Vec<f64> = sums[i].iter().map(|sum| *sum / layer_voxels).collect();
            let mid = per_z_mean.len() / 2;
            SpeciesProfile {
                species: *species,
                name: species_name(*species).to_string(),
                total_concentration: totals[i],
                max_value: if max_values[i].is_finite() {
                    max_values[i]
                } else {
                    0.0
                },
                surface_mean: per_z_mean.first().copied().unwrap_or(0.0),
                middle_mean: per_z_mean.get(mid).copied().unwrap_or(0.0),
                deep_mean: per_z_mean.last().copied().unwrap_or(0.0),
                per_z_mean,
            }
        })
        .collect();

    let oxidant = profiles.iter().find(|profile| profile.species == 1);
    let reductant = profiles.iter().find(|profile| profile.species == 2);
    let oxidant_penetration_z =
        oxidant.and_then(|profile| profile.per_z_mean.iter().position(|value| *value < 0.01));
    let reductant_penetration_z =
        reductant.and_then(|profile| profile.per_z_mean.iter().rposition(|value| *value < 0.01));
    let redox_overlap_layers = match (oxidant, reductant) {
        (Some(ox), Some(red)) => ox
            .per_z_mean
            .iter()
            .zip(red.per_z_mean.iter())
            .filter(|(ox, red)| **ox > 0.01 && **red > 0.01)
            .count(),
        _ => 0,
    };

    Ok(SnapshotChemistry {
        tick,
        nonfinite_values,
        negative_values,
        oxidant_penetration_z,
        reductant_penetration_z,
        redox_overlap_layers,
        species_profiles: profiles,
    })
}

fn full_rulesets_enabled(run_dir: &Path) -> AnalysisResult<bool> {
    let meta_json: serde_json::Value =
        serde_json::from_slice(&fs::read(run_dir.join("run_meta.json"))?)?;
    let mode = meta_json
        .get("ruleset_output_mode")
        .and_then(|value| value.as_str())
        .unwrap_or("off");
    Ok(mode == "full" || mode == "both")
}

fn summarize_rulesets(run_dir: &Path, tick: u64, meta: &RunMeta) -> AnalysisResult<RulesetSummary> {
    let meta_json: serde_json::Value =
        serde_json::from_slice(&fs::read(run_dir.join("run_meta.json"))?)?;
    let mode = meta_json
        .get("ruleset_output_mode")
        .and_then(|value| value.as_str())
        .unwrap_or("off");
    if mode != "full" && mode != "both" {
        return Err("full ruleset dumps were not enabled for this run".into());
    }
    let pattern = meta_json
        .get("ruleset_full_file_pattern")
        .and_then(|value| value.as_str())
        .ok_or("run_meta.json is missing ruleset_full_file_pattern")?;
    let path = snapshot_path(run_dir, pattern, tick)?;
    let capacity = u64::from(meta.grid_x)
        .checked_mul(u64::from(meta.grid_y))
        .and_then(|v| v.checked_mul(u64::from(meta.grid_z)))
        .ok_or("ruleset byte limit overflow")?;
    let byte_limit = u64::from(RULESET_FULL_HEADER_SIZE)
        .checked_add(capacity * u64::from(RULESET_FULL_CANONICAL_SIZE))
        .and_then(|v| v.checked_add(capacity * u64::from(RULESET_FULL_CELL_REF_STRIDE)))
        .ok_or("ruleset byte limit overflow")?;
    let bytes = read_binary_payload_bounded(&path, byte_limit)?;
    parse_ruleset_summary(tick, &bytes)
}

fn parse_ruleset_summary(tick: u64, bytes: &[u8]) -> AnalysisResult<RulesetSummary> {
    if bytes.len() < RULESET_FULL_HEADER_SIZE as usize {
        return Err("ruleset file is too small for header".into());
    }
    if bytes[0..4] != RULESET_FULL_MAGIC {
        return Err("ruleset magic mismatch".into());
    }
    let version = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
    let flags = u32::from_le_bytes(bytes[8..12].try_into().unwrap());
    let unique_count = u32::from_le_bytes(bytes[12..16].try_into().unwrap());
    let cell_count = u32::from_le_bytes(bytes[16..20].try_into().unwrap());
    let ruleset_size = u32::from_le_bytes(bytes[20..24].try_into().unwrap());
    let transport_stride = match (version, ruleset_size) {
        (RULESET_FULL_FORMAT_VERSION, RULESET_FULL_CANONICAL_SIZE) => TRANSPORT_V2_STRIDE,
        (RULESET_FULL_FORMAT_VERSION_V1, RULESET_FULL_CANONICAL_SIZE_V1) => TRANSPORT_V1_STRIDE,
        _ => {
            return Err(format!(
                "ruleset version/size mismatch: version {version}, size {ruleset_size}"
            )
            .into());
        }
    };
    if usize::try_from(ruleset_size)
        .ok()
        .is_none_or(|size| size < RECEPTOR_PAYLOAD_BYTES + TRANSPORTER_COUNT * transport_stride)
    {
        return Err(
            format!("ruleset size too small for transporter payload: {ruleset_size}").into(),
        );
    }
    if flags != 0 {
        return Err(format!("ruleset flags must be zero, got {flags}").into());
    }

    let refs_off =
        RULESET_FULL_HEADER_SIZE as usize + unique_count as usize * ruleset_size as usize;
    let expected_len = refs_off + cell_count as usize * RULESET_FULL_CELL_REF_STRIDE as usize;
    if bytes.len() != expected_len {
        return Err(format!(
            "ruleset file length {}, expected {expected_len}",
            bytes.len()
        )
        .into());
    }

    let mut counts: HashMap<u32, u32> = HashMap::new();
    for index in 0..cell_count as usize {
        let off = refs_off + index * RULESET_FULL_CELL_REF_STRIDE as usize;
        let dict_id = u32::from_le_bytes(bytes[off + 6..off + 10].try_into().unwrap());
        if dict_id >= unique_count {
            return Err(format!("cell ref {index} has out-of-range dict id {dict_id}").into());
        }
        *counts.entry(dict_id).or_insert(0) += 1;
    }
    let (dominant_dict_id, dominant_count) = counts
        .iter()
        .max_by_key(|(_, count)| **count)
        .map(|(dict_id, count)| (Some(*dict_id), *count))
        .unwrap_or((None, 0));
    let shannon_diversity = if cell_count > 0 {
        counts
            .values()
            .map(|count| {
                let p = *count as f64 / cell_count as f64;
                -p * p.log(E)
            })
            .sum()
    } else {
        0.0
    };
    let transporters = summarize_transporters_from_rulesets(
        bytes,
        ruleset_size as usize,
        unique_count,
        &counts,
        dominant_dict_id,
        transport_stride,
    )?;
    Ok(RulesetSummary {
        tick,
        unique_count,
        cell_count,
        dominant_dict_id,
        dominant_count,
        dominant_fraction: fraction(u64::from(dominant_count), u64::from(cell_count)),
        shannon_diversity,
        transporters,
    })
}

fn summarize_transporters_from_rulesets(
    bytes: &[u8],
    ruleset_size: usize,
    unique_count: u32,
    cell_counts: &HashMap<u32, u32>,
    dominant_dict_id: Option<u32>,
    transport_stride: usize,
) -> AnalysisResult<TransporterSummary> {
    let dict_off = RULESET_FULL_HEADER_SIZE as usize;
    let mut active_slots = 0u64;
    let mut uptake_dominant_slots = 0u64;
    let mut secretion_dominant_slots = 0u64;
    let mut bidirectional_slots = 0u64;
    let mut gated_slots = 0u64;
    let mut abs_gate_weight_sum = 0.0f64;
    let mut uptake_rate_sum = 0.0f64;
    let mut secretion_rate_sum = 0.0f64;
    let mut pair_acc: HashMap<(u8, u8), TransportPairAccumulator> = HashMap::new();
    let mut dominant_genotype = None;

    for dict_id in 0..unique_count {
        let count = u64::from(*cell_counts.get(&dict_id).unwrap_or(&0));
        let off = dict_off + dict_id as usize * ruleset_size;
        let payload = &bytes[off..off + ruleset_size];
        let (genotype, transporters) = summarize_genotype_transport(payload, transport_stride)?;
        if dominant_dict_id == Some(dict_id) {
            dominant_genotype = Some(GenotypeTransportSummary {
                dict_id,
                active_slots: genotype.active_slots,
                gated_slots: genotype.gated_slots,
                uptake_dominant_slots: genotype.uptake_dominant_slots,
                secretion_dominant_slots: genotype.secretion_dominant_slots,
            });
        }
        if count == 0 {
            continue;
        }
        for tp in transporters {
            let uptake = tp.uptake_rate.max(0.0);
            let secretion = tp.secrete_rate.max(0.0);
            let is_active =
                uptake > ACTIVE_TRANSPORT_THRESHOLD || secretion > ACTIVE_TRANSPORT_THRESHOLD;
            if !is_active {
                continue;
            }
            active_slots += count;
            uptake_rate_sum += f64::from(uptake) * count as f64;
            secretion_rate_sum += f64::from(secretion) * count as f64;
            if uptake > ACTIVE_TRANSPORT_THRESHOLD && secretion > ACTIVE_TRANSPORT_THRESHOLD {
                bidirectional_slots += count;
            }
            if uptake > secretion + ACTIVE_TRANSPORT_THRESHOLD {
                uptake_dominant_slots += count;
            } else if secretion > uptake + ACTIVE_TRANSPORT_THRESHOLD {
                secretion_dominant_slots += count;
            }
            if tp.gate_weight.abs() > ACTIVE_TRANSPORT_THRESHOLD {
                gated_slots += count;
                abs_gate_weight_sum += f64::from(tp.gate_weight.abs()) * count as f64;
            }
            let pair = pair_acc
                .entry((tp.ext_species, tp.int_species))
                .or_default();
            pair.active_slots += count;
            pair.uptake_rate_sum += f64::from(uptake) * count as f64;
            pair.secrete_rate_sum += f64::from(secretion) * count as f64;
            if tp.gate_weight.abs() > ACTIVE_TRANSPORT_THRESHOLD {
                pair.gated_slots += count;
            }
        }
    }

    let mut common_pairs: Vec<_> = pair_acc
        .into_iter()
        .map(|((ext_species, int_species), acc)| TransportPairSummary {
            ext_species,
            int_species,
            active_slots: acc.active_slots,
            avg_uptake_rate: if acc.active_slots > 0 {
                acc.uptake_rate_sum / acc.active_slots as f64
            } else {
                0.0
            },
            avg_secrete_rate: if acc.active_slots > 0 {
                acc.secrete_rate_sum / acc.active_slots as f64
            } else {
                0.0
            },
            gated_slots: acc.gated_slots,
        })
        .collect();
    common_pairs.sort_by(|a, b| {
        b.active_slots
            .cmp(&a.active_slots)
            .then_with(|| a.ext_species.cmp(&b.ext_species))
            .then_with(|| a.int_species.cmp(&b.int_species))
    });
    common_pairs.truncate(MAX_COMMON_TRANSPORT_PAIRS);

    let cell_count: u64 = cell_counts.values().map(|count| u64::from(*count)).sum();
    Ok(TransporterSummary {
        active_slots,
        avg_active_per_cell: fraction(active_slots, cell_count),
        uptake_dominant_slots,
        secretion_dominant_slots,
        bidirectional_slots,
        gated_slots,
        avg_abs_gate_weight: if gated_slots > 0 {
            abs_gate_weight_sum / gated_slots as f64
        } else {
            0.0
        },
        uptake_rate_sum,
        secretion_rate_sum,
        common_pairs,
        dominant_genotype,
    })
}

fn summarize_transporter_pair_timeline(
    ruleset_timeline: &[RulesetSummary],
) -> Vec<TransportPairTimeline> {
    if ruleset_timeline.is_empty() {
        return Vec::new();
    }

    let mut totals: HashMap<(u8, u8), u64> = HashMap::new();
    for rulesets in ruleset_timeline {
        for pair in &rulesets.transporters.common_pairs {
            *totals
                .entry((pair.ext_species, pair.int_species))
                .or_insert(0) += pair.active_slots;
        }
    }

    let mut selected = Vec::new();
    if let Some(latest) = ruleset_timeline.last() {
        for pair in &latest.transporters.common_pairs {
            push_unique_pair(&mut selected, (pair.ext_species, pair.int_species));
        }
    }
    for rulesets in ruleset_timeline {
        if let Some(pair) = rulesets.transporters.common_pairs.first() {
            push_unique_pair(&mut selected, (pair.ext_species, pair.int_species));
        }
    }

    let mut by_total = totals
        .iter()
        .map(|(pair, total)| (*pair, *total))
        .collect::<Vec<_>>();
    by_total.sort_by(|((a_ext, a_int), a_total), ((b_ext, b_int), b_total)| {
        b_total
            .cmp(a_total)
            .then_with(|| a_ext.cmp(b_ext))
            .then_with(|| a_int.cmp(b_int))
    });
    for (pair, _) in by_total {
        push_unique_pair(&mut selected, pair);
        if selected.len() >= MAX_COMMON_TRANSPORT_PAIRS {
            break;
        }
    }

    selected
        .into_iter()
        .take(MAX_COMMON_TRANSPORT_PAIRS)
        .map(|(ext_species, int_species)| {
            let points = ruleset_timeline
                .iter()
                .map(|rulesets| {
                    rulesets
                        .transporters
                        .common_pairs
                        .iter()
                        .find(|pair| {
                            pair.ext_species == ext_species && pair.int_species == int_species
                        })
                        .map(|pair| TransportPairTimelinePoint {
                            tick: rulesets.tick,
                            active_slots: pair.active_slots,
                            avg_uptake_rate: pair.avg_uptake_rate,
                            avg_secrete_rate: pair.avg_secrete_rate,
                            gated_slots: pair.gated_slots,
                        })
                        .unwrap_or(TransportPairTimelinePoint {
                            tick: rulesets.tick,
                            active_slots: 0,
                            avg_uptake_rate: 0.0,
                            avg_secrete_rate: 0.0,
                            gated_slots: 0,
                        })
                })
                .collect::<Vec<_>>();
            TransportPairTimeline {
                ext_species,
                int_species,
                total_active_slots: *totals.get(&(ext_species, int_species)).unwrap_or(&0),
                first_tick: points
                    .iter()
                    .find(|point| point.active_slots > 0)
                    .map(|point| point.tick),
                latest_tick: points
                    .iter()
                    .rev()
                    .find(|point| point.active_slots > 0)
                    .map(|point| point.tick),
                points,
            }
        })
        .collect()
}

fn push_unique_pair(pairs: &mut Vec<(u8, u8)>, pair: (u8, u8)) {
    if pairs.len() < MAX_COMMON_TRANSPORT_PAIRS && !pairs.contains(&pair) {
        pairs.push(pair);
    }
}

fn summarize_genotype_transport(
    payload: &[u8],
    transport_stride: usize,
) -> AnalysisResult<(
    GenotypeTransportStats,
    [ParsedTransporter; TRANSPORTER_COUNT],
)> {
    let transporters = parse_transporters(payload, transport_stride)?;
    let mut stats = GenotypeTransportStats::default();
    for tp in &transporters {
        let uptake = tp.uptake_rate.max(0.0);
        let secretion = tp.secrete_rate.max(0.0);
        if uptake > ACTIVE_TRANSPORT_THRESHOLD || secretion > ACTIVE_TRANSPORT_THRESHOLD {
            stats.active_slots += 1;
            if tp.gate_weight.abs() > ACTIVE_TRANSPORT_THRESHOLD {
                stats.gated_slots += 1;
            }
            if uptake > secretion + ACTIVE_TRANSPORT_THRESHOLD {
                stats.uptake_dominant_slots += 1;
            } else if secretion > uptake + ACTIVE_TRANSPORT_THRESHOLD {
                stats.secretion_dominant_slots += 1;
            }
        }
    }
    Ok((stats, transporters))
}

fn parse_transporters(
    payload: &[u8],
    transport_stride: usize,
) -> AnalysisResult<[ParsedTransporter; TRANSPORTER_COUNT]> {
    let mut out = [ParsedTransporter {
        uptake_rate: 0.0,
        secrete_rate: 0.0,
        ext_species: 0,
        int_species: 0,
        gate_weight: 0.0,
    }; TRANSPORTER_COUNT];
    let transport_off = RECEPTOR_PAYLOAD_BYTES;
    for (i, slot) in out.iter_mut().enumerate() {
        let off = transport_off + i * transport_stride;
        let end = off + transport_stride;
        if end > payload.len() {
            return Err("ruleset transporter payload truncated".into());
        }
        let uptake_rate = f32::from_le_bytes(payload[off..off + 4].try_into().unwrap());
        let secrete_rate = f32::from_le_bytes(payload[off + 4..off + 8].try_into().unwrap());
        let gate_weight = if transport_stride >= TRANSPORT_V2_STRIDE {
            f32::from_le_bytes(payload[off + 11..off + 15].try_into().unwrap())
        } else {
            0.0
        };
        *slot = ParsedTransporter {
            uptake_rate: finite_rate_for_analysis(uptake_rate),
            secrete_rate: finite_rate_for_analysis(secrete_rate),
            ext_species: payload[off + 8],
            int_species: payload[off + 9],
            gate_weight: if gate_weight.is_finite() {
                gate_weight
            } else {
                0.0
            },
        };
    }
    Ok(out)
}

fn finite_rate_for_analysis(value: f32) -> f32 {
    if value.is_finite() { value } else { 0.0 }
}

fn classify_run_findings(analysis: &RunAnalysis) -> Vec<Finding> {
    let mut findings = Vec::new();
    if let Some(traj) = &analysis.trajectory {
        if traj.growth_factor >= 5.0 {
            findings.push(finding(
                "surface_bloom_candidate",
                FindingLevel::Interesting,
                "Strong population expansion",
                format!(
                    "Population grew {:.1}x from {} to {} cells.",
                    traj.growth_factor, traj.initial_population, traj.final_population
                ),
                vec![format!("max population {}", traj.max_population)],
            ));
        } else if traj.final_population * 2 < traj.initial_population {
            findings.push(finding(
                "population_collapse",
                FindingLevel::Warning,
                "Population collapse or remnant state",
                format!(
                    "Population fell from {} to {} cells.",
                    traj.initial_population, traj.final_population
                ),
                vec![format!("total deaths {}", traj.total_deaths)],
            ));
        }
        if traj.total_divisions > 0 && traj.total_deaths as f64 / traj.total_divisions as f64 > 1.0
        {
            findings.push(finding(
                "high_turnover",
                FindingLevel::Interesting,
                "High turnover ecology",
                "Deaths exceeded divisions across the tick log.".to_string(),
                vec![format!(
                    "deaths/divisions = {:.2}",
                    traj.total_deaths as f64 / traj.total_divisions as f64
                )],
            ));
        }
    }
    if let Some(zonation) = &analysis.zonation {
        if zonation.middle_count == 0 && zonation.top_count + zonation.deep_count > 0 {
            findings.push(finding(
                "empty_middle",
                FindingLevel::Interesting,
                "Empty middle zone",
                "Living cells avoid the middle third while occupying surface and/or deep layers."
                    .to_string(),
                vec![
                    format!("top {}", zonation.top_count),
                    format!("middle {}", zonation.middle_count),
                    format!("deep {}", zonation.deep_count),
                ],
            ));
        }
        if zonation.top_fraction >= 0.8 {
            findings.push(finding(
                "surface_dominance",
                FindingLevel::Interesting,
                "Surface-dominated population",
                format!(
                    "{:.1}% of living cells are in the surface third.",
                    zonation.top_fraction * 100.0
                ),
                Vec::new(),
            ));
        }
        if zonation.deep_fraction >= 0.5 {
            findings.push(finding(
                "deep_dominance",
                FindingLevel::Interesting,
                "Deep-zone population",
                format!(
                    "{:.1}% of living cells are in the deep third.",
                    zonation.deep_fraction * 100.0
                ),
                Vec::new(),
            ));
        }
    }
    if let Some(latest_chemistry) = analysis.chemistry.last() {
        if latest_chemistry.nonfinite_values > 0 {
            findings.push(finding(
                "nonfinite_field_values",
                FindingLevel::Warning,
                "Non-finite chemistry values",
                "The latest inspected field contains non-finite floats.".to_string(),
                vec![format!("count {}", latest_chemistry.nonfinite_values)],
            ));
        }
        if latest_chemistry.redox_overlap_layers > 0 {
            findings.push(finding(
                "redox_overlap",
                FindingLevel::Info,
                "Oxidant and reductant overlap",
                format!(
                    "The latest inspected field has {} z-layers where both oxidant and reductant exceed 0.01.",
                    latest_chemistry.redox_overlap_layers
                ),
                Vec::new(),
            ));
        }
        if let Some(species0) = latest_chemistry
            .species_profiles
            .iter()
            .find(|profile| profile.species == 0)
            && species0.total_concentration > 1.0
        {
            findings.push(finding(
                "external_species0_pool",
                FindingLevel::Info,
                "External species 0 accumulated",
                "The latest inspected field contains a nontrivial pool of external species 0."
                    .to_string(),
                vec![
                    format!("total {:.3}", species0.total_concentration),
                    format!("max {:.3}", species0.max_value),
                ],
            ));
        }
    }
    if let Some(rulesets) = &analysis.rulesets {
        if rulesets.cell_count > 0 && rulesets.dominant_fraction < 0.25 {
            findings.push(finding(
                "high_genotype_diversity",
                FindingLevel::Interesting,
                "Diffuse genotype landscape",
                "No single ruleset dominates the living population.".to_string(),
                vec![format!(
                    "dominant genotype {:.1}%",
                    rulesets.dominant_fraction * 100.0
                )],
            ));
        } else if rulesets.dominant_fraction > 0.75 && rulesets.cell_count > 5 {
            findings.push(finding(
                "dominant_clone",
                FindingLevel::Interesting,
                "Dominant genotype",
                "A single ruleset accounts for most living cells.".to_string(),
                vec![format!(
                    "dominant genotype {:.1}%",
                    rulesets.dominant_fraction * 100.0
                )],
            ));
        }
    }
    findings
}

fn classify_comparison_findings(runs: &[RunComparisonEntry]) -> Vec<Finding> {
    let mut findings = Vec::new();
    let best = runs
        .iter()
        .filter_map(|run| run.final_population.map(|pop| (run, pop)))
        .max_by_key(|(_, pop)| *pop);
    let worst = runs
        .iter()
        .filter_map(|run| run.final_population.map(|pop| (run, pop)))
        .min_by_key(|(_, pop)| *pop);
    if let (Some((best, best_pop)), Some((worst, worst_pop))) = (best, worst)
        && best.name != worst.name
        && best_pop > worst_pop.saturating_mul(2)
    {
        findings.push(finding(
            "population_outcome_split",
            FindingLevel::Interesting,
            "Runs separate strongly by final population",
            format!(
                "{} ended at {} cells while {} ended at {} cells.",
                best.name, best_pop, worst.name, worst_pop
            ),
            Vec::new(),
        ));
    }
    findings
}

impl RunComparisonEntry {
    fn from_analysis(analysis: &RunAnalysis) -> Self {
        Self {
            name: run_name(&analysis.run_dir),
            run_dir: analysis.run_dir.clone(),
            final_population: analysis
                .trajectory
                .as_ref()
                .map(|trajectory| trajectory.final_population),
            max_population: analysis
                .trajectory
                .as_ref()
                .map(|trajectory| trajectory.max_population),
            growth_factor: analysis
                .trajectory
                .as_ref()
                .map(|trajectory| trajectory.growth_factor),
            final_avg_energy: analysis
                .trajectory
                .as_ref()
                .map(|trajectory| trajectory.final_avg_energy),
            top_fraction: analysis
                .zonation
                .as_ref()
                .map(|zonation| zonation.top_fraction),
            middle_fraction: analysis
                .zonation
                .as_ref()
                .map(|zonation| zonation.middle_fraction),
            deep_fraction: analysis
                .zonation
                .as_ref()
                .map(|zonation| zonation.deep_fraction),
            genotype_unique_count: analysis
                .rulesets
                .as_ref()
                .map(|rulesets| rulesets.unique_count),
            genotype_dominant_fraction: analysis
                .rulesets
                .as_ref()
                .map(|rulesets| rulesets.dominant_fraction),
            transporter_active_per_cell: analysis
                .rulesets
                .as_ref()
                .map(|rulesets| rulesets.transporters.avg_active_per_cell),
            transporter_gated_slots: analysis
                .rulesets
                .as_ref()
                .map(|rulesets| rulesets.transporters.gated_slots),
        }
    }
}

fn finding(
    id: &str,
    level: FindingLevel,
    title: &str,
    narrative: String,
    evidence: Vec<String>,
) -> Finding {
    Finding {
        id: id.to_string(),
        level,
        title: title.to_string(),
        narrative,
        evidence,
    }
}

fn fraction(part: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        part as f64 / total as f64
    }
}

fn species_name(species: usize) -> &'static str {
    match species {
        0 => "free_energy",
        1 => "oxidant",
        2 => "reductant",
        3 => "carbon",
        4 => "organic",
        _ => "unknown",
    }
}

fn starter_label(starter_type: u8) -> &'static str {
    match starter_type {
        0 => "phototroph",
        1 => "chemolithotroph",
        2 => "anaerobe",
        _ => "other",
    }
}

fn run_name(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_else(|| path.to_str().unwrap_or("run"))
        .to_string()
}

fn display_opt_usize(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "n/a".to_string())
}

fn display_opt_u64(value: Option<u64>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "n/a".to_string())
}

fn display_opt_f64(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.3}"))
        .unwrap_or_else(|| "n/a".to_string())
}

fn display_opt_pct(value: Option<f64>) -> String {
    value
        .map(|value| format!("{:.1}", value * 100.0))
        .unwrap_or_else(|| "n/a".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sampled_ticks_choose_first_middle_last() {
        assert_eq!(sample_ticks(&[]), Vec::<u64>::new());
        assert_eq!(sample_ticks(&[10]), vec![10]);
        assert_eq!(sample_ticks(&[0, 80, 159]), vec![0, 80, 159]);
        assert_eq!(sample_ticks(&[0, 10, 20, 30]), vec![0, 20, 30]);
    }

    #[test]
    fn ticks_csv_parser_reads_trajectory_and_z_counts() {
        let rows = parse_ticks_csv(
            "tick,population,avg_energy,avg_enzyme_a,avg_enzyme_b,avg_active_rxns,divisions_this_tick,deaths_this_tick,z0_cells,z1_cells\n\
             0,2,1.0,0,0,5,0,0,1,1\n\
             1,3,1.5,0,0,5,1,0,2,1\n",
        )
        .unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[1].population, 3);
        assert_eq!(rows[1].z_counts, vec![2, 1]);
        let summary = summarize_trajectory(&rows).unwrap();
        assert_eq!(summary.total_divisions, 1);
        assert_eq!(summary.growth_factor, 1.5);
    }

    #[test]
    fn field_summary_aggregates_depth_profiles() {
        let meta = RunMeta::new(1, 1, 2, 5, 0, true, false);
        let mut floats = Vec::new();
        for z in 0..2 {
            for species in 0..5 {
                floats.push((z * 10 + species) as f32);
            }
        }
        let bytes: Vec<u8> = floats
            .iter()
            .flat_map(|value| value.to_le_bytes())
            .collect();
        let summary = summarize_field(7, &meta, &bytes).unwrap();
        let oxidant = summary
            .species_profiles
            .iter()
            .find(|profile| profile.species == 1)
            .unwrap();
        let species0 = summary
            .species_profiles
            .iter()
            .find(|profile| profile.species == 0)
            .unwrap();
        assert_eq!(summary.tick, 7);
        assert_eq!(species0.per_z_mean, vec![0.0, 10.0]);
        assert_eq!(species0.total_concentration, 10.0);
        assert_eq!(species0.max_value, 10.0);
        assert_eq!(oxidant.per_z_mean, vec![1.0, 11.0]);
    }

    #[test]
    fn cell_summary_reports_starter_ancestry_by_zone_and_energy() {
        let cells = vec![
            LoadedCell {
                pos: [0, 0, 0],
                lineage_id: 1,
                starter_type: 0,
                energy: 1.0,
            },
            LoadedCell {
                pos: [0, 0, 5],
                lineage_id: 2,
                starter_type: 0,
                energy: 3.0,
            },
            LoadedCell {
                pos: [0, 0, 8],
                lineage_id: 3,
                starter_type: 2,
                energy: 2.0,
            },
        ];

        let summary = summarize_cells(42, 9, &cells);

        assert_eq!(summary.starter_counts, [2, 0, 1, 0]);
        let photo = summary
            .starter_ancestry
            .iter()
            .find(|starter| starter.starter_type == 0)
            .unwrap();
        assert_eq!(photo.label, "phototroph");
        assert_eq!(photo.count, 2);
        assert!((photo.avg_energy - 2.0).abs() < 1e-9);
        assert_eq!(
            (photo.top_count, photo.middle_count, photo.deep_count),
            (1, 1, 0)
        );
        let anaerobe = summary
            .starter_ancestry
            .iter()
            .find(|starter| starter.starter_type == 2)
            .unwrap();
        assert_eq!(
            (
                anaerobe.top_count,
                anaerobe.middle_count,
                anaerobe.deep_count
            ),
            (0, 0, 1)
        );
    }

    #[test]
    fn ruleset_summary_parses_dictionary_histogram() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&RULESET_FULL_MAGIC);
        bytes.extend_from_slice(&RULESET_FULL_FORMAT_VERSION.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&2u32.to_le_bytes());
        bytes.extend_from_slice(&3u32.to_le_bytes());
        bytes.extend_from_slice(&RULESET_FULL_CANONICAL_SIZE.to_le_bytes());
        bytes.extend(vec![0; RULESET_FULL_CANONICAL_SIZE as usize * 2]);
        for (z, dict_id) in [(0u16, 0u32), (1u16, 1u32), (2u16, 1u32)] {
            bytes.extend_from_slice(&0u16.to_le_bytes());
            bytes.extend_from_slice(&0u16.to_le_bytes());
            bytes.extend_from_slice(&z.to_le_bytes());
            bytes.extend_from_slice(&dict_id.to_le_bytes());
        }
        let summary = parse_ruleset_summary(5, &bytes).unwrap();
        assert_eq!(summary.unique_count, 2);
        assert_eq!(summary.cell_count, 3);
        assert_eq!(summary.dominant_dict_id, Some(1));
        assert_eq!(summary.dominant_count, 2);
        assert_eq!(summary.transporters.active_slots, 0);
    }

    fn set_transport(
        payload: &mut [u8],
        slot: usize,
        uptake_rate: f32,
        secrete_rate: f32,
        ext_species: u8,
        int_species: u8,
        gate_receptor: u8,
        gate_weight: f32,
    ) {
        let off = RECEPTOR_PAYLOAD_BYTES + slot * TRANSPORT_V2_STRIDE;
        payload[off..off + 4].copy_from_slice(&uptake_rate.to_le_bytes());
        payload[off + 4..off + 8].copy_from_slice(&secrete_rate.to_le_bytes());
        payload[off + 8] = ext_species;
        payload[off + 9] = int_species;
        payload[off + 10] = gate_receptor;
        payload[off + 11..off + 15].copy_from_slice(&gate_weight.to_le_bytes());
    }

    fn set_v1_transport(
        payload: &mut [u8],
        slot: usize,
        uptake_rate: f32,
        secrete_rate: f32,
        ext_species: u8,
        int_species: u8,
    ) {
        let off = RECEPTOR_PAYLOAD_BYTES + slot * TRANSPORT_V1_STRIDE;
        payload[off..off + 4].copy_from_slice(&uptake_rate.to_le_bytes());
        payload[off + 4..off + 8].copy_from_slice(&secrete_rate.to_le_bytes());
        payload[off + 8] = ext_species;
        payload[off + 9] = int_species;
    }

    #[test]
    fn ruleset_summary_parses_v1_transporters_as_ungated() {
        let mut dict = vec![0u8; RULESET_FULL_CANONICAL_SIZE_V1 as usize];
        set_v1_transport(&mut dict, 0, 0.4, 0.1, 2, 3);

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&RULESET_FULL_MAGIC);
        bytes.extend_from_slice(&RULESET_FULL_FORMAT_VERSION_V1.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&1u32.to_le_bytes());
        bytes.extend_from_slice(&1u32.to_le_bytes());
        bytes.extend_from_slice(&RULESET_FULL_CANONICAL_SIZE_V1.to_le_bytes());
        bytes.extend_from_slice(&dict);
        bytes.extend_from_slice(&4u16.to_le_bytes());
        bytes.extend_from_slice(&5u16.to_le_bytes());
        bytes.extend_from_slice(&6u16.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());

        let summary = parse_ruleset_summary(9, &bytes).unwrap();

        assert_eq!(summary.unique_count, 1);
        assert_eq!(summary.cell_count, 1);
        assert_eq!(summary.transporters.active_slots, 1);
        assert_eq!(summary.transporters.gated_slots, 0);
        assert_eq!(summary.transporters.avg_abs_gate_weight, 0.0);
        assert_eq!(summary.transporters.uptake_dominant_slots, 1);
        let common = &summary.transporters.common_pairs[0];
        assert_eq!((common.ext_species, common.int_species), (2, 3));
        assert_eq!(common.gated_slots, 0);
        assert!((common.avg_uptake_rate - 0.4).abs() < 1e-6);
        assert!((common.avg_secrete_rate - 0.1).abs() < 1e-6);
    }

    #[test]
    fn ruleset_summary_reports_transporter_ecology() {
        let mut dict0 = vec![0u8; RULESET_FULL_CANONICAL_SIZE as usize];
        set_transport(&mut dict0, 0, 1.0, 0.0, 1, 2, 1, 0.5);
        set_transport(&mut dict0, 1, 0.0, 2.0, 3, 4, 0, 0.0);
        let mut dict1 = vec![0u8; RULESET_FULL_CANONICAL_SIZE as usize];
        set_transport(&mut dict1, 0, 0.2, 0.3, 1, 2, 2, -0.25);

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&RULESET_FULL_MAGIC);
        bytes.extend_from_slice(&RULESET_FULL_FORMAT_VERSION.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&2u32.to_le_bytes());
        bytes.extend_from_slice(&3u32.to_le_bytes());
        bytes.extend_from_slice(&RULESET_FULL_CANONICAL_SIZE.to_le_bytes());
        bytes.extend_from_slice(&dict0);
        bytes.extend_from_slice(&dict1);
        for (z, dict_id) in [(0u16, 0u32), (1u16, 1u32), (2u16, 1u32)] {
            bytes.extend_from_slice(&0u16.to_le_bytes());
            bytes.extend_from_slice(&0u16.to_le_bytes());
            bytes.extend_from_slice(&z.to_le_bytes());
            bytes.extend_from_slice(&dict_id.to_le_bytes());
        }

        let summary = parse_ruleset_summary(5, &bytes).unwrap();
        let transporters = &summary.transporters;

        assert_eq!(transporters.active_slots, 4);
        assert!((transporters.avg_active_per_cell - 4.0 / 3.0).abs() < 1e-9);
        assert_eq!(transporters.uptake_dominant_slots, 1);
        assert_eq!(transporters.secretion_dominant_slots, 3);
        assert_eq!(transporters.bidirectional_slots, 2);
        assert_eq!(transporters.gated_slots, 3);
        assert!((transporters.avg_abs_gate_weight - 1.0 / 3.0).abs() < 1e-6);
        let dominant = transporters.dominant_genotype.as_ref().unwrap();
        assert_eq!(dominant.dict_id, 1);
        assert_eq!(dominant.active_slots, 1);
        assert_eq!(dominant.gated_slots, 1);
        assert_eq!(dominant.secretion_dominant_slots, 1);

        let common = &transporters.common_pairs[0];
        assert_eq!((common.ext_species, common.int_species), (1, 2));
        assert_eq!(common.active_slots, 3);
        assert!((common.avg_uptake_rate - (1.4 / 3.0)).abs() < 1e-6);
        assert!((common.avg_secrete_rate - 0.2).abs() < 1e-6);
        assert_eq!(common.gated_slots, 3);
    }

    fn pair(ext_species: u8, int_species: u8, active_slots: u64) -> TransportPairSummary {
        TransportPairSummary {
            ext_species,
            int_species,
            active_slots,
            avg_uptake_rate: active_slots as f64 / 10.0,
            avg_secrete_rate: active_slots as f64 / 20.0,
            gated_slots: active_slots / 2,
        }
    }

    fn rulesets_with_pairs(tick: u64, common_pairs: Vec<TransportPairSummary>) -> RulesetSummary {
        let active_slots = common_pairs.iter().map(|pair| pair.active_slots).sum();
        let gated_slots = common_pairs.iter().map(|pair| pair.gated_slots).sum();
        RulesetSummary {
            tick,
            unique_count: 1,
            cell_count: 10,
            dominant_dict_id: Some(0),
            dominant_count: 10,
            dominant_fraction: 1.0,
            shannon_diversity: 0.0,
            transporters: TransporterSummary {
                active_slots,
                avg_active_per_cell: active_slots as f64 / 10.0,
                uptake_dominant_slots: active_slots,
                secretion_dominant_slots: 0,
                bidirectional_slots: 0,
                gated_slots,
                avg_abs_gate_weight: 0.5,
                uptake_rate_sum: 0.0,
                secretion_rate_sum: 0.0,
                common_pairs,
                dominant_genotype: None,
            },
        }
    }

    #[test]
    fn transporter_pair_timeline_keeps_latest_and_changing_leaders() {
        let timeline = vec![
            rulesets_with_pairs(0, vec![pair(1, 2, 9), pair(5, 6, 2)]),
            rulesets_with_pairs(10, vec![pair(3, 4, 8), pair(1, 2, 4)]),
            rulesets_with_pairs(20, vec![pair(7, 8, 6), pair(3, 4, 5)]),
        ];

        let pair_timeline = summarize_transporter_pair_timeline(&timeline);

        assert_eq!(
            pair_timeline
                .iter()
                .map(|pair| (pair.ext_species, pair.int_species))
                .collect::<Vec<_>>(),
            vec![(7, 8), (3, 4), (1, 2), (5, 6)]
        );
        assert_eq!(
            pair_timeline[0]
                .points
                .iter()
                .map(|point| (point.tick, point.active_slots))
                .collect::<Vec<_>>(),
            vec![(0, 0), (10, 0), (20, 6)]
        );
        assert_eq!(pair_timeline[2].first_tick, Some(0));
        assert_eq!(pair_timeline[2].latest_tick, Some(10));
        assert_eq!(pair_timeline[2].total_active_slots, 13);
    }

    #[test]
    fn markdown_renders_transporter_pair_timeline() {
        let ruleset_timeline = vec![
            rulesets_with_pairs(0, vec![pair(1, 2, 4)]),
            rulesets_with_pairs(10, vec![pair(3, 4, 6), pair(1, 2, 2)]),
        ];
        let analysis = RunAnalysis {
            run_dir: PathBuf::from("/tmp/example_run"),
            grid: [1, 1, 1],
            available_ticks: vec![0, 10],
            sampled_ticks: vec![0, 10],
            trajectory: None,
            zonation: None,
            cells: None,
            cell_timeline: Vec::new(),
            chemistry: Vec::new(),
            rulesets: ruleset_timeline.last().cloned(),
            transporter_pair_timeline: summarize_transporter_pair_timeline(&ruleset_timeline),
            ruleset_timeline,
            findings: Vec::new(),
            warnings: Vec::new(),
        };

        let markdown = render_run_markdown(&analysis);

        assert!(markdown.contains("### Transporter Pair Timeline"));
        assert!(markdown.contains("| tick | ext3->int4 | ext1->int2 |"));
        assert!(markdown.contains("| 0 | - | 4 (g2) |"));
        assert!(markdown.contains("| 10 | 6 (g3) | 2 (g1) |"));
    }

    #[test]
    fn report_writer_with_no_outputs_does_not_create_directory() {
        let out_dir = PathBuf::from("/tmp/marl_analysis_no_output_test");
        let _ = fs::remove_dir_all(&out_dir);
        let analysis = RunAnalysis {
            run_dir: PathBuf::from("/tmp/example_run"),
            grid: [1, 1, 1],
            available_ticks: Vec::new(),
            sampled_ticks: Vec::new(),
            trajectory: None,
            zonation: None,
            cells: None,
            cell_timeline: Vec::new(),
            chemistry: Vec::new(),
            rulesets: None,
            ruleset_timeline: Vec::new(),
            transporter_pair_timeline: Vec::new(),
            findings: Vec::new(),
            warnings: Vec::new(),
        };

        write_run_reports(&analysis, &out_dir, false, false).unwrap();
        assert!(!out_dir.exists());
    }
}
