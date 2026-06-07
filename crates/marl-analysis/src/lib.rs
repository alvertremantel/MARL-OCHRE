use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::f64::consts::E;
use std::fs;
use std::path::{Path, PathBuf};

use marl_config::stoich::STOICH_STARTER_TYPE_COUNT;
use marl_config::{ChemicalComposition, S_EXT, external_species_descriptor};
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

const STOICH_COMPACT_MERGE_REL_TOLERANCE: f64 = 1.0e-6;
const STOICH_EVENT_CSV_AMOUNT_DECIMALS: f64 = 1.0e-6;

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
    pub rng_seed: Option<u64>,
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
    pub stoich_v2: Option<StoichV2Analysis>,
    pub byproduct_calibration: Option<ByproductCalibrationSummary>,
    pub findings: Vec<Finding>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ComparisonAnalysis {
    pub runs: Vec<RunComparisonEntry>,
    pub paired_byproduct_population: Vec<ByproductPairedComparison>,
    pub byproduct_strength_aggregates: Vec<ByproductStrengthAggregate>,
    pub findings: Vec<Finding>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunComparisonEntry {
    pub name: String,
    pub run_dir: PathBuf,
    pub rng_seed: Option<u64>,
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
    pub stoich_total_events: Option<u64>,
    pub stoich_imbalanced_events: Option<u64>,
    pub reaction_byproduct_events: Option<u64>,
    pub reaction_byproduct_amount: Option<f64>,
    pub reaction_leakage_events: Option<u64>,
    pub reaction_leakage_energy_to_heat: Option<f64>,
    pub byproduct_final_pool: Option<f64>,
    pub byproduct_retained_fraction: Option<f64>,
    pub byproduct_signed_excess_final_pool: Option<f64>,
    pub byproduct_signed_excess_retained_fraction: Option<f64>,
    pub byproduct_excess_final_pool: Option<f64>,
    pub byproduct_excess_retained_fraction: Option<f64>,
    pub byproduct_cross_feeding_candidates: Option<u64>,
    pub byproduct_public_pool_candidates: Option<u64>,
    pub byproduct_adjusted_species: Vec<ByproductAdjustedSpeciesComparison>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ByproductAdjustedSpeciesComparison {
    pub species_index: i16,
    pub name: String,
    pub role: String,
    pub produced_amount: f64,
    pub baseline_field_total: f64,
    pub final_field_total: f64,
    pub signed_excess_final_pool: f64,
    pub clipped_excess_final_pool: f64,
    pub clipped_excess_retained_fraction: Option<f64>,
    pub adjusted_interpretation: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ByproductPairedComparison {
    pub rng_seed: Option<u64>,
    pub baseline_name: String,
    pub run_name: String,
    pub byproduct_amount: Option<f64>,
    pub delta_final_population: Option<i64>,
    pub final_population_ratio: Option<f64>,
    pub delta_growth_factor: Option<f64>,
    pub delta_final_avg_energy: Option<f64>,
    pub byproduct_signed_excess_retained_fraction: Option<f64>,
    pub byproduct_excess_retained_fraction: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ByproductStrengthAggregate {
    pub label: String,
    pub run_count: usize,
    pub mean_byproduct_amount: Option<f64>,
    pub mean_delta_final_population: Option<f64>,
    pub sd_delta_final_population: Option<f64>,
    pub min_delta_final_population: Option<f64>,
    pub max_delta_final_population: Option<f64>,
    pub mean_final_population_ratio: Option<f64>,
    pub mean_delta_growth_factor: Option<f64>,
    pub mean_delta_final_avg_energy: Option<f64>,
    pub mean_signed_excess_retained_fraction: Option<f64>,
    pub mean_clipped_excess_retained_fraction: Option<f64>,
    pub sd_clipped_excess_retained_fraction: Option<f64>,
    pub min_clipped_excess_retained_fraction: Option<f64>,
    pub max_clipped_excess_retained_fraction: Option<f64>,
    pub positive_population_delta_count: usize,
    pub positive_signed_excess_count: usize,
    pub clipped_excess_warning_count: usize,
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
    pub descriptor: SpeciesDescriptorProfile,
    pub total_concentration: f64,
    pub max_value: f64,
    pub surface_mean: f64,
    pub middle_mean: f64,
    pub deep_mean: f64,
    pub per_z_mean: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SpeciesDescriptorProfile {
    pub bond_energy: f32,
    pub extracellular_stability: f32,
    pub membrane_permeability: f32,
    pub storage_density: f32,
    pub work_coupling: f32,
    pub composition: ChemicalComposition,
}

impl From<marl_config::ExternalSpeciesDescriptor> for SpeciesDescriptorProfile {
    fn from(descriptor: marl_config::ExternalSpeciesDescriptor) -> Self {
        Self {
            bond_energy: descriptor.bond_energy,
            extracellular_stability: descriptor.extracellular_stability,
            membrane_permeability: descriptor.membrane_permeability,
            storage_density: descriptor.storage_density,
            work_coupling: descriptor.work_coupling,
            composition: descriptor.composition,
        }
    }
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
    pub external_species: Vec<TransportSpeciesSummary>,
    pub common_pairs: Vec<TransportPairSummary>,
    pub dominant_genotype: Option<GenotypeTransportSummary>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TransportSpeciesSummary {
    pub ext_species: u8,
    pub active_slots: u64,
    pub uptake_active_slots: u64,
    pub secretion_active_slots: u64,
    pub uptake_dominant_slots: u64,
    pub secretion_dominant_slots: u64,
    pub gated_slots: u64,
    pub uptake_rate_sum: f64,
    pub secretion_rate_sum: f64,
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
pub struct StoichV2Analysis {
    pub summary_present: bool,
    pub events_present: bool,
    pub total_ticks: Option<u64>,
    pub enforcement: Option<String>,
    pub gross_residual_abs_sum: Option<f64>,
    pub total_events: u64,
    pub imbalanced_events: u64,
    pub reaction_byproduct_events: u64,
    pub reaction_byproduct_amount: f64,
    pub reaction_byproduct_model_abs: f64,
    pub reaction_byproduct_residual_abs: f64,
    pub reaction_byproduct_by_species: Vec<StoichSpeciesEventSummary>,
    pub reaction_byproduct_by_species_starter: Vec<StoichSpeciesStarterEventSummary>,
    pub reaction_leakage_events: u64,
    pub reaction_leakage_amount: f64,
    pub reaction_leakage_energy_to_heat: f64,
    pub transport_flux_by_species: Vec<StoichTransportFluxSummary>,
    pub transport_flux_by_species_starter: Vec<StoichTransportStarterFluxSummary>,
}

#[derive(Debug, Clone, Serialize)]
pub struct StoichSpeciesEventSummary {
    pub species_index: i16,
    pub amount: f64,
    pub events: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct StoichSpeciesStarterEventSummary {
    pub species_index: i16,
    pub starter_type: u8,
    pub amount: f64,
    pub events: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct StoichTransportFluxSummary {
    pub species_index: i16,
    pub uptake_amount: f64,
    pub uptake_events: u64,
    pub secretion_amount: f64,
    pub secretion_events: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct StoichTransportStarterFluxSummary {
    pub species_index: i16,
    pub starter_type: u8,
    pub uptake_amount: f64,
    pub uptake_events: u64,
    pub secretion_amount: f64,
    pub secretion_events: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ByproductCalibrationSummary {
    pub total_byproduct_amount: f64,
    pub final_byproduct_pool: Option<f64>,
    pub retained_fraction: Option<f64>,
    pub field_tick: Option<u64>,
    pub ruleset_tick: Option<u64>,
    pub missing_field_species: u64,
    pub cross_feeding_candidates: u64,
    pub public_pool_candidates: u64,
    pub species: Vec<ByproductSpeciesCalibration>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ByproductSpeciesCalibration {
    pub species_index: i16,
    pub name: String,
    pub role: String,
    pub produced_amount: f64,
    pub final_field_total: Option<f64>,
    pub retained_fraction: Option<f64>,
    pub uptake_active_slots: u64,
    pub secretion_active_slots: u64,
    pub uptake_rate_sum: f64,
    pub secretion_rate_sum: f64,
    pub interpretation: String,
}

#[derive(Debug, Clone, Default)]
struct StoichEventKindAccumulator {
    events: u64,
    amount: f64,
    model_abs: f64,
    residual_abs: f64,
    reservoir_energy: f64,
}

#[derive(Debug, Clone, Default)]
struct StoichV2EventAccumulator {
    total_events: u64,
    imbalanced_events: u64,
    byproduct: StoichEventKindAccumulator,
    leakage: StoichEventKindAccumulator,
    byproduct_by_species: HashMap<i16, StoichSpeciesEventSummary>,
    byproduct_by_species_starter: HashMap<(i16, u8), StoichSpeciesStarterEventSummary>,
    transport_flux_by_species: HashMap<i16, StoichTransportFluxSummary>,
    transport_flux_by_species_starter: HashMap<(i16, u8), StoichTransportStarterFluxSummary>,
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
const BYPRODUCT_CANDIDATE_MIN_AMOUNT: f64 = 1e-6;
const BYPRODUCT_CANDIDATE_MIN_TOTAL_FRACTION: f64 = 1e-4;
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

#[derive(Debug, Clone, Default)]
struct TransportSpeciesAccumulator {
    active_slots: u64,
    uptake_active_slots: u64,
    secretion_active_slots: u64,
    uptake_dominant_slots: u64,
    secretion_dominant_slots: u64,
    gated_slots: u64,
    uptake_rate_sum: f64,
    secretion_rate_sum: f64,
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

    let full_rulesets_available = if cfg.include_rulesets && latest_tick.is_some() {
        full_rulesets_enabled(run_dir)?
    } else {
        false
    };
    let mut ruleset_timeline = Vec::new();
    if cfg.include_rulesets && latest_tick.is_some() {
        if full_rulesets_available {
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
    let stoich_v2 = match read_stoich_v2_analysis(run_dir) {
        Ok(summary) => summary,
        Err(err) => {
            warnings.push(format!("stoich v2 analysis skipped: {err}"));
            None
        }
    };
    let calibration_field = available_ticks.last().copied().and_then(|tick| {
        if let Some(summary) = chemistry.iter().find(|summary| summary.tick == tick) {
            Some(summary.clone())
        } else {
            match load_field_bytes(run_dir, tick, &meta)
                .and_then(|bytes| summarize_field(tick, &meta, &bytes))
            {
                Ok(summary) => Some(summary),
                Err(err) => {
                    warnings.push(format!(
                        "tick {tick}: byproduct calibration field analysis skipped: {err}"
                    ));
                    None
                }
            }
        }
    });
    let calibration_rulesets = if cfg.include_rulesets && full_rulesets_available {
        available_ticks.last().copied().and_then(|tick| {
            if let Some(summary) = ruleset_timeline.iter().find(|summary| summary.tick == tick) {
                Some(summary.clone())
            } else {
                match summarize_rulesets(run_dir, tick, &meta) {
                    Ok(summary) => Some(summary),
                    Err(err) => {
                        warnings.push(format!(
                            "tick {tick}: byproduct calibration ruleset analysis skipped: {err}"
                        ));
                        None
                    }
                }
            }
        })
    } else {
        None
    };
    let byproduct_calibration = summarize_byproduct_calibration(
        stoich_v2.as_ref(),
        calibration_field.as_ref(),
        calibration_rulesets.as_ref(),
    );

    let mut analysis = RunAnalysis {
        run_dir: run_dir.to_path_buf(),
        grid: [meta.grid_x, meta.grid_y, meta.grid_z],
        rng_seed: read_run_rng_seed(run_dir),
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
        stoich_v2,
        byproduct_calibration,
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
    let mut analyses = Vec::new();
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
                analyses.push(analysis);
            }
            Err(err) => warnings.push(format!("{}: analysis failed: {err}", run_name(run_dir))),
        }
    }
    let mut runs = analyses
        .iter()
        .map(RunComparisonEntry::from_analysis)
        .collect::<Vec<_>>();
    warnings.extend(apply_byproduct_baseline_adjustment(&analyses, &mut runs));
    let paired_byproduct_population = summarize_paired_byproduct_population(&runs);
    let byproduct_strength_aggregates =
        summarize_byproduct_strength_aggregates(&paired_byproduct_population);
    let findings = classify_comparison_findings(&runs);
    Ok(ComparisonAnalysis {
        runs,
        paired_byproduct_population,
        byproduct_strength_aggregates,
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
    if let Some(chemistry) = analysis.chemistry.last() {
        let roles = chemistry
            .species_profiles
            .iter()
            .map(|profile| {
                format!(
                    "ext{} {}:{} total={:.2}",
                    profile.species,
                    profile.name,
                    chemical_role(profile),
                    profile.total_concentration
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        out.push_str(&format!(
            "  chemistry roles at tick {}: {roles}\n",
            chemistry.tick
        ));
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
    if let Some(stoich) = &analysis.stoich_v2 {
        if stoich.events_present {
            out.push_str(&format!(
                "  stoich v2: events={}, imbalanced={}, byproduct {:.4} in {} events, leakage {:.4} to heat\n",
                stoich.total_events,
                stoich.imbalanced_events,
                stoich.reaction_byproduct_amount,
                stoich.reaction_byproduct_events,
                stoich.reaction_leakage_energy_to_heat
            ));
            if !stoich.reaction_byproduct_by_species.is_empty() {
                let species = stoich
                    .reaction_byproduct_by_species
                    .iter()
                    .map(|summary| {
                        format!(
                            "ext{}:{:.4}/{}",
                            summary.species_index, summary.amount, summary.events
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(", ");
                out.push_str(&format!("  reaction byproducts by species: {species}\n"));
            }
            if let Some(producers) =
                format_byproduct_starter_rows(&stoich.reaction_byproduct_by_species_starter, 5)
            {
                out.push_str(&format!("  reaction byproduct producers: {producers}\n"));
            }
            if let Some(fluxes) = format_transport_fluxes(&stoich.transport_flux_by_species, 5) {
                out.push_str(&format!("  transport flux by species: {fluxes}\n"));
            }
            if let Some(fluxes) =
                format_transport_starter_fluxes(&stoich.transport_flux_by_species_starter, 5)
            {
                out.push_str(&format!("  transport flux by starter: {fluxes}\n"));
            }
        } else {
            out.push_str(&format!(
                "  stoich v2: summary present={}, event-level reaction totals unavailable\n",
                stoich.summary_present
            ));
        }
    }
    if let Some(calibration) = &analysis.byproduct_calibration {
        out.push_str(&format!(
            "  byproduct calibration: produced={:.4}, final_pool={}, retained={}, field_tick={:?}, ruleset_tick={:?}, missing_field_species={}, uptake_pressure_candidates={}, public_pool_candidates={}\n",
            calibration.total_byproduct_amount,
            calibration
                .final_byproduct_pool
                .map(|value| format!("{value:.4}"))
                .unwrap_or_else(|| "n/a".to_string()),
            calibration
                .retained_fraction
                .map(|value| format!("{value:.3}"))
                .unwrap_or_else(|| "n/a".to_string()),
            calibration.field_tick,
            calibration.ruleset_tick,
            calibration.missing_field_species,
            calibration.cross_feeding_candidates,
            calibration.public_pool_candidates
        ));
        let species = calibration
            .species
            .iter()
            .take(4)
            .map(|species| {
                format!(
                    "ext{} {}:{} prod={:.4} pool={} uptake_slots={}",
                    species.species_index,
                    species.name,
                    species.interpretation,
                    species.produced_amount,
                    species
                        .final_field_total
                        .map(|value| format!("{value:.4}"))
                        .unwrap_or_else(|| "n/a".to_string()),
                    species.uptake_active_slots
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        out.push_str(&format!("  byproduct calibration species: {species}\n"));
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
            "  {}: rng_seed={}, final_pop={:?}, growth={:?}, top={:?}, dominant_genotype={:?}, active_transporters={:?}, byproduct={}, gross_retained={}, signed_excess_retained={}, clipped_excess_retained={}, leakage_heat={}\n",
            run.name,
            run.rng_seed
                .map(|seed| seed.to_string())
                .unwrap_or_else(|| "n/a".to_string()),
            run.final_population,
            run.growth_factor.map(|v| format!("{v:.2}x")),
            run.top_fraction.map(|v| format!("{:.1}%", v * 100.0)),
            run.genotype_dominant_fraction
                .map(|v| format!("{:.1}%", v * 100.0)),
            run.transporter_active_per_cell
                .map(|v| format!("{v:.2}/cell")),
            run.reaction_byproduct_amount
                .map(|value| format!("{value:.4}"))
                .unwrap_or_else(|| "n/a".to_string()),
            run.byproduct_retained_fraction
                .map(|value| format!("{value:.3}"))
                .unwrap_or_else(|| "n/a".to_string()),
            run.byproduct_signed_excess_retained_fraction
                .map(|value| format!("{value:.3}"))
                .unwrap_or_else(|| "n/a".to_string()),
            run.byproduct_excess_retained_fraction
                .map(|value| format!("{value:.3}"))
                .unwrap_or_else(|| "n/a".to_string()),
            run.reaction_leakage_energy_to_heat
                .map(|value| format!("{value:.4}"))
                .unwrap_or_else(|| "n/a".to_string())
        ));
    }
    if !analysis.paired_byproduct_population.is_empty() {
        out.push_str("  paired byproduct population effects:\n");
        for pair in &analysis.paired_byproduct_population {
            out.push_str(&format!(
                "    {} vs {}: rng_seed={}, byproduct={}, delta_pop={}, pop_ratio={}, delta_growth={}, delta_energy={}, signed_excess_retained={}, clipped_excess_retained={}\n",
                pair.run_name,
                pair.baseline_name,
                pair.rng_seed
                    .map(|seed| seed.to_string())
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.byproduct_amount
                    .map(|value| format!("{value:.4}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.delta_final_population
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.final_population_ratio
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.delta_growth_factor
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.delta_final_avg_energy
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.byproduct_signed_excess_retained_fraction
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.byproduct_excess_retained_fraction
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string())
            ));
        }
    }
    if !analysis.byproduct_strength_aggregates.is_empty() {
        out.push_str("  byproduct strength aggregates:\n");
        for aggregate in &analysis.byproduct_strength_aggregates {
            out.push_str(&format!(
                "    {}: n={}, mean_delta_pop={} sd={} range={}..{}, positive_pop={}/{}, mean_pop_ratio={}, mean_signed_excess_retained={}, positive_signed_excess={}/{}, mean_clipped_excess_retained={} sd={} range={}..{}, excess_warnings={}/{}\n",
                aggregate.label,
                aggregate.run_count,
                aggregate
                    .mean_delta_final_population
                    .map(|value| format!("{value:.1}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate
                    .sd_delta_final_population
                    .map(|value| format!("{value:.1}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate
                    .min_delta_final_population
                    .map(|value| format!("{value:.1}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate
                    .max_delta_final_population
                    .map(|value| format!("{value:.1}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate.positive_population_delta_count,
                aggregate.run_count,
                aggregate
                    .mean_final_population_ratio
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate
                    .mean_signed_excess_retained_fraction
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate.positive_signed_excess_count,
                aggregate.run_count,
                aggregate
                    .mean_clipped_excess_retained_fraction
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate
                    .sd_clipped_excess_retained_fraction
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate
                    .min_clipped_excess_retained_fraction
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate
                    .max_clipped_excess_retained_fraction
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                aggregate.clipped_excess_warning_count,
                aggregate.run_count
            ));
        }
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

fn format_transport_fluxes(fluxes: &[StoichTransportFluxSummary], limit: usize) -> Option<String> {
    let mut nonzero = fluxes
        .iter()
        .filter(|flux| {
            flux.uptake_amount > 0.0
                || flux.secretion_amount > 0.0
                || flux.uptake_events > 0
                || flux.secretion_events > 0
        })
        .collect::<Vec<_>>();
    nonzero.sort_by(|left, right| {
        let left_total = left.uptake_amount + left.secretion_amount;
        let right_total = right.uptake_amount + right.secretion_amount;
        right_total
            .partial_cmp(&left_total)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.species_index.cmp(&right.species_index))
    });
    if nonzero.is_empty() {
        return None;
    }
    Some(
        nonzero
            .into_iter()
            .take(limit)
            .map(|flux| {
                let name = usize::try_from(flux.species_index)
                    .ok()
                    .map(external_species_descriptor)
                    .map(|descriptor| descriptor.name)
                    .unwrap_or("unknown");
                format!(
                    "ext{}({}): uptake {:.4}/{} secretion {:.4}/{}",
                    flux.species_index,
                    name,
                    flux.uptake_amount,
                    flux.uptake_events,
                    flux.secretion_amount,
                    flux.secretion_events
                )
            })
            .collect::<Vec<_>>()
            .join(", "),
    )
}

fn format_byproduct_starter_rows(
    rows: &[StoichSpeciesStarterEventSummary],
    limit: usize,
) -> Option<String> {
    let mut nonzero = rows
        .iter()
        .filter(|row| row.amount > 0.0 || row.events > 0)
        .collect::<Vec<_>>();
    nonzero.sort_by(|left, right| {
        right
            .amount
            .partial_cmp(&left.amount)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.species_index.cmp(&right.species_index))
            .then_with(|| left.starter_type.cmp(&right.starter_type))
    });
    if nonzero.is_empty() {
        return None;
    }
    Some(
        nonzero
            .into_iter()
            .take(limit)
            .map(|row| {
                let name = usize::try_from(row.species_index)
                    .ok()
                    .map(external_species_descriptor)
                    .map(|descriptor| descriptor.name)
                    .unwrap_or("unknown");
                format!(
                    "ext{}({}) {}:{:.4}/{}",
                    row.species_index,
                    name,
                    starter_label(row.starter_type),
                    row.amount,
                    row.events
                )
            })
            .collect::<Vec<_>>()
            .join(", "),
    )
}

fn format_transport_starter_fluxes(
    fluxes: &[StoichTransportStarterFluxSummary],
    limit: usize,
) -> Option<String> {
    let mut nonzero = fluxes
        .iter()
        .filter(|flux| {
            flux.uptake_amount > 0.0
                || flux.secretion_amount > 0.0
                || flux.uptake_events > 0
                || flux.secretion_events > 0
        })
        .collect::<Vec<_>>();
    nonzero.sort_by(|left, right| {
        let left_total = left.uptake_amount + left.secretion_amount;
        let right_total = right.uptake_amount + right.secretion_amount;
        right_total
            .partial_cmp(&left_total)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.species_index.cmp(&right.species_index))
            .then_with(|| left.starter_type.cmp(&right.starter_type))
    });
    if nonzero.is_empty() {
        return None;
    }
    Some(
        nonzero
            .into_iter()
            .take(limit)
            .map(|flux| {
                let name = usize::try_from(flux.species_index)
                    .ok()
                    .map(external_species_descriptor)
                    .map(|descriptor| descriptor.name)
                    .unwrap_or("unknown");
                format!(
                    "ext{}({}) {}: uptake {:.4}/{} secretion {:.4}/{}",
                    flux.species_index,
                    name,
                    starter_label(flux.starter_type),
                    flux.uptake_amount,
                    flux.uptake_events,
                    flux.secretion_amount,
                    flux.secretion_events
                )
            })
            .collect::<Vec<_>>()
            .join(", "),
    )
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
        if let Some(latest) = analysis.chemistry.last()
            && !latest.species_profiles.is_empty()
        {
            out.push_str("\n### Latest Chemical Roles\n\n");
            out.push_str("| species | role | total | max | bond_energy | stability | permeability | storage | work |\n");
            out.push_str("|---|---|---:|---:|---:|---:|---:|---:|---:|\n");
            for profile in &latest.species_profiles {
                out.push_str(&format!(
                    "| ext{} {} | {} | {:.3} | {:.3} | {:.2} | {:.2} | {:.2} | {:.2} | {:.2} |\n",
                    profile.species,
                    profile.name,
                    chemical_role(profile),
                    profile.total_concentration,
                    profile.max_value,
                    profile.descriptor.bond_energy,
                    profile.descriptor.extracellular_stability,
                    profile.descriptor.membrane_permeability,
                    profile.descriptor.storage_density,
                    profile.descriptor.work_coupling
                ));
            }
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
    if let Some(stoich) = &analysis.stoich_v2 {
        out.push_str("\n## Stoichiometry V2\n\n");
        out.push_str(&format!(
            "- Summary present: {}\n- Events present: {}\n- Total ticks: {}\n- Enforcement: {}\n- Total events: {}\n- Imbalanced events: {}\n- Gross residual abs sum: {}\n",
            stoich.summary_present,
            stoich.events_present,
            stoich
                .total_ticks
                .map(|value| value.to_string())
                .unwrap_or_else(|| "n/a".to_string()),
            stoich.enforcement.as_deref().unwrap_or("n/a"),
            stoich.total_events,
            stoich.imbalanced_events,
            stoich
                .gross_residual_abs_sum
                .map(|value| format!("{value:.6}"))
                .unwrap_or_else(|| "n/a".to_string())
        ));
        if stoich.events_present {
            out.push_str(&format!(
                "- Reaction byproducts: {} events, amount {:.6}, model_abs {:.6}, residual_abs {:.6}\n- Reaction leakage: {} events, amount {:.6}, heat energy {:.6}\n",
                stoich.reaction_byproduct_events,
                stoich.reaction_byproduct_amount,
                stoich.reaction_byproduct_model_abs,
                stoich.reaction_byproduct_residual_abs,
                stoich.reaction_leakage_events,
                stoich.reaction_leakage_amount,
                stoich.reaction_leakage_energy_to_heat
            ));
            if !stoich.reaction_byproduct_by_species.is_empty() {
                out.push_str("\n### Reaction Byproducts By Species\n\n");
                out.push_str("| species | events | amount |\n");
                out.push_str("|---:|---:|---:|\n");
                for species in &stoich.reaction_byproduct_by_species {
                    out.push_str(&format!(
                        "| {} | {} | {:.6} |\n",
                        species.species_index, species.events, species.amount
                    ));
                }
            }
            let byproduct_starters = stoich
                .reaction_byproduct_by_species_starter
                .iter()
                .filter(|row| row.amount > 0.0 || row.events > 0)
                .collect::<Vec<_>>();
            if !byproduct_starters.is_empty() {
                out.push_str("\n### Reaction Byproduct Producers By Starter\n\n");
                out.push_str("| species | name | starter | events | amount |\n");
                out.push_str("|---:|---|---|---:|---:|\n");
                for row in byproduct_starters {
                    let name = usize::try_from(row.species_index)
                        .ok()
                        .map(external_species_descriptor)
                        .map(|descriptor| descriptor.name)
                        .unwrap_or("unknown");
                    out.push_str(&format!(
                        "| {} | {} | {} | {} | {:.6} |\n",
                        row.species_index,
                        name,
                        starter_label(row.starter_type),
                        row.events,
                        row.amount
                    ));
                }
            }
            let transport_fluxes = stoich
                .transport_flux_by_species
                .iter()
                .filter(|flux| {
                    flux.uptake_amount > 0.0
                        || flux.secretion_amount > 0.0
                        || flux.uptake_events > 0
                        || flux.secretion_events > 0
                })
                .collect::<Vec<_>>();
            if !transport_fluxes.is_empty() {
                out.push_str("\n### Transport Flux By Species\n\n");
                out.push_str("| species | name | uptake_events | uptake_amount | secretion_events | secretion_amount | net_uptake |\n");
                out.push_str("|---:|---|---:|---:|---:|---:|---:|\n");
                for flux in transport_fluxes {
                    let name = usize::try_from(flux.species_index)
                        .ok()
                        .map(external_species_descriptor)
                        .map(|descriptor| descriptor.name)
                        .unwrap_or("unknown");
                    out.push_str(&format!(
                        "| {} | {} | {} | {:.6} | {} | {:.6} | {:.6} |\n",
                        flux.species_index,
                        name,
                        flux.uptake_events,
                        flux.uptake_amount,
                        flux.secretion_events,
                        flux.secretion_amount,
                        flux.uptake_amount - flux.secretion_amount
                    ));
                }
            }
            let transport_starters = stoich
                .transport_flux_by_species_starter
                .iter()
                .filter(|flux| {
                    flux.uptake_amount > 0.0
                        || flux.secretion_amount > 0.0
                        || flux.uptake_events > 0
                        || flux.secretion_events > 0
                })
                .collect::<Vec<_>>();
            if !transport_starters.is_empty() {
                out.push_str("\n### Transport Flux By Species And Starter\n\n");
                out.push_str("| species | name | starter | uptake_events | uptake_amount | secretion_events | secretion_amount | net_uptake |\n");
                out.push_str("|---:|---|---|---:|---:|---:|---:|---:|\n");
                for flux in transport_starters {
                    let name = usize::try_from(flux.species_index)
                        .ok()
                        .map(external_species_descriptor)
                        .map(|descriptor| descriptor.name)
                        .unwrap_or("unknown");
                    out.push_str(&format!(
                        "| {} | {} | {} | {} | {:.6} | {} | {:.6} | {:.6} |\n",
                        flux.species_index,
                        name,
                        starter_label(flux.starter_type),
                        flux.uptake_events,
                        flux.uptake_amount,
                        flux.secretion_events,
                        flux.secretion_amount,
                        flux.uptake_amount - flux.secretion_amount
                    ));
                }
            }
        } else {
            out.push_str(
                "- Event-level reaction byproduct/leakage totals: unavailable (`stoich_v2_events.csv` not present)\n",
            );
        }
    }
    if let Some(calibration) = &analysis.byproduct_calibration {
        out.push_str("\n## Byproduct Calibration\n\n");
        out.push_str(&format!(
            "- Produced amount: {:.6}\n- Final byproduct pool: {}\n- Retained fraction: {}\n- Field tick: {}\n- Ruleset tick: {}\n- Missing field species: {}\n- Cross-feeding candidates: {}\n- Public-pool candidates: {}\n",
            calibration.total_byproduct_amount,
            calibration
                .final_byproduct_pool
                .map(|value| format!("{value:.6}"))
                .unwrap_or_else(|| "n/a".to_string()),
            calibration
                .retained_fraction
                .map(|value| format!("{value:.6}"))
                .unwrap_or_else(|| "n/a".to_string()),
            calibration
                .field_tick
                .map(|tick| tick.to_string())
                .unwrap_or_else(|| "n/a".to_string()),
            calibration
                .ruleset_tick
                .map(|tick| tick.to_string())
                .unwrap_or_else(|| "n/a".to_string()),
            calibration.missing_field_species,
            calibration.cross_feeding_candidates,
            calibration.public_pool_candidates
        ));
        out.push_str("\n| species | role | produced | final_pool | retained | uptake_slots | secretion_slots | interpretation |\n");
        out.push_str("|---:|---|---:|---:|---:|---:|---:|---|\n");
        for species in &calibration.species {
            out.push_str(&format!(
                "| {} {} | {} | {:.6} | {} | {} | {} | {} | {} |\n",
                species.species_index,
                species.name,
                species.role,
                species.produced_amount,
                species
                    .final_field_total
                    .map(|value| format!("{value:.6}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                species
                    .retained_fraction
                    .map(|value| format!("{value:.6}"))
                    .unwrap_or_else(|| "n/a".to_string()),
                species.uptake_active_slots,
                species.secretion_active_slots,
                species.interpretation
            ));
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
    out.push_str("| run | rng_seed | final_pop | max_pop | growth | final_energy | top% | mid% | deep% | unique_genotypes | dominant_genotype% | active_transporters/cell | gated_slots |\n");
    out.push_str("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for run in &analysis.runs {
        out.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            run.name,
            run.rng_seed
                .map(|seed| seed.to_string())
                .unwrap_or_else(|| "n/a".to_string()),
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
    if analysis
        .runs
        .iter()
        .any(|run| run.stoich_total_events.is_some())
    {
        out.push_str("\n## Stoichiometry V2 Comparison\n\n");
        out.push_str("| run | events | imbalanced | byproduct_events | byproduct_amount | gross_final_byproduct_pool | gross_retained_fraction | signed_excess_final_pool | signed_excess_retained_fraction | clipped_excess_final_pool | clipped_excess_retained_fraction | adjusted_uptake_pressure_candidates | adjusted_public_pool_candidates | leakage_events | leakage_heat |\n");
        out.push_str(
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
        );
        for run in &analysis.runs {
            out.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
                run.name,
                display_opt_u64(run.stoich_total_events),
                display_opt_u64(run.stoich_imbalanced_events),
                display_opt_u64(run.reaction_byproduct_events),
                display_opt_f64(run.reaction_byproduct_amount),
                display_opt_f64(run.byproduct_final_pool),
                display_opt_f64(run.byproduct_retained_fraction),
                display_opt_f64(run.byproduct_signed_excess_final_pool),
                display_opt_f64(run.byproduct_signed_excess_retained_fraction),
                display_opt_f64(run.byproduct_excess_final_pool),
                display_opt_f64(run.byproduct_excess_retained_fraction),
                display_opt_u64(run.byproduct_cross_feeding_candidates),
                display_opt_u64(run.byproduct_public_pool_candidates),
                display_opt_u64(run.reaction_leakage_events),
                display_opt_f64(run.reaction_leakage_energy_to_heat)
            ));
        }
    }
    if analysis
        .runs
        .iter()
        .any(|run| !run.byproduct_adjusted_species.is_empty())
    {
        out.push_str("\n## Baseline-Adjusted Byproduct Species\n\n");
        out.push_str("| run | species | role | produced | baseline_pool | final_pool | signed_excess | clipped_excess | clipped_retained | adjusted_interpretation |\n");
        out.push_str("|---|---:|---|---:|---:|---:|---:|---:|---:|---|\n");
        for run in &analysis.runs {
            for species in &run.byproduct_adjusted_species {
                out.push_str(&format!(
                    "| {} | {} {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
                    run.name,
                    species.species_index,
                    species.name,
                    species.role,
                    display_f64(species.produced_amount),
                    display_f64(species.baseline_field_total),
                    display_f64(species.final_field_total),
                    display_f64(species.signed_excess_final_pool),
                    display_f64(species.clipped_excess_final_pool),
                    display_opt_f64(species.clipped_excess_retained_fraction),
                    species.adjusted_interpretation
                ));
            }
        }
    }
    if !analysis.paired_byproduct_population.is_empty() {
        out.push_str("\n## Paired Byproduct Population Effects\n\n");
        out.push_str("| rng_seed | run | baseline | byproduct_amount | delta_final_pop | final_pop_ratio | delta_growth | delta_final_energy | signed_excess_retained_fraction | clipped_excess_retained_fraction |\n");
        out.push_str("|---:|---|---|---:|---:|---:|---:|---:|---:|---:|\n");
        for pair in &analysis.paired_byproduct_population {
            out.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
                pair.rng_seed
                    .map(|seed| seed.to_string())
                    .unwrap_or_else(|| "n/a".to_string()),
                pair.run_name,
                pair.baseline_name,
                display_opt_f64(pair.byproduct_amount),
                pair.delta_final_population
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "n/a".to_string()),
                display_opt_f64(pair.final_population_ratio),
                display_opt_f64(pair.delta_growth_factor),
                display_opt_f64(pair.delta_final_avg_energy),
                display_opt_f64(pair.byproduct_signed_excess_retained_fraction),
                display_opt_f64(pair.byproduct_excess_retained_fraction)
            ));
        }
    }
    if !analysis.byproduct_strength_aggregates.is_empty() {
        out.push_str("\n## Byproduct Strength Aggregates\n\n");
        out.push_str("| label | n | mean_byproduct_amount | mean_delta_final_pop | sd_delta_final_pop | delta_final_pop_range | positive_pop_delta | mean_pop_ratio | mean_delta_growth | mean_delta_final_energy | mean_signed_excess_retained | positive_signed_excess | mean_clipped_excess_retained | sd_clipped_excess_retained | clipped_excess_range | clipped_excess_warnings |\n");
        out.push_str(
            "|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|\n",
        );
        for aggregate in &analysis.byproduct_strength_aggregates {
            out.push_str(&format!(
                "| {} | {} | {} | {} | {} | {}..{} | {}/{} | {} | {} | {} | {} | {}/{} | {} | {} | {}..{} | {}/{} |\n",
                aggregate.label,
                aggregate.run_count,
                display_opt_f64(aggregate.mean_byproduct_amount),
                display_opt_f64(aggregate.mean_delta_final_population),
                display_opt_f64(aggregate.sd_delta_final_population),
                display_opt_f64(aggregate.min_delta_final_population),
                display_opt_f64(aggregate.max_delta_final_population),
                aggregate.positive_population_delta_count,
                aggregate.run_count,
                display_opt_f64(aggregate.mean_final_population_ratio),
                display_opt_f64(aggregate.mean_delta_growth_factor),
                display_opt_f64(aggregate.mean_delta_final_avg_energy),
                display_opt_f64(aggregate.mean_signed_excess_retained_fraction),
                aggregate.positive_signed_excess_count,
                aggregate.run_count,
                display_opt_f64(aggregate.mean_clipped_excess_retained_fraction),
                display_opt_f64(aggregate.sd_clipped_excess_retained_fraction),
                display_opt_f64(aggregate.min_clipped_excess_retained_fraction),
                display_opt_f64(aggregate.max_clipped_excess_retained_fraction),
                aggregate.clipped_excess_warning_count,
                aggregate.run_count
            ));
        }
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

fn read_stoich_v2_analysis(run_dir: &Path) -> AnalysisResult<Option<StoichV2Analysis>> {
    let summary_path = run_dir.join("stoich_v2_summary.json");
    let events_path = run_dir.join("stoich_v2_events.csv");
    let summary_present = summary_path.exists();
    let csv_events_present = events_path.exists();
    if !summary_present && !csv_events_present {
        return Ok(None);
    }

    let (total_ticks, enforcement, gross_residual_abs_sum, summary_accumulator) = if summary_present
    {
        let text = fs::read_to_string(&summary_path)?;
        let value: serde_json::Value = serde_json::from_str(&text)?;
        (
            value.get("total_ticks").and_then(|value| value.as_u64()),
            value
                .get("enforcement")
                .and_then(|value| value.as_str())
                .map(str::to_owned),
            value
                .get("gross_residual_abs_sum")
                .and_then(|value| value.as_f64()),
            parse_stoich_v2_summary_accumulator(&value),
        )
    } else {
        (None, None, None, None)
    };
    let summary_events_present = summary_accumulator.is_some();
    let events_present = csv_events_present || summary_events_present;

    let accumulator = if csv_events_present {
        let mut accumulator = parse_stoich_v2_events(&fs::read_to_string(&events_path)?)?;
        if let Some(summary_accumulator) = summary_accumulator {
            if stoich_amounts_match(
                accumulator.byproduct.amount,
                summary_accumulator.byproduct.amount,
                accumulator.byproduct.events,
            ) {
                accumulator.byproduct_by_species_starter =
                    summary_accumulator.byproduct_by_species_starter;
            }
            accumulator.transport_flux_by_species = summary_accumulator.transport_flux_by_species;
            accumulator.transport_flux_by_species_starter =
                summary_accumulator.transport_flux_by_species_starter;
        }
        accumulator
    } else {
        summary_accumulator.unwrap_or_default()
    };

    let mut byproduct_by_species = accumulator
        .byproduct_by_species
        .into_values()
        .collect::<Vec<_>>();
    byproduct_by_species.sort_by_key(|summary| summary.species_index);
    let mut byproduct_by_species_starter = accumulator
        .byproduct_by_species_starter
        .into_values()
        .collect::<Vec<_>>();
    byproduct_by_species_starter
        .sort_by_key(|summary| (summary.species_index, summary.starter_type));
    let mut transport_flux_by_species = accumulator
        .transport_flux_by_species
        .into_values()
        .collect::<Vec<_>>();
    transport_flux_by_species.sort_by_key(|summary| summary.species_index);
    let mut transport_flux_by_species_starter = accumulator
        .transport_flux_by_species_starter
        .into_values()
        .collect::<Vec<_>>();
    transport_flux_by_species_starter
        .sort_by_key(|summary| (summary.species_index, summary.starter_type));

    Ok(Some(StoichV2Analysis {
        summary_present,
        events_present,
        total_ticks,
        enforcement,
        gross_residual_abs_sum,
        total_events: accumulator.total_events,
        imbalanced_events: accumulator.imbalanced_events,
        reaction_byproduct_events: accumulator.byproduct.events,
        reaction_byproduct_amount: accumulator.byproduct.amount,
        reaction_byproduct_model_abs: accumulator.byproduct.model_abs,
        reaction_byproduct_residual_abs: accumulator.byproduct.residual_abs,
        reaction_byproduct_by_species: byproduct_by_species,
        reaction_byproduct_by_species_starter: byproduct_by_species_starter,
        reaction_leakage_events: accumulator.leakage.events,
        reaction_leakage_amount: accumulator.leakage.amount,
        reaction_leakage_energy_to_heat: accumulator.leakage.reservoir_energy,
        transport_flux_by_species,
        transport_flux_by_species_starter,
    }))
}

fn stoich_amounts_match(left: f64, right: f64, rounded_event_count: u64) -> bool {
    let scale = left.abs().max(right.abs()).max(1.0);
    let relative_slack = STOICH_COMPACT_MERGE_REL_TOLERANCE * scale;
    let csv_rounding_slack = 0.5 * STOICH_EVENT_CSV_AMOUNT_DECIMALS * rounded_event_count as f64;
    (left - right).abs() <= relative_slack + csv_rounding_slack
}

fn parse_summary_species_index(entry: &serde_json::Value) -> Option<i16> {
    let raw = entry.get("species_index")?.as_i64()?;
    let species_index = i16::try_from(raw).ok()?;
    let species = usize::try_from(species_index).ok()?;
    (species < S_EXT).then_some(species_index)
}

fn parse_summary_starter_type(entry: &serde_json::Value) -> Option<u8> {
    let raw = entry.get("starter_type")?.as_u64()?;
    let starter_type = u8::try_from(raw).ok()?;
    (usize::from(starter_type) < STOICH_STARTER_TYPE_COUNT).then_some(starter_type)
}

fn parse_stoich_v2_summary_accumulator(
    value: &serde_json::Value,
) -> Option<StoichV2EventAccumulator> {
    let ledger = value.get("ledger")?;
    let byproduct = ledger.get("reaction_byproduct")?;
    let leakage = ledger.get("reaction_leakage")?;
    let mut accumulator = StoichV2EventAccumulator::default();
    accumulator.total_events = value
        .get("stages")
        .and_then(|stages| stages.as_array())
        .map(|stages| {
            stages
                .iter()
                .filter_map(|stage| stage.get("summary"))
                .filter_map(|summary| summary.get("event_count"))
                .filter_map(|count| count.as_u64())
                .sum()
        })
        .unwrap_or(0);
    accumulator.imbalanced_events = value
        .get("stages")
        .and_then(|stages| stages.as_array())
        .map(|stages| {
            stages
                .iter()
                .filter_map(|stage| stage.get("summary"))
                .filter_map(|summary| summary.get("imbalanced_event_count"))
                .filter_map(|count| count.as_u64())
                .sum()
        })
        .unwrap_or(0);
    accumulator.byproduct.events = byproduct
        .get("events")
        .and_then(|value| value.as_u64())
        .unwrap_or(0);
    accumulator.byproduct.amount = byproduct
        .get("amount")
        .and_then(|value| value.as_f64())
        .unwrap_or(0.0);
    accumulator.byproduct.model_abs = byproduct
        .get("model_abs")
        .and_then(|value| value.as_f64())
        .unwrap_or(0.0);
    accumulator.byproduct.residual_abs = byproduct
        .get("residual_abs")
        .and_then(|value| value.as_f64())
        .unwrap_or(0.0);
    accumulator.leakage.events = leakage
        .get("events")
        .and_then(|value| value.as_u64())
        .unwrap_or(0);
    accumulator.leakage.amount = leakage
        .get("amount")
        .and_then(|value| value.as_f64())
        .unwrap_or(0.0);
    accumulator.leakage.model_abs = leakage
        .get("model_abs")
        .and_then(|value| value.as_f64())
        .unwrap_or(0.0);
    accumulator.leakage.residual_abs = leakage
        .get("residual_abs")
        .and_then(|value| value.as_f64())
        .unwrap_or(0.0);
    accumulator.leakage.reservoir_energy = leakage
        .get("reservoir_energy")
        .and_then(|value| value.as_f64())
        .unwrap_or(0.0);

    if let Some(species) = ledger
        .get("reaction_byproduct_by_species")
        .and_then(|value| value.as_array())
    {
        for entry in species {
            let Some(events) = entry.get("events").and_then(|value| value.as_u64()) else {
                continue;
            };
            if events == 0 {
                continue;
            }
            let Some(species_index) = parse_summary_species_index(entry) else {
                continue;
            };
            accumulator.byproduct_by_species.insert(
                species_index,
                StoichSpeciesEventSummary {
                    species_index,
                    amount: entry
                        .get("amount")
                        .and_then(|value| value.as_f64())
                        .unwrap_or(0.0),
                    events,
                },
            );
        }
    }

    if let Some(species) = ledger
        .get("reaction_byproduct_by_species_starter")
        .and_then(|value| value.as_array())
    {
        for entry in species {
            let Some(events) = entry.get("events").and_then(|value| value.as_u64()) else {
                continue;
            };
            let amount = entry
                .get("amount")
                .and_then(|value| value.as_f64())
                .unwrap_or(0.0);
            if events == 0 && amount == 0.0 {
                continue;
            }
            let Some(species_index) = parse_summary_species_index(entry) else {
                continue;
            };
            let Some(starter_type) = parse_summary_starter_type(entry) else {
                continue;
            };
            accumulator.byproduct_by_species_starter.insert(
                (species_index, starter_type),
                StoichSpeciesStarterEventSummary {
                    species_index,
                    starter_type,
                    amount,
                    events,
                },
            );
        }
    }

    if let Some(species) = ledger
        .get("transport_flux_by_species")
        .and_then(|value| value.as_array())
    {
        for entry in species {
            let uptake_events = entry
                .get("uptake_events")
                .and_then(|value| value.as_u64())
                .unwrap_or(0);
            let secretion_events = entry
                .get("secretion_events")
                .and_then(|value| value.as_u64())
                .unwrap_or(0);
            let uptake_amount = entry
                .get("uptake_amount")
                .and_then(|value| value.as_f64())
                .unwrap_or(0.0);
            let secretion_amount = entry
                .get("secretion_amount")
                .and_then(|value| value.as_f64())
                .unwrap_or(0.0);
            if uptake_events == 0
                && secretion_events == 0
                && uptake_amount == 0.0
                && secretion_amount == 0.0
            {
                continue;
            }
            let Some(species_index) = parse_summary_species_index(entry) else {
                continue;
            };
            accumulator.transport_flux_by_species.insert(
                species_index,
                StoichTransportFluxSummary {
                    species_index,
                    uptake_amount,
                    uptake_events,
                    secretion_amount,
                    secretion_events,
                },
            );
        }
    }

    if let Some(species) = ledger
        .get("transport_flux_by_species_starter")
        .and_then(|value| value.as_array())
    {
        for entry in species {
            let uptake_events = entry
                .get("uptake_events")
                .and_then(|value| value.as_u64())
                .unwrap_or(0);
            let secretion_events = entry
                .get("secretion_events")
                .and_then(|value| value.as_u64())
                .unwrap_or(0);
            let uptake_amount = entry
                .get("uptake_amount")
                .and_then(|value| value.as_f64())
                .unwrap_or(0.0);
            let secretion_amount = entry
                .get("secretion_amount")
                .and_then(|value| value.as_f64())
                .unwrap_or(0.0);
            if uptake_events == 0
                && secretion_events == 0
                && uptake_amount == 0.0
                && secretion_amount == 0.0
            {
                continue;
            }
            let Some(species_index) = parse_summary_species_index(entry) else {
                continue;
            };
            let Some(starter_type) = parse_summary_starter_type(entry) else {
                continue;
            };
            accumulator.transport_flux_by_species_starter.insert(
                (species_index, starter_type),
                StoichTransportStarterFluxSummary {
                    species_index,
                    starter_type,
                    uptake_amount,
                    uptake_events,
                    secretion_amount,
                    secretion_events,
                },
            );
        }
    }

    Some(accumulator)
}

fn parse_stoich_v2_events(text: &str) -> AnalysisResult<StoichV2EventAccumulator> {
    let mut lines = text.lines();
    let header = lines.next().ok_or("stoich_v2_events.csv is empty")?;
    let columns: Vec<&str> = header.split(',').collect();
    let index = |name: &str| -> AnalysisResult<usize> {
        columns
            .iter()
            .position(|column| *column == name)
            .ok_or_else(|| format!("stoich_v2_events.csv missing column {name}").into())
    };
    let kind_i = index("kind")?;
    let species_i = index("species_index")?;
    let amount_i = index("amount")?;
    let model_c_i = index("model_c")?;
    let model_h_i = index("model_h")?;
    let model_o_i = index("model_o")?;
    let model_s_i = index("model_s")?;
    let model_redox_i = index("model_redox")?;
    let model_energy_i = index("model_energy")?;
    let reservoir_energy_i = index("reservoir_energy")?;
    let residual_abs_i = index("residual_abs")?;
    let balanced_i = index("balanced")?;

    let mut accumulator = StoichV2EventAccumulator::default();
    for (line_no, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let values: Vec<&str> = line.split(',').collect();
        let get = |i: usize| -> AnalysisResult<&str> {
            values.get(i).copied().ok_or_else(|| {
                format!(
                    "stoich_v2_events.csv line {} missing column {}",
                    line_no + 2,
                    i
                )
                .into()
            })
        };
        let kind = get(kind_i)?;
        let species_index = get(species_i)?.parse::<i16>()?;
        let amount = get(amount_i)?.parse::<f64>()?;
        let model_abs = [
            model_c_i,
            model_h_i,
            model_o_i,
            model_s_i,
            model_redox_i,
            model_energy_i,
        ]
        .into_iter()
        .map(|i| get(i).and_then(|value| value.parse::<f64>().map_err(|e| e.into())))
        .collect::<AnalysisResult<Vec<_>>>()?
        .into_iter()
        .map(f64::abs)
        .sum::<f64>();
        let reservoir_energy = get(reservoir_energy_i)?.parse::<f64>()?;
        let residual_abs = get(residual_abs_i)?.parse::<f64>()?;
        let balanced = get(balanced_i)?.parse::<u8>()? != 0;

        accumulator.total_events += 1;
        if !balanced {
            accumulator.imbalanced_events += 1;
        }
        match kind {
            "reaction_byproduct" => {
                accumulator.byproduct.events += 1;
                accumulator.byproduct.amount += amount;
                accumulator.byproduct.model_abs += model_abs;
                accumulator.byproduct.residual_abs += residual_abs;
                let species_summary = accumulator
                    .byproduct_by_species
                    .entry(species_index)
                    .or_insert(StoichSpeciesEventSummary {
                        species_index,
                        amount: 0.0,
                        events: 0,
                    });
                species_summary.amount += amount;
                species_summary.events += 1;
            }
            "reaction_leakage" => {
                accumulator.leakage.events += 1;
                accumulator.leakage.amount += amount;
                accumulator.leakage.model_abs += model_abs;
                accumulator.leakage.residual_abs += residual_abs;
                accumulator.leakage.reservoir_energy += reservoir_energy;
            }
            _ => {}
        }
    }
    Ok(accumulator)
}

fn summarize_byproduct_calibration(
    stoich: Option<&StoichV2Analysis>,
    latest_chemistry: Option<&SnapshotChemistry>,
    latest_rulesets: Option<&RulesetSummary>,
) -> Option<ByproductCalibrationSummary> {
    let stoich = stoich.filter(|stoich| stoich.events_present)?;
    if stoich.reaction_byproduct_amount <= 0.0 || stoich.reaction_byproduct_by_species.is_empty() {
        return None;
    }

    let mut species = Vec::new();
    let mut final_pool_sum = 0.0;
    let mut final_pool_seen = false;
    let mut missing_field_species = 0;
    let mut cross_feeding_candidates = 0;
    let mut public_pool_candidates = 0;

    for byproduct in &stoich.reaction_byproduct_by_species {
        let species_index = byproduct.species_index;
        let species_usize = usize::try_from(species_index).ok();
        let profile = species_usize.and_then(|species| {
            latest_chemistry.and_then(|chemistry| {
                chemistry
                    .species_profiles
                    .iter()
                    .find(|profile| profile.species == species)
            })
        });
        let transport = species_usize.and_then(|species| {
            latest_rulesets.and_then(|rulesets| {
                rulesets
                    .transporters
                    .external_species
                    .iter()
                    .find(|summary| summary.ext_species as usize == species)
            })
        });
        let descriptor = species_usize
            .map(external_species_descriptor)
            .unwrap_or_else(|| external_species_descriptor(usize::MAX));
        let descriptor_profile: SpeciesDescriptorProfile = descriptor.into();
        let final_field_total = profile.map(|profile| profile.total_concentration);
        if let Some(total) = final_field_total {
            final_pool_seen = true;
            final_pool_sum += total;
        } else {
            missing_field_species += 1;
        }
        let retained_fraction = final_field_total.and_then(|total| {
            (byproduct.amount > 0.0).then_some((total / byproduct.amount).max(0.0))
        });
        let uptake_active_slots = transport
            .map(|summary| summary.uptake_active_slots)
            .unwrap_or(0);
        let secretion_active_slots = transport
            .map(|summary| summary.secretion_active_slots)
            .unwrap_or(0);
        let interpretation = classify_byproduct_calibration(
            retained_fraction,
            transport
                .map(|summary| summary.uptake_rate_sum)
                .unwrap_or(0.0),
            transport
                .map(|summary| summary.secretion_rate_sum)
                .unwrap_or(0.0),
        );
        if is_meaningful_byproduct_candidate(byproduct.amount, stoich.reaction_byproduct_amount) {
            if interpretation == "cross_feeding_candidate" {
                cross_feeding_candidates += 1;
            } else if interpretation == "public_pool_candidate" {
                public_pool_candidates += 1;
            }
        }
        species.push(ByproductSpeciesCalibration {
            species_index,
            name: descriptor.name.to_string(),
            role: chemical_role_from_descriptor(&descriptor_profile).to_string(),
            produced_amount: byproduct.amount,
            final_field_total,
            retained_fraction,
            uptake_active_slots,
            secretion_active_slots,
            uptake_rate_sum: transport
                .map(|summary| summary.uptake_rate_sum)
                .unwrap_or(0.0),
            secretion_rate_sum: transport
                .map(|summary| summary.secretion_rate_sum)
                .unwrap_or(0.0),
            interpretation: interpretation.to_string(),
        });
    }

    let complete_field_coverage = final_pool_seen && missing_field_species == 0;
    let final_byproduct_pool = complete_field_coverage.then_some(final_pool_sum);
    let retained_fraction = final_byproduct_pool.and_then(|pool| {
        (stoich.reaction_byproduct_amount > 0.0).then_some(pool / stoich.reaction_byproduct_amount)
    });

    Some(ByproductCalibrationSummary {
        total_byproduct_amount: stoich.reaction_byproduct_amount,
        final_byproduct_pool,
        retained_fraction,
        field_tick: latest_chemistry.map(|chemistry| chemistry.tick),
        ruleset_tick: latest_rulesets.map(|rulesets| rulesets.tick),
        missing_field_species,
        cross_feeding_candidates,
        public_pool_candidates,
        species,
    })
}

fn classify_byproduct_calibration(
    retained_fraction: Option<f64>,
    uptake_rate_sum: f64,
    secretion_rate_sum: f64,
) -> &'static str {
    let has_uptake_pressure = uptake_rate_sum > ACTIVE_TRANSPORT_THRESHOLD as f64;
    if let Some(retained) = retained_fraction {
        if uptake_rate_sum > secretion_rate_sum && retained < 0.5 {
            "cross_feeding_candidate"
        } else if retained >= 0.75 && !has_uptake_pressure {
            "public_pool_candidate"
        } else if retained >= 0.75 {
            "accumulating_pool"
        } else if has_uptake_pressure {
            "uptake_pressure"
        } else {
            "low_retention"
        }
    } else if has_uptake_pressure {
        "uptake_pressure_no_field_snapshot"
    } else {
        "production_only"
    }
}

fn is_meaningful_byproduct_candidate(produced_amount: f64, total_byproduct_amount: f64) -> bool {
    produced_amount >= BYPRODUCT_CANDIDATE_MIN_AMOUNT
        && produced_amount
            >= total_byproduct_amount.max(0.0) * BYPRODUCT_CANDIDATE_MIN_TOTAL_FRACTION
}

fn apply_byproduct_baseline_adjustment(
    analyses: &[RunAnalysis],
    runs: &mut [RunComparisonEntry],
) -> Vec<String> {
    let mut warnings = Vec::new();
    let baseline_candidates = analyses
        .iter()
        .filter(|analysis| {
            analysis
                .stoich_v2
                .as_ref()
                .filter(|stoich| stoich.events_present)
                .is_some_and(|stoich| stoich.reaction_byproduct_amount == 0.0)
                && analysis.chemistry.last().is_some()
        })
        .collect::<Vec<_>>();
    if baseline_candidates.is_empty() {
        return warnings;
    }

    let mut baselines_by_seed: HashMap<u64, &RunAnalysis> = HashMap::new();
    let mut ambiguous_baseline_seeds = HashSet::new();
    for baseline_analysis in &baseline_candidates {
        let Some(seed) = baseline_analysis.rng_seed else {
            if baseline_candidates.len() > 1 {
                warnings.push(format!(
                    "byproduct excess baseline: zero-byproduct baseline {} has no rng_seed and cannot be used for seed-paired adjustment",
                    run_name(&baseline_analysis.run_dir)
                ));
            }
            continue;
        };
        if ambiguous_baseline_seeds.contains(&seed) {
            warnings.push(format!(
                "byproduct excess baseline: additional duplicate zero-byproduct baseline {} for rng_seed {seed} is ignored because that seed is already ambiguous",
                run_name(&baseline_analysis.run_dir)
            ));
        } else if let Some(previous) = baselines_by_seed.remove(&seed) {
            ambiguous_baseline_seeds.insert(seed);
            warnings.push(format!(
                "byproduct excess baseline: duplicate zero-byproduct baselines for rng_seed {seed}; skipping seed-paired adjustment for this seed because {} and {} are ambiguous",
                run_name(&previous.run_dir),
                run_name(&baseline_analysis.run_dir)
            ));
        } else {
            baselines_by_seed.insert(seed, *baseline_analysis);
        }
    }

    for (analysis, run) in analyses.iter().zip(runs.iter_mut()) {
        let Some(calibration) = &analysis.byproduct_calibration else {
            continue;
        };
        let Some(produced) = run.reaction_byproduct_amount.filter(|amount| *amount > 0.0) else {
            continue;
        };
        let baseline_analysis = if let Some(seed) = analysis.rng_seed {
            if ambiguous_baseline_seeds.contains(&seed) {
                warnings.push(format!(
                    "{}: byproduct excess baseline skipped because rng_seed {seed} has multiple zero-byproduct baselines",
                    run.name
                ));
                continue;
            }
            match baselines_by_seed.get(&seed).copied() {
                Some(baseline) => baseline,
                None if baseline_candidates.len() == 1 => {
                    let baseline = baseline_candidates[0];
                    if let Some(baseline_seed) = baseline.rng_seed {
                        warnings.push(format!(
                            "{}: byproduct excess baseline skipped because run rng_seed {seed} does not match the only zero-byproduct baseline rng_seed {baseline_seed}",
                            run.name
                        ));
                        continue;
                    }
                    baseline
                }
                None => {
                    warnings.push(format!(
                        "{}: byproduct excess baseline skipped because no zero-byproduct baseline with rng_seed {seed} was found",
                        run.name
                    ));
                    continue;
                }
            }
        } else if baseline_candidates.len() == 1 {
            baseline_candidates[0]
        } else {
            warnings.push(format!(
                "{}: byproduct excess baseline skipped because run has no rng_seed and multiple zero-byproduct baselines are available",
                run.name
            ));
            continue;
        };
        let Some(baseline_chemistry) = baseline_analysis.chemistry.last() else {
            continue;
        };
        if analysis.grid != baseline_analysis.grid {
            warnings.push(format!(
                "{}: byproduct excess baseline skipped because grid {:?} differs from baseline {:?}",
                run.name, analysis.grid, baseline_analysis.grid
            ));
            continue;
        }
        if calibration.field_tick != Some(baseline_chemistry.tick) {
            warnings.push(format!(
                "{}: byproduct excess baseline skipped because field tick {:?} differs from baseline tick {}",
                run.name, calibration.field_tick, baseline_chemistry.tick
            ));
            continue;
        }
        let mut signed_excess_pool = 0.0;
        let mut clipped_excess_pool = 0.0;
        let mut adjusted_uptake_pressure_candidates = 0;
        let mut adjusted_public_pool_candidates = 0;
        let mut adjusted_species = Vec::new();
        let mut complete = true;
        for species in &calibration.species {
            let Some(run_total) = species.final_field_total else {
                complete = false;
                break;
            };
            let Some(baseline_total) = baseline_chemistry
                .species_profiles
                .iter()
                .find(|profile| profile.species as i16 == species.species_index)
                .map(|profile| profile.total_concentration)
            else {
                complete = false;
                break;
            };
            let signed_excess = run_total - baseline_total;
            let clipped_excess = signed_excess.max(0.0);
            signed_excess_pool += signed_excess;
            clipped_excess_pool += clipped_excess;
            let adjusted_retained_fraction =
                (species.produced_amount > 0.0).then_some(clipped_excess / species.produced_amount);
            let adjusted_interpretation = classify_byproduct_calibration(
                adjusted_retained_fraction,
                species.uptake_rate_sum,
                species.secretion_rate_sum,
            );
            let meaningful_candidate =
                is_meaningful_byproduct_candidate(species.produced_amount, produced);
            if meaningful_candidate {
                if adjusted_interpretation == "cross_feeding_candidate" {
                    adjusted_uptake_pressure_candidates += 1;
                } else if adjusted_interpretation == "public_pool_candidate" {
                    adjusted_public_pool_candidates += 1;
                }
            }
            let adjusted_interpretation = if meaningful_candidate {
                adjusted_interpretation
            } else if clipped_excess > 0.0 {
                "trace_byproduct_field_shift"
            } else {
                "trace_byproduct"
            };
            adjusted_species.push(ByproductAdjustedSpeciesComparison {
                species_index: species.species_index,
                name: species.name.clone(),
                role: species.role.clone(),
                produced_amount: species.produced_amount,
                baseline_field_total: baseline_total,
                final_field_total: run_total,
                signed_excess_final_pool: signed_excess,
                clipped_excess_final_pool: clipped_excess,
                clipped_excess_retained_fraction: meaningful_candidate
                    .then_some(adjusted_retained_fraction)
                    .flatten(),
                adjusted_interpretation: adjusted_interpretation.to_string(),
            });
        }
        if complete {
            run.byproduct_signed_excess_final_pool = Some(signed_excess_pool);
            run.byproduct_signed_excess_retained_fraction = Some(signed_excess_pool / produced);
            run.byproduct_excess_final_pool = Some(clipped_excess_pool);
            run.byproduct_excess_retained_fraction = Some(clipped_excess_pool / produced);
            run.byproduct_cross_feeding_candidates = Some(adjusted_uptake_pressure_candidates);
            run.byproduct_public_pool_candidates = Some(adjusted_public_pool_candidates);
            run.byproduct_adjusted_species = adjusted_species;
        }
    }
    warnings
}

fn summarize_paired_byproduct_population(
    runs: &[RunComparisonEntry],
) -> Vec<ByproductPairedComparison> {
    let baseline_candidates = runs
        .iter()
        .filter(|run| {
            run.reaction_byproduct_amount
                .is_some_and(|amount| amount == 0.0)
        })
        .collect::<Vec<_>>();
    if baseline_candidates.is_empty() {
        return Vec::new();
    }

    let mut seed_baselines: HashMap<u64, &RunComparisonEntry> = HashMap::new();
    let mut duplicate_seed_baselines = HashSet::new();
    for baseline in &baseline_candidates {
        if let Some(seed) = baseline.rng_seed
            && seed_baselines.insert(seed, baseline).is_some()
        {
            duplicate_seed_baselines.insert(seed);
        }
    }

    let mut pairs = Vec::new();
    for run in runs {
        let Some(byproduct_amount) = run.reaction_byproduct_amount.filter(|amount| *amount > 0.0)
        else {
            continue;
        };
        let baseline = if let Some(seed) = run.rng_seed {
            if duplicate_seed_baselines.contains(&seed) {
                continue;
            }
            seed_baselines.get(&seed).copied()
        } else if baseline_candidates.len() == 1 {
            baseline_candidates.first().copied()
        } else {
            None
        };
        let Some(baseline) = baseline else {
            continue;
        };

        let delta_final_population = match (run.final_population, baseline.final_population) {
            (Some(run_pop), Some(baseline_pop)) => Some(run_pop as i64 - baseline_pop as i64),
            _ => None,
        };
        let final_population_ratio = match (run.final_population, baseline.final_population) {
            (Some(run_pop), Some(baseline_pop)) if baseline_pop > 0 => {
                Some(run_pop as f64 / baseline_pop as f64)
            }
            _ => None,
        };
        let delta_growth_factor = match (run.growth_factor, baseline.growth_factor) {
            (Some(run_growth), Some(baseline_growth)) => Some(run_growth - baseline_growth),
            _ => None,
        };
        let delta_final_avg_energy = match (run.final_avg_energy, baseline.final_avg_energy) {
            (Some(run_energy), Some(baseline_energy)) => Some(run_energy - baseline_energy),
            _ => None,
        };

        pairs.push(ByproductPairedComparison {
            rng_seed: run.rng_seed,
            baseline_name: baseline.name.clone(),
            run_name: run.name.clone(),
            byproduct_amount: Some(byproduct_amount),
            delta_final_population,
            final_population_ratio,
            delta_growth_factor,
            delta_final_avg_energy,
            byproduct_signed_excess_retained_fraction: run
                .byproduct_signed_excess_retained_fraction,
            byproduct_excess_retained_fraction: run.byproduct_excess_retained_fraction,
        });
    }

    pairs
}

fn summarize_byproduct_strength_aggregates(
    pairs: &[ByproductPairedComparison],
) -> Vec<ByproductStrengthAggregate> {
    let mut groups: BTreeMap<String, Vec<&ByproductPairedComparison>> = BTreeMap::new();
    for pair in pairs {
        let label = byproduct_strength_label(&pair.run_name);
        groups.entry(label).or_default().push(pair);
    }

    groups
        .into_iter()
        .map(|(label, pairs)| {
            let delta_population_stats = summarize_numbers(
                pairs
                    .iter()
                    .filter_map(|pair| pair.delta_final_population.map(|value| value as f64)),
            );
            let clipped_retention_stats = summarize_numbers(
                pairs
                    .iter()
                    .filter_map(|pair| pair.byproduct_excess_retained_fraction),
            );
            let positive_population_delta_count = pairs
                .iter()
                .filter(|pair| pair.delta_final_population.is_some_and(|delta| delta > 0))
                .count();
            let positive_signed_excess_count = pairs
                .iter()
                .filter(|pair| {
                    pair.byproduct_signed_excess_retained_fraction
                        .is_some_and(|fraction| fraction > 0.0)
                })
                .count();
            let clipped_excess_warning_count = pairs
                .iter()
                .filter(|pair| {
                    pair.byproduct_excess_retained_fraction
                        .is_some_and(|fraction| fraction >= 0.75)
                })
                .count();
            ByproductStrengthAggregate {
                label,
                run_count: pairs.len(),
                mean_byproduct_amount: mean_f64(
                    pairs.iter().filter_map(|pair| pair.byproduct_amount),
                ),
                mean_delta_final_population: delta_population_stats.mean,
                sd_delta_final_population: delta_population_stats.stddev,
                min_delta_final_population: delta_population_stats.min,
                max_delta_final_population: delta_population_stats.max,
                mean_final_population_ratio: mean_f64(
                    pairs.iter().filter_map(|pair| pair.final_population_ratio),
                ),
                mean_delta_growth_factor: mean_f64(
                    pairs.iter().filter_map(|pair| pair.delta_growth_factor),
                ),
                mean_delta_final_avg_energy: mean_f64(
                    pairs.iter().filter_map(|pair| pair.delta_final_avg_energy),
                ),
                mean_signed_excess_retained_fraction: mean_f64(
                    pairs
                        .iter()
                        .filter_map(|pair| pair.byproduct_signed_excess_retained_fraction),
                ),
                mean_clipped_excess_retained_fraction: clipped_retention_stats.mean,
                sd_clipped_excess_retained_fraction: clipped_retention_stats.stddev,
                min_clipped_excess_retained_fraction: clipped_retention_stats.min,
                max_clipped_excess_retained_fraction: clipped_retention_stats.max,
                positive_population_delta_count,
                positive_signed_excess_count,
                clipped_excess_warning_count,
            }
        })
        .collect()
}

#[derive(Debug, Clone, Copy)]
struct NumberSummary {
    mean: Option<f64>,
    stddev: Option<f64>,
    min: Option<f64>,
    max: Option<f64>,
}

fn summarize_numbers(values: impl Iterator<Item = f64>) -> NumberSummary {
    let values = values.filter(|value| value.is_finite()).collect::<Vec<_>>();
    if values.is_empty() {
        return NumberSummary {
            mean: None,
            stddev: None,
            min: None,
            max: None,
        };
    }
    let mean = values.iter().copied().sum::<f64>() / values.len() as f64;
    let variance = if values.len() > 1 {
        values
            .iter()
            .map(|value| {
                let delta = value - mean;
                delta * delta
            })
            .sum::<f64>()
            / (values.len() - 1) as f64
    } else {
        0.0
    };
    NumberSummary {
        mean: Some(mean),
        stddev: Some(variance.sqrt()),
        min: values.iter().copied().reduce(f64::min),
        max: values.iter().copied().reduce(f64::max),
    }
}

fn byproduct_strength_label(run_name: &str) -> String {
    run_name
        .rsplit_once('_')
        .map(|(_, suffix)| suffix.to_string())
        .filter(|suffix| suffix.starts_with("byp"))
        .unwrap_or_else(|| run_name.to_string())
}

fn mean_f64(values: impl Iterator<Item = f64>) -> Option<f64> {
    let mut count = 0usize;
    let mut sum = 0.0;
    for value in values {
        if value.is_finite() {
            count += 1;
            sum += value;
        }
    }
    (count > 0).then_some(sum / count as f64)
}

fn read_run_rng_seed(run_dir: &Path) -> Option<u64> {
    let summary = fs::read_to_string(run_dir.join("summary.md")).ok()?;
    summary
        .lines()
        .find_map(|line| line.strip_prefix("- RNG seed: "))
        .and_then(|value| value.trim().parse().ok())
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

    let tracked_species = (0..s_ext).collect::<Vec<_>>();
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
            let descriptor = external_species_descriptor(*species);
            SpeciesProfile {
                species: *species,
                name: descriptor.name.to_string(),
                descriptor: descriptor.into(),
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
    let mut species_acc: HashMap<u8, TransportSpeciesAccumulator> = HashMap::new();
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
            let species = species_acc.entry(tp.ext_species).or_default();
            species.active_slots += count;
            species.uptake_rate_sum += f64::from(uptake) * count as f64;
            species.secretion_rate_sum += f64::from(secretion) * count as f64;
            if uptake > ACTIVE_TRANSPORT_THRESHOLD {
                species.uptake_active_slots += count;
            }
            if secretion > ACTIVE_TRANSPORT_THRESHOLD {
                species.secretion_active_slots += count;
            }
            if uptake > secretion + ACTIVE_TRANSPORT_THRESHOLD {
                species.uptake_dominant_slots += count;
            } else if secretion > uptake + ACTIVE_TRANSPORT_THRESHOLD {
                species.secretion_dominant_slots += count;
            }
            if tp.gate_weight.abs() > ACTIVE_TRANSPORT_THRESHOLD {
                species.gated_slots += count;
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

    let mut external_species = species_acc
        .into_iter()
        .map(|(ext_species, acc)| TransportSpeciesSummary {
            ext_species,
            active_slots: acc.active_slots,
            uptake_active_slots: acc.uptake_active_slots,
            secretion_active_slots: acc.secretion_active_slots,
            uptake_dominant_slots: acc.uptake_dominant_slots,
            secretion_dominant_slots: acc.secretion_dominant_slots,
            gated_slots: acc.gated_slots,
            uptake_rate_sum: acc.uptake_rate_sum,
            secretion_rate_sum: acc.secretion_rate_sum,
        })
        .collect::<Vec<_>>();
    external_species.sort_by(|a, b| {
        b.active_slots
            .cmp(&a.active_slots)
            .then_with(|| a.ext_species.cmp(&b.ext_species))
    });

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
        external_species,
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
    if let Some(calibration) = &analysis.byproduct_calibration {
        if calibration.cross_feeding_candidates > 0 {
            findings.push(finding(
                "byproduct_cross_feeding_candidate",
                FindingLevel::Interesting,
                "Byproduct uptake-pressure candidate",
                format!(
                    "{} byproduct species show uptake pressure with low final retention.",
                    calibration.cross_feeding_candidates
                ),
                calibration
                    .species
                    .iter()
                    .filter(|species| species.interpretation == "cross_feeding_candidate")
                    .map(|species| {
                        format!(
                            "ext{} {} produced {:.3}, retained {}",
                            species.species_index,
                            species.name,
                            species.produced_amount,
                            species
                                .retained_fraction
                                .map(|value| format!("{value:.3}"))
                                .unwrap_or_else(|| "n/a".to_string())
                        )
                    })
                    .collect(),
            ));
        }
        if calibration.public_pool_candidates > 0 {
            findings.push(finding(
                "byproduct_public_pool_candidate",
                FindingLevel::Warning,
                "Byproduct public-pool candidate",
                format!(
                    "{} byproduct species accumulated without detected uptake pressure.",
                    calibration.public_pool_candidates
                ),
                calibration
                    .species
                    .iter()
                    .filter(|species| species.interpretation == "public_pool_candidate")
                    .map(|species| {
                        format!(
                            "ext{} {} produced {:.3}, retained {}",
                            species.species_index,
                            species.name,
                            species.produced_amount,
                            species
                                .retained_fraction
                                .map(|value| format!("{value:.3}"))
                                .unwrap_or_else(|| "n/a".to_string())
                        )
                    })
                    .collect(),
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
    let cross_feeding = runs
        .iter()
        .filter_map(|run| {
            run.byproduct_cross_feeding_candidates
                .map(|count| (run, count))
        })
        .max_by_key(|(_, count)| *count);
    if let Some((run, count)) = cross_feeding
        && count > 0
    {
        findings.push(finding(
            "comparison_byproduct_cross_feeding_candidate",
            FindingLevel::Interesting,
            "Byproduct uptake-pressure candidate in comparison",
            format!(
                "{} had {} byproduct species with baseline-adjusted low retention and uptake pressure.",
                run.name, count
            ),
            Vec::new(),
        ));
    }
    let public_pool = runs
        .iter()
        .filter_map(|run| {
            run.byproduct_public_pool_candidates
                .map(|count| (run, count))
        })
        .max_by_key(|(_, count)| *count);
    if let Some((run, count)) = public_pool
        && count > 0
    {
        findings.push(finding(
            "comparison_byproduct_public_pool_candidate",
            FindingLevel::Warning,
            "Byproduct public-pool candidate in comparison",
            format!(
                "{} had {} byproduct species accumulating without detected uptake pressure.",
                run.name, count
            ),
            Vec::new(),
        ));
    }
    let excess_retained = runs
        .iter()
        .filter_map(|run| {
            run.byproduct_excess_retained_fraction
                .map(|fraction| (run, fraction))
        })
        .max_by(|(_, a), (_, b)| a.total_cmp(b));
    if let Some((run, fraction)) = excess_retained
        && fraction >= 0.75
    {
        findings.push(finding(
            "comparison_byproduct_excess_pool_candidate",
            FindingLevel::Warning,
            "Control-adjusted byproduct pool accumulation",
            format!(
                "{} retained {:.1}% of produced byproduct as field excess over the zero-byproduct baseline.",
                run.name,
                fraction * 100.0
            ),
            Vec::new(),
        ));
    }
    findings
}

impl RunComparisonEntry {
    fn from_analysis(analysis: &RunAnalysis) -> Self {
        let event_stoich = analysis
            .stoich_v2
            .as_ref()
            .filter(|stoich| stoich.events_present);
        let byproduct_calibration = analysis.byproduct_calibration.as_ref();
        Self {
            name: run_name(&analysis.run_dir),
            run_dir: analysis.run_dir.clone(),
            rng_seed: analysis.rng_seed,
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
            stoich_total_events: event_stoich.map(|stoich| stoich.total_events),
            stoich_imbalanced_events: event_stoich.map(|stoich| stoich.imbalanced_events),
            reaction_byproduct_events: event_stoich.map(|stoich| stoich.reaction_byproduct_events),
            reaction_byproduct_amount: event_stoich.map(|stoich| stoich.reaction_byproduct_amount),
            reaction_leakage_events: event_stoich.map(|stoich| stoich.reaction_leakage_events),
            reaction_leakage_energy_to_heat: event_stoich
                .map(|stoich| stoich.reaction_leakage_energy_to_heat),
            byproduct_final_pool: byproduct_calibration
                .and_then(|calibration| calibration.final_byproduct_pool),
            byproduct_retained_fraction: byproduct_calibration
                .and_then(|calibration| calibration.retained_fraction),
            byproduct_signed_excess_final_pool: None,
            byproduct_signed_excess_retained_fraction: None,
            byproduct_excess_final_pool: None,
            byproduct_excess_retained_fraction: None,
            byproduct_cross_feeding_candidates: byproduct_calibration
                .map(|calibration| calibration.cross_feeding_candidates),
            byproduct_public_pool_candidates: byproduct_calibration
                .map(|calibration| calibration.public_pool_candidates),
            byproduct_adjusted_species: Vec::new(),
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

fn chemical_role(profile: &SpeciesProfile) -> &'static str {
    chemical_role_from_descriptor(&profile.descriptor)
}

fn chemical_role_from_descriptor(descriptor: &SpeciesDescriptorProfile) -> &'static str {
    let composition = &descriptor.composition;
    if descriptor.work_coupling >= 0.75 && descriptor.bond_energy >= 0.5 {
        "work_currency"
    } else if composition.signal_group >= 0.5 {
        "signal"
    } else if composition.structural_group >= 0.5 {
        "structural"
    } else if descriptor.storage_density >= 0.4 {
        "storage"
    } else if composition.oxidizing_power >= 0.5 {
        "oxidant"
    } else if composition.reducing_power >= 0.5 {
        "reductant"
    } else if composition.toxin_group >= 0.5 {
        "toxin"
    } else if composition.carbon_backbone >= 0.5 {
        "carbon_source"
    } else {
        "inert"
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

fn display_f64(value: f64) -> String {
    format!("{value:.3}")
}

fn display_opt_pct(value: Option<f64>) -> String {
    value
        .map(|value| format!("{:.1}", value * 100.0))
        .unwrap_or_else(|| "n/a".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_dir(name: &str) -> PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!("{name}_{}_{}", std::process::id(), nanos));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn sampled_ticks_choose_first_middle_last() {
        assert_eq!(sample_ticks(&[]), Vec::<u64>::new());
        assert_eq!(sample_ticks(&[10]), vec![10]);
        assert_eq!(sample_ticks(&[0, 80, 159]), vec![0, 80, 159]);
        assert_eq!(sample_ticks(&[0, 10, 20, 30]), vec![0, 20, 30]);
    }

    #[test]
    fn run_rng_seed_parser_reads_summary_metadata() {
        let valid = test_dir("marl_analysis_seed_valid");
        fs::write(
            valid.join("summary.md"),
            "# Run Summary\n\n- RNG algorithm: chacha12\n- RNG seed: 41001\n",
        )
        .unwrap();
        assert_eq!(read_run_rng_seed(&valid), Some(41001));

        let malformed = test_dir("marl_analysis_seed_malformed");
        fs::write(malformed.join("summary.md"), "- RNG seed: nope\n").unwrap();
        assert_eq!(read_run_rng_seed(&malformed), None);

        let missing = test_dir("marl_analysis_seed_missing");
        fs::write(missing.join("summary.md"), "- Ticks: 10\n").unwrap();
        assert_eq!(read_run_rng_seed(&missing), None);

        let _ = fs::remove_dir_all(&valid);
        let _ = fs::remove_dir_all(&malformed);
        let _ = fs::remove_dir_all(&missing);
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
        assert_eq!(summary.species_profiles.len(), 5);
        assert_eq!(species0.per_z_mean, vec![0.0, 10.0]);
        assert_eq!(species0.total_concentration, 10.0);
        assert_eq!(species0.max_value, 10.0);
        assert_eq!(species0.name, "free_energy");
        assert_eq!(species0.descriptor.work_coupling, 1.0);
        assert_eq!(chemical_role(species0), "work_currency");
        assert_eq!(oxidant.per_z_mean, vec![1.0, 11.0]);
        assert_eq!(oxidant.descriptor.composition.oxidizing_power, 1.0);
        assert_eq!(chemical_role(oxidant), "oxidant");
    }

    #[test]
    fn field_summary_tracks_all_external_species_for_byproduct_pools() {
        let meta = RunMeta::new(1, 1, 1, 12, 0, true, false);
        let floats = (0..12).map(|species| species as f32).collect::<Vec<_>>();
        let bytes = floats
            .iter()
            .flat_map(|value| value.to_le_bytes())
            .collect::<Vec<_>>();

        let summary = summarize_field(3, &meta, &bytes).unwrap();

        assert_eq!(summary.species_profiles.len(), 12);
        let structural = summary
            .species_profiles
            .iter()
            .find(|profile| profile.species == 7)
            .unwrap();
        assert_eq!(structural.name, "structural");
        assert_eq!(structural.total_concentration, 7.0);
        assert_eq!(chemical_role(structural), "structural");
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
        assert_eq!(transporters.external_species.len(), 2);
        assert_eq!(transporters.external_species[0].ext_species, 1);
        assert_eq!(transporters.external_species[0].active_slots, 3);
        assert_eq!(transporters.external_species[0].uptake_active_slots, 3);
        assert_eq!(transporters.external_species[0].secretion_active_slots, 2);
        assert_eq!(transporters.external_species[0].uptake_dominant_slots, 1);
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
                external_species: common_pairs
                    .iter()
                    .map(|pair| TransportSpeciesSummary {
                        ext_species: pair.ext_species,
                        active_slots: pair.active_slots,
                        uptake_active_slots: pair.active_slots,
                        secretion_active_slots: (pair.avg_secrete_rate > 0.0)
                            .then_some(pair.active_slots)
                            .unwrap_or(0),
                        uptake_dominant_slots: pair.active_slots,
                        secretion_dominant_slots: 0,
                        gated_slots: pair.gated_slots,
                        uptake_rate_sum: pair.avg_uptake_rate * pair.active_slots as f64,
                        secretion_rate_sum: pair.avg_secrete_rate * pair.active_slots as f64,
                    })
                    .collect(),
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
            rng_seed: None,
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
            stoich_v2: None,
            byproduct_calibration: None,
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
    fn comparison_reports_include_stoich_v2_metrics() {
        let run = |name: &str, byproduct: Option<f64>| RunComparisonEntry {
            name: name.to_string(),
            run_dir: PathBuf::from(format!("/tmp/{name}")),
            rng_seed: None,
            final_population: Some(10),
            max_population: Some(12),
            growth_factor: Some(1.2),
            final_avg_energy: Some(0.5),
            top_fraction: Some(0.4),
            middle_fraction: Some(0.5),
            deep_fraction: Some(0.1),
            genotype_unique_count: None,
            genotype_dominant_fraction: None,
            transporter_active_per_cell: None,
            transporter_gated_slots: None,
            stoich_total_events: byproduct.map(|_| 100),
            stoich_imbalanced_events: byproduct.map(|_| 3),
            reaction_byproduct_events: byproduct.map(|_| 7),
            reaction_byproduct_amount: byproduct,
            reaction_leakage_events: byproduct.map(|_| 11),
            reaction_leakage_energy_to_heat: byproduct.map(|value| value / 2.0),
            byproduct_final_pool: byproduct.map(|value| value / 4.0),
            byproduct_retained_fraction: byproduct.map(|_| 0.25),
            byproduct_signed_excess_final_pool: byproduct.map(|value| value / 10.0),
            byproduct_signed_excess_retained_fraction: byproduct.map(|_| 0.1),
            byproduct_excess_final_pool: byproduct.map(|value| value / 5.0),
            byproduct_excess_retained_fraction: byproduct.map(|_| 0.2),
            byproduct_cross_feeding_candidates: byproduct.map(|_| 1),
            byproduct_public_pool_candidates: byproduct.map(|_| 0),
            byproduct_adjusted_species: Vec::new(),
        };
        let analysis = ComparisonAnalysis {
            runs: vec![run("with_events", Some(2.5)), run("without_events", None)],
            paired_byproduct_population: Vec::new(),
            byproduct_strength_aggregates: vec![ByproductStrengthAggregate {
                label: "byp008".to_string(),
                run_count: 2,
                mean_byproduct_amount: Some(10.0),
                mean_delta_final_population: Some(25.0),
                sd_delta_final_population: Some(5.0),
                min_delta_final_population: Some(20.0),
                max_delta_final_population: Some(30.0),
                mean_final_population_ratio: Some(1.2),
                mean_delta_growth_factor: Some(2.0),
                mean_delta_final_avg_energy: Some(-0.1),
                mean_signed_excess_retained_fraction: Some(-0.2),
                mean_clipped_excess_retained_fraction: Some(0.4),
                sd_clipped_excess_retained_fraction: Some(0.1),
                min_clipped_excess_retained_fraction: Some(0.3),
                max_clipped_excess_retained_fraction: Some(0.5),
                positive_population_delta_count: 1,
                positive_signed_excess_count: 0,
                clipped_excess_warning_count: 0,
            }],
            findings: Vec::new(),
            warnings: Vec::new(),
        };

        let terminal = render_comparison_terminal(&analysis);
        let markdown = render_comparison_markdown(&analysis);

        assert!(terminal.contains("byproduct=2.5000"));
        assert!(terminal.contains("gross_retained=0.250"));
        assert!(terminal.contains("signed_excess_retained=0.100"));
        assert!(terminal.contains("clipped_excess_retained=0.200"));
        assert!(terminal.contains("leakage_heat=1.2500"));
        assert!(terminal.contains("without_events:"));
        assert!(terminal.contains("byproduct=n/a"));
        assert!(markdown.contains("## Stoichiometry V2 Comparison"));
        assert!(markdown.contains(
            "| with_events | n/a | 10 | 12 | 1.200 | 0.500 | 40.0 | 50.0 | 10.0 | n/a | n/a | n/a | n/a |"
        ));
        assert!(markdown.contains(
            "| with_events | 100 | 3 | 7 | 2.500 | 0.625 | 0.250 | 0.250 | 0.100 | 0.500 | 0.200 | 1 | 0 | 11 | 1.250 |"
        ));
        assert!(markdown.contains(
            "| without_events | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |"
        ));
        assert!(markdown.contains("## Byproduct Strength Aggregates"));
        assert!(markdown.contains(
            "| byp008 | 2 | 10.000 | 25.000 | 5.000 | 20.000..30.000 | 1/2 | 1.200 | 2.000 | -0.100 | -0.200 | 0/2 | 0.400 | 0.100 | 0.300..0.500 | 0/2 |"
        ));
    }

    fn species_profile(species: usize, total_concentration: f64) -> SpeciesProfile {
        let descriptor = external_species_descriptor(species);
        SpeciesProfile {
            species,
            name: descriptor.name.to_string(),
            descriptor: descriptor.into(),
            total_concentration,
            max_value: total_concentration,
            surface_mean: total_concentration,
            middle_mean: total_concentration,
            deep_mean: total_concentration,
            per_z_mean: vec![total_concentration],
        }
    }

    fn minimal_analysis(
        name: &str,
        chemistry: Vec<SnapshotChemistry>,
        stoich_v2: Option<StoichV2Analysis>,
        byproduct_calibration: Option<ByproductCalibrationSummary>,
    ) -> RunAnalysis {
        RunAnalysis {
            run_dir: PathBuf::from(format!("/tmp/{name}")),
            grid: [1, 1, 1],
            rng_seed: None,
            available_ticks: Vec::new(),
            sampled_ticks: Vec::new(),
            trajectory: None,
            zonation: None,
            cells: None,
            cell_timeline: Vec::new(),
            chemistry,
            rulesets: None,
            ruleset_timeline: Vec::new(),
            transporter_pair_timeline: Vec::new(),
            stoich_v2,
            byproduct_calibration,
            findings: Vec::new(),
            warnings: Vec::new(),
        }
    }

    fn seeded(mut analysis: RunAnalysis, seed: u64) -> RunAnalysis {
        analysis.rng_seed = Some(seed);
        analysis
    }

    fn zero_byproduct_stoich() -> StoichV2Analysis {
        StoichV2Analysis {
            summary_present: false,
            events_present: true,
            total_ticks: None,
            enforcement: None,
            gross_residual_abs_sum: None,
            total_events: 1,
            imbalanced_events: 0,
            reaction_byproduct_events: 0,
            reaction_byproduct_amount: 0.0,
            reaction_byproduct_model_abs: 0.0,
            reaction_byproduct_residual_abs: 0.0,
            reaction_byproduct_by_species: Vec::new(),
            reaction_leakage_events: 0,
            reaction_leakage_amount: 0.0,
            reaction_leakage_energy_to_heat: 0.0,
            transport_flux_by_species: Vec::new(),
            reaction_byproduct_by_species_starter: Vec::new(),
            transport_flux_by_species_starter: Vec::new(),
        }
    }

    fn byproduct_stoich(species_index: i16, amount: f64) -> StoichV2Analysis {
        StoichV2Analysis {
            summary_present: false,
            events_present: true,
            total_ticks: None,
            enforcement: None,
            gross_residual_abs_sum: None,
            total_events: 1,
            imbalanced_events: 0,
            reaction_byproduct_events: 1,
            reaction_byproduct_amount: amount,
            reaction_byproduct_model_abs: 0.0,
            reaction_byproduct_residual_abs: 0.0,
            reaction_byproduct_by_species: vec![StoichSpeciesEventSummary {
                species_index,
                amount,
                events: 1,
            }],
            reaction_leakage_events: 0,
            reaction_leakage_amount: 0.0,
            reaction_leakage_energy_to_heat: 0.0,
            transport_flux_by_species: Vec::new(),
            reaction_byproduct_by_species_starter: Vec::new(),
            transport_flux_by_species_starter: Vec::new(),
        }
    }

    fn byproduct_calibration(
        species_index: i16,
        final_field_total: f64,
        produced_amount: f64,
    ) -> ByproductCalibrationSummary {
        ByproductCalibrationSummary {
            total_byproduct_amount: produced_amount,
            final_byproduct_pool: Some(final_field_total),
            retained_fraction: Some(final_field_total / produced_amount),
            field_tick: Some(10),
            ruleset_tick: None,
            missing_field_species: 0,
            cross_feeding_candidates: 0,
            public_pool_candidates: 0,
            species: vec![ByproductSpeciesCalibration {
                species_index,
                name: external_species_descriptor(species_index as usize)
                    .name
                    .to_string(),
                role: "carbon_source".to_string(),
                produced_amount,
                final_field_total: Some(final_field_total),
                retained_fraction: Some(final_field_total / produced_amount),
                uptake_active_slots: 0,
                secretion_active_slots: 0,
                uptake_rate_sum: 0.0,
                secretion_rate_sum: 0.0,
                interpretation: "accumulating_pool".to_string(),
            }],
        }
    }

    #[test]
    fn comparison_baseline_adjusts_byproduct_field_pool() {
        let baseline = minimal_analysis(
            "baseline",
            vec![SnapshotChemistry {
                tick: 10,
                nonfinite_values: 0,
                negative_values: 0,
                oxidant_penetration_z: None,
                reductant_penetration_z: None,
                redox_overlap_layers: 0,
                species_profiles: vec![species_profile(4, 100.0)],
            }],
            Some(StoichV2Analysis {
                summary_present: false,
                events_present: true,
                total_ticks: None,
                enforcement: None,
                gross_residual_abs_sum: None,
                total_events: 1,
                imbalanced_events: 0,
                reaction_byproduct_events: 0,
                reaction_byproduct_amount: 0.0,
                reaction_byproduct_model_abs: 0.0,
                reaction_byproduct_residual_abs: 0.0,
                reaction_byproduct_by_species: Vec::new(),
                reaction_leakage_events: 0,
                reaction_leakage_amount: 0.0,
                reaction_leakage_energy_to_heat: 0.0,
                transport_flux_by_species: Vec::new(),
                reaction_byproduct_by_species_starter: Vec::new(),
                transport_flux_by_species_starter: Vec::new(),
            }),
            None,
        );
        let byproduct = minimal_analysis(
            "byproduct",
            vec![SnapshotChemistry {
                tick: 10,
                nonfinite_values: 0,
                negative_values: 0,
                oxidant_penetration_z: None,
                reductant_penetration_z: None,
                redox_overlap_layers: 0,
                species_profiles: vec![species_profile(4, 112.0)],
            }],
            Some(StoichV2Analysis {
                summary_present: false,
                events_present: true,
                total_ticks: None,
                enforcement: None,
                gross_residual_abs_sum: None,
                total_events: 1,
                imbalanced_events: 0,
                reaction_byproduct_events: 1,
                reaction_byproduct_amount: 20.0,
                reaction_byproduct_model_abs: 0.0,
                reaction_byproduct_residual_abs: 0.0,
                reaction_byproduct_by_species: vec![StoichSpeciesEventSummary {
                    species_index: 4,
                    amount: 20.0,
                    events: 1,
                }],
                reaction_leakage_events: 0,
                reaction_leakage_amount: 0.0,
                reaction_leakage_energy_to_heat: 0.0,
                transport_flux_by_species: Vec::new(),
                reaction_byproduct_by_species_starter: Vec::new(),
                transport_flux_by_species_starter: Vec::new(),
            }),
            Some(ByproductCalibrationSummary {
                total_byproduct_amount: 20.0,
                final_byproduct_pool: Some(112.0),
                retained_fraction: Some(5.6),
                field_tick: Some(10),
                ruleset_tick: None,
                missing_field_species: 0,
                cross_feeding_candidates: 0,
                public_pool_candidates: 0,
                species: vec![ByproductSpeciesCalibration {
                    species_index: 4,
                    name: "organic".to_string(),
                    role: "carbon_source".to_string(),
                    produced_amount: 20.0,
                    final_field_total: Some(112.0),
                    retained_fraction: Some(5.6),
                    uptake_active_slots: 0,
                    secretion_active_slots: 0,
                    uptake_rate_sum: 0.0,
                    secretion_rate_sum: 0.0,
                    interpretation: "accumulating_pool".to_string(),
                }],
            }),
        );
        let analyses = vec![baseline, byproduct];
        let mut runs = analyses
            .iter()
            .map(RunComparisonEntry::from_analysis)
            .collect::<Vec<_>>();

        let warnings = apply_byproduct_baseline_adjustment(&analyses, &mut runs);

        assert!(warnings.is_empty());
        assert_eq!(runs[0].byproduct_excess_final_pool, None);
        assert_eq!(runs[1].byproduct_excess_final_pool, Some(12.0));
        assert_eq!(runs[1].byproduct_excess_retained_fraction, Some(0.6));
        assert_eq!(runs[1].byproduct_signed_excess_final_pool, Some(12.0));
        assert_eq!(runs[1].byproduct_signed_excess_retained_fraction, Some(0.6));
        let findings = classify_comparison_findings(&runs);
        assert!(
            !findings
                .iter()
                .any(|finding| finding.id == "comparison_byproduct_excess_pool_candidate")
        );
    }

    #[test]
    fn comparison_candidates_use_baseline_adjusted_species_excess() {
        let baseline = minimal_analysis(
            "baseline",
            vec![SnapshotChemistry {
                tick: 10,
                nonfinite_values: 0,
                negative_values: 0,
                oxidant_penetration_z: None,
                reductant_penetration_z: None,
                redox_overlap_layers: 0,
                species_profiles: vec![species_profile(4, 100.0), species_profile(7, 10.0)],
            }],
            Some(zero_byproduct_stoich()),
            None,
        );
        let byproduct = minimal_analysis(
            "byproduct",
            vec![SnapshotChemistry {
                tick: 10,
                nonfinite_values: 0,
                negative_values: 0,
                oxidant_penetration_z: None,
                reductant_penetration_z: None,
                redox_overlap_layers: 0,
                species_profiles: vec![species_profile(4, 100.5), species_profile(7, 9.0)],
            }],
            Some(StoichV2Analysis {
                summary_present: false,
                events_present: true,
                total_ticks: None,
                enforcement: None,
                gross_residual_abs_sum: None,
                total_events: 2,
                imbalanced_events: 0,
                reaction_byproduct_events: 2,
                reaction_byproduct_amount: 1.000000001,
                reaction_byproduct_model_abs: 0.0,
                reaction_byproduct_residual_abs: 0.0,
                reaction_byproduct_by_species: vec![
                    StoichSpeciesEventSummary {
                        species_index: 4,
                        amount: 1.0,
                        events: 1,
                    },
                    StoichSpeciesEventSummary {
                        species_index: 7,
                        amount: 1e-9,
                        events: 1,
                    },
                ],
                reaction_leakage_events: 0,
                reaction_leakage_amount: 0.0,
                reaction_leakage_energy_to_heat: 0.0,
                transport_flux_by_species: Vec::new(),
                reaction_byproduct_by_species_starter: Vec::new(),
                transport_flux_by_species_starter: Vec::new(),
            }),
            Some(ByproductCalibrationSummary {
                total_byproduct_amount: 1.000000001,
                final_byproduct_pool: Some(109.5),
                retained_fraction: Some(109.5),
                field_tick: Some(10),
                ruleset_tick: None,
                missing_field_species: 0,
                cross_feeding_candidates: 1,
                public_pool_candidates: 1,
                species: vec![
                    ByproductSpeciesCalibration {
                        species_index: 4,
                        name: "organic".to_string(),
                        role: "carbon_source".to_string(),
                        produced_amount: 1.0,
                        final_field_total: Some(100.5),
                        retained_fraction: Some(100.5),
                        uptake_active_slots: 0,
                        secretion_active_slots: 0,
                        uptake_rate_sum: 0.0,
                        secretion_rate_sum: 0.0,
                        interpretation: "public_pool_candidate".to_string(),
                    },
                    ByproductSpeciesCalibration {
                        species_index: 7,
                        name: "structural".to_string(),
                        role: "structural".to_string(),
                        produced_amount: 1e-9,
                        final_field_total: Some(9.0),
                        retained_fraction: Some(9.0e9),
                        uptake_active_slots: 1,
                        secretion_active_slots: 0,
                        uptake_rate_sum: 1.0,
                        secretion_rate_sum: 0.0,
                        interpretation: "cross_feeding_candidate".to_string(),
                    },
                ],
            }),
        );
        let analyses = vec![baseline, byproduct];
        let mut runs = analyses
            .iter()
            .map(RunComparisonEntry::from_analysis)
            .collect::<Vec<_>>();

        let warnings = apply_byproduct_baseline_adjustment(&analyses, &mut runs);

        assert!(warnings.is_empty());
        assert_eq!(runs[1].byproduct_retained_fraction, Some(109.5));
        assert_eq!(runs[1].byproduct_excess_final_pool, Some(0.5));
        assert!((runs[1].byproduct_excess_retained_fraction.unwrap() - 0.5).abs() < 1e-9);
        assert_eq!(runs[1].byproduct_cross_feeding_candidates, Some(0));
        assert_eq!(runs[1].byproduct_public_pool_candidates, Some(0));
        assert_eq!(
            runs[1].byproduct_adjusted_species[1].adjusted_interpretation,
            "trace_byproduct"
        );
    }

    #[test]
    fn comparison_baseline_adjustment_pairs_by_rng_seed() {
        let baseline_seed_1 = seeded(
            minimal_analysis(
                "baseline_s1",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 100.0)],
                }],
                Some(zero_byproduct_stoich()),
                None,
            ),
            41001,
        );
        let baseline_seed_2 = seeded(
            minimal_analysis(
                "baseline_s2",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 200.0)],
                }],
                Some(zero_byproduct_stoich()),
                None,
            ),
            41002,
        );
        let byproduct_seed_1 = seeded(
            minimal_analysis(
                "byproduct_s1",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 112.0)],
                }],
                Some(byproduct_stoich(4, 20.0)),
                Some(byproduct_calibration(4, 112.0, 20.0)),
            ),
            41001,
        );
        let byproduct_seed_2 = seeded(
            minimal_analysis(
                "byproduct_s2",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 215.0)],
                }],
                Some(byproduct_stoich(4, 20.0)),
                Some(byproduct_calibration(4, 215.0, 20.0)),
            ),
            41002,
        );
        let analyses = vec![
            baseline_seed_1,
            baseline_seed_2,
            byproduct_seed_1,
            byproduct_seed_2,
        ];
        let mut runs = analyses
            .iter()
            .map(RunComparisonEntry::from_analysis)
            .collect::<Vec<_>>();

        let warnings = apply_byproduct_baseline_adjustment(&analyses, &mut runs);

        assert!(warnings.is_empty());
        assert_eq!(runs[2].byproduct_excess_final_pool, Some(12.0));
        assert_eq!(runs[2].byproduct_excess_retained_fraction, Some(0.6));
        assert_eq!(runs[3].byproduct_excess_final_pool, Some(15.0));
        assert_eq!(runs[3].byproduct_excess_retained_fraction, Some(0.75));
    }

    #[test]
    fn comparison_summarizes_paired_byproduct_population_effects() {
        let run = |name: &str,
                   seed: u64,
                   byproduct_amount: f64,
                   final_population: u64,
                   growth_factor: f64,
                   final_avg_energy: f64|
         -> RunComparisonEntry {
            RunComparisonEntry {
                name: name.to_string(),
                run_dir: PathBuf::from(format!("/tmp/{name}")),
                rng_seed: Some(seed),
                final_population: Some(final_population),
                max_population: Some(final_population),
                growth_factor: Some(growth_factor),
                final_avg_energy: Some(final_avg_energy),
                top_fraction: None,
                middle_fraction: None,
                deep_fraction: None,
                genotype_unique_count: None,
                genotype_dominant_fraction: None,
                transporter_active_per_cell: None,
                transporter_gated_slots: None,
                stoich_total_events: Some(1),
                stoich_imbalanced_events: Some(0),
                reaction_byproduct_events: Some((byproduct_amount > 0.0) as u64),
                reaction_byproduct_amount: Some(byproduct_amount),
                reaction_leakage_events: Some(0),
                reaction_leakage_energy_to_heat: Some(0.0),
                byproduct_final_pool: None,
                byproduct_retained_fraction: None,
                byproduct_signed_excess_final_pool: None,
                byproduct_signed_excess_retained_fraction: None,
                byproduct_excess_final_pool: None,
                byproduct_excess_retained_fraction: Some(0.25),
                byproduct_cross_feeding_candidates: Some(0),
                byproduct_public_pool_candidates: Some(0),
                byproduct_adjusted_species: Vec::new(),
            }
        };
        let runs = vec![
            run("s1_byp000", 41001, 0.0, 100, 10.0, 0.5),
            run("s1_byp008", 41001, 8.0, 125, 12.5, 0.6),
        ];

        let pairs = summarize_paired_byproduct_population(&runs);

        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].baseline_name, "s1_byp000");
        assert_eq!(pairs[0].run_name, "s1_byp008");
        assert_eq!(pairs[0].delta_final_population, Some(25));
        assert_eq!(pairs[0].final_population_ratio, Some(1.25));
        assert_eq!(pairs[0].delta_growth_factor, Some(2.5));
        assert_eq!(pairs[0].delta_final_avg_energy, Some(0.09999999999999998));

        let aggregates = summarize_byproduct_strength_aggregates(&pairs);
        assert_eq!(aggregates.len(), 1);
        assert_eq!(aggregates[0].label, "byp008");
        assert_eq!(aggregates[0].run_count, 1);
        assert_eq!(aggregates[0].mean_delta_final_population, Some(25.0));
        assert_eq!(aggregates[0].positive_population_delta_count, 1);
    }

    #[test]
    fn comparison_baseline_adjustment_rejects_mismatched_explicit_seed() {
        let baseline = seeded(
            minimal_analysis(
                "baseline_s1",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 100.0)],
                }],
                Some(zero_byproduct_stoich()),
                None,
            ),
            41001,
        );
        let byproduct = seeded(
            minimal_analysis(
                "byproduct_s2",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 112.0)],
                }],
                Some(byproduct_stoich(4, 20.0)),
                Some(byproduct_calibration(4, 112.0, 20.0)),
            ),
            41002,
        );
        let analyses = vec![baseline, byproduct];
        let mut runs = analyses
            .iter()
            .map(RunComparisonEntry::from_analysis)
            .collect::<Vec<_>>();

        let warnings = apply_byproduct_baseline_adjustment(&analyses, &mut runs);

        assert!(warnings.iter().any(|warning| {
            warning.contains("run rng_seed 41002") && warning.contains("baseline rng_seed 41001")
        }));
        assert_eq!(runs[1].byproduct_excess_final_pool, None);
        assert_eq!(runs[1].byproduct_excess_retained_fraction, None);
    }

    #[test]
    fn comparison_baseline_adjustment_skips_duplicate_seed_baselines() {
        let baseline_a = seeded(
            minimal_analysis(
                "baseline_a",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 100.0)],
                }],
                Some(zero_byproduct_stoich()),
                None,
            ),
            41001,
        );
        let baseline_b = seeded(
            minimal_analysis(
                "baseline_b",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 105.0)],
                }],
                Some(zero_byproduct_stoich()),
                None,
            ),
            41001,
        );
        let byproduct = seeded(
            minimal_analysis(
                "byproduct_s1",
                vec![SnapshotChemistry {
                    tick: 10,
                    nonfinite_values: 0,
                    negative_values: 0,
                    oxidant_penetration_z: None,
                    reductant_penetration_z: None,
                    redox_overlap_layers: 0,
                    species_profiles: vec![species_profile(4, 112.0)],
                }],
                Some(byproduct_stoich(4, 20.0)),
                Some(byproduct_calibration(4, 112.0, 20.0)),
            ),
            41001,
        );
        let analyses = vec![baseline_a, baseline_b, byproduct];
        let mut runs = analyses
            .iter()
            .map(RunComparisonEntry::from_analysis)
            .collect::<Vec<_>>();

        let warnings = apply_byproduct_baseline_adjustment(&analyses, &mut runs);

        assert!(warnings.iter().any(|warning| {
            warning.contains("duplicate zero-byproduct baselines for rng_seed 41001")
        }));
        assert!(warnings.iter().any(|warning| {
            warning.contains("byproduct_s1")
                && warning.contains("multiple zero-byproduct baselines")
        }));
        assert_eq!(runs[2].byproduct_excess_final_pool, None);
        assert_eq!(runs[2].byproduct_excess_retained_fraction, None);
    }

    #[test]
    fn comparison_baseline_adjustment_requires_matching_field_tick() {
        let baseline = minimal_analysis(
            "baseline",
            vec![SnapshotChemistry {
                tick: 5,
                nonfinite_values: 0,
                negative_values: 0,
                oxidant_penetration_z: None,
                reductant_penetration_z: None,
                redox_overlap_layers: 0,
                species_profiles: vec![species_profile(4, 100.0)],
            }],
            Some(StoichV2Analysis {
                summary_present: false,
                events_present: true,
                total_ticks: None,
                enforcement: None,
                gross_residual_abs_sum: None,
                total_events: 1,
                imbalanced_events: 0,
                reaction_byproduct_events: 0,
                reaction_byproduct_amount: 0.0,
                reaction_byproduct_model_abs: 0.0,
                reaction_byproduct_residual_abs: 0.0,
                reaction_byproduct_by_species: Vec::new(),
                reaction_leakage_events: 0,
                reaction_leakage_amount: 0.0,
                reaction_leakage_energy_to_heat: 0.0,
                transport_flux_by_species: Vec::new(),
                reaction_byproduct_by_species_starter: Vec::new(),
                transport_flux_by_species_starter: Vec::new(),
            }),
            None,
        );
        let byproduct = minimal_analysis(
            "byproduct",
            vec![SnapshotChemistry {
                tick: 10,
                nonfinite_values: 0,
                negative_values: 0,
                oxidant_penetration_z: None,
                reductant_penetration_z: None,
                redox_overlap_layers: 0,
                species_profiles: vec![species_profile(4, 112.0)],
            }],
            Some(StoichV2Analysis {
                summary_present: false,
                events_present: true,
                total_ticks: None,
                enforcement: None,
                gross_residual_abs_sum: None,
                total_events: 1,
                imbalanced_events: 0,
                reaction_byproduct_events: 1,
                reaction_byproduct_amount: 20.0,
                reaction_byproduct_model_abs: 0.0,
                reaction_byproduct_residual_abs: 0.0,
                reaction_byproduct_by_species: vec![StoichSpeciesEventSummary {
                    species_index: 4,
                    amount: 20.0,
                    events: 1,
                }],
                reaction_leakage_events: 0,
                reaction_leakage_amount: 0.0,
                reaction_leakage_energy_to_heat: 0.0,
                transport_flux_by_species: Vec::new(),
                reaction_byproduct_by_species_starter: Vec::new(),
                transport_flux_by_species_starter: Vec::new(),
            }),
            Some(ByproductCalibrationSummary {
                total_byproduct_amount: 20.0,
                final_byproduct_pool: Some(112.0),
                retained_fraction: Some(5.6),
                field_tick: Some(10),
                ruleset_tick: None,
                missing_field_species: 0,
                cross_feeding_candidates: 0,
                public_pool_candidates: 0,
                species: vec![ByproductSpeciesCalibration {
                    species_index: 4,
                    name: "organic".to_string(),
                    role: "carbon_source".to_string(),
                    produced_amount: 20.0,
                    final_field_total: Some(112.0),
                    retained_fraction: Some(5.6),
                    uptake_active_slots: 0,
                    secretion_active_slots: 0,
                    uptake_rate_sum: 0.0,
                    secretion_rate_sum: 0.0,
                    interpretation: "accumulating_pool".to_string(),
                }],
            }),
        );
        let analyses = vec![baseline, byproduct];
        let mut runs = analyses
            .iter()
            .map(RunComparisonEntry::from_analysis)
            .collect::<Vec<_>>();

        let warnings = apply_byproduct_baseline_adjustment(&analyses, &mut runs);

        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("field tick"))
        );
        assert_eq!(runs[1].byproduct_excess_final_pool, None);
        assert_eq!(runs[1].byproduct_excess_retained_fraction, None);
    }

    #[test]
    fn byproduct_calibration_classifies_cross_feeding_and_public_pools() {
        let stoich = StoichV2Analysis {
            summary_present: false,
            events_present: true,
            total_ticks: None,
            enforcement: None,
            gross_residual_abs_sum: None,
            total_events: 2,
            imbalanced_events: 0,
            reaction_byproduct_events: 2,
            reaction_byproduct_amount: 11.0,
            reaction_byproduct_model_abs: 0.0,
            reaction_byproduct_residual_abs: 0.0,
            reaction_byproduct_by_species: vec![
                StoichSpeciesEventSummary {
                    species_index: 4,
                    amount: 10.0,
                    events: 1,
                },
                StoichSpeciesEventSummary {
                    species_index: 7,
                    amount: 1.0,
                    events: 1,
                },
            ],
            reaction_leakage_events: 0,
            reaction_leakage_amount: 0.0,
            reaction_leakage_energy_to_heat: 0.0,
            transport_flux_by_species: Vec::new(),
            reaction_byproduct_by_species_starter: Vec::new(),
            transport_flux_by_species_starter: Vec::new(),
        };
        let chemistry = SnapshotChemistry {
            tick: 20,
            nonfinite_values: 0,
            negative_values: 0,
            oxidant_penetration_z: None,
            reductant_penetration_z: None,
            redox_overlap_layers: 0,
            species_profiles: vec![species_profile(4, 2.0), species_profile(7, 1.0)],
        };
        let rulesets = rulesets_with_pairs(20, vec![pair(4, 4, 12)]);

        let calibration =
            summarize_byproduct_calibration(Some(&stoich), Some(&chemistry), Some(&rulesets))
                .unwrap();

        assert_eq!(calibration.cross_feeding_candidates, 1);
        assert_eq!(calibration.public_pool_candidates, 1);
        assert_eq!(calibration.field_tick, Some(20));
        assert_eq!(calibration.ruleset_tick, Some(20));
        assert_eq!(calibration.missing_field_species, 0);
        assert!((calibration.final_byproduct_pool.unwrap() - 3.0).abs() < 1e-9);
        assert!((calibration.retained_fraction.unwrap() - (3.0 / 11.0)).abs() < 1e-9);
        let organic = calibration
            .species
            .iter()
            .find(|species| species.species_index == 4)
            .unwrap();
        assert_eq!(organic.interpretation, "cross_feeding_candidate");
        assert_eq!(organic.uptake_active_slots, 12);
        let structural = calibration
            .species
            .iter()
            .find(|species| species.species_index == 7)
            .unwrap();
        assert_eq!(structural.interpretation, "public_pool_candidate");

        let analysis = RunAnalysis {
            run_dir: PathBuf::from("/tmp/example_run"),
            grid: [1, 1, 1],
            rng_seed: None,
            available_ticks: Vec::new(),
            sampled_ticks: Vec::new(),
            trajectory: None,
            zonation: None,
            cells: None,
            cell_timeline: Vec::new(),
            chemistry: vec![chemistry],
            rulesets: Some(rulesets),
            ruleset_timeline: Vec::new(),
            transporter_pair_timeline: Vec::new(),
            stoich_v2: Some(stoich),
            byproduct_calibration: Some(calibration),
            findings: Vec::new(),
            warnings: Vec::new(),
        };
        let findings = classify_run_findings(&analysis);
        assert!(
            findings
                .iter()
                .any(|finding| { finding.id == "byproduct_cross_feeding_candidate" })
        );
        assert!(
            findings
                .iter()
                .any(|finding| { finding.id == "byproduct_public_pool_candidate" })
        );
        let terminal = render_run_terminal(&analysis);
        let markdown = render_run_markdown(&analysis);
        assert!(terminal.contains("byproduct calibration: produced=11.0000"));
        assert!(markdown.contains("## Byproduct Calibration"));
        assert!(markdown.contains("cross_feeding_candidate"));
        assert!(markdown.contains("public_pool_candidate"));
    }

    #[test]
    fn byproduct_calibration_requires_complete_field_coverage_for_aggregate_retention() {
        let stoich = StoichV2Analysis {
            summary_present: false,
            events_present: true,
            total_ticks: None,
            enforcement: None,
            gross_residual_abs_sum: None,
            total_events: 2,
            imbalanced_events: 0,
            reaction_byproduct_events: 2,
            reaction_byproduct_amount: 2.0,
            reaction_byproduct_model_abs: 0.0,
            reaction_byproduct_residual_abs: 0.0,
            reaction_byproduct_by_species: vec![
                StoichSpeciesEventSummary {
                    species_index: 4,
                    amount: 1.0,
                    events: 1,
                },
                StoichSpeciesEventSummary {
                    species_index: 7,
                    amount: 1.0,
                    events: 1,
                },
            ],
            reaction_leakage_events: 0,
            reaction_leakage_amount: 0.0,
            reaction_leakage_energy_to_heat: 0.0,
            transport_flux_by_species: Vec::new(),
            reaction_byproduct_by_species_starter: Vec::new(),
            transport_flux_by_species_starter: Vec::new(),
        };
        let chemistry = SnapshotChemistry {
            tick: 9,
            nonfinite_values: 0,
            negative_values: 0,
            oxidant_penetration_z: None,
            reductant_penetration_z: None,
            redox_overlap_layers: 0,
            species_profiles: vec![species_profile(4, 0.25)],
        };

        let calibration = summarize_byproduct_calibration(Some(&stoich), Some(&chemistry), None)
            .expect("byproduct production should produce calibration summary");

        assert_eq!(calibration.field_tick, Some(9));
        assert_eq!(calibration.ruleset_tick, None);
        assert_eq!(calibration.missing_field_species, 1);
        assert_eq!(calibration.final_byproduct_pool, None);
        assert_eq!(calibration.retained_fraction, None);
        assert_eq!(
            calibration
                .species
                .iter()
                .find(|species| species.species_index == 4)
                .unwrap()
                .final_field_total,
            Some(0.25)
        );
        assert_eq!(
            calibration
                .species
                .iter()
                .find(|species| species.species_index == 7)
                .unwrap()
                .final_field_total,
            None
        );
    }

    #[test]
    fn stoich_v2_events_report_reaction_byproducts_and_leakage() {
        let events = "tick,stage,kind,actor_id,species_index,template_id,amount,model_c,model_h,model_o,model_s,model_redox,model_energy,reservoir_c,reservoir_h,reservoir_o,reservoir_s,reservoir_redox,reservoir_energy,residual_abs,balanced\n\
0,reactions,reaction_byproduct,7,4,255,0.500000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1\n\
1,reactions,reaction_byproduct,8,7,255,0.250000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1\n\
1,reactions,reaction_leakage,8,0,255,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,-0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,1\n";

        let summary = parse_stoich_v2_events(events).unwrap();

        assert_eq!(summary.total_events, 3);
        assert_eq!(summary.imbalanced_events, 0);
        assert_eq!(summary.byproduct.events, 2);
        assert!((summary.byproduct.amount - 0.75).abs() < 1e-9);
        assert_eq!(summary.leakage.events, 1);
        assert!((summary.leakage.reservoir_energy - 0.125).abs() < 1e-9);

        let analysis = RunAnalysis {
            run_dir: PathBuf::from("/tmp/example_run"),
            grid: [1, 1, 1],
            rng_seed: None,
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
            stoich_v2: Some(StoichV2Analysis {
                summary_present: false,
                events_present: true,
                total_ticks: None,
                enforcement: None,
                gross_residual_abs_sum: None,
                total_events: summary.total_events,
                imbalanced_events: summary.imbalanced_events,
                reaction_byproduct_events: summary.byproduct.events,
                reaction_byproduct_amount: summary.byproduct.amount,
                reaction_byproduct_model_abs: summary.byproduct.model_abs,
                reaction_byproduct_residual_abs: summary.byproduct.residual_abs,
                reaction_byproduct_by_species: summary.byproduct_by_species.into_values().collect(),
                reaction_leakage_events: summary.leakage.events,
                reaction_leakage_amount: summary.leakage.amount,
                reaction_leakage_energy_to_heat: summary.leakage.reservoir_energy,
                transport_flux_by_species: Vec::new(),
                reaction_byproduct_by_species_starter: Vec::new(),
                transport_flux_by_species_starter: Vec::new(),
            }),
            byproduct_calibration: None,
            findings: Vec::new(),
            warnings: Vec::new(),
        };
        let terminal = render_run_terminal(&analysis);
        let markdown = render_run_markdown(&analysis);

        assert!(terminal.contains("byproduct 0.7500 in 2 events"));
        assert!(terminal.contains("leakage 0.1250 to heat"));
        assert!(markdown.contains("## Stoichiometry V2"));
        assert!(markdown.contains("| 4 | 1 | 0.500000 |"));
    }

    #[test]
    fn stoich_v2_analysis_handles_missing_and_summary_only_files() {
        let no_stoich_dir = PathBuf::from("/tmp/marl_analysis_no_stoich_v2_test");
        let _ = fs::remove_dir_all(&no_stoich_dir);
        fs::create_dir_all(&no_stoich_dir).unwrap();
        assert!(read_stoich_v2_analysis(&no_stoich_dir).unwrap().is_none());

        let summary_only_dir = PathBuf::from("/tmp/marl_analysis_summary_only_stoich_v2_test");
        let _ = fs::remove_dir_all(&summary_only_dir);
        fs::create_dir_all(&summary_only_dir).unwrap();
        fs::write(
            summary_only_dir.join("stoich_v2_summary.json"),
            r#"{
                "schema_version": 2,
                "total_ticks": 5,
                "enforcement": "audit",
                "gross_residual_abs_sum": 1.25
            }"#,
        )
        .unwrap();

        let summary = read_stoich_v2_analysis(&summary_only_dir)
            .unwrap()
            .expect("summary-only stoich analysis should be present");
        assert!(summary.summary_present);
        assert!(!summary.events_present);
        assert_eq!(summary.total_ticks, Some(5));
        assert_eq!(summary.enforcement.as_deref(), Some("audit"));
        assert_eq!(summary.gross_residual_abs_sum, Some(1.25));

        let analysis = RunAnalysis {
            run_dir: summary_only_dir.clone(),
            grid: [1, 1, 1],
            rng_seed: None,
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
            stoich_v2: Some(summary),
            byproduct_calibration: None,
            findings: Vec::new(),
            warnings: Vec::new(),
        };
        let terminal = render_run_terminal(&analysis);
        let markdown = render_run_markdown(&analysis);

        assert!(terminal.contains("event-level reaction totals unavailable"));
        assert!(!terminal.contains("byproduct 0.0000"));
        assert!(markdown.contains("Event-level reaction byproduct/leakage totals: unavailable"));
        assert!(!markdown.contains("Reaction byproducts: 0 events"));

        let _ = fs::remove_dir_all(&no_stoich_dir);
        let _ = fs::remove_dir_all(&summary_only_dir);
    }

    #[test]
    fn stoich_v2_analysis_uses_summary_calibration_totals_without_event_csv() {
        let dir = test_dir("marl_analysis_summary_stoich_totals");
        fs::write(
            dir.join("stoich_v2_summary.json"),
            r#"{
                "schema_version": 2,
                "total_ticks": 8,
                "enforcement": "audit",
                "ledger": {
                    "reaction_byproduct": {
                        "events": 3,
                        "amount": 1.25,
                        "model_abs": 0.5,
                        "residual_abs": 0.0,
                        "reservoir_energy": 0.0
                    },
                    "reaction_byproduct_by_species": [
                        { "species_index": 0, "amount": 0.0, "events": 0 },
                        { "species_index": 4, "amount": 1.25, "events": 3 },
                        { "species_index": 70000, "amount": 99.0, "events": 1 }
                    ],
                    "reaction_byproduct_by_species_starter": [
                        { "species_index": 4, "starter_type": 1, "amount": 1.25, "events": 3 },
                        { "species_index": 4, "starter_type": 99, "amount": 99.0, "events": 1 },
                        { "species_index": 70000, "starter_type": 1, "amount": 99.0, "events": 1 }
                    ],
                    "reaction_leakage": {
                        "events": 2,
                        "amount": 0.75,
                        "model_abs": 0.75,
                        "residual_abs": 0.0,
                        "reservoir_energy": 0.75
                    },
                    "transport_flux_by_species": [
                        {
                            "species_index": 0,
                            "uptake_amount": 2.5,
                            "uptake_events": 4,
                            "secretion_amount": 0.25,
                            "secretion_events": 1
                        },
                        {
                            "species_index": 4,
                            "uptake_amount": 0.0,
                            "uptake_events": 0,
                            "secretion_amount": 1.5,
                            "secretion_events": 3
                        },
                        {
                            "species_index": 70000,
                            "uptake_amount": 99.0,
                            "uptake_events": 1,
                            "secretion_amount": 0.0,
                            "secretion_events": 0
                        }
                    ],
                    "transport_flux_by_species_starter": [
                        {
                            "species_index": 0,
                            "starter_type": 2,
                            "uptake_amount": 2.5,
                            "uptake_events": 4,
                            "secretion_amount": 0.25,
                            "secretion_events": 1
                        },
                        {
                            "species_index": 4,
                            "starter_type": 1,
                            "uptake_amount": 0.0,
                            "uptake_events": 0,
                            "secretion_amount": 1.5,
                            "secretion_events": 3
                        },
                        {
                            "species_index": 0,
                            "starter_type": 99,
                            "uptake_amount": 99.0,
                            "uptake_events": 1,
                            "secretion_amount": 0.0,
                            "secretion_events": 0
                        },
                        {
                            "species_index": 70000,
                            "starter_type": 1,
                            "uptake_amount": 99.0,
                            "uptake_events": 1,
                            "secretion_amount": 0.0,
                            "secretion_events": 0
                        }
                    ]
                },
                "stages": [
                    {
                        "stage": "reactions",
                        "summary": {
                            "event_count": 9,
                            "strict_rejection_count": 0,
                            "imbalanced_event_count": 1,
                            "gross_model_abs": 0.0,
                            "gross_reservoir_abs": 0.0,
                            "gross_residual_abs": 0.0,
                            "net_model_delta": { "c": 0.0, "h": 0.0, "o": 0.0, "s": 0.0, "redox": 0.0, "energy": 0.0 },
                            "net_reservoir_delta": { "c": 0.0, "h": 0.0, "o": 0.0, "s": 0.0, "redox": 0.0, "energy": 0.0 }
                        }
                    }
                ],
                "gross_residual_abs_sum": 0.0
            }"#,
        )
        .unwrap();

        let analysis = read_stoich_v2_analysis(&dir)
            .unwrap()
            .expect("summary-only compact stoich totals should be available");

        assert!(analysis.summary_present);
        assert!(analysis.events_present);
        assert_eq!(analysis.total_events, 9);
        assert_eq!(analysis.imbalanced_events, 1);
        assert_eq!(analysis.reaction_byproduct_events, 3);
        assert!((analysis.reaction_byproduct_amount - 1.25).abs() < 1e-9);
        assert_eq!(analysis.reaction_byproduct_by_species.len(), 1);
        assert_eq!(analysis.reaction_byproduct_by_species[0].species_index, 4);
        assert_eq!(analysis.reaction_byproduct_by_species_starter.len(), 1);
        assert_eq!(
            analysis.reaction_byproduct_by_species_starter[0].starter_type,
            1
        );
        assert!((analysis.reaction_byproduct_by_species_starter[0].amount - 1.25).abs() < 1e-9);
        assert_eq!(analysis.reaction_leakage_events, 2);
        assert!((analysis.reaction_leakage_energy_to_heat - 0.75).abs() < 1e-9);
        assert_eq!(analysis.transport_flux_by_species.len(), 2);
        assert_eq!(analysis.transport_flux_by_species[0].species_index, 0);
        assert!((analysis.transport_flux_by_species[0].uptake_amount - 2.5).abs() < 1e-9);
        assert_eq!(analysis.transport_flux_by_species[0].uptake_events, 4);
        assert!((analysis.transport_flux_by_species[0].secretion_amount - 0.25).abs() < 1e-9);
        assert_eq!(analysis.transport_flux_by_species[1].species_index, 4);
        assert!((analysis.transport_flux_by_species[1].secretion_amount - 1.5).abs() < 1e-9);
        assert_eq!(analysis.transport_flux_by_species_starter.len(), 2);
        assert_eq!(
            analysis.transport_flux_by_species_starter[0].starter_type,
            2
        );
        assert!((analysis.transport_flux_by_species_starter[0].uptake_amount - 2.5).abs() < 1e-9);
        assert_eq!(
            analysis.transport_flux_by_species_starter[1].species_index,
            4
        );
        assert_eq!(
            analysis.transport_flux_by_species_starter[1].starter_type,
            1
        );

        let run = minimal_analysis("summary_stoich_render", Vec::new(), Some(analysis), None);
        let terminal = render_run_terminal(&run);
        let markdown = render_run_markdown(&run);
        assert!(
            terminal
                .contains("reaction byproduct producers: ext4(organic) chemolithotroph:1.2500/3")
        );
        assert!(terminal.contains(
            "transport flux by starter: ext0(free_energy) anaerobe: uptake 2.5000/4 secretion 0.2500/1"
        ));
        assert!(markdown.contains("### Reaction Byproduct Producers By Starter"));
        assert!(markdown.contains("### Transport Flux By Species And Starter"));
        assert!(markdown.contains("| 4 | organic | chemolithotroph | 3 | 1.250000 |"));
        assert!(
            markdown.contains(
                "| 0 | free_energy | anaerobe | 4 | 2.500000 | 1 | 0.250000 | 2.250000 |"
            )
        );

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn stoich_v2_analysis_merges_csv_events_with_summary_transport_flux() {
        let dir = test_dir("marl_analysis_csv_summary_transport_flux");
        fs::write(
            dir.join("stoich_v2_summary.json"),
            r#"{
                "schema_version": 2,
                "total_ticks": 3,
                "enforcement": "audit",
                "ledger": {
                    "reaction_byproduct": {
                        "events": 99,
                        "amount": 99.0,
                        "model_abs": 99.0,
                        "residual_abs": 99.0,
                        "reservoir_energy": 0.0
                    },
                    "reaction_byproduct_by_species": [
                        { "species_index": 4, "amount": 99.0, "events": 99 }
                    ],
                    "reaction_byproduct_by_species_starter": [
                        { "species_index": 4, "starter_type": 0, "amount": 99.0, "events": 99 }
                    ],
                    "reaction_leakage": {
                        "events": 0,
                        "amount": 0.0,
                        "model_abs": 0.0,
                        "residual_abs": 0.0,
                        "reservoir_energy": 0.0
                    },
                    "transport_flux_by_species": [
                        {
                            "species_index": 0,
                            "uptake_amount": 0.125,
                            "uptake_events": 2,
                            "secretion_amount": 4.5,
                            "secretion_events": 7
                        }
                    ],
                    "transport_flux_by_species_starter": [
                        {
                            "species_index": 0,
                            "starter_type": 0,
                            "uptake_amount": 0.125,
                            "uptake_events": 2,
                            "secretion_amount": 4.5,
                            "secretion_events": 7
                        }
                    ]
                },
                "stages": [
                    {
                        "stage": "reactions",
                        "summary": {
                            "event_count": 99,
                            "strict_rejection_count": 0,
                            "imbalanced_event_count": 0,
                            "gross_model_abs": 0.0,
                            "gross_reservoir_abs": 0.0,
                            "gross_residual_abs": 0.0,
                            "net_model_delta": { "c": 0.0, "h": 0.0, "o": 0.0, "s": 0.0, "redox": 0.0, "energy": 0.0 },
                            "net_reservoir_delta": { "c": 0.0, "h": 0.0, "o": 0.0, "s": 0.0, "redox": 0.0, "energy": 0.0 }
                        }
                    }
                ],
                "gross_residual_abs_sum": 0.0
            }"#,
        )
        .unwrap();
        fs::write(
            dir.join("stoich_v2_events.csv"),
            "tick,stage,kind,actor_id,species_index,template_id,amount,model_c,model_h,model_o,model_s,model_redox,model_energy,reservoir_c,reservoir_h,reservoir_o,reservoir_s,reservoir_redox,reservoir_energy,residual_abs,balanced\n\
0,reactions,reaction_byproduct,7,4,255,0.500000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1\n",
        )
        .unwrap();

        let analysis = read_stoich_v2_analysis(&dir)
            .unwrap()
            .expect("csv plus summary stoich analysis should be present");

        assert!(analysis.summary_present);
        assert!(analysis.events_present);
        assert_eq!(analysis.total_ticks, Some(3));
        assert_eq!(analysis.total_events, 1);
        assert_eq!(analysis.reaction_byproduct_events, 1);
        assert!((analysis.reaction_byproduct_amount - 0.5).abs() < 1e-9);
        assert_eq!(analysis.reaction_byproduct_by_species_starter.len(), 0);
        assert_eq!(analysis.transport_flux_by_species.len(), 1);
        let flux = &analysis.transport_flux_by_species[0];
        assert_eq!(flux.species_index, 0);
        assert!((flux.uptake_amount - 0.125).abs() < 1e-9);
        assert_eq!(flux.uptake_events, 2);
        assert!((flux.secretion_amount - 4.5).abs() < 1e-9);
        assert_eq!(flux.secretion_events, 7);
        assert_eq!(analysis.transport_flux_by_species_starter.len(), 1);
        let starter_flux = &analysis.transport_flux_by_species_starter[0];
        assert_eq!(starter_flux.species_index, 0);
        assert_eq!(starter_flux.starter_type, 0);
        assert!((starter_flux.secretion_amount - 4.5).abs() < 1e-9);

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn stoich_v2_analysis_allows_csv_rounding_when_merging_summary_byproduct_starters() {
        let dir = test_dir("marl_analysis_csv_summary_rounding_byproduct_starters");
        fs::write(
            dir.join("stoich_v2_summary.json"),
            r#"{
                "schema_version": 2,
                "total_ticks": 3,
                "enforcement": "audit",
                "ledger": {
                    "reaction_byproduct": {
                        "events": 3,
                        "amount": 0.0000015,
                        "model_abs": 0.0,
                        "residual_abs": 0.0,
                        "reservoir_energy": 0.0
                    },
                    "reaction_byproduct_by_species": [
                        { "species_index": 4, "amount": 0.0000015, "events": 3 }
                    ],
                    "reaction_byproduct_by_species_starter": [
                        { "species_index": 4, "starter_type": 2, "amount": 0.0000015, "events": 3 }
                    ],
                    "reaction_leakage": {
                        "events": 0,
                        "amount": 0.0,
                        "model_abs": 0.0,
                        "residual_abs": 0.0,
                        "reservoir_energy": 0.0
                    },
                    "transport_flux_by_species": [],
                    "transport_flux_by_species_starter": []
                },
                "stages": [
                    {
                        "stage": "reactions",
                        "summary": {
                            "event_count": 3,
                            "strict_rejection_count": 0,
                            "imbalanced_event_count": 0,
                            "gross_model_abs": 0.0,
                            "gross_reservoir_abs": 0.0,
                            "gross_residual_abs": 0.0,
                            "net_model_delta": { "c": 0.0, "h": 0.0, "o": 0.0, "s": 0.0, "redox": 0.0, "energy": 0.0 },
                            "net_reservoir_delta": { "c": 0.0, "h": 0.0, "o": 0.0, "s": 0.0, "redox": 0.0, "energy": 0.0 }
                        }
                    }
                ],
                "gross_residual_abs_sum": 0.0
            }"#,
        )
        .unwrap();
        fs::write(
            dir.join("stoich_v2_events.csv"),
            "tick,stage,kind,actor_id,species_index,template_id,amount,model_c,model_h,model_o,model_s,model_redox,model_energy,reservoir_c,reservoir_h,reservoir_o,reservoir_s,reservoir_redox,reservoir_energy,residual_abs,balanced\n\
0,reactions,reaction_byproduct,7,4,255,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1\n\
1,reactions,reaction_byproduct,7,4,255,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1\n\
2,reactions,reaction_byproduct,7,4,255,0.000001,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1\n",
        )
        .unwrap();

        let analysis = read_stoich_v2_analysis(&dir)
            .unwrap()
            .expect("csv plus same-run rounded summary should be present");

        assert_eq!(analysis.reaction_byproduct_events, 3);
        assert!((analysis.reaction_byproduct_amount - 0.000003).abs() < 1e-12);
        assert_eq!(analysis.reaction_byproduct_by_species_starter.len(), 1);
        assert_eq!(
            analysis.reaction_byproduct_by_species_starter[0].starter_type,
            2
        );
        assert!(
            (analysis.reaction_byproduct_by_species_starter[0].amount - 0.0000015).abs() < 1e-12
        );

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn stoich_v2_analysis_handles_events_without_summary() {
        let events_only_dir = PathBuf::from("/tmp/marl_analysis_events_only_stoich_v2_test");
        let _ = fs::remove_dir_all(&events_only_dir);
        fs::create_dir_all(&events_only_dir).unwrap();
        fs::write(
            events_only_dir.join("stoich_v2_events.csv"),
            "tick,stage,kind,actor_id,species_index,template_id,amount,model_c,model_h,model_o,model_s,model_redox,model_energy,reservoir_c,reservoir_h,reservoir_o,reservoir_s,reservoir_redox,reservoir_energy,residual_abs,balanced\n\
0,reactions,reaction_byproduct,7,4,255,0.500000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1\n",
        )
        .unwrap();

        let summary = read_stoich_v2_analysis(&events_only_dir)
            .unwrap()
            .expect("events-only stoich analysis should be present");
        assert!(!summary.summary_present);
        assert!(summary.events_present);
        assert_eq!(summary.total_ticks, None);
        assert_eq!(summary.reaction_byproduct_events, 1);
        assert!((summary.reaction_byproduct_amount - 0.5).abs() < 1e-9);

        let _ = fs::remove_dir_all(&events_only_dir);
    }

    #[test]
    fn report_writer_with_no_outputs_does_not_create_directory() {
        let out_dir = PathBuf::from("/tmp/marl_analysis_no_output_test");
        let _ = fs::remove_dir_all(&out_dir);
        let analysis = RunAnalysis {
            run_dir: PathBuf::from("/tmp/example_run"),
            grid: [1, 1, 1],
            rng_seed: None,
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
            stoich_v2: None,
            byproduct_calibration: None,
            findings: Vec::new(),
            warnings: Vec::new(),
        };

        write_run_reports(&analysis, &out_dir, false, false).unwrap();
        assert!(!out_dir.exists());
    }
}
