use std::error::Error;
use std::path::PathBuf;

use marl_analysis::{
    AnalysisConfig, ScanMode, analyze_run, compare_runs, render_comparison_terminal,
    render_run_terminal, write_comparison_reports, write_run_reports,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputMode {
    All,
    JsonOnly,
    MarkdownOnly,
    TerminalOnly,
}

impl OutputMode {
    fn writes_json(self) -> bool {
        matches!(self, Self::All | Self::JsonOnly)
    }

    fn writes_markdown(self) -> bool {
        matches!(self, Self::All | Self::MarkdownOnly)
    }

    fn prints_terminal(self) -> bool {
        matches!(self, Self::All | Self::TerminalOnly)
    }
}

#[derive(Debug)]
struct CliOptions {
    cfg: AnalysisConfig,
    out_dir: Option<PathBuf>,
    output_mode: OutputMode,
    paths: Vec<PathBuf>,
}

fn main() {
    if let Err(err) = run() {
        eprintln!("error: {err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let Some(command) = args.next() else {
        print_help();
        return Ok(());
    };
    let rest: Vec<String> = args.collect();
    match command.as_str() {
        "run" => run_one(parse_options(rest, true)?)?,
        "compare" => run_compare(parse_options(rest, false)?)?,
        "-h" | "--help" | "help" => print_help(),
        other => return Err(format!("unknown command {other:?}").into()),
    }
    Ok(())
}

fn run_one(options: CliOptions) -> Result<(), Box<dyn Error>> {
    if options.paths.len() != 1 {
        return Err("run expects exactly one run directory".into());
    }
    let run_dir = &options.paths[0];
    let analysis = analyze_run(run_dir, &options.cfg)?;
    let out_dir = options.out_dir.unwrap_or_else(|| run_dir.join("analysis"));
    write_run_reports(
        &analysis,
        &out_dir,
        options.output_mode.writes_json(),
        options.output_mode.writes_markdown(),
    )?;
    if options.output_mode.prints_terminal() {
        print!("{}", render_run_terminal(&analysis));
        if options.output_mode != OutputMode::TerminalOnly {
            println!("reports written to {}", out_dir.display());
        }
    }
    Ok(())
}

fn run_compare(options: CliOptions) -> Result<(), Box<dyn Error>> {
    if options.paths.len() < 2 {
        return Err("compare expects at least two run directories".into());
    }
    let analysis = compare_runs(&options.paths, &options.cfg)?;
    let out_dir = options.out_dir.unwrap_or_else(|| PathBuf::from("analysis"));
    write_comparison_reports(
        &analysis,
        &out_dir,
        options.output_mode.writes_json(),
        options.output_mode.writes_markdown(),
    )?;
    if options.output_mode.prints_terminal() {
        print!("{}", render_comparison_terminal(&analysis));
        if options.output_mode != OutputMode::TerminalOnly {
            println!("reports written to {}", out_dir.display());
        }
    }
    Ok(())
}

fn parse_options(args: Vec<String>, allow_single_path: bool) -> Result<CliOptions, Box<dyn Error>> {
    let mut cfg = AnalysisConfig::default();
    let mut out_dir = None;
    let mut output_mode = OutputMode::All;
    let mut paths = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--out-dir" => {
                let value = args
                    .get(i + 1)
                    .ok_or("--out-dir requires a directory argument")?;
                out_dir = Some(PathBuf::from(value));
                i += 2;
            }
            "--all-snapshots" => {
                cfg.scan_mode = ScanMode::All;
                i += 1;
            }
            "--latest-only" => {
                cfg.scan_mode = ScanMode::LatestOnly;
                i += 1;
            }
            "--ticks" => {
                let value = args.get(i + 1).ok_or("--ticks requires a comma list")?;
                cfg.scan_mode = ScanMode::Explicit(parse_tick_list(value)?);
                i += 2;
            }
            "--no-rulesets" => {
                cfg.include_rulesets = false;
                i += 1;
            }
            "--json-only" => {
                output_mode = set_output_mode(output_mode, OutputMode::JsonOnly)?;
                i += 1;
            }
            "--markdown-only" => {
                output_mode = set_output_mode(output_mode, OutputMode::MarkdownOnly)?;
                i += 1;
            }
            "--terminal-only" => {
                output_mode = set_output_mode(output_mode, OutputMode::TerminalOnly)?;
                i += 1;
            }
            "-h" | "--help" => {
                print_help();
                std::process::exit(0);
            }
            value if value.starts_with('-') => {
                return Err(format!("unknown flag {value:?}").into());
            }
            value => {
                paths.push(PathBuf::from(value));
                i += 1;
            }
        }
    }

    if allow_single_path && paths.is_empty() {
        return Err("missing run directory".into());
    }

    Ok(CliOptions {
        cfg,
        out_dir,
        output_mode,
        paths,
    })
}

fn set_output_mode(current: OutputMode, next: OutputMode) -> Result<OutputMode, Box<dyn Error>> {
    if current != OutputMode::All && current != next {
        return Err("choose only one of --json-only, --markdown-only, or --terminal-only".into());
    }
    Ok(next)
}

fn parse_tick_list(raw: &str) -> Result<Vec<u64>, Box<dyn Error>> {
    let ticks: Result<Vec<_>, _> = raw
        .split(',')
        .filter(|part| !part.trim().is_empty())
        .map(|part| part.trim().parse::<u64>())
        .collect();
    let mut ticks = ticks?;
    if ticks.is_empty() {
        return Err("--ticks must name at least one tick".into());
    }
    ticks.sort_unstable();
    ticks.dedup();
    Ok(ticks)
}

fn print_help() {
    println!(
        "marl-analyze\n\n\
         Usage:\n  \
           marl-analyze run <RUN_DIR> [options]\n  \
           marl-analyze compare <RUN_DIR> <RUN_DIR>... [options]\n\n\
         Options:\n  \
           --out-dir <DIR>       Write reports to DIR\n  \
           --all-snapshots       Inspect every available binary field snapshot\n  \
           --latest-only         Inspect only the latest binary field snapshot\n  \
           --ticks <A,B,C>       Inspect explicit snapshot ticks\n  \
           --no-rulesets         Skip full ruleset sidecar analysis\n  \
           --json-only           Write JSON only\n  \
           --markdown-only       Write Markdown only\n  \
           --terminal-only       Print terminal report only\n"
    );
}
