use anyhow::Result;
use clap::{Parser, Subcommand};
use ruspladder::annotation::AnnotationFilters;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "ruspladder",
    version,
    about = "SplAdder Rust migration (in development)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Build, merge, quantify graphs and report alternative splicing events.
    Build(Box<ruspladder::build_cli::BuildArgs>),
    /// Differentially test counted events between two conditions.
    Test(Box<ruspladder::test_cli::TestArgs>),
    /// Prepare annotation graphs. Alignment preparation is still being migrated.
    Prep {
        #[arg(short, long)]
        annotation: PathBuf,
        #[arg(long)]
        filter_overlap_genes: bool,
        #[arg(long)]
        filter_overlap_exons: bool,
        #[arg(long)]
        filter_overlap_transcripts: bool,
        #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..=64))]
        parallel: u16,
        #[arg(short, long)]
        verbose: bool,
    },
}

fn run() -> Result<()> {
    match Cli::parse().command {
        Command::Build(options) => ruspladder::build_cli::run(&options)?,
        Command::Test(options) => ruspladder::test_cli::run(&options)?,
        Command::Prep {
            annotation,
            filter_overlap_genes,
            filter_overlap_exons,
            filter_overlap_transcripts,
            parallel,
            verbose,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(parallel as usize)
                .build()?;
            let result = pool.install(|| {
                ruspladder::build_cli::prepare_annotation(
                    &annotation,
                    AnnotationFilters {
                        overlap_genes: filter_overlap_genes,
                        overlap_exons: filter_overlap_exons,
                        overlap_transcripts: filter_overlap_transcripts,
                    },
                )
            })?;
            if verbose {
                eprintln!(
                    "Prepared {} genes from {}",
                    result.len(),
                    annotation.display()
                );
            }
        }
    }
    Ok(())
}

fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("ruspladder: {error:#}");
            std::process::ExitCode::FAILURE
        }
    }
}
