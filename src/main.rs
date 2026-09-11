use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use ruspladder::{
    annotation::{AnnotationFilters, read_annotation},
    cache,
};
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
        #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..))]
        parallel: u16,
        #[arg(short, long)]
        verbose: bool,
    },
}

fn run() -> Result<()> {
    match Cli::parse().command {
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
                read_annotation(
                    &annotation,
                    AnnotationFilters {
                        overlap_genes: filter_overlap_genes,
                        overlap_exons: filter_overlap_exons,
                        overlap_transcripts: filter_overlap_transcripts,
                    },
                )
            })?;
            let mut cache_name = annotation.as_os_str().to_owned();
            cache_name.push(".ruspladder.hdf5");
            let cache_path = PathBuf::from(cache_name);
            cache::write_genes(&cache_path, &result.genes)?;
            for (suffix, names) in result.excluded {
                let mut report = annotation.as_os_str().to_owned();
                report.push(format!(".genes_excluded_{suffix}"));
                std::fs::write(&report, names.join("\n") + "\n")
                    .context("write annotation exclusion report")?;
            }
            if verbose {
                eprintln!(
                    "Prepared {} genes in {}",
                    result.genes.len(),
                    cache_path.display()
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
