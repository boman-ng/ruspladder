use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "ruspladder",
    version,
    about = "Rust implementation of SplAdder's non-visual alternative splicing analysis"
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
    /// Prepare annotation and sparse alignment summaries.
    Prep(Box<ruspladder::prep::PrepArgs>),
}

fn run() -> Result<()> {
    match Cli::parse().command {
        Command::Build(options) => ruspladder::build_cli::run(&options)?,
        Command::Test(options) => ruspladder::test_cli::run(&options)?,
        Command::Prep(options) => ruspladder::prep::run(&options)?,
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
