use anyhow::Result;
use ruspladder::{reads::ReadOptions, sparse};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};

#[derive(Deserialize)]
struct Request {
    bam: PathBuf,
    output: PathBuf,
    chromosomes: Vec<String>,
    reference: Option<PathBuf>,
    options: ReadOptions,
    parallel: usize,
    window: usize,
    unstranded: bool,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        sparse::write_summary(
            &r.bam,
            &r.output,
            &r.chromosomes,
            r.reference.as_deref(),
            &r.options,
            sparse::SummaryOptions {
                parallel: r.parallel,
                window: r.window,
                unstranded: r.unstranded,
            },
        )?;
        println!("{{\"result\":\"ok\"}}");
    }
    Ok(())
}
