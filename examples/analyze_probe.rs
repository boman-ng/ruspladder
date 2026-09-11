use anyhow::Result;
use ruspladder::{
    analyze::{self, AnalysisOptions},
    annotation::Gene,
    events::{Event, EventType},
    verify::VerifyOptions,
};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    genes: Vec<Gene>,
    events: Vec<Event>,
    kind: EventType,
    counts: PathBuf,
    output: PathBuf,
    samples: Vec<String>,
    sample_idx: Vec<usize>,
    options: VerifyOptions,
    min_reads: f64,
    parallel: usize,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let samples: Vec<_> = r.samples.iter().map(String::as_str).collect();
        let options = AnalysisOptions {
            samples: &samples,
            sample_idx: &r.sample_idx,
            verify: r.options,
            psi_min_reads: r.min_reads,
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(r.parallel)
            .build()?;
        let confirmed = pool.install(|| {
            analyze::analyze(&r.genes, &r.events, r.kind, &r.counts, &r.output, &options)
        })?;
        println!("{}", serde_json::json!({"confirmed": confirmed}));
    }
    Ok(())
}
