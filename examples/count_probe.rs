use anyhow::Result;
use ruspladder::{annotation::Gene, count, reads::ReadOptions};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};

#[derive(Deserialize)]
struct Request {
    genes: Vec<Gene>,
    bams: Vec<PathBuf>,
    reference: Option<PathBuf>,
    options: ReadOptions,
    parallel: usize,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut r: Request = serde_json::from_str(&line?)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(r.parallel)
            .build()?;
        let mut counts = Vec::new();
        for bam in &r.bams {
            counts.push(pool.install(|| {
                count::count_sample(&mut r.genes, bam, r.reference.as_deref(), &r.options)
            })?);
        }
        println!(
            "{}",
            serde_json::json!({"genes": r.genes, "counts": counts})
        );
    }
    Ok(())
}
