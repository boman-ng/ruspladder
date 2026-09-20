use anyhow::Result;
use ruspladder::{annotation::Gene, count, count_io::CountWriter, expression, reads::ReadOptions};
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
    samples: Vec<String>,
    counts: PathBuf,
    expression: Option<PathBuf>,
    sample_idx: Option<Vec<usize>>,
    readlen: f64,
    parallel: usize,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut r: Request = serde_json::from_str(&line?)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(r.parallel)
            .build()?;
        let options = ReadOptions::default();
        let first = pool.install(|| {
            count::count_sample(&mut r.genes, &r.bams[0], r.reference.as_deref(), &options)
        })?;
        let mut writer = CountWriter::create(
            &r.counts,
            &r.genes,
            &r.samples.iter().map(String::as_str).collect::<Vec<_>>(),
            r.bams.len(),
        )?;
        writer.append(&first)?;
        drop(first);
        for bam in &r.bams[1..] {
            writer.append(&pool.install(|| {
                count::count_sample(&mut r.genes, bam, r.reference.as_deref(), &options)
            })?)?;
        }
        writer.finish()?;
        if let Some(path) = &r.expression {
            expression::compute(
                &r.genes,
                &r.counts,
                r.readlen,
                r.sample_idx.as_deref(),
                Some(path),
            )?;
        }
        println!("{{\"ok\":true}}");
    }
    Ok(())
}
