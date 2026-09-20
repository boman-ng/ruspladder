use anyhow::{Result, ensure};
use ruspladder::{annotation::Gene, cache, merge};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};

#[derive(Deserialize)]
struct Request {
    samples: Vec<Vec<Gene>>,
    chunksize: Option<usize>,
    min_count: Option<u64>,
    cache: PathBuf,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let samples = if let Some(size) = r.chunksize {
            let mut source = r.samples.into_iter();
            let mut merged = Vec::new();
            loop {
                let chunk: Vec<_> = source.by_ref().take(size).map(Ok).collect();
                if chunk.is_empty() {
                    break;
                }
                merged.push(merge::merge_samples(chunk)?);
            }
            merged
        } else {
            r.samples
        };
        let mut genes = merge::merge_samples(samples.into_iter().map(Ok))?;
        if let Some(count) = r.min_count {
            merge::filter_edge_support(&mut genes, count)?;
        }
        cache::write_genes(&r.cache, &genes)?;
        ensure!(
            cache::read_genes(&r.cache)? == genes,
            "merged graph cache roundtrip mismatch"
        );
        println!("{}", serde_json::to_string(&genes)?);
    }
    Ok(())
}
