//! Reference-test interface for new intron edges and coverage query side effects.
use anyhow::Result;
use ruspladder::{
    annotation::Gene,
    graph::Interval,
    intron_edges::{self, IntronOptions},
};
use serde::Deserialize;
use std::io::{self, BufRead};

#[derive(Deserialize)]
struct Request {
    genes: Vec<Gene>,
    introns: Vec<Vec<Interval>>,
    options: IntronOptions,
    coverage_depth: u64,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut r: Request = serde_json::from_str(&line?)?;
        let inserted = intron_edges::insert_edges(&mut r.genes, &r.introns, r.options, |gene| {
            Ok(vec![r.coverage_depth; (gene.stop - gene.start) as usize])
        })?;
        println!(
            "{}",
            serde_json::json!({"genes": r.genes, "inserted": inserted})
        );
    }
    Ok(())
}
