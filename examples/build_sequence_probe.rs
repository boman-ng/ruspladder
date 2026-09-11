use anyhow::Result;
use ruspladder::{
    annotation::Gene,
    build_graph::{self, GraphOptions},
};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    genes: Vec<Gene>,
    samples: Vec<Vec<PathBuf>>,
    options: GraphOptions,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut r: Request = serde_json::from_str(&line?)?;
        let mut result = Vec::new();
        for bams in r.samples {
            let (genes, inserted) =
                build_graph::generate(r.genes.clone(), &bams, &mut r.options, None)?;
            result.push(serde_json::json!({"genes":genes,"inserted":inserted,"read_filter":r.options.reads.filter}));
        }
        println!("{}", serde_json::to_string(&result)?);
    }
    Ok(())
}
