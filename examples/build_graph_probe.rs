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
    bams: Vec<PathBuf>,
    options: GraphOptions,
    reference: Option<PathBuf>,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut r: Request = serde_json::from_str(&line?)?;
        let (genes, inserted) =
            build_graph::generate(r.genes, &r.bams, &mut r.options, r.reference.as_deref())?;
        println!(
            "{}",
            serde_json::json!({"genes": genes, "inserted": inserted, "read_filter": r.options.reads.filter})
        );
    }
    Ok(())
}
