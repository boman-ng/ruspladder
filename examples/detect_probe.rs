//! Test-only raw detector interface, before coordinate collection or filtering.
use anyhow::Result;
use ruspladder::{detect, graph::SpliceGraph};
use serde::Deserialize;
use std::io::{self, BufRead};

#[derive(Deserialize)]
struct Request {
    graph: SpliceGraph,
    strand: char,
    edge_limit: usize,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let request: Request = serde_json::from_str(&line?)?;
        println!(
            "{}",
            serde_json::to_string(&detect::detect(
                &request.graph,
                request.strand,
                request.edge_limit
            ))?
        );
    }
    Ok(())
}
