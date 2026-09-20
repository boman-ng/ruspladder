//! Test-only event collection interface, including public coordinate properties.
use anyhow::Result;
use ruspladder::{annotation::Gene, events};
use serde::Deserialize;
use std::{
    collections::BTreeMap,
    io::{self, BufRead},
};

#[derive(Deserialize)]
struct Request {
    genes: Vec<Gene>,
    chromosomes: BTreeMap<String, usize>,
    edge_limit: usize,
    curate: bool,
    parallel: usize,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let request: Request = serde_json::from_str(&line?)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(request.parallel)
            .build()?;
        let found = pool.install(|| {
            events::collect(
                &request.genes,
                &request.chromosomes,
                request.edge_limit,
                request.curate,
            )
        });
        let result: BTreeMap<_, Vec<_>> = found.into_iter().map(|(kind, events)| (kind.as_str(), events.into_iter().map(|event| {
            serde_json::json!({"coords": event.coords(), "inner_coords": event.inner_coords(), "introns": event.introns(),
                "span": event.span(), "coordinate_strings": event.coordinate_strings(), "event": event})
        }).collect())).collect();
        println!("{}", serde_json::to_string(&result)?);
    }
    Ok(())
}
