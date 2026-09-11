use anyhow::Result;
use ruspladder::{
    quantify::Quantified,
    testing::{self, InputOptions},
};
use serde::Deserialize;
use std::io::{self, BufRead};
#[derive(Deserialize)]
struct Request {
    counts: Vec<Vec<f64>>,
    samples: usize,
    cap_expression: bool,
    events: EventInput,
    group_a: usize,
    options: InputOptions,
}
#[derive(Deserialize)]
struct EventInput {
    coverage: [Vec<Vec<f64>>; 2],
    psi: Vec<Vec<Option<f64>>>,
    gene_idx: Vec<usize>,
    event_idx: Vec<usize>,
    samples: Vec<String>,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let expression = testing::prepare_expression(r.counts, r.samples, r.cap_expression)?;
        let events = Quantified {
            coverage: r.events.coverage,
            psi: r
                .events
                .psi
                .into_iter()
                .map(|row| row.into_iter().map(|x| x.unwrap_or(f64::NAN)).collect())
                .collect(),
            gene_idx: r.events.gene_idx,
            event_idx: r.events.event_idx,
            samples: r.events.samples,
        };
        let prepared = testing::prepare_events(events, &expression, r.group_a, r.options)?;
        let means = prepared.as_ref().map(|p| {
            testing::means_and_fold_changes(&p.counts, &p.size_factors, r.group_a)
                .into_iter()
                .map(|row| row.map(|v| v.to_string()))
                .collect::<Vec<_>>()
        });
        println!(
            "{}",
            serde_json::json!({"expression":expression,"prepared":prepared,"means":means})
        );
    }
    Ok(())
}
