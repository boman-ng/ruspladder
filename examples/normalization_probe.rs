use anyhow::Result;
use ruspladder::expression::{self, Normalization};
use serde::Deserialize;
use std::io::{self, BufRead};
#[derive(Deserialize)]
struct Request {
    counts: Vec<Vec<f64>>,
    samples: usize,
    kind: Normalization,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let value = match expression::size_factors(&r.counts, r.samples, r.kind) {
            Ok(values) => serde_json::json!({"values":values}),
            Err(error) => serde_json::json!({"error":error.to_string()}),
        };
        println!("{value}");
    }
    Ok(())
}
