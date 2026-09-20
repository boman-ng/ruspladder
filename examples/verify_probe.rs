use anyhow::Result;
use ruspladder::{
    annotation::Gene,
    count::Counts,
    events::Event,
    verify::{self, VerifyOptions},
};
use serde::Deserialize;
use std::io::{self, BufRead};

#[derive(Deserialize)]
struct Request {
    gene: Gene,
    event: Event,
    counts: Counts,
    options: VerifyOptions,
    min_reads: f64,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let result = verify::verify_event(&r.event, &r.gene, &r.counts, r.options)?;
        let psi = verify::psi(&result.info, r.event.event_type, r.min_reads);
        println!(
            "{}",
            serde_json::json!({"verified": result.verified, "info": result.info, "psi": psi})
        );
    }
    Ok(())
}
