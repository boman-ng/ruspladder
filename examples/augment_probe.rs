//! Differential-test interface for coverage-dependent graph augmentation.
use anyhow::Result;
use ruspladder::{
    annotation::Gene,
    augment::{self, CassetteOptions, RetentionOptions},
    graph::Interval,
};
use serde::Deserialize;
use std::io::{self, BufRead};

#[derive(Deserialize)]
struct Request {
    gene: Gene,
    track: Vec<u64>,
    introns: Vec<Interval>,
    retention: Option<RetentionOptions>,
    cassette: Option<CassetteOptions>,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut request: Request = serde_json::from_str(&line?)?;
        let mut inserted = 0;
        if let Some(options) = request.cassette {
            inserted += augment::insert_cassettes(
                &mut request.gene,
                &request.introns,
                &request.track,
                options,
            )?;
        }
        if let Some(options) = request.retention {
            inserted += augment::insert_retentions(&mut request.gene, &request.track, options)?;
        }
        println!(
            "{}",
            serde_json::json!({"gene": request.gene, "inserted": inserted})
        );
    }
    Ok(())
}
