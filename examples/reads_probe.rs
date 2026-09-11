//! Test-only interface for comparing collapsed read evidence with SplAdder.
use anyhow::Result;
use ruspladder::reads::{AlignmentReader, ReadOptions};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};

#[derive(Deserialize)]
struct Request {
    bam: PathBuf,
    reference: Option<PathBuf>,
    chromosome: String,
    start: i64,
    stop: i64,
    #[serde(default)]
    options: ReadOptions,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let request: Request = serde_json::from_str(&line?)?;
        let mut reader = AlignmentReader::open(&request.bam, request.reference.as_deref())?;
        let result = reader.region(
            &request.chromosome,
            request.start,
            request.stop,
            &request.options,
        )?;
        println!("{}", serde_json::to_string(&result)?);
    }
    Ok(())
}
