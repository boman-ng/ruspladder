use anyhow::Result;
use ruspladder::{
    events::Event,
    statistics::TestingResult,
    test_output::{self, Metadata, Options},
    testing::Prepared,
};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    output: PathBuf,
    options: Options,
    metadata: Metadata,
    events: Vec<Event>,
    prepared: Prepared,
    result: TestingResult,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        test_output::write_results(
            &r.output,
            &r.options,
            &r.metadata,
            &r.events,
            &r.prepared,
            &r.result,
        )?;
        println!("{{\"ok\":true}}");
    }
    Ok(())
}
