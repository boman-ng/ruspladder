use anyhow::Result;
use ruspladder::{
    events::Event,
    output::{self, Format},
};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    events: Vec<Event>,
    format: Format,
    counts: PathBuf,
    output: PathBuf,
    samples: Vec<String>,
    indices: Vec<usize>,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        output::write(
            &r.output,
            r.format,
            &r.events,
            &r.counts,
            &r.samples.iter().map(String::as_str).collect::<Vec<_>>(),
            &r.indices,
        )?;
        println!("{{\"ok\":true}}");
    }
    Ok(())
}
