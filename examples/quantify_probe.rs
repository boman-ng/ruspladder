use anyhow::Result;
use ruspladder::{events::EventType, quantify};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    path: PathBuf,
    group_a: Vec<usize>,
    group_b: Vec<usize>,
    kind: EventType,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        println!(
            "{}",
            serde_json::to_string(&quantify::from_counted_events(
                &r.path, &r.group_a, &r.group_b, r.kind
            )?)?
        );
    }
    Ok(())
}
