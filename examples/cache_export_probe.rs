use anyhow::Result;
use ruspladder::cache;
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    graph: Option<PathBuf>,
    events: Option<PathBuf>,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let value = if let Some(path) = r.graph {
            serde_json::to_value(cache::read_genes(&path)?)?
        } else {
            serde_json::to_value(cache::read_events(&r.events.unwrap())?)?
        };
        println!("{}", value);
    }
    Ok(())
}
