use anyhow::Result;
use ruspladder::count_io;
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    paths: Vec<PathBuf>,
    output: PathBuf,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        count_io::collect(&r.paths, &r.output)?;
        println!("{{\"ok\":true}}");
    }
    Ok(())
}
