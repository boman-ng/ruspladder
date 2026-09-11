use anyhow::Result;
use ruspladder::{reads::ReadOptions, sparse::SparseReader};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    path: PathBuf,
    chromosome: String,
    start: i64,
    stop: i64,
    options: ReadOptions,
    unstranded: bool,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let v = SparseReader::open(&r.path)?.query(
            &r.chromosome,
            r.start,
            r.stop,
            &r.options,
            true,
            r.unstranded,
        )?;
        println!(
            "{}",
            serde_json::json!({"coverage":v.coverage,"introns_plus":v.introns_plus,"introns_minus":v.introns_minus})
        );
    }
    Ok(())
}
