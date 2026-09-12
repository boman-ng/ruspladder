//! Isolate the existing graph-cache read/write costs on a real graph.
use anyhow::{Result, ensure};
use ruspladder::cache;
use std::{path::Path, time::Instant};

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    ensure!(args.len() == 3, "usage: cache_bench INPUT OUTPUT");
    let start = Instant::now();
    let genes = cache::read_genes(Path::new(&args[1]))?;
    let read = start.elapsed().as_secs_f64();
    let start = Instant::now();
    cache::write_genes(Path::new(&args[2]), &genes)?;
    let write = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let roundtrip = cache::read_genes(Path::new(&args[2]))?;
    ensure!(genes == roundtrip, "graph cache roundtrip mismatch");
    println!(
        "{}",
        serde_json::json!({"genes":genes.len(), "read_seconds":read,
        "write_seconds":write,"reread_seconds":start.elapsed().as_secs_f64()})
    );
    Ok(())
}
