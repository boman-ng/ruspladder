//! Reference-only adapter: import already decoded Python fixtures into native caches.
use anyhow::{Result, ensure};
use ruspladder::{annotation::Gene, cache, events::Event};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};
#[derive(Deserialize)]
struct Request {
    genes: Vec<Gene>,
    graph: PathBuf,
    events: Vec<EventFile>,
}
#[derive(Deserialize)]
struct EventFile {
    path: PathBuf,
    values: Vec<Event>,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        cache::write_genes(&r.graph, &r.genes)?;
        ensure!(
            cache::read_genes(&r.graph)? == r.genes,
            "graph cache roundtrip mismatch"
        );
        for event in r.events {
            cache::write_events(&event.path, &event.values)?;
            ensure!(
                cache::read_events(&event.path)? == event.values,
                "event cache roundtrip mismatch"
            );
            let before = std::fs::read(&event.path)?;
            let failed: Result<()> = cache::atomic_write(&event.path, |temporary| {
                std::fs::write(temporary, b"incomplete")?;
                anyhow::bail!("injected cache write failure")
            });
            ensure!(
                failed.is_err() && std::fs::read(&event.path)? == before,
                "failed cache write replaced completed data"
            );
        }
        println!("{{\"ok\":true}}");
    }
    Ok(())
}
