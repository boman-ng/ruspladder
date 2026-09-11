//! Test-only annotation reader and full HDF5 roundtrip comparison.
use anyhow::{Result, ensure};
use ruspladder::{
    annotation::{AnnotationFilters, read_annotation},
    cache,
};
use std::path::Path;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    ensure!(
        args.len() == 4,
        "usage: annotation_probe ANNOTATION CACHE FILTER_BITS"
    );
    let flags: u8 = args[3].parse()?;
    let annotation = read_annotation(
        Path::new(&args[1]),
        AnnotationFilters {
            overlap_genes: flags & 1 != 0,
            overlap_exons: flags & 2 != 0,
            overlap_transcripts: flags & 4 != 0,
        },
    )?;
    cache::write_genes(Path::new(&args[2]), &annotation.genes)?;
    let genes = cache::read_genes(Path::new(&args[2]))?;
    ensure!(
        genes == annotation.genes,
        "HDF5 graph cache roundtrip mismatch"
    );
    println!(
        "{}",
        serde_json::json!({"genes": genes, "excluded": annotation.excluded})
    );
    Ok(())
}
