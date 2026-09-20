use anyhow::Result;
use ruspladder::{
    annotation::Gene,
    introns::{self, IntronLists},
    reference::Reference,
};
use serde::Deserialize;
use std::{
    io::{self, BufRead},
    path::PathBuf,
};

#[derive(Deserialize)]
struct Request {
    genes: Vec<Gene>,
    introns: IntronLists,
    offset: Option<i64>,
    reference: Option<PathBuf>,
    lenient: bool,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut r: Request = serde_json::from_str(&line?)?;
        if let Some(path) = r.reference {
            introns::filter_consensus(
                &mut r.introns,
                &r.genes,
                &Reference::open(&path)?,
                r.lenient,
            )?;
        }
        if let Some(offset) = r.offset {
            introns::filter_ambiguous(&mut r.introns, &r.genes, offset)?;
        }
        println!("{}", serde_json::to_string(&r.introns)?);
    }
    Ok(())
}
