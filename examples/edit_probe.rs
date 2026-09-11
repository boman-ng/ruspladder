//! Test-only graph editing and gene labeling interface.
use anyhow::Result;
use ruspladder::{
    annotation::Gene,
    editgraph::{self, RemoveExons},
};
use serde::Deserialize;
use std::io::{self, BufRead};

#[derive(Deserialize)]
struct Request {
    gene: Gene,
    remove: Option<RemoveExons>,
    merge_duplicates: bool,
    label: bool,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let mut request: Request = serde_json::from_str(&line?)?;
        if let Some(options) = request.remove {
            editgraph::remove_short_exons(&mut request.gene, options);
        }
        if request.merge_duplicates {
            editgraph::merge_duplicate_exons(&mut request.gene.splicegraph);
        }
        if request.label {
            request.gene.label_alt();
        }
        println!("{}", serde_json::to_string(&request.gene)?);
    }
    Ok(())
}
