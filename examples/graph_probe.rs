//! Test-only JSON transport for differential checks against the Python classes.
use anyhow::Result;
use ruspladder::graph::{Interval, SegmentGraph, SpliceGraph};
use serde::Deserialize;
use std::io::{self, BufRead};

#[derive(Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Operation {
    Sort,
    Uniquify,
    UpdateTerminals,
    Subset {
        indices: Vec<usize>,
    },
    Intron {
        left: Vec<usize>,
        right: Vec<usize>,
        keep_end: bool,
        keep_start: bool,
    },
    Cassette {
        exon: Interval,
        before: Vec<usize>,
        after: Vec<usize>,
    },
    Retention {
        left: usize,
        right: usize,
    },
}

#[derive(Deserialize)]
struct Request {
    transcripts: Vec<Vec<Interval>>,
    #[serde(default)]
    operations: Vec<Operation>,
}

fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let request: Request = serde_json::from_str(&line?)?;
        let mut graph = SpliceGraph::from_transcripts(&request.transcripts);
        for op in request.operations {
            match op {
                Operation::Sort => graph.sort(),
                Operation::Uniquify => graph.uniquify(),
                Operation::UpdateTerminals => graph.update_terminals(),
                Operation::Subset { indices } => graph.reorder(&indices),
                Operation::Intron {
                    left,
                    right,
                    keep_end,
                    keep_start,
                } => graph.add_intron(&left, keep_end, &right, keep_start),
                Operation::Cassette {
                    exon,
                    before,
                    after,
                } => graph.add_cassette_exon(exon, &before, &after),
                Operation::Retention { left, right } => graph.add_intron_retention(left, right),
            }
        }
        let segments = SegmentGraph::from_splice_graph(&graph);
        println!(
            "{}",
            serde_json::json!({"graph": graph, "segments": segments,
            "non_alt": segments.non_alternative_segments()})
        );
    }
    Ok(())
}
