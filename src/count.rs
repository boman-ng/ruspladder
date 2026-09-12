// Adapted from SplAdder v3.1.1 count.count_graph_coverage (BSD-3-Clause).
use crate::{
    annotation::Gene,
    graph::SegmentGraph,
    reads::{EvidenceReader, ReadOptions},
};
use anyhow::{Result, ensure};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, path::Path};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Counts {
    pub segments: Vec<f64>,
    pub seg_pos: Vec<f64>,
    pub edges: Vec<[u64; 2]>,
    /// NumPy promotes mixed int64 edge IDs / uint64 sums to float64.
    pub edges_float: bool,
}

fn count_gene(
    gene: &mut Gene,
    reader: &mut EvidenceReader,
    options: &ReadOptions,
) -> Result<Counts> {
    if gene.segmentgraph.segments.is_empty() {
        gene.segmentgraph = SegmentGraph::from_splice_graph(&gene.splicegraph);
    }
    let segments = &gene.segmentgraph.segments;
    ensure!(
        !segments.is_empty(),
        "cannot count a graph without segments"
    );
    gene.start = *segments.iter().flatten().min().unwrap();
    gene.stop = *segments.iter().flatten().max().unwrap();
    let mut options = options.clone();
    options.filter = None;
    options.strand = Some(gene.strand);
    let evidence = reader.region(&gene.chr, gene.start, gene.stop, &options)?;
    ensure!(
        evidence.coverage.len() == (gene.stop - gene.start) as usize,
        "counted segments extend beyond sparse coverage"
    );
    let mut result = Counts {
        segments: Vec::with_capacity(segments.len()),
        seg_pos: Vec::with_capacity(segments.len()),
        edges: Vec::new(),
        edges_float: false,
    };
    for &[start, stop] in segments {
        let coverage =
            &evidence.coverage[(start - gene.start) as usize..(stop - gene.start) as usize];
        result
            .segments
            .push(coverage.iter().map(|&x| x as f64).sum::<f64>() / coverage.len() as f64);
        result
            .seg_pos
            .push(coverage.iter().filter(|&&x| x > 0).count() as f64);
    }
    let mut introns = BTreeMap::<[i64; 2], u64>::new();
    for [a, b, n] in evidence
        .introns_plus
        .into_iter()
        .chain(evidence.introns_minus)
    {
        if a > gene.start && b < gene.stop {
            *introns.entry([a, b]).or_default() += n as u64;
        }
    }
    for &[a, b] in &gene.segmentgraph.edges {
        let matched = introns.get(&[segments[a][1], segments[b][0]]);
        result.edges_float |= matched.is_some();
        result.edges.push([
            (a * segments.len() + b) as u64,
            matched.copied().unwrap_or(0),
        ]);
    }
    result.edges_float |= result.edges.is_empty();
    Ok(result)
}

pub fn count_sample(
    genes: &mut [Gene],
    bam: &Path,
    reference: Option<&Path>,
    options: &ReadOptions,
) -> Result<Vec<Counts>> {
    let mut options = options.clone();
    options.filter = None;
    genes
        .par_iter_mut()
        .with_max_len(64)
        .map_init(
            || EvidenceReader::open(bam, reference, &options),
            |reader, gene| {
                let reader = reader
                    .as_mut()
                    .map_err(|error| anyhow::anyhow!("{error:#}"))?;
                count_gene(gene, reader, &options)
            },
        )
        .collect()
}
