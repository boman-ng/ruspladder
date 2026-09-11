// Adapted from SplAdder v3.1.1 merge.py and editgraph.filter_by_edgecount.
// BSD-3-Clause; see licenses/SplAdder-BSD.txt.
use crate::{
    annotation::Gene,
    graph::{SegmentGraph, SpliceGraph},
};
use anyhow::{Result, ensure};
use std::collections::{BTreeMap, BTreeSet};

fn initialize_counts(gene: &mut Gene) {
    gene.edge_count = Some(
        gene.splicegraph
            .edges
            .iter()
            .map(|row| row.iter().map(|&j| (j, 1)).collect())
            .collect(),
    );
}

fn merge_one(base: &mut Gene, mut incoming: Gene) -> Result<()> {
    if incoming.edge_count.is_none() {
        initialize_counts(&mut incoming);
    }
    let old_vertices = &base.splicegraph.vertices;
    let ignore_new = old_vertices.len() > 10000;
    let vertices: Vec<_> = if ignore_new {
        old_vertices.clone()
    } else {
        old_vertices
            .iter()
            .chain(&incoming.splicegraph.vertices)
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    };
    let index: BTreeMap<_, _> = vertices.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let mut counts = vec![BTreeMap::<usize, u64>::new(); vertices.len()];
    let mut graph = SpliceGraph {
        terminals: vec![[false; 2]; vertices.len()],
        edges: vec![Vec::new(); vertices.len()],
        vertices,
    };
    for gene in [&*base, &incoming] {
        let source_counts = gene.edge_count.as_ref().unwrap();
        ensure!(
            source_counts.len() == gene.splicegraph.vertices.len(),
            "edge-count shape disagrees with graph during merge"
        );
        let mapping: Vec<_> = gene
            .splicegraph
            .vertices
            .iter()
            .map(|v| index.get(v).copied())
            .collect();
        for (i, row) in source_counts.iter().enumerate() {
            let Some(a) = mapping[i] else {
                continue;
            };
            for &(j, value) in row {
                if let Some(b) = mapping[j] {
                    *counts[a].entry(b).or_default() += value;
                }
            }
        }
        if !ignore_new {
            for (i, row) in gene.splicegraph.edges.iter().enumerate() {
                for &j in row {
                    graph.connect(mapping[i].unwrap(), mapping[j].unwrap());
                }
            }
        }
    }
    base.edge_count = Some(
        counts
            .into_iter()
            .map(|row| row.into_iter().filter(|&(_, n)| n > 0).collect())
            .collect(),
    );
    if !ignore_new {
        graph.update_terminals();
        base.splicegraph = graph;
    }
    Ok(())
}

/// Merge sample graphs from the same annotation, in the supplied sample order.
/// The first sample owns gene metadata/bounds; support counts sum across samples.
pub fn merge_samples(samples: impl IntoIterator<Item = Result<Vec<Gene>>>) -> Result<Vec<Gene>> {
    let mut result: Vec<Gene> = Vec::new();
    for (sample_index, sample) in samples.into_iter().enumerate() {
        let mut sample = sample?;
        ensure!(!sample.is_empty(), "cannot merge an empty sample graph");
        sample.sort_by(|a, b| a.name.cmp(&b.name));
        for gene in &mut sample {
            gene.splicegraph.uniquify();
        }
        if sample_index == 0 {
            if sample[0].edge_count.is_none() {
                for gene in &mut sample {
                    initialize_counts(gene);
                }
            }
            result = sample;
            continue;
        }
        ensure!(
            result.len() == sample.len()
                && result.iter().zip(&sample).all(|(a, b)| a.name == b.name),
            "sample graphs have different gene sets; SplAdder v3.1.1's append-gene path is undefined for these inputs"
        );
        for (base, incoming) in result.iter_mut().zip(sample) {
            merge_one(base, incoming)?;
        }
    }
    ensure!(!result.is_empty(), "no graphs to merge");
    for gene in &mut result {
        gene.label_alt();
        if !gene.segmentgraph.segments.is_empty() {
            gene.segmentgraph = SegmentGraph::from_splice_graph(&gene.splicegraph);
        }
    }
    let order = crate::sort::argsort_i64(&result.iter().map(|g| g.start).collect::<Vec<_>>());
    let mut slots: Vec<_> = result.into_iter().map(Some).collect();
    result = order
        .into_iter()
        .map(|i| slots[i].take().unwrap())
        .collect();
    result.sort_by(|a, b| a.chr.cmp(&b.chr));
    Ok(result)
}

pub fn filter_edge_support(genes: &mut Vec<Gene>, min_count: u64) -> Result<()> {
    let mut retained = Vec::new();
    for mut gene in std::mem::take(genes) {
        let counts = gene
            .edge_count
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("graph validation requires edge support counts"))?;
        let original_isolated: BTreeSet<_> = gene
            .splicegraph
            .edges
            .iter()
            .enumerate()
            .filter_map(|(i, row)| row.is_empty().then_some(i))
            .collect();
        let mut edges = vec![Vec::new(); counts.len()];
        for (i, row) in counts.iter().enumerate() {
            for &(j, count) in row {
                if count >= min_count {
                    edges[i].push(j);
                }
            }
            if min_count == 0 {
                edges[i] = (0..counts.len()).collect();
            }
        }
        let keep: Vec<_> = (0..edges.len())
            .filter(|&i| original_isolated.contains(&i) || !edges[i].is_empty())
            .collect();
        gene.splicegraph.edges = edges;
        if !keep.is_empty() {
            gene.splicegraph.reorder(&keep);
            // Source leaves the original edge_count matrix and terminals.
            retained.push(gene);
        }
    }
    *genes = retained;
    Ok(())
}
