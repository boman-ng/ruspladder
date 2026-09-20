// Adapted from SplAdder v3.1.1 editgraph.py and merge.py (BSD-3-Clause).
// See licenses/SplAdder-BSD.txt.
use crate::{annotation::Gene, graph::SpliceGraph};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct RemoveExons {
    pub terminal_short_len: i64,
    pub terminal_short_extend: i64,
    pub min_exon_len: i64,
    pub min_exon_len_remove: i64,
}

/// Preserve upstream's j-2 predecessor search and first/last vertex extension.
pub fn remove_short_exons(gene: &mut Gene, options: RemoveExons) {
    let graph = &mut gene.splicegraph;
    if graph.vertices.is_empty() {
        return;
    }
    if graph.vertices[0][1] - graph.vertices[0][0] < options.terminal_short_len {
        graph.vertices[0][0] -= options.terminal_short_extend;
        gene.start = gene.start.min(graph.vertices[0][0]);
    }
    let last = graph.vertices.len() - 1;
    if graph.vertices[last][1] - graph.vertices[last][0] < options.terminal_short_len {
        graph.vertices[last][1] += options.terminal_short_extend;
        gene.stop = gene.stop.max(graph.vertices[last][1]);
    }
    let mut remove = vec![false; graph.vertices.len()];
    for (j, removed) in remove.iter_mut().enumerate().skip(1) {
        if graph.vertices[j][1] - graph.vertices[j][0] >= options.min_exon_len {
            continue;
        }
        let next = (j + 1..graph.vertices.len()).find(|&k| {
            graph.vertices[k][1] - graph.vertices[k][0] >= options.min_exon_len_remove
                && graph.connected(j, k)
        });
        let previous = (0..j.saturating_sub(1)).rev().find(|&k| {
            graph.vertices[k][1] - graph.vertices[k][0] >= options.min_exon_len_remove
                && graph.connected(k, j)
        });
        if let (Some(a), Some(b)) = (previous, next) {
            graph.connect(a, b);
            *removed = graph.vertices[j][1] - graph.vertices[j][0] < options.min_exon_len_remove;
        }
    }
    graph.reorder(
        &(0..remove.len())
            .filter(|&i| !remove[i])
            .collect::<Vec<_>>(),
    );
}

/// Unlike SpliceGraph::uniquify, this routine ORs duplicate terminal flags and
/// leaves already-unique graphs in their original order.
pub fn merge_duplicate_exons(graph: &mut SpliceGraph) {
    use std::collections::BTreeSet;
    if graph.vertices.iter().collect::<BTreeSet<_>>().len() == graph.vertices.len() {
        return;
    }
    graph.sort();
    let mut merged = SpliceGraph::default();
    let mut mapping = Vec::with_capacity(graph.vertices.len());
    for (i, &vertex) in graph.vertices.iter().enumerate() {
        let j = if merged.vertices.last() == Some(&vertex) {
            let j = merged.vertices.len() - 1;
            for t in 0..2 {
                merged.terminals[j][t] |= graph.terminals[i][t];
            }
            j
        } else {
            merged.push_vertex(vertex, graph.terminals[i])
        };
        mapping.push(j);
    }
    for (i, neighbors) in graph.edges.iter().enumerate() {
        for &j in neighbors {
            merged.connect(mapping[i], mapping[j]);
        }
    }
    *graph = merged;
}
