// Adapted from SplAdder v3.1.1 alt_splice/detect.py (BSD-3-Clause).
// See licenses/SplAdder-BSD.txt. Enumeration and tie order follow the reference.
use crate::graph::SpliceGraph;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlternativeSite {
    pub common: usize,
    pub alternatives: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MultipleSkip {
    pub first: usize,
    pub skipped: Vec<usize>,
    pub last: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Detections {
    pub exon_skip: Vec<[usize; 3]>,
    pub intron_retention: Vec<[usize; 3]>,
    pub alt_5prime: Vec<AlternativeSite>,
    pub alt_3prime: Vec<AlternativeSite>,
    pub mutex_exons: Vec<[usize; 4]>,
    pub mult_exon_skip: Vec<MultipleSkip>,
}

pub fn detect(graph: &SpliceGraph, strand: char, edge_limit: usize) -> Detections {
    let mut result = Detections {
        mult_exon_skip: multiple_skips(graph, edge_limit),
        ..Default::default()
    };
    // Despite the option name, five upstream detectors limit vertex count.
    // The multiple-skip detector instead limits candidate donor/acceptor pairs.
    if graph.vertices.len() > edge_limit {
        return result;
    }
    for a in 0..graph.vertices.len() {
        for b in graph.successors(a) {
            for c in graph.successors(b) {
                if graph.connected(a, c) {
                    result.exon_skip.push([a, b, c]);
                }
            }
        }
    }
    let mut introns = BTreeSet::new();
    for a in 0..graph.vertices.len() {
        for b in graph.successors(a) {
            let intron = [graph.vertices[a][1], graph.vertices[b][0]];
            if introns.contains(&intron) {
                continue;
            }
            let retained = graph
                .vertices
                .iter()
                .enumerate()
                .filter(|(_, e)| e[0] < intron[0] && e[1] > intron[1])
                .min_by_key(|(_, e)| e[1] - e[0]);
            if let Some((retained, _)) = retained {
                result.intron_retention.push([a, b, retained]);
                introns.insert(intron);
            }
        }
    }
    for a in 0..graph.vertices.len() {
        let successors: Vec<_> = graph.successors(a).collect();
        for (offset, &b) in successors.iter().enumerate() {
            for &c in &successors[offset + 1..] {
                if graph.connected(b, c) || graph.vertices[c][0] < graph.vertices[b][1] {
                    continue;
                }
                for d in graph.successors(c) {
                    if graph.connected(b, d) {
                        result.mutex_exons.push([a, b, c, d]);
                    }
                }
            }
        }
    }
    for common in 0..graph.vertices.len().saturating_sub(2) {
        let neighbors: Vec<_> = graph.successors(common).collect();
        let alternatives = alternative_sites(graph, &neighbors, 0);
        if alternatives.len() >= 2 {
            let group = AlternativeSite {
                common,
                alternatives,
            };
            match strand {
                '+' => result.alt_3prime.push(group),
                '-' => result.alt_5prime.push(group),
                _ => {}
            }
        }
    }
    for common in 2..graph.vertices.len() {
        let neighbors: Vec<_> = graph.predecessors(common).collect();
        let alternatives = alternative_sites(graph, &neighbors, 1);
        if alternatives.len() >= 2 {
            let group = AlternativeSite {
                common,
                alternatives,
            };
            match strand {
                '+' => result.alt_5prime.push(group),
                '-' => result.alt_3prime.push(group),
                _ => {}
            }
        }
    }
    result
}

fn alternative_sites(graph: &SpliceGraph, neighbors: &[usize], side: usize) -> Vec<usize> {
    const MIN_OVERLAP: i64 = 11;
    let mut sites = BTreeSet::new();
    let mut indices = Vec::new();
    for (i, &a) in neighbors.iter().enumerate() {
        for &b in &neighbors[i + 1..] {
            let x = graph.vertices[a];
            let y = graph.vertices[b];
            if x[side] == y[side] || x[1].min(y[1]) - x[0].max(y[0]) < MIN_OVERLAP {
                continue;
            }
            for vertex in [a, b] {
                if sites.insert(graph.vertices[vertex][side]) {
                    indices.push(vertex);
                }
            }
        }
    }
    indices
}

pub fn multiple_skips(graph: &SpliceGraph, edge_limit: usize) -> Vec<MultipleSkip> {
    let n = graph.vertices.len();
    let mut pairs = Vec::new();
    // The reference tests A^3, A^4, ... for a path with >=3 edges beside a
    // direct edge. A topological longest-path scan on the same DAG expresses
    // exactly that reachability test without dense matrix powers/count overflow.
    for first in 0..n {
        let mut longest: Vec<Option<usize>> = vec![None; n];
        longest[first] = Some(0);
        for at in first..n {
            if let Some(length) = longest[at] {
                for next in graph.successors(at) {
                    longest[next] = Some(longest[next].unwrap_or(0).max(length + 1));
                }
            }
        }
        for last in graph.successors(first) {
            if longest[last].is_some_and(|length| length >= 3) {
                pairs.push((first, last));
            }
        }
    }
    if pairs.len() > edge_limit {
        return Vec::new();
    }
    let mut result = Vec::new();
    for (first, last) in pairs {
        let (short, short_path) = floyd_paths(graph, (first, last), 1.0);
        let (long, long_path) = floyd_paths(graph, (first, last), -1.0);
        let index = first * n + last;
        let path = if short[index].is_finite() && short[index] > 2.0 {
            &short_path
        } else if long[index].is_finite() && -long[index] > 2.0 {
            &long_path
        } else {
            continue;
        };
        let mut backtrace = vec![path[index]];
        while *backtrace.last().unwrap() > first {
            backtrace.push(path[first * n + backtrace.last().unwrap()]);
        }
        backtrace.pop();
        backtrace.reverse();
        result.push(MultipleSkip {
            first,
            skipped: backtrace,
            last,
        });
    }
    result
}

/// Floyd-Warshall distances and predecessor selection, including strict `<`
/// updates and ascending intermediate indices, follow the upstream fastcore.
fn floyd_paths(
    graph: &SpliceGraph,
    removed: (usize, usize),
    weight: f64,
) -> (Vec<f64>, Vec<usize>) {
    let n = graph.vertices.len();
    let mut distance = vec![f64::INFINITY; n * n];
    let mut predecessor = vec![0; n * n];
    for i in 0..n {
        for j in graph.successors(i) {
            if (i, j) != removed {
                distance[i * n + j] = weight;
                predecessor[i * n + j] = i;
            }
        }
        distance[i * n + i] = 0.0;
        predecessor[i * n + i] = i;
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let candidate = distance[i * n + k] + distance[k * n + j];
                if candidate < distance[i * n + j] {
                    distance[i * n + j] = candidate;
                    predecessor[i * n + j] = predecessor[k * n + j];
                }
            }
        }
    }
    (distance, predecessor)
}
