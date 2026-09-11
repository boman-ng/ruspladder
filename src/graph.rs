// Algorithms adapted from SplAdder v3.1.1 classes/splicegraph.py,
// classes/segmentgraph.py, and classes/gene.py (BSD-3-Clause).
// Copyright (c) 2008-2012 Cheng Soon Ong, Gunnar Raetsch, Andre Kahles;
// 2012-2019 Gunnar Raetsch, Andre Kahles. See licenses/SplAdder-BSD.txt.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// Genomic coordinates use SplAdder's zero-based, half-open convention.
pub type Interval = [i64; 2];

/// Symmetric splice graph. Neighbor lists are sorted, unique vertex indices.
/// Terminal flags are stored independently: upstream mutations do not always
/// recompute them, and that distinction is observable by later algorithms.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SpliceGraph {
    pub vertices: Vec<Interval>,
    pub edges: Vec<Vec<usize>>,
    pub terminals: Vec<[bool; 2]>,
}

impl SpliceGraph {
    pub fn from_transcripts(transcripts: &[Vec<Interval>]) -> Self {
        let mut graph = Self::default();
        for exons in transcripts {
            if exons.len() == 1 {
                // Upstream deliberately appends even a duplicate single exon.
                graph.push_vertex(exons[0], [false; 2]);
            } else {
                for pair in exons.windows(2) {
                    // Look up BOTH endpoints before adding either one, matching
                    // from_gene even for repeated coordinates within a transcript.
                    let a = graph.vertices.iter().rposition(|e| *e == pair[0]);
                    let b = graph.vertices.iter().rposition(|e| *e == pair[1]);
                    let a = a.unwrap_or_else(|| graph.push_vertex(pair[0], [false; 2]));
                    let b = b.unwrap_or_else(|| graph.push_vertex(pair[1], [false; 2]));
                    graph.connect(a, b);
                }
            }
        }
        let order =
            crate::sort::argsort_i64(&graph.vertices.iter().map(|e| e[0]).collect::<Vec<_>>());
        graph.reorder(&order);
        graph.update_terminals();
        graph
    }

    pub fn push_vertex(&mut self, exon: Interval, terminals: [bool; 2]) -> usize {
        let index = self.vertices.len();
        self.vertices.push(exon);
        self.edges.push(Vec::new());
        self.terminals.push(terminals);
        index
    }

    pub fn connected(&self, a: usize, b: usize) -> bool {
        self.edges[a].binary_search(&b).is_ok()
    }

    pub fn connect(&mut self, a: usize, b: usize) {
        for (from, to) in [(a, b), (b, a)] {
            if let Err(at) = self.edges[from].binary_search(&to) {
                self.edges[from].insert(at, to);
            }
        }
    }

    pub fn disconnect(&mut self, a: usize, b: usize) {
        for (from, to) in [(a, b), (b, a)] {
            if let Ok(at) = self.edges[from].binary_search(&to) {
                self.edges[from].remove(at);
            }
        }
    }

    pub fn predecessors(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        self.edges[i].iter().copied().take_while(move |&j| j < i)
    }

    pub fn successors(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        self.edges[i].iter().copied().filter(move |&j| j > i)
    }

    pub fn update_terminals(&mut self) {
        for (i, neighbors) in self.edges.iter().enumerate() {
            self.terminals[i] = [
                !neighbors.iter().any(|&j| j <= i),
                !neighbors.iter().any(|&j| j >= i),
            ];
        }
    }

    /// Also implements upstream subset(): omitted vertices and edges disappear.
    pub fn reorder(&mut self, order: &[usize]) {
        let mut mapping = vec![None; self.vertices.len()];
        for (new, &old) in order.iter().enumerate() {
            mapping[old] = Some(new);
        }
        let vertices = order.iter().map(|&i| self.vertices[i]).collect();
        let terminals = order.iter().map(|&i| self.terminals[i]).collect();
        let edges = order
            .iter()
            .map(|&i| {
                let mut neighbors: Vec<_> =
                    self.edges[i].iter().filter_map(|&j| mapping[j]).collect();
                neighbors.sort_unstable();
                neighbors
            })
            .collect();
        *self = Self {
            vertices,
            edges,
            terminals,
        };
    }

    pub fn sort(&mut self) {
        let mut order: Vec<_> = (0..self.vertices.len()).collect();
        order.sort_by_key(|&i| self.vertices[i]);
        self.reorder(&order);
    }

    pub fn uniquify(&mut self) {
        self.sort();
        let mut result = Self::default();
        let mut mapping = Vec::with_capacity(self.vertices.len());
        for (i, &exon) in self.vertices.iter().enumerate() {
            let index = if result.vertices.last() == Some(&exon) {
                let index = result.vertices.len() - 1;
                // Upstream retains the last duplicate's terminal flags.
                result.terminals[index] = self.terminals[i];
                index
            } else {
                result.push_vertex(exon, self.terminals[i])
            };
            mapping.push(index);
        }
        for (i, neighbors) in self.edges.iter().enumerate() {
            for &j in neighbors {
                result.connect(mapping[i], mapping[j]);
            }
        }
        *self = result;
    }

    pub fn add_intron(
        &mut self,
        left: &[usize],
        keep_end: bool,
        right: &[usize],
        keep_start: bool,
    ) {
        if !right.is_empty() {
            let mut duplicate = Vec::new();
            if keep_end {
                duplicate.extend(
                    left.iter()
                        .copied()
                        .filter(|&i| !self.edges[i].iter().any(|&j| j >= i)),
                );
            }
            if keep_start {
                duplicate.extend(
                    right
                        .iter()
                        .copied()
                        .filter(|&i| !self.edges[i].iter().any(|&j| j <= i)),
                );
            }
            for old in duplicate {
                let neighbors = self.edges[old].clone();
                let new = self.push_vertex(self.vertices[old], self.terminals[old]);
                for j in neighbors {
                    self.connect(new, j);
                }
            }
        }
        for &a in left {
            for &b in right {
                self.connect(a, b);
            }
        }
        self.uniquify();
    }

    pub fn add_cassette_exon(&mut self, exon: Interval, before: &[usize], after: &[usize]) {
        let new = self.push_vertex(exon, [false; 2]);
        for &i in before.iter().chain(after) {
            self.connect(i, new);
        }
    }

    pub fn add_intron_retention(&mut self, left: usize, right: usize) {
        let neighbors: Vec<_> = self.edges[left]
            .iter()
            .copied()
            .filter(|&i| i <= left)
            .chain(self.edges[right].iter().copied().filter(|&i| i >= right))
            .collect();
        let new = self.push_vertex(
            [self.vertices[left][0], self.vertices[right][1]],
            [self.terminals[left][0], self.terminals[right][1]],
        );
        for i in neighbors {
            self.connect(i, new);
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SegmentGraph {
    pub segments: Vec<Interval>,
    /// For each splice vertex, indices of the segments it contains.
    pub matches: Vec<Vec<usize>>,
    /// Directed donor -> acceptor edges in genomic coordinate order.
    pub edges: Vec<[usize; 2]>,
}

impl SegmentGraph {
    pub fn from_splice_graph(graph: &SpliceGraph) -> Self {
        let mut delta = BTreeMap::<i64, i64>::new();
        for &[start, end] in &graph.vertices {
            *delta.entry(start).or_default() += 1;
            *delta.entry(end).or_default() -= 1;
        }
        let mut segments = Vec::new();
        let mut previous = None;
        let mut active = 0;
        for (position, change) in delta {
            if let Some(start) = previous
                && active > 0
            {
                segments.push([start, position]);
            }
            active += change;
            previous = Some(position);
        }
        let matches: Vec<Vec<usize>> = graph
            .vertices
            .iter()
            .map(|exon| {
                segments
                    .iter()
                    .enumerate()
                    .filter_map(|(j, segment)| {
                        (exon[0] <= segment[0] && exon[1] >= segment[1]).then_some(j)
                    })
                    .collect()
            })
            .collect();
        let mut edges = BTreeSet::new();
        for (i, neighbors) in graph.edges.iter().enumerate() {
            for &j in neighbors.iter().filter(|&&j| j >= i) {
                if let (Some(&donor), Some(&acceptor)) = (matches[i].last(), matches[j].first()) {
                    edges.insert([donor, acceptor]);
                }
            }
        }
        Self {
            segments,
            matches,
            edges: edges.into_iter().collect(),
        }
    }

    pub fn non_alternative_segments(&self) -> Vec<usize> {
        let mut keep = vec![true; self.segments.len()];
        for &[donor, acceptor] in &self.edges {
            if acceptor > donor {
                keep[donor + 1..acceptor].fill(false);
            }
        }
        keep.into_iter()
            .enumerate()
            .filter_map(|(i, keep)| keep.then_some(i))
            .collect()
    }
}
