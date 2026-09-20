// Adapted from SplAdder v3.1.1 editgraph.py (BSD-3-Clause).
// See licenses/SplAdder-BSD.txt.
use crate::{annotation::Gene, graph::Interval};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct RetentionOptions {
    pub min_retention_cov: f64,
    pub min_retention_region: f64,
    pub min_retention_max_exon_fold_diff: f64,
    pub min_retention_rel_cov: f64,
    pub max_retention_rel_cov: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct CassetteOptions {
    pub min_cassette_region: f64,
    pub min_cassette_cov: f64,
    pub min_cassette_rel_diff: f64,
}

fn mean(values: &[u64]) -> f64 {
    values.iter().map(|&v| v as f64).sum::<f64>() / values.len() as f64
}

fn median(values: &[u64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    (sorted[(sorted.len() - 1) / 2] as f64 + sorted[sorted.len() / 2] as f64) / 2.0
}

fn coverage(track: &[u64], start: i64, interval: Interval) -> Result<&[u64]> {
    let [a, b] = interval.map(|p| p - start);
    ensure!(
        a >= 0 && b >= 0 && a as usize <= track.len() && b as usize <= track.len(),
        "coverage interval {interval:?} is outside gene track"
    );
    Ok(if a < b {
        &track[a as usize..b as usize]
    } else {
        &[]
    })
}

pub fn insert_retentions(
    gene: &mut Gene,
    track: &[u64],
    options: RetentionOptions,
) -> Result<usize> {
    let graph = &mut gene.splicegraph;
    let exon_coverage: Vec<_> = graph
        .vertices
        .iter()
        .map(|&v| coverage(track, gene.start, v).map(median))
        .collect::<Result<_>>()?;
    let mut pending = BTreeSet::new();
    for (k, neighbors) in graph.edges.iter().enumerate() {
        for &l in neighbors.iter().filter(|&&l| l >= k) {
            if graph
                .vertices
                .iter()
                .any(|v| v[0] > graph.vertices[k][1] && v[1] < graph.vertices[l][0])
            {
                continue;
            }
            let cov = coverage(
                track,
                gene.start,
                [graph.vertices[k][1], graph.vertices[l][0]],
            )?;
            let average = mean(cov);
            let fraction =
                cov.iter().filter(|&&v| v as f64 > 0.5 * average).count() as f64 / cov.len() as f64;
            let exon_average = (exon_coverage[k] + exon_coverage[l]) / 2.0;
            if median(cov) > options.min_retention_cov
                && fraction > options.min_retention_region
                && exon_coverage[k].max(exon_coverage[l])
                    / (1e-6 + exon_coverage[k].min(exon_coverage[l]))
                    <= options.min_retention_max_exon_fold_diff
                && average >= options.min_retention_rel_cov * exon_average
                && average <= options.max_retention_rel_cov * exon_average
            {
                pending.insert([k, l]);
            }
        }
    }
    let inserted = pending.len();
    if inserted == 0 {
        return Ok(0);
    }
    // For a nonnegative DAG adjacency matrix, exp(A)[i,j] > 0 exactly
    // when j is reachable from i. The source only inspects this predicate.
    // Build the same candidate paths without computing an unused exponential.
    for k in (0..graph.vertices.len()).rev() {
        let successors: Vec<_> = pending.iter().filter(|p| p[0] == k).map(|p| p[1]).collect();
        for j in successors {
            let ends: Vec<_> = pending.iter().filter(|p| p[0] == j).map(|p| p[1]).collect();
            pending.extend(ends.into_iter().map(|l| [k, l]));
        }
    }
    loop {
        let next = pending.iter().find(|p| p[0] < p[1]).copied();
        if let Some([k, l]) = next {
            pending.remove(&[k, l]);
            graph.add_intron_retention(k, l);
        }
        let order =
            crate::sort::argsort_i64(&graph.vertices.iter().map(|v| v[0]).collect::<Vec<_>>());
        let mut reverse = vec![0; order.len()];
        for (i, &old) in order.iter().enumerate() {
            reverse[old] = i;
        }
        pending = pending
            .into_iter()
            .map(|[a, b]| [reverse[a], reverse[b]])
            .collect();
        graph.reorder(&order);
        if next.is_none() {
            break;
        }
    }
    Ok(inserted)
}

pub fn insert_cassettes(
    gene: &mut Gene,
    introns: &[Interval],
    track: &[u64],
    options: CassetteOptions,
) -> Result<usize> {
    let graph = &mut gene.splicegraph;
    let mut all: BTreeSet<_> = introns.iter().copied().collect();
    for (k, neighbors) in graph.edges.iter().enumerate() {
        for &l in neighbors.iter().filter(|&&l| l > k) {
            all.insert([graph.vertices[k][1], graph.vertices[l][0]]);
        }
    }
    let all: Vec<_> = all
        .into_iter()
        .filter(|p| p[1] > gene.start && p[0] < gene.stop)
        .collect();
    let starts: BTreeSet<_> = all.iter().map(|p| p[0]).collect();
    let ends: BTreeSet<_> = all.iter().map(|p| p[1]).collect();
    let mut pending = Vec::new();
    for (k, left) in all.iter().enumerate() {
        for right in &all[k + 1..] {
            if left[1] >= right[0]
                || !graph.vertices.iter().any(|v| v[1] == left[0])
                || !graph.vertices.iter().any(|v| v[0] == right[1])
            {
                continue;
            }
            let exon = [left[1], right[0]];
            if graph
                .vertices
                .iter()
                .any(|v| v[0] < exon[1] && v[1] > exon[0])
            {
                continue;
            }
            let cov = coverage(track, gene.start, exon)?;
            let pre_start = ends
                .range(..exon[0])
                .next_back()
                .copied()
                .unwrap_or(gene.start);
            let post_stop = starts
                .range((
                    std::ops::Bound::Excluded(exon[1]),
                    std::ops::Bound::Unbounded,
                ))
                .next()
                .copied()
                .unwrap_or(gene.stop);
            let before = coverage(track, gene.start, [pre_start, exon[0]])?;
            let after = coverage(track, gene.start, [exon[1], post_stop])?;
            let n_before = before.len().min(cov.len());
            let n_after = after.len().min(cov.len());
            let average = mean(cov);
            let fraction =
                cov.iter().filter(|&&v| v as f64 > 0.2 * average).count() as f64 / cov.len() as f64;
            if fraction > options.min_cassette_region && median(cov) > options.min_cassette_cov {
                // Python [-0:] denotes the entire slice; NaN medians remain
                // NaN through max(nan, 1), hence comparisons fail as upstream.
                let suffix = |values: &[u64], count: usize| {
                    if count == 0 { 0 } else { values.len() - count }
                };
                let med_after = median(&after[..n_after]);
                let med_before = median(&before[suffix(before, n_before)..]);
                let floor = |x: f64| if x.is_nan() { x } else { x.max(1.0) };
                if median(&cov[suffix(cov, n_after)..]) / floor(med_after) - 1.0
                    >= options.min_cassette_rel_diff
                    && median(&cov[..n_before]) / floor(med_before) - 1.0
                        >= options.min_cassette_rel_diff
                {
                    pending.push((exon, left[0], right[1]));
                }
            }
        }
    }
    let original = graph.vertices.len();
    for &(exon, donor, acceptor) in &pending {
        let before: Vec<_> = (0..original)
            .filter(|&i| graph.vertices[i][1] == donor)
            .collect();
        let after: Vec<_> = (0..original)
            .filter(|&i| graph.vertices[i][0] == acceptor)
            .collect();
        graph.add_cassette_exon(exon, &before, &after);
    }
    if !pending.is_empty() {
        // Both source reorder calls are observable for equal starts.
        for _ in 0..2 {
            let order =
                crate::sort::argsort_i64(&graph.vertices.iter().map(|v| v[0]).collect::<Vec<_>>());
            graph.reorder(&order);
        }
    }
    Ok(pending.len())
}
