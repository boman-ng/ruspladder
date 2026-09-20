// Adapted from SplAdder v3.1.1 alt_splice/verify.py and helpers.compute_psi.
// BSD-3-Clause; see licenses/SplAdder-BSD.txt.
use crate::{
    annotation::Gene,
    count::Counts,
    events::{Event, EventType},
    graph::Interval,
    numeric,
};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct SkipOptions {
    pub min_skip_rel_cov: f64,
    pub min_non_skip_count: f64,
    pub min_skip_count: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct RetentionOptions {
    pub min_retention_cov: f64,
    pub min_retention_region: f64,
    pub min_retention_rel_cov: f64,
    pub min_non_retention_count: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct AltPrimeOptions {
    pub min_diff_rel_cov: f64,
    pub min_intron_count: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct MutexOptions {
    pub min_skip_rel_cov: f64,
    pub min_conf_count: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct VerifyOptions {
    pub use_anno_support: bool,
    pub exon_skip: SkipOptions,
    pub mult_exon_skip: SkipOptions,
    pub intron_retention: RetentionOptions,
    pub alt_prime: AltPrimeOptions,
    pub mutex_exons: MutexOptions,
}

#[derive(Debug, Serialize)]
pub struct Verified {
    pub verified: Vec<bool>,
    pub info: Vec<f64>,
}

impl Verified {
    pub fn empty(kind: EventType) -> Self {
        let (flags, features) = match kind {
            EventType::ExonSkip => (4, 7),
            EventType::MultExonSkip => (5, 10),
            EventType::MutexExons => (4, 9),
            _ => (2, 6),
        };
        let mut info = vec![0.0; features];
        info[0] = 1.0;
        Self {
            verified: vec![false; flags],
            info,
        }
    }
}

struct Evidence<'a> {
    gene: &'a Gene,
    counts: &'a Counts,
}

impl Evidence<'_> {
    fn segments(&self, exon: Interval, fallback: bool) -> Result<Vec<usize>> {
        let graph = &self.gene.segmentgraph;
        let mut result = Vec::new();
        for (i, &v) in self.gene.splicegraph.vertices.iter().enumerate() {
            if v == exon {
                result.extend(&graph.matches[i]);
            }
        }
        if result.is_empty() && fallback {
            result.extend(
                graph
                    .segments
                    .iter()
                    .enumerate()
                    .filter_map(|(i, v)| (v[0] >= exon[0] && v[1] <= exon[1]).then_some(i)),
            );
        }
        ensure!(
            !result.is_empty(),
            "event exon {exon:?} has no matching segments"
        );
        if !fallback {
            result.sort_unstable();
        }
        Ok(result)
    }

    fn length(&self, indices: &[usize]) -> f64 {
        indices
            .iter()
            .map(|&i| {
                let v = self.gene.segmentgraph.segments[i];
                v[1] - v[0]
            })
            .sum::<i64>() as f64
    }

    fn mean(&self, indices: &[usize]) -> f64 {
        let values: Vec<_> = indices
            .iter()
            .map(|&i| {
                let v = self.gene.segmentgraph.segments[i];
                self.counts.segments[i] * (v[1] - v[0]) as f64
            })
            .collect();
        numeric::sum(&values) / self.length(indices)
    }

    fn edge(&self, from: &[usize], to: &[usize], required: bool) -> Result<f64> {
        let id = (*from.last().unwrap() * self.gene.segmentgraph.segments.len() + to[0]) as u64;
        let mut values = self.counts.edges.iter().filter(|v| v[0] == id);
        let first = values.next();
        if required {
            ensure!(
                first.is_some() && values.next().is_none(),
                "event requires exactly one count for segment edge {id}"
            );
        }
        Ok(first.map_or(0.0, |v| v[1] as f64))
    }
}

fn difference(a: &[usize], b: &[usize]) -> Vec<usize> {
    a.iter()
        .copied()
        .filter(|x| !b.contains(x))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn intersection(a: &[usize], b: &[usize]) -> Vec<usize> {
    a.iter()
        .copied()
        .filter(|x| b.contains(x))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

pub fn verify_event(
    event: &Event,
    gene: &Gene,
    counts: &Counts,
    options: VerifyOptions,
) -> Result<Verified> {
    let kind = event.event_type;
    let mut result = Verified::empty(kind);
    // verify_wrapper takes this path before inspecting event coordinates.
    if counts.edges.is_empty() {
        return Ok(result);
    }
    let negative = event
        .exons1
        .iter()
        .chain(&event.exons2)
        .flatten()
        .any(|&p| p < 0);
    let alt_prime = matches!(kind, EventType::Alt3prime | EventType::Alt5prime);
    let invalid_length = !alt_prime
        && event
            .exons1
            .iter()
            .chain(&event.exons2)
            .any(|v| v[1] - v[0] < 1);
    let invalid_alt = alt_prime
        && event.exons1[0][1] != event.exons2[0][1]
        && event.exons1[1][0] != event.exons2[1][0];
    let invalid_mutex = kind == EventType::MutexExons && {
        let a = event.exons1[1];
        let b = event.exons2[1];
        (a[1] > b[0] && a[0] < b[0]) || (b[1] > a[0] && b[0] < a[0])
    };
    if negative || invalid_length || invalid_alt || invalid_mutex {
        result.info[0] = 0.0;
        return Ok(result);
    }
    let evidence = Evidence { gene, counts };
    let annotated = |a: Interval, b: Interval| {
        options.use_anno_support && gene.introns_anno.contains(&[a[1], b[0]])
    };
    let exon_annotated = |exon: Interval| {
        options.use_anno_support && gene.exons.iter().flatten().any(|&v| v == exon)
    };
    let info = &mut result.info;
    let verified = &mut result.verified;
    match kind {
        EventType::ExonSkip | EventType::MultExonSkip => {
            let multiple = kind == EventType::MultExonSkip;
            let config = if multiple {
                options.mult_exon_skip
            } else {
                options.exon_skip
            };
            let exons = &event.exons2;
            let pre = evidence.segments(exons[0], false)?;
            let last_exon = *exons.last().unwrap();
            let after = evidence.segments(last_exon, false)?;
            let middle = exons[1..exons.len() - 1]
                .iter()
                .map(|&v| evidence.segments(v, false))
                .collect::<Result<Vec<_>>>()?;
            let all: Vec<_> = if multiple {
                middle
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect()
            } else {
                middle[0].clone()
            };
            info[1] = evidence.mean(&pre);
            info[2] = evidence.mean(&all);
            info[3] = evidence.mean(&after);
            let anno_pre = annotated(exons[0], exons[1]);
            let anno_after = annotated(exons[exons.len() - 2], last_exon);
            verified[0] = info[2] >= config.min_skip_rel_cov * (info[1] + info[3]) / 2.0
                || (anno_pre && anno_after);
            info[4] = evidence.edge(&pre, &middle[0], !multiple)?;
            info[5] = evidence.edge(middle.last().unwrap(), &after, !multiple)?;
            info[6] = evidence.edge(&pre, &after, !multiple)?;
            verified[1] = info[4] >= config.min_non_skip_count || anno_pre;
            verified[2] = info[5] >= config.min_non_skip_count || anno_after;
            if multiple {
                for pair in middle.windows(2) {
                    info[7] += evidence.edge(&pair[0], &pair[1], false)?;
                }
                info[8] = (exons.len() - 2) as f64;
                info[9] = exons[1..exons.len() - 1]
                    .iter()
                    .map(|v| v[1] - v[0])
                    .sum::<i64>() as f64;
                let inner_annotated = options.use_anno_support
                    && exons[1..exons.len() - 1]
                        .windows(2)
                        .all(|p| gene.introns_anno.contains(&[p[0][1], p[1][0]]));
                verified[3] = info[7] / info[8] >= config.min_non_skip_count || inner_annotated;
                verified[4] = info[6] >= config.min_skip_count || annotated(exons[0], last_exon);
            } else {
                verified[3] = info[6] >= config.min_skip_count || annotated(exons[0], last_exon);
            }
        }
        EventType::IntronRetention => {
            let a = event.exons1[0];
            let b = event.exons1[1];
            let left = evidence.segments(a, false)?;
            let right = evidence.segments(b, false)?;
            let all: Vec<_> = (left[0]..*right.last().unwrap()).collect();
            let intron = difference(&difference(&all, &left), &right);
            ensure!(
                !intron.is_empty(),
                "retention event has no intronic segments"
            );
            info[1] = evidence.mean(&left);
            info[2] = evidence.mean(&intron);
            info[3] = evidence.mean(&right);
            info[5] = numeric::sum(
                &intron
                    .iter()
                    .map(|&i| counts.seg_pos[i])
                    .collect::<Vec<_>>(),
            ) / evidence.length(&intron);
            let config = options.intron_retention;
            verified[0] = (info[2] > config.min_retention_cov
                && info[5] > config.min_retention_region
                && info[2] >= config.min_retention_rel_cov * (info[1] + info[3]) / 2.0)
                || exon_annotated([a[0], b[1]]);
            info[4] = evidence.edge(&left, &right, true)?;
            verified[1] = info[4] >= config.min_non_retention_count || annotated(a, b);
        }
        EventType::Alt3prime | EventType::Alt5prime => {
            let a = evidence.segments(event.exons1[0], true)?;
            let b = evidence.segments(event.exons1[1], true)?;
            let c = evidence.segments(event.exons2[0], true)?;
            let d = evidence.segments(event.exons2[1], true)?;
            let (constant1, constant2, diff) = if a == c {
                let mut diff = difference(&b, &d);
                if diff.is_empty() {
                    diff = difference(&d, &b);
                }
                (a.clone(), intersection(&b, &d), diff)
            } else {
                ensure!(
                    b == d,
                    "both exons differ in an alternative splice-site event"
                );
                let mut diff = difference(&a, &c);
                if diff.is_empty() {
                    diff = difference(&c, &a);
                }
                (intersection(&c, &a), b.clone(), diff)
            };
            info[1] = evidence.mean(&constant1);
            info[2] = evidence.mean(&diff);
            info[3] = evidence.mean(&constant2);
            let length1 = evidence.length(&constant1);
            let length2 = evidence.length(&constant2);
            let anno = annotated(event.exons1[0], event.exons1[1])
                && annotated(event.exons2[0], event.exons2[1]);
            verified[0] = info[2]
                >= options.alt_prime.min_diff_rel_cov * (info[1] * length1 + info[3] * length2)
                    / (length1 + length2)
                || anno;
            info[4] = evidence.edge(&a, &b, true)?;
            info[5] = evidence.edge(&c, &d, true)?;
            verified[1] = info[4].min(info[5]) >= options.alt_prime.min_intron_count || anno;
        }
        EventType::MutexExons => {
            let exons = [
                event.exons1[0],
                event.exons1[1],
                event.exons2[1],
                *event.exons1.last().unwrap(),
            ];
            let segments = exons
                .iter()
                .map(|&v| evidence.segments(v, false))
                .collect::<Result<Vec<_>>>()?;
            for i in 0..4 {
                info[i + 1] = evidence.mean(&segments[i]);
            }
            let config = options.mutex_exons;
            verified[0] = info[2] >= config.min_skip_rel_cov * (info[1] + info[4]) / 2.0
                || exon_annotated(exons[1]);
            verified[1] = info[3] >= config.min_skip_rel_cov * (info[1] + info[4]) / 2.0
                || exon_annotated(exons[2]);
            for (feature, a, b) in [(5, 0, 1), (6, 0, 2), (7, 1, 3), (8, 2, 3)] {
                info[feature] = evidence.edge(&segments[a], &segments[b], false)?;
            }
            verified[2] = info[5].min(info[7]) >= config.min_conf_count
                || (annotated(exons[0], exons[1]) && annotated(exons[1], exons[3]));
            verified[3] = info[6].min(info[8]) >= config.min_conf_count
                || (annotated(exons[0], exons[2]) && annotated(exons[2], exons[3]));
        }
    }
    Ok(result)
}

pub fn psi(info: &[f64], kind: EventType, min_reads: f64) -> (f64, f64, f64) {
    let (a, b) = match kind {
        EventType::ExonSkip => (info[4].min(info[5]), info[6]),
        EventType::IntronRetention => (info[2], info[4]),
        EventType::Alt3prime | EventType::Alt5prime => (info[5], info[4]),
        EventType::MutexExons => (info[6] + info[8], info[5] + info[7]),
        EventType::MultExonSkip => (info[4] + info[5] + info[7], (info[8] + 1.0) * info[6]),
    };
    (
        if a + b < min_reads {
            f64::NAN
        } else {
            a / (a + b)
        },
        a,
        b,
    )
}
