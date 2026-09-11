// Adapted from SplAdder v3.1.1 classes/event.py and alt_splice/{collect,events}.py.
// BSD-3-Clause; see licenses/SplAdder-BSD.txt.
use crate::{
    annotation::Gene,
    detect::{self, AlternativeSite},
    graph::Interval,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EventType {
    ExonSkip,
    IntronRetention,
    #[serde(rename = "alt_3prime")]
    Alt3prime,
    #[serde(rename = "alt_5prime")]
    Alt5prime,
    MultExonSkip,
    MutexExons,
}

impl EventType {
    pub const ALL: [Self; 6] = [
        Self::ExonSkip,
        Self::IntronRetention,
        Self::Alt3prime,
        Self::Alt5prime,
        Self::MultExonSkip,
        Self::MutexExons,
    ];
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExonSkip => "exon_skip",
            Self::IntronRetention => "intron_retention",
            Self::Alt3prime => "alt_3prime",
            Self::Alt5prime => "alt_5prime",
            Self::MultExonSkip => "mult_exon_skip",
            Self::MutexExons => "mutex_exons",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Event {
    pub event_type: EventType,
    pub chr: String,
    pub strand: char,
    pub exons1: Vec<Interval>,
    pub exons2: Vec<Interval>,
    pub gene_name: Vec<String>,
    pub gene_idx: usize,
    pub id: usize,
    pub annotated: u8,
}

impl Event {
    fn new(
        kind: EventType,
        gene: &Gene,
        gene_idx: usize,
        mut exons1: Vec<Interval>,
        mut exons2: Vec<Interval>,
    ) -> Self {
        exons1.sort();
        exons2.sort();
        let mut annotated = 3;
        if exons1
            .windows(2)
            .any(|p| !gene.introns_anno.contains(&[p[0][1], p[1][0]]))
        {
            annotated -= 1;
        }
        if exons2
            .windows(2)
            .any(|p| !gene.introns_anno.contains(&[p[0][1], p[1][0]]))
        {
            annotated -= 2;
        }
        Self {
            event_type: kind,
            chr: gene.chr.clone(),
            strand: gene.strand,
            exons1,
            exons2,
            gene_name: vec![gene.name.clone()],
            gene_idx,
            id: 0,
            annotated,
        }
    }

    fn shorter_first(&mut self) {
        if self.exons1.iter().map(|e| e[1] - e[0]).sum::<i64>()
            > self.exons2.iter().map(|e| e[1] - e[0]).sum::<i64>()
        {
            // Annotation is assigned BEFORE the swap in upstream; don't
            // silently reinterpret the bit field when preserving compatibility.
            std::mem::swap(&mut self.exons1, &mut self.exons2);
        }
    }

    pub fn span(&self) -> i64 {
        let min = self
            .exons1
            .iter()
            .chain(&self.exons2)
            .flatten()
            .min()
            .unwrap();
        let max = self
            .exons1
            .iter()
            .chain(&self.exons2)
            .flatten()
            .max()
            .unwrap();
        max - min
    }

    pub fn coords(&self) -> Vec<i64> {
        let mut coords: Vec<_> = if self.event_type == EventType::MultExonSkip {
            self.exons1
                .iter()
                .flatten()
                .take(4)
                .chain(self.exons2.iter().rev().take(2).flatten())
                .copied()
                .collect()
        } else {
            self.exons1
                .iter()
                .chain(&self.exons2)
                .flatten()
                .copied()
                .collect()
        };
        coords.sort_unstable();
        coords
    }

    pub fn inner_coords(&self) -> Vec<i64> {
        let mut one: Vec<_> = self.exons1.iter().flatten().copied().collect();
        let mut two: Vec<_> = self.exons2.iter().flatten().copied().collect();
        let mut coords = if self.event_type == EventType::MultExonSkip {
            two.sort_unstable();
            [
                two[1..4].to_vec(),
                two[two.len() - 4..two.len() - 1].to_vec(),
            ]
            .concat()
        } else if self.event_type == EventType::MutexExons {
            [
                one[1..4].to_vec(),
                self.exons2[1].to_vec(),
                vec![self.exons1[2][0]],
            ]
            .concat()
        } else {
            one.sort_unstable();
            two.sort_unstable();
            [
                one[1..one.len() - 1].to_vec(),
                two[1..two.len() - 1].to_vec(),
            ]
            .concat()
        };
        coords.sort_unstable();
        if self.event_type != EventType::MutexExons {
            coords.dedup();
        }
        coords
    }

    pub fn introns(&self) -> Vec<Interval> {
        self.exons1
            .windows(2)
            .chain(self.exons2.windows(2))
            .map(|p| [p[0][1], p[1][0]])
            .collect()
    }

    pub fn coordinate_strings(&self) -> (String, String) {
        let one: BTreeSet<_> = self
            .exons1
            .iter()
            .map(|e| format!("{}-{}", e[0], e[1]))
            .collect();
        let two: BTreeSet<_> = self
            .exons2
            .iter()
            .map(|e| format!("{}-{}", e[0], e[1]))
            .collect();
        let both: Vec<_> = one.union(&two).cloned().collect();
        let starts: Vec<i64> = both
            .iter()
            .map(|s| s.split('-').next().unwrap().parse().unwrap())
            .collect();
        let order = crate::sort::argsort_i64(&starts);
        let labels: Vec<_> = order.iter().map(|&i| both[i].as_str()).collect();
        let usage: Vec<_> = labels
            .iter()
            .map(|s| {
                if one.contains(*s) && two.contains(*s) {
                    "0"
                } else {
                    "1"
                }
            })
            .collect();
        (labels.join(":"), usage.join(":"))
    }
}

pub fn collect_gene(
    gene: &Gene,
    gene_idx: usize,
    edge_limit: usize,
) -> BTreeMap<EventType, Vec<Event>> {
    let mut result: BTreeMap<_, _> = EventType::ALL
        .into_iter()
        .map(|kind| (kind, Vec::new()))
        .collect();
    if gene.is_alt != Some(true) {
        return result;
    }
    let found = detect::detect(&gene.splicegraph, gene.strand, edge_limit);
    let vertices = &gene.splicegraph.vertices;
    let mut add = |kind, one, two, shorter| {
        let mut event = Event::new(kind, gene, gene_idx, one, two);
        if shorter {
            event.shorter_first();
        }
        result.get_mut(&kind).unwrap().push(event);
    };
    for [a, b, _] in found.intron_retention {
        add(
            EventType::IntronRetention,
            vec![vertices[a], vertices[b]],
            vec![[vertices[a][0], vertices[b][1]]],
            false,
        );
    }
    for [a, b, c] in found.exon_skip {
        add(
            EventType::ExonSkip,
            vec![vertices[a], vertices[c]],
            vec![vertices[a], vertices[b], vertices[c]],
            false,
        );
    }
    for skip in found.mult_exon_skip {
        let mut indices = vec![skip.first];
        indices.extend(skip.skipped);
        indices.push(skip.last);
        add(
            EventType::MultExonSkip,
            vec![vertices[skip.first], vertices[skip.last]],
            indices.iter().map(|&i| vertices[i]).collect(),
            false,
        );
    }
    for [a, b, c, d] in found.mutex_exons {
        add(
            EventType::MutexExons,
            vec![vertices[a], vertices[b], vertices[d]],
            vec![vertices[a], vertices[c], vertices[d]],
            true,
        );
    }
    for (kind, groups) in [
        (EventType::Alt5prime, found.alt_5prime),
        (EventType::Alt3prime, found.alt_3prime),
    ] {
        for AlternativeSite {
            common,
            alternatives,
        } in groups
        {
            for (i, &a) in alternatives.iter().enumerate() {
                for &b in &alternatives[i + 1..] {
                    if vertices[a][0] >= vertices[b][1] || vertices[a][1] <= vertices[b][0] {
                        continue;
                    }
                    add(
                        kind,
                        vec![vertices[a], vertices[common]],
                        vec![vertices[b], vertices[common]],
                        true,
                    );
                }
            }
        }
    }
    result
}

pub fn post_process(mut events: Vec<Event>, chromosomes: &BTreeMap<String, usize>) -> Vec<Event> {
    events.retain(|event| event.coords().iter().all(|&position| position > 0));
    for event in &mut events {
        event.exons1.sort();
        event.exons2.sort();
    }
    events.retain(|event| event.introns().iter().all(|intron| intron[1] > intron[0]));
    events
        .sort_by_cached_key(|event| (chromosomes[&event.chr], event.strand == '-', event.coords()));
    events.sort_by_cached_key(|event| {
        (
            chromosomes[&event.chr],
            event.strand == '-',
            event.inner_coords(),
        )
    });
    let mut unique: Vec<Event> = Vec::new();
    for event in events {
        if let Some(last) = unique.last_mut()
            && last.chr == event.chr
            && last.strand == event.strand
            && last.inner_coords() == event.inner_coords()
        {
            if event.span() <= last.span() {
                *last = event;
            }
        } else {
            unique.push(event);
        }
    }
    for (i, event) in unique.iter_mut().enumerate() {
        event.id = i + 1;
    }
    unique
}

pub fn curate_alt_prime(events: &mut Vec<Event>) {
    events.retain_mut(|event| {
        if event.exons1[1][0] - event.exons1[0][1] < 1
            || event.exons2[1][0] - event.exons2[0][1] < 1
        {
            return false;
        }
        let common_left = event.exons1[0] == event.exons2[0];
        let common_right = event.exons1[1] == event.exons2[1];
        let disjoint = |a: Interval, b: Interval| a[1] <= b[0] || a[0] >= b[1];
        if (common_left && disjoint(event.exons1[1], event.exons2[1]))
            || (common_right && disjoint(event.exons1[0], event.exons2[0]))
        {
            return true;
        }
        if common_left {
            let end = event.exons1[1][1].min(event.exons2[1][1]);
            event.exons1[1][1] = end;
            event.exons2[1][1] = end;
            event.shorter_first();
        } else if common_right {
            let start = event.exons1[0][0].max(event.exons2[0][0]);
            event.exons1[0][0] = start;
            event.exons2[0][0] = start;
            event.shorter_first();
        }
        true
    });
}

pub fn collect(
    genes: &[Gene],
    chromosomes: &BTreeMap<String, usize>,
    edge_limit: usize,
    curate: bool,
) -> BTreeMap<EventType, Vec<Event>> {
    use rayon::prelude::*;
    let per_gene: Vec<_> = genes
        .par_iter()
        .enumerate()
        .map(|(i, gene)| collect_gene(gene, i, edge_limit))
        .collect();
    let mut result: BTreeMap<_, Vec<_>> = EventType::ALL
        .into_iter()
        .map(|kind| (kind, Vec::new()))
        .collect();
    for found in per_gene {
        for (kind, events) in found {
            result.get_mut(&kind).unwrap().extend(events);
        }
    }
    for (&kind, events) in &mut result {
        *events = post_process(std::mem::take(events), chromosomes);
        if curate && matches!(kind, EventType::Alt3prime | EventType::Alt5prime) {
            curate_alt_prime(events);
        }
    }
    result
}
