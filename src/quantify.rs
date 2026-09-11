// Adapted from SplAdder v3.1.1 alt_splice.quantify_from_counted_events.
// BSD-3-Clause. spladder_test always sets use_exon_counts=False and gen_event_ids=False.
use crate::{events::EventType, hdf5io};
use anyhow::{Result, ensure};
use hdf5::File;
use serde::Serialize;
use std::{collections::BTreeSet, path::Path};

#[derive(Debug, Serialize)]
pub struct Quantified {
    pub coverage: [Vec<Vec<f64>>; 2],
    pub psi: Vec<Vec<f64>>,
    pub gene_idx: Vec<usize>,
    pub event_idx: Vec<usize>,
    pub samples: Vec<String>,
}

pub fn from_counted_events(
    path: &Path,
    group_a: &[usize],
    group_b: &[usize],
    kind: EventType,
) -> Result<Quantified> {
    let file = File::open(path)?;
    let confirmed = file.dataset("conf_idx")?.read_raw::<u64>()?;
    let event_idx = file
        .dataset(if file.link_exists("filter_idx") {
            "filter_idx"
        } else {
            "conf_idx"
        })?
        .read_raw::<u64>()?
        .into_iter()
        .map(|i| i as usize)
        .collect();
    let names = hdf5io::read_strings(&file.dataset("event_features")?)?;
    let labels = hdf5io::read_strings(&file.dataset("samples")?)?;
    let counts = file.dataset("event_counts")?;
    let psi = file.dataset("psi")?;
    let genes = file.dataset("gene_idx")?.read_raw::<u64>()?;
    ensure!(
        !group_a.is_empty() || !group_b.is_empty(),
        "quantification requires samples"
    );
    for group in [group_a, group_b] {
        ensure!(
            group
                .iter()
                .all(|&i| i < labels.len() && i < counts.shape()[0])
                && group.iter().copied().collect::<BTreeSet<_>>().len() == group.len(),
            "invalid or duplicate sample selection"
        );
    }
    let (first, second): (&[&str], &[&str]) = match kind {
        EventType::ExonSkip => (&["e1e3_conf"], &["e1e2_conf", "e2e3_conf"]),
        EventType::MultExonSkip => (&["e1e3_conf"], &["e1e2_conf", "e2e3_conf", "sum_e2_conf"]),
        EventType::IntronRetention => (&["e1e3_conf"], &[]),
        EventType::MutexExons => (&["e1e2_conf", "e2e4_conf"], &["e1e3_conf", "e3e4_conf"]),
        _ => (&["e1e3_conf"], &["e2_conf"]),
    };
    let fields = [first, second].map(|fields| {
        fields
            .iter()
            .map(|&name| {
                names
                    .iter()
                    .position(|n| n == name)
                    .ok_or_else(|| anyhow::anyhow!("missing event feature {name}"))
            })
            .collect::<Result<Vec<_>>>()
    });
    let [first, second] = fields;
    let fields = [first?, second?];
    let selected: Vec<_> = group_a.iter().chain(group_b).copied().collect();
    let mut sorted_a = group_a.to_vec();
    let mut sorted_b = group_b.to_vec();
    sorted_a.sort_unstable();
    sorted_b.sort_unstable();
    let psi_order: Vec<_> = sorted_a.into_iter().chain(sorted_b).collect();
    let mut result = Quantified {
        coverage: [Vec::new(), Vec::new()],
        psi: Vec::new(),
        gene_idx: Vec::new(),
        event_idx,
        samples: selected.iter().map(|&i| labels[i].clone()).collect(),
    };
    if result.samples[0].ends_with("npz") {
        for sample in &mut result.samples {
            // The upstream regex has an unescaped dot, consuming any separator.
            if sample.len() >= 4 && sample[sample.len() - 3..].eq_ignore_ascii_case("npz") {
                sample.truncate(sample.len() - 4);
            }
        }
    }
    for event in confirmed {
        let event = event as usize;
        ensure!(
            event < genes.len() && event < counts.shape()[2],
            "confirmed event index out of range"
        );
        let evidence = counts.read_slice_2d::<f64, _>((.., .., event))?;
        for (isoform, fields) in fields.iter().enumerate() {
            result.coverage[isoform].push(
                selected
                    .iter()
                    .map(|&sample| {
                        if fields.is_empty() {
                            0.0
                        } else if fields.iter().any(|&f| evidence[[sample, f]].is_nan()) {
                            f64::NAN
                        } else {
                            fields
                                .iter()
                                .map(|&f| evidence[[sample, f]])
                                .fold(f64::INFINITY, f64::min)
                                .floor()
                        }
                    })
                    .collect(),
            );
        }
        let values = psi.read_slice_1d::<f64, _>((.., event))?;
        // Source restores coverage and labels after sorting each group, but
        // leaves PSI in sorted group order. Retain this observable behavior.
        result
            .psi
            .push(psi_order.iter().map(|&i| values[i]).collect());
        result.gene_idx.push(genes[event] as usize);
    }
    Ok(result)
}
