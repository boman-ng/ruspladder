// Adapted from SplAdder v3.1.1 alt_splice/{analyze,verify}.py (BSD-3-Clause).
use crate::{
    annotation::Gene,
    count_io::CountReader,
    events::{Event, EventType},
    hdf5io,
    verify::{self, Verified, VerifyOptions},
};
use anyhow::{Result, ensure};
use hdf5::File;
use ndarray::{ArrayView2, ArrayView3};
use rayon::prelude::*;
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

pub fn features(kind: EventType) -> &'static [&'static str] {
    match kind {
        EventType::ExonSkip => &[
            "valid",
            "e1_cov",
            "e2_cov",
            "e3_cov",
            "e1e2_conf",
            "e2e3_conf",
            "e1e3_conf",
        ],
        EventType::MultExonSkip => &[
            "valid",
            "e1_cov",
            "e2_cov",
            "e3_cov",
            "e1e2_conf",
            "e2e3_conf",
            "e1e3_conf",
            "sum_e2_conf",
            "num_e2",
            "len_e2",
        ],
        EventType::IntronRetention => &[
            "valid",
            "e1_cov",
            "e2_cov",
            "e3_cov",
            "e1e3_conf",
            "e2_cov_region",
        ],
        EventType::MutexExons => &[
            "valid",
            "e1_cov",
            "e2_cov",
            "e3_cov",
            "e4_cov",
            "e1e2_conf",
            "e1e3_conf",
            "e2e4_conf",
            "e3e4_conf",
        ],
        _ => &[
            "valid",
            "e1_cov",
            "e2_cov",
            "e3_cov",
            "e1e3_conf",
            "e2_conf",
        ],
    }
}

fn metadata(file: &File, genes: &[Gene], samples: &[&str], kind: EventType) -> Result<()> {
    hdf5io::strings(file, "samples", &[samples.len()], samples, false)?;
    file.link_soft("/samples", "strains")?;
    hdf5io::strings(
        file,
        "event_features",
        &[features(kind).len()],
        features(kind),
        false,
    )?;
    for (name, values) in [
        (
            "gene_names",
            genes.iter().map(|g| g.name.as_str()).collect::<Vec<_>>(),
        ),
        ("gene_chr", genes.iter().map(|g| g.chr.as_str()).collect()),
        (
            "gene_strand",
            genes
                .iter()
                .map(|g| if g.strand == '+' { "+" } else { "-" })
                .collect(),
        ),
    ] {
        hdf5io::strings(file, name, &[genes.len()], &values, false)?;
    }
    hdf5io::array(
        file,
        "gene_pos",
        &[genes.len(), 2],
        &genes
            .iter()
            .flat_map(|g| [g.start, g.stop])
            .collect::<Vec<_>>(),
        false,
    )?;
    Ok(())
}

fn event_positions(event: &Event) -> Vec<i64> {
    let a = &event.exons1;
    let b = &event.exons2;
    match event.event_type {
        EventType::ExonSkip => b.iter().flatten().copied().collect(),
        EventType::IntronRetention => a.iter().flatten().copied().collect(),
        EventType::MultExonSkip => [b[0], b[1], b[b.len() - 2], b[b.len() - 1]]
            .into_iter()
            .flatten()
            .collect(),
        EventType::MutexExons => [
            a[0][0], a[1][0], b[1][0], b[2][0], a[0][1], a[1][1], b[1][1], b[2][1],
        ]
        .to_vec(),
        _ => a
            .iter()
            .zip(b)
            .map(|(a, b)| [a[0], a[1], b[0], b[1]])
            .collect::<BTreeSet<_>>()
            .into_iter()
            .flatten()
            .collect(),
    }
}

pub struct AnalysisOptions<'a> {
    pub samples: &'a [&'a str],
    pub sample_idx: &'a [usize],
    pub verify: VerifyOptions,
    pub psi_min_reads: f64,
}

/// Verify bounded event batches in parallel and write them in public event order.
pub fn analyze(
    genes: &[Gene],
    events: &[Event],
    kind: EventType,
    counts_path: &Path,
    output: &Path,
    options: &AnalysisOptions<'_>,
) -> Result<Vec<usize>> {
    let file = File::create(output)?;
    metadata(&file, genes, options.samples, kind)?;
    if events.is_empty() {
        hdf5io::array(&file, "event_counts", &[1], &[0i64], false)?;
        file.close()?;
        return Ok(Vec::new());
    }
    ensure!(
        events
            .iter()
            .all(|e| e.event_type == kind && e.gene_idx < genes.len()),
        "event type or gene index mismatch"
    );
    let ns = options.sample_idx.len();
    let nf = features(kind).len();
    let nv = Verified::empty(kind).verified.len();
    let ne = events.len();
    let counts_out = file
        .new_dataset::<f64>()
        .shape((ns, nf, ne))
        .deflate(4)
        .chunk((1, nf, ne.min(128)))
        .create("event_counts")?;
    let verified_out = file
        .new_dataset::<bool>()
        .shape((ns, nv, ne))
        .deflate(4)
        .chunk((1, nv, ne.min(128)))
        .create("verified")?;
    let psi_out = file
        .new_dataset::<f64>()
        .shape((ns, ne))
        .deflate(4)
        .chunk((1, ne.min(128)))
        .create("psi")?;
    let iso1_out = file
        .new_dataset::<i32>()
        .shape((ns, ne))
        .deflate(4)
        .chunk((1, ne.min(128)))
        .create("iso1")?;
    let iso2_out = file
        .new_dataset::<i32>()
        .shape((ns, ne))
        .deflate(4)
        .chunk((1, ne.min(128)))
        .create("iso2")?;
    let reader = CountReader::open(counts_path, genes.len())?;
    let mut num_verified = vec![0i64; nv * ne];
    for (chunk, events) in events.chunks(128).enumerate() {
        let start = chunk * 128;
        let n = events.len();
        let mut evidence = BTreeMap::new();
        for event in events {
            if let std::collections::btree_map::Entry::Vacant(slot) = evidence.entry(event.gene_idx)
            {
                slot.insert(reader.gene(event.gene_idx, options.sample_idx)?);
            }
        }
        let result: Vec<Vec<Verified>> = events
            .par_iter()
            .map(|event| {
                evidence[&event.gene_idx]
                    .iter()
                    .map(|counts| {
                        verify::verify_event(event, &genes[event.gene_idx], counts, options.verify)
                    })
                    .collect::<Result<_>>()
            })
            .collect::<Result<_>>()?;
        let mut counts = vec![0.0; ns * nf * n];
        let mut flags = vec![false; ns * nv * n];
        let mut psi = vec![0.0; ns * n];
        let mut iso1 = vec![0i32; ns * n];
        let mut iso2 = vec![0i32; ns * n];
        for (e, result) in result.iter().enumerate() {
            for (s, result) in result.iter().enumerate() {
                for (f, &value) in result.info.iter().enumerate() {
                    counts[(s * nf + f) * n + e] = value;
                }
                for (v, &value) in result.verified.iter().enumerate() {
                    flags[(s * nv + v) * n + e] = value;
                    num_verified[v * ne + start + e] += i64::from(value);
                }
                let (p, a, b) = verify::psi(&result.info, kind, options.psi_min_reads);
                psi[s * n + e] = p;
                // NumPy's float->int32 conversion yields INT_MIN for NaN/overflow.
                let cast = |x: f64| {
                    if x.is_finite() && x.trunc() >= i32::MIN as f64 && x.trunc() <= i32::MAX as f64
                    {
                        x as i32
                    } else {
                        i32::MIN
                    }
                };
                iso1[s * n + e] = cast(b);
                iso2[s * n + e] = cast(a);
            }
        }
        counts_out.write_slice(
            ArrayView3::from_shape((ns, nf, n), &counts)?,
            (.., .., start..start + n),
        )?;
        verified_out.write_slice(
            ArrayView3::from_shape((ns, nv, n), &flags)?,
            (.., .., start..start + n),
        )?;
        psi_out.write_slice(
            ArrayView2::from_shape((ns, n), &psi)?,
            (.., start..start + n),
        )?;
        iso1_out.write_slice(
            ArrayView2::from_shape((ns, n), &iso1)?,
            (.., start..start + n),
        )?;
        iso2_out.write_slice(
            ArrayView2::from_shape((ns, n), &iso2)?,
            (.., start..start + n),
        )?;
    }
    hdf5io::array(
        &file,
        "gene_idx",
        &[ne],
        &events.iter().map(|e| e.gene_idx as i64).collect::<Vec<_>>(),
        true,
    )?;
    let positions: Vec<_> = events.iter().map(event_positions).collect();
    ensure!(
        positions.iter().all(|p| p.len() == positions[0].len()),
        "event positions have inconsistent lengths"
    );
    let shape = if kind == EventType::MutexExons {
        vec![ne, 2, 4]
    } else {
        vec![ne, positions[0].len()]
    };
    hdf5io::array(
        &file,
        "event_pos",
        &shape,
        &positions.into_iter().flatten().collect::<Vec<_>>(),
        false,
    )?;
    let confirmed: Vec<_> = (0..ne)
        .map(|e| (0..nv).map(|v| num_verified[v * ne + e]).min().unwrap())
        .collect();
    hdf5io::array(&file, "num_verified", &[nv, ne], &num_verified, false)?;
    hdf5io::array(&file, "confirmed", &[ne], &confirmed, false)?;
    let indices: Vec<_> = confirmed
        .iter()
        .enumerate()
        .filter_map(|(i, &n)| (n >= 1).then_some(i as i64))
        .collect();
    if !indices.is_empty() {
        hdf5io::array(&file, "conf_idx", &[indices.len()], &indices, false)?;
    }
    file.close()?;
    Ok(indices.into_iter().map(|i| i as usize).collect())
}
