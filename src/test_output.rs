// Adapted from SplAdder v3.1.1 spladder_test.py output section (BSD-3-Clause).
use crate::{
    correction::{self, Correction},
    events::{Event, EventType},
    hdf5io, sort,
    statistics::TestingResult,
    testing::{self, Prepared},
};
use anyhow::{Result, ensure};
use serde::Deserialize;
use std::{
    collections::BTreeSet,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Deserialize)]
pub struct Metadata {
    pub samples: Vec<String>,
    pub gene_ids: Vec<String>,
    pub gene_symbols: Vec<String>,
}
#[derive(Deserialize)]
pub struct Options {
    pub confidence: u8,
    pub kind: EventType,
    pub correction: Correction,
    pub group_a: usize,
    pub labels: [String; 2],
}

fn float(value: f64) -> String {
    if value.is_nan() {
        return "nan".into();
    }
    let repr = format!("{value:?}");
    if let Some((mantissa, exponent)) = repr.split_once('e') {
        let exponent: i32 = exponent.parse().unwrap();
        format!("{mantissa}e{exponent:+03}")
    } else {
        repr
    }
}

pub fn write_results(
    outdir: &Path,
    options: &Options,
    metadata: &Metadata,
    events: &[Event],
    prepared: &Prepared,
    result: &TestingResult,
) -> Result<()> {
    let n = result.pvalues.len();
    ensure!(
        n == prepared.event_idx.len()
            && n == prepared.gene_idx.len()
            && result.coverage.len() == n
            && result.dispersion_raw.len() == n
            && result.dispersion_adjusted.len() == n
            && prepared.event_idx.iter().all(|&i| i < events.len())
            && prepared
                .gene_idx
                .iter()
                .all(|&i| i < metadata.gene_ids.len() && i < metadata.gene_symbols.len()),
        "test result shape mismatch"
    );
    std::fs::create_dir_all(outdir)?;
    let corrected = correction::adjust(&result.pvalues, options.correction)?;
    let means =
        testing::means_and_fold_changes(&result.coverage, &prepared.size_factors, options.group_a);
    let order = sort::argsort_f64(&result.pvalues);
    let mut header = vec![
        "event_id",
        "chrm",
        "exon_pos",
        "alt_usage",
        "gene_id",
        "gene_name",
        "p_val",
        "p_val_adj",
        "dPSI",
        "mean_event_count_A",
        "mean_event_count_B",
        "log2FC_event_count",
        "mean_gene_exp_A",
        "mean_gene_exp_B",
        "log2FC_gene_exp",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect::<Vec<_>>();
    let tag = format!("C{}_{}", options.confidence, options.kind.as_str());
    let mut ordinary = BufWriter::new(File::create(
        outdir.join(format!("test_results_{tag}.tsv")),
    )?);
    let mut unique = BufWriter::new(File::create(
        outdir.join(format!("test_results_{tag}.gene_unique.tsv")),
    )?);
    let mut extended = BufWriter::new(File::create(
        outdir.join(format!("test_results_extended_{tag}.tsv")),
    )?);
    writeln!(ordinary, "{}", header.join("\t"))?;
    writeln!(unique, "{}", header.join("\t"))?;
    for prefix in ["event_count", "gene_exp"] {
        header.extend(metadata.samples.iter().map(|s| format!("{prefix}:{s}")));
    }
    header.extend(["disp_raw".into(), "disp_adj".into()]);
    writeln!(extended, "{}", header.join("\t"))?;
    let mut seen = BTreeSet::new();
    for i in order {
        let event = &events[prepared.event_idx[i]];
        let gene = prepared.gene_idx[i];
        let (positions, usage) = event.coordinate_strings();
        let mut row = vec![
            format!("{}.{}", options.kind.as_str(), prepared.event_idx[i] + 1),
            event.chr.clone(),
            positions,
            usage,
            metadata.gene_ids[gene].clone(),
            metadata.gene_symbols[gene].clone(),
            float(result.pvalues[i]),
            float(corrected[i]),
            float(prepared.delta_psi[i]),
        ];
        row.extend(means[i].iter().map(|&x| float(x)));
        let line = row.join("\t");
        writeln!(ordinary, "{line}")?;
        if seen.insert(&metadata.gene_ids[gene]) {
            writeln!(unique, "{line}")?;
        }
        row.extend(
            result.coverage[i]
                .iter()
                .zip(&prepared.size_factors)
                .map(|(&x, &s)| float(x / s)),
        );
        row.extend([
            float(result.dispersion_raw[i]),
            float(result.dispersion_adjusted[i]),
        ]);
        writeln!(extended, "{}", row.join("\t"))?;
    }
    ordinary.flush()?;
    unique.flush()?;
    extended.flush()?;
    // Replace the internal Python pickle setup cache with native HDF5.
    let setup = hdf5::File::create(outdir.join(format!("test_setup_{tag}.hdf5")))?;
    for name in ["gene_samples", "event_samples"] {
        hdf5io::strings(
            &setup,
            name,
            &[metadata.samples.len()],
            &metadata
                .samples
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            false,
        )?;
    }
    for (name, matrix) in [
        ("dmatrix0", &prepared.null),
        ("dmatrix1", &prepared.alternative),
    ] {
        hdf5io::array(
            &setup,
            name,
            &[matrix.len(), matrix[0].len()],
            &matrix
                .iter()
                .flatten()
                .map(|&v| v as i64)
                .collect::<Vec<_>>(),
            false,
        )?;
    }
    hdf5io::strings(&setup, "event_type", &[], &[options.kind.as_str()], false)?;
    hdf5io::strings(
        &setup,
        "labels",
        &[2],
        &options
            .labels
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        false,
    )?;
    setup.flush()?;
    Ok(())
}
