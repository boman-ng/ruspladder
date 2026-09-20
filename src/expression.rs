// Adapted from SplAdder v3.1.1 count.py (BSD-3-Clause).
use crate::{annotation::Gene, hdf5io, numeric};
use anyhow::{Result, ensure};
use hdf5::File;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Normalization {
    Geomean,
    Tc,
    Uq,
}

fn percentile(values: &mut [f64], fraction: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    if values.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    values.sort_unstable_by(f64::total_cmp);
    let index = fraction * (values.len() - 1) as f64;
    let lower = index.floor() as usize;
    let upper = index.ceil() as usize;
    if lower == upper {
        values[lower]
    } else {
        values[lower] * (upper as f64 - index) + values[upper] * (index - lower as f64)
    }
}

pub fn size_factors(counts: &[Vec<f64>], samples: usize, kind: Normalization) -> Result<Vec<f64>> {
    size_factors_impl(counts, samples, kind, numeric::sum)
}

/// `_get_gene_expression` selects columns by advanced indexing, creating a
/// Fortran-order array. NumPy reduces these noncontiguous rows sequentially.
pub(crate) fn size_factors_column_subset(counts: &[Vec<f64>], samples: usize) -> Result<Vec<f64>> {
    size_factors_impl(
        counts,
        samples,
        Normalization::Geomean,
        if counts.len() > 1 {
            numeric::sum_strided
        } else {
            numeric::sum
        },
    )
}

fn size_factors_impl(
    counts: &[Vec<f64>],
    samples: usize,
    kind: Normalization,
    sum: fn(&[f64]) -> f64,
) -> Result<Vec<f64>> {
    ensure!(
        !matches!(kind, Normalization::Uq),
        "SplAdder v3.1.1 uq normalization fails: scoreatpercentile is not defined"
    );
    ensure!(
        counts.iter().all(|row| row.len() == samples),
        "normalization sample shape mismatch"
    );
    let means: Vec<_> = counts
        .iter()
        .map(|row| {
            let logs: Vec<_> = row.iter().map(|&x| (x + 1.0).ln()).collect();
            (sum(&logs) / samples as f64).exp()
        })
        .collect();
    Ok((0..samples)
        .map(|s| {
            let mut values: Vec<_> = counts
                .iter()
                .enumerate()
                .filter_map(|(g, row)| match kind {
                    Normalization::Geomean => (row[s] > 0.0).then_some(row[s] / means[g]),
                    _ => Some(row[s]),
                })
                .collect();
            match kind {
                Normalization::Geomean => percentile(&mut values, 0.5),
                Normalization::Tc => numeric::sum(&values) / 100_000_000.0,
                Normalization::Uq => unreachable!(),
            }
        })
        .collect())
}

#[derive(Debug, Serialize)]
pub struct Expression {
    pub raw_count: Vec<Vec<f64>>,
    pub raw_count_non_alt: Vec<Vec<f64>>,
    pub size_factors: Vec<f64>,
    pub size_factors_non_alt: Vec<f64>,
}

pub fn compute(
    genes: &[Gene],
    counts_path: &Path,
    readlen: f64,
    sample_idx: Option<&[usize]>,
    output: Option<&Path>,
) -> Result<Expression> {
    let input = File::open(counts_path)?;
    let samples = hdf5io::read_strings(&input.dataset("samples")?)?;
    let indices: Vec<_> =
        sample_idx.map_or_else(|| (0..samples.len()).collect(), <[usize]>::to_vec);
    let segments = input.dataset("segments")?;
    ensure!(
        readlen > 0.0
            && indices
                .iter()
                .all(|&s| s < segments.shape()[1] && s < samples.len()),
        "invalid expression sample index or read length"
    );
    let lengths = input.dataset("seg_len")?.read_raw::<i64>()?;
    let gene_ids = input.dataset("gene_ids_segs")?.read_raw::<i64>()?;
    let names = hdf5io::read_strings(&input.dataset("gene_names")?)?;
    let mut first = vec![None; genes.len()];
    for (row, &id) in gene_ids.iter().enumerate() {
        ensure!(
            id >= 0 && (id as usize) < first.len(),
            "invalid segment gene index"
        );
        first[id as usize].get_or_insert(row);
    }
    let mut result = Expression {
        raw_count: Vec::new(),
        raw_count_non_alt: Vec::new(),
        size_factors: Vec::new(),
        size_factors_non_alt: Vec::new(),
    };
    for (g, gene) in genes.iter().enumerate() {
        ensure!(
            names.get(g) == Some(&gene.name),
            "expression gene order mismatch"
        );
        let start = first[g].ok_or_else(|| anyhow::anyhow!("gene has no counted segments"))?;
        let end = start + gene.segmentgraph.segments.len();
        ensure!(end <= gene_ids.len(), "expression segment range mismatch");
        let values = segments.read_slice_2d::<f64, _>((start..end, ..))?;
        let non_alt = gene.segmentgraph.non_alternative_segments();
        ensure!(
            !non_alt.is_empty(),
            "SplAdder v3.1.1 requires non-alternative segments for expression"
        );
        let quantify = |selected: &[usize]| -> Vec<f64> {
            indices
                .iter()
                .map(|&s| {
                    selected
                        .iter()
                        .map(|&i| values[[i, s]] * lengths[start + i] as f64)
                        .sum::<f64>()
                        / readlen
                })
                .collect()
        };
        result
            .raw_count
            .push(quantify(&(0..end - start).collect::<Vec<_>>()));
        result.raw_count_non_alt.push(quantify(&non_alt));
    }
    result.size_factors = size_factors(&result.raw_count, indices.len(), Normalization::Geomean)?;
    result.size_factors_non_alt = size_factors(
        &result.raw_count_non_alt,
        indices.len(),
        Normalization::Geomean,
    )?;
    if let Some(path) = output {
        let out = File::create(path)?;
        // NumPy slicing preserves the input byte width, including the S255
        // labels produced by per-sample count collection.
        let sample_dataset = out
            .new_dataset_builder()
            .empty_as(&input.dataset("samples")?.dtype()?.to_descriptor()?)
            .shape(indices.len())
            .deflate(4)
            .chunk_min_kb(64)
            .create("samples")?;
        hdf5io::write_strings(
            &sample_dataset,
            &indices
                .iter()
                .map(|&i| samples[i].as_str())
                .collect::<Vec<_>>(),
        )?;
        out.link_soft("/samples", "strains")?;
        hdf5io::strings(
            &out,
            "gene_ids",
            &[genes.len()],
            &genes.iter().map(|g| g.name.as_str()).collect::<Vec<_>>(),
            true,
        )?;
        hdf5io::strings(
            &out,
            "gene_symbols",
            &[genes.len()],
            &genes
                .iter()
                .map(|g| g.symbol.as_deref().unwrap_or("None"))
                .collect::<Vec<_>>(),
            true,
        )?;
        for (name, values) in [
            ("raw_count", &result.raw_count),
            ("raw_count_non_alt", &result.raw_count_non_alt),
        ] {
            hdf5io::array(
                &out,
                name,
                &[genes.len(), indices.len()],
                &values.iter().flatten().copied().collect::<Vec<_>>(),
                true,
            )?;
        }
        hdf5io::array(
            &out,
            "size_factors",
            &[indices.len()],
            &result.size_factors,
            true,
        )?;
        hdf5io::array(
            &out,
            "size_factors_non_alt",
            &[indices.len()],
            &result.size_factors_non_alt,
            true,
        )?;
        out.close()?;
    }
    Ok(result)
}
