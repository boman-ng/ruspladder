// Adapted from SplAdder v3.1.1 spladder_test.py (BSD-3-Clause).
use crate::{
    expression::{self, Normalization},
    numeric,
    quantify::Quantified,
};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Deserialize)]
pub struct InputOptions {
    pub cap_outliers: bool,
    pub max_zero_fraction: f64,
    pub min_dpsi: f64,
}

#[derive(Debug, Serialize)]
pub struct ExpressionInput {
    pub counts: Vec<Vec<f64>>,
    pub size_factors: Vec<f64>,
    pub capped: usize,
}

// scipy.stats.scoreatpercentile uses a weighted sum, unlike np.percentile's lerp.
fn percentile(values: &[f64], fraction: f64) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);
    let index = fraction * (sorted.len() - 1) as f64;
    let lo = index.floor() as usize;
    let hi = index.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        sorted[lo] * (hi as f64 - index) + sorted[hi] * (index - lo as f64)
    }
}

fn cap(row: &mut [f64], factors: Option<&[f64]>, multiplier: f64) -> usize {
    let logs: Vec<_> = row
        .iter()
        .enumerate()
        .map(|(i, &x)| (x / factors.map_or(1.0, |f| f[i]) + 1.0).log2())
        .collect();
    let lower = percentile(&logs, 0.25);
    let upper = percentile(&logs, 0.75);
    let spread = upper - lower;
    if spread <= 0.0 || spread.is_nan() {
        return 0;
    }
    let threshold = upper + multiplier * spread;
    let value = 2.0_f64.powf(threshold) - 1.0;
    let mut capped = 0;
    for (i, x) in row.iter_mut().enumerate() {
        if logs[i] > threshold {
            *x = value * factors.map_or(1.0, |f| f[i]);
            capped += 1;
        }
    }
    capped
}

pub fn prepare_expression(
    mut counts: Vec<Vec<f64>>,
    samples: usize,
    cap_outliers: bool,
) -> Result<ExpressionInput> {
    ensure!(samples > 0, "testing requires samples");
    let size_factors = expression::size_factors(&counts, samples, Normalization::Geomean)?;
    let capped = if cap_outliers {
        counts
            .iter_mut()
            .map(|row| cap(row, Some(&size_factors), 1.5))
            .sum()
    } else {
        0
    };
    Ok(ExpressionInput {
        counts,
        size_factors,
        capped,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Prepared {
    pub counts: Vec<Vec<f64>>,
    pub size_factors: Vec<f64>,
    pub event_size_factors: Vec<f64>,
    pub null: Vec<Vec<f64>>,
    pub alternative: Vec<Vec<f64>>,
    pub selected: Vec<bool>,
    pub delta_psi: Vec<f64>,
    pub event_idx: Vec<usize>,
    pub gene_idx: Vec<usize>,
}

fn nanmean(values: &[f64]) -> f64 {
    let count = values.iter().filter(|x| !x.is_nan()).count();
    numeric::sum(
        &values
            .iter()
            .map(|&x| if x.is_nan() { 0.0 } else { x })
            .collect::<Vec<_>>(),
    ) / count as f64
}

pub fn prepare_events(
    mut events: Quantified,
    expression: &ExpressionInput,
    group_a: usize,
    options: InputOptions,
) -> Result<Option<Prepared>> {
    let samples = expression.size_factors.len();
    ensure!(
        group_a > 0 && group_a < samples,
        "testing requires two nonempty conditions"
    );
    let n = events.gene_idx.len();
    ensure!(
        events
            .coverage
            .iter()
            .all(|c| c.len() == n && c.iter().all(|r| r.len() == samples))
            && events.psi.len() == n
            && events.psi.iter().all(|r| r.len() == samples)
            && events.event_idx.len() == n
            && events.gene_idx.iter().all(|&g| g < expression.counts.len()),
        "test quantification shape mismatch"
    );
    if options.cap_outliers {
        for rows in &mut events.coverage {
            for row in rows {
                cap(row, None, 3.0);
            }
        }
    }
    let combined: Vec<_> = events.coverage.iter().flatten().cloned().collect();
    let event_size_factors = expression::size_factors(&combined, samples, Normalization::Geomean)?;
    let eligible = |row: &[f64]| {
        let low = |values: &[f64]| {
            values.iter().filter(|&&x| x <= 1.0).count() as f64 / values.len() as f64
        };
        low(&row[..group_a]) <= options.max_zero_fraction
            || low(&row[group_a..]) <= options.max_zero_fraction
    };
    let masks = events
        .coverage
        .each_ref()
        .map(|rows| rows.iter().map(|r| eligible(r)).collect::<Vec<_>>());
    let retained: Vec<_> = (0..n).filter(|&i| masks[0][i] || masks[1][i]).collect();
    if retained.is_empty() {
        return Ok(None);
    }
    let delta_psi: Vec<_> = retained
        .iter()
        .map(|&i| {
            let row = &events.psi[i];
            if row[..group_a].iter().all(|x| x.is_nan())
                || row[group_a..].iter().all(|x| x.is_nan())
            {
                0.0
            } else {
                nanmean(&row[..group_a]) - nanmean(&row[group_a..])
            }
        })
        .collect();
    let mut counts = Vec::with_capacity(retained.len() * 2);
    let mut selected = Vec::with_capacity(retained.len() * 2);
    for (isoform, mask) in masks.iter().enumerate() {
        for (j, &i) in retained.iter().enumerate() {
            let row = events.coverage[isoform][i]
                .iter()
                .copied()
                .chain(
                    expression.counts[events.gene_idx[i]]
                        .iter()
                        .map(|&x| x + 0.5),
                )
                .map(f64::round_ties_even)
                .collect();
            counts.push(row);
            selected.push(mask[i] && delta_psi[j].abs() >= options.min_dpsi);
        }
    }
    let alternative: Vec<_> = (0..samples * 2)
        .map(|i| {
            vec![
                1.0,
                if i < group_a { 1.0 } else { 0.0 },
                if i % samples < group_a { 1.0 } else { 0.0 },
                if i >= samples { 1.0 } else { 0.0 },
            ]
        })
        .collect();
    let null = alternative
        .iter()
        .map(|row| vec![row[0], row[2], row[3]])
        .collect();
    let size_factors = event_size_factors
        .iter()
        .chain(&expression.size_factors)
        .copied()
        .collect();
    Ok(Some(Prepared {
        counts,
        size_factors,
        event_size_factors,
        null,
        alternative,
        selected,
        delta_psi,
        event_idx: retained.iter().map(|&i| events.event_idx[i]).collect(),
        gene_idx: retained.iter().map(|&i| events.gene_idx[i]).collect(),
    }))
}

pub fn means_and_fold_changes(
    coverage: &[Vec<f64>],
    factors: &[f64],
    group_a: usize,
) -> Vec<[f64; 6]> {
    let samples = factors.len() / 2;
    coverage
        .iter()
        .map(|row| {
            let normalized: Vec<_> = row.iter().zip(factors).map(|(&x, &f)| x / f).collect();
            let a = nanmean(&normalized[..group_a]);
            let b = nanmean(&normalized[group_a..samples]);
            let ga = nanmean(&normalized[samples..samples + group_a]);
            let gb = nanmean(&normalized[samples + group_a..]);
            [a, b, a.log2() - b.log2(), ga, gb, ga.log2() - gb.log2()]
        })
        .collect()
}
