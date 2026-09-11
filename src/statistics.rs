// Adapted from SplAdder v3.1.1 spladder_test.py (BSD-3-Clause).
use crate::{
    glm::{self, Family},
    likelihood, numeric, optimize,
};
use anyhow::{Result, ensure};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Deserialize)]
pub struct TestOptions {
    pub min_count: f64,
    pub max_zero_fraction: f64,
}

#[derive(Debug, Serialize)]
pub struct Dispersions {
    pub values: Vec<f64>,
    pub converged: Vec<bool>,
}

fn estimate_one(
    counts: &[f64],
    design: &[Vec<f64>],
    sf: &[f64],
    min_count: f64,
    prior: Option<(f64, f64)>,
) -> Result<(f64, bool)> {
    let response: Vec<_> = counts.iter().map(|&x| x.trunc()).collect();
    if prior.is_none()
        && (response.iter().zip(sf).map(|(&y, &sf)| y / sf).sum::<f64>() < min_count
            || response.iter().filter(|&&y| y == 0.0).count() as f64 / response.len() as f64 > 0.6)
    {
        return Ok((f64::NAN, false));
    }
    let offset: Vec<_> = sf.iter().map(|s| s.ln()).collect();
    let mut dispersion = 0.1;
    for _ in 0..10 {
        let fit = glm::fit(
            &response,
            design,
            &offset,
            Family::NegativeBinomial { alpha: dispersion },
        )?;
        let previous = dispersion;
        dispersion = optimize::bounded(
            |disp| {
                Ok(-match prior {
                    None => likelihood::adjusted(disp, design, &response, &fit.mu),
                    Some((fitted, prior)) => {
                        likelihood::adjusted_shrink(disp, design, &response, &fit.mu, fitted, prior)
                    }
                })
            },
            [0.0, 10.0],
            1e-5,
        )?
        .x;
        if (dispersion.ln() - previous.ln()).abs() < 1e-4 {
            return Ok((dispersion, true));
        }
    }
    Ok((dispersion, false))
}

fn pack(values: Vec<(f64, bool)>) -> Dispersions {
    let (values, converged) = values.into_iter().unzip();
    Dispersions { values, converged }
}

pub fn estimate(
    counts: &[Vec<f64>],
    design: &[Vec<f64>],
    sf: &[f64],
    min_count: f64,
) -> Result<Dispersions> {
    ensure!(
        counts.iter().all(|c| c.len() == sf.len()),
        "dispersion sample shape mismatch"
    );
    let result = pack(
        counts
            .par_iter()
            .map(|counts| estimate_one(counts, design, sf, min_count, None))
            .collect::<Result<_>>()?,
    );
    ensure!(
        result.converged.iter().any(|&c| c),
        "none of the dispersion estimates converged"
    );
    Ok(result)
}

fn percentile(values: &[f64], fraction: f64) -> f64 {
    if values.is_empty() || values.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);
    let index = (sorted.len() - 1) as f64 * fraction;
    let lower = index.floor() as usize;
    let upper = index.ceil() as usize;
    let weight = index - lower as f64;
    let delta = sorted[upper] - sorted[lower];
    if weight < 0.5 {
        sorted[lower] + delta * weight
    } else {
        sorted[upper] - delta * (1.0 - weight)
    }
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() || values.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);
    let mid = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

#[derive(Debug, Serialize)]
pub struct Trend {
    pub values: Vec<f64>,
    pub coefficients: Vec<f64>,
    pub indices: Vec<usize>,
}

/// Source replaces raw NaNs with the lower fitted boundary in place.
pub fn fit_trend(counts: &[Vec<f64>], raw: &mut Dispersions, sf: &[f64]) -> Result<Trend> {
    let means: Vec<_> = counts
        .iter()
        .map(|row| {
            numeric::sum(&row.iter().zip(sf).map(|(&y, &s)| y / s).collect::<Vec<_>>())
                / sf.len() as f64
        })
        .collect();
    let mut unique: Vec<_> = raw
        .values
        .iter()
        .zip(&raw.converged)
        .filter_map(|(&x, &c)| c.then_some(x))
        .collect();
    unique.sort_unstable_by(f64::total_cmp);
    unique.dedup();
    let lower = percentile(&unique, 0.01);
    let upper = percentile(&unique, 0.99);
    for x in &mut raw.values {
        if x.is_nan() {
            *x = lower;
        }
    }
    let indices: Vec<_> = raw
        .values
        .iter()
        .enumerate()
        .filter_map(|(i, &x)| (x > lower && x < upper).then_some(i))
        .collect();
    let design: Vec<_> = indices.iter().map(|&i| vec![1.0 / means[i], 1.0]).collect();
    let response: Vec<_> = indices.iter().map(|&i| raw.values[i]).collect();
    let fit = glm::fit(
        &response,
        &design,
        &vec![0.0; indices.len()],
        Family::GammaIdentity,
    )?;
    let values = raw
        .values
        .iter()
        .zip(means)
        .map(|(&raw, mean)| {
            if raw.is_nan() {
                raw
            } else {
                fit.params[0] / mean + fit.params[1]
            }
        })
        .collect();
    Ok(Trend {
        values,
        coefficients: fit.params,
        indices,
    })
}

pub fn variance_prior(raw: &[f64], fitted: &[f64], indices: &[usize], sample_variance: f64) -> f64 {
    let residuals: Vec<_> = indices
        .iter()
        .map(|&i| raw[i].ln() - fitted[i].ln())
        .collect();
    let middle = median(&residuals);
    let deviation = median(
        &residuals
            .iter()
            .map(|x| (x - middle).abs())
            .collect::<Vec<_>>(),
    ) * 1.4826;
    let prior = deviation.powi(2) - sample_variance;
    if prior < 0.1 { 0.1 } else { prior }
}

pub fn adjust(
    counts: &[Vec<f64>],
    design: &[Vec<f64>],
    raw: &Dispersions,
    trend: &Trend,
    sf: &[f64],
) -> Result<(Dispersions, f64)> {
    let sample_variance =
        likelihood::trigamma((design.len() as f64 - design[0].len() as f64) / 2.0);
    let prior = variance_prior(&raw.values, &trend.values, &trend.indices, sample_variance);
    let estimates = counts
        .par_iter()
        .enumerate()
        .map(|(i, counts)| {
            if raw.values[i].is_nan() {
                Ok((f64::NAN, false))
            } else {
                estimate_one(counts, design, sf, 0.0, Some((trend.values[i], prior)))
            }
        })
        .collect::<Result<_>>()?;
    Ok((pack(estimates), prior))
}

pub fn test(
    counts: &[Vec<f64>],
    adjusted: &[f64],
    sf: &[f64],
    null: &[Vec<f64>],
    alternative: &[Vec<f64>],
    selected: &[bool],
    max_zero_fraction: f64,
) -> Result<Vec<f64>> {
    ensure!(
        !null.is_empty() && alternative[0].len() == null[0].len() + 1,
        "SplAdder's LRT requires one additional coefficient"
    );
    let offset: Vec<_> = sf.iter().map(|s| s.ln()).collect();
    Ok(counts
        .par_iter()
        .enumerate()
        .map(|(i, counts)| {
            if adjusted[i].is_nan() || !selected[i] {
                return f64::NAN;
            }
            let response: Vec<_> = counts.iter().map(|y| y.trunc()).collect();
            if response[..response.len() / 2]
                .iter()
                .filter(|&&y| y == 0.0)
                .count() as f64
                > max_zero_fraction * response.len() as f64 / 2.0
            {
                return f64::NAN;
            }
            let family = Family::NegativeBinomial { alpha: adjusted[i] };
            let (Ok(fit0), Ok(fit1)) = (
                glm::fit(&response, null, &offset, family),
                glm::fit(&response, alternative, &offset, family),
            ) else {
                return f64::NAN;
            };
            likelihood::chi2_pvalue(fit0.deviance - fit1.deviance)
        })
        .collect())
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TestingResult {
    pub pvalues: Vec<f64>,
    pub coverage: Vec<Vec<f64>>,
    pub dispersion_raw: Vec<f64>,
    pub dispersion_adjusted: Vec<f64>,
}

pub fn run(
    counts: &[Vec<f64>],
    null: &[Vec<f64>],
    alternative: &[Vec<f64>],
    sf: &[f64],
    selected: &[bool],
    options: TestOptions,
) -> Result<TestingResult> {
    ensure!(
        !counts.is_empty() && counts.len().is_multiple_of(2),
        "testing requires paired isoform count rows"
    );
    let mut raw = estimate(counts, alternative, sf, options.min_count)?;
    let trend = fit_trend(counts, &mut raw, sf)?;
    let (adjusted, _) = adjust(counts, alternative, &raw, &trend, sf)?;
    let pvalues = test(
        counts,
        &adjusted.values,
        sf,
        null,
        alternative,
        selected,
        options.max_zero_fraction,
    )?;
    let offset = counts.len() / 2;
    let mut result = TestingResult {
        pvalues: Vec::new(),
        coverage: Vec::new(),
        dispersion_raw: Vec::new(),
        dispersion_adjusted: Vec::new(),
    };
    for i in 0..offset {
        let (a, b) = (pvalues[i], pvalues[i + offset]);
        // The reference chooses the larger of the two finite isoform p-values;
        // a tie keeps the first. One available isoform is used on its own.
        let second = (!b.is_nan() && a.is_nan()) || (!a.is_nan() && !b.is_nan() && b > a);
        let row = if second { i + offset } else { i };
        let p = pvalues[row];
        result
            .pvalues
            .push(if p.is_nan() || p > 1.0 { 1.0 } else { p });
        result.coverage.push(counts[row].clone());
        result.dispersion_raw.push(raw.values[row]);
        result.dispersion_adjusted.push(adjusted.values[row]);
    }
    Ok(result)
}
