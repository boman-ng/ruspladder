// Adapted from SplAdder v3.1.1 spladder_test.adj_pval and statsmodels 0.14.4
// stats/multitest.py. BSD-3-Clause; see licenses/{SplAdder,statsmodels}-BSD.txt.
use crate::numeric;
use anyhow::{Result, ensure};
use serde::Deserialize;

#[derive(Clone, Copy, Deserialize)]
pub enum Correction {
    BH,
    Bonferroni,
    Holm,
    Hochberg,
    Hommel,
    BY,
    TSBH,
}

pub fn adjust(pvalues: &[f64], method: Correction) -> Result<Vec<f64>> {
    let mut order: Vec<_> = (0..pvalues.len())
        .filter(|&i| !pvalues[i].is_nan())
        .collect();
    ensure!(
        !order.is_empty(),
        "SplAdder's multiple-testing correction requires non-NaN p-values"
    );
    ensure!(
        !matches!(method, Correction::TSBH),
        "SplAdder v3.1.1 passes unrecognized statsmodels method tsbh"
    );
    order.sort_unstable_by(|&a, &b| pvalues[a].total_cmp(&pvalues[b]));
    let p: Vec<_> = order.iter().map(|&i| pvalues[i]).collect();
    let n = p.len();
    let mut adjusted = p.clone();
    match method {
        Correction::Bonferroni => {
            for x in &mut adjusted {
                *x *= n as f64;
            }
        }
        Correction::Holm | Correction::Hochberg => {
            for (i, x) in adjusted.iter_mut().enumerate() {
                *x *= (n - i) as f64;
            }
            if matches!(method, Correction::Holm) {
                for i in 1..n {
                    adjusted[i] = adjusted[i].max(adjusted[i - 1]);
                }
            } else {
                for i in (0..n - 1).rev() {
                    adjusted[i] = adjusted[i].min(adjusted[i + 1]);
                }
            }
        }
        Correction::Hommel => {
            for m in (2..=n).rev() {
                let bound = p[n - m..]
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| m as f64 * x / (i + 1) as f64)
                    .fold(f64::INFINITY, f64::min);
                for x in &mut adjusted[n - m..] {
                    *x = x.max(bound);
                }
                for i in 0..n - m {
                    adjusted[i] = adjusted[i].max((m as f64 * p[i]).min(bound));
                }
            }
        }
        Correction::BH | Correction::BY => {
            let harmonic = if matches!(method, Correction::BY) {
                numeric::sum(&(1..=n).map(|i| 1.0 / i as f64).collect::<Vec<_>>())
            } else {
                1.0
            };
            for (i, x) in adjusted.iter_mut().enumerate() {
                *x /= ((i + 1) as f64 / n as f64) / harmonic;
            }
            for i in (0..n - 1).rev() {
                adjusted[i] = adjusted[i].min(adjusted[i + 1]);
            }
        }
        Correction::TSBH => unreachable!(),
    }
    let mut result = pvalues.to_vec();
    for (i, x) in order.into_iter().zip(adjusted) {
        result[i] = x.min(1.0);
    }
    Ok(result)
}
