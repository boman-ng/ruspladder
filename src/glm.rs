// Adapted from statsmodels 0.14.4 GLM._fit_irls, _MinimalWLS, families and links.
// BSD-3-Clause; see licenses/statsmodels-BSD.txt. Linear solves reuse the reference OpenBLAS/LAPACK kernels.
use crate::{linalg, numeric};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Family {
    NegativeBinomial { alpha: f64 },
    GammaIdentity,
}

fn clean(value: f64) -> f64 {
    if value < f64::EPSILON {
        f64::EPSILON
    } else {
        value
    }
}

impl Family {
    fn link(self, mu: f64) -> f64 {
        match self {
            Self::NegativeBinomial { .. } => clean(mu).ln(),
            Self::GammaIdentity => mu,
        }
    }
    fn inverse(self, eta: f64) -> f64 {
        match self {
            Self::NegativeBinomial { .. } => eta.exp(),
            Self::GammaIdentity => eta,
        }
    }
    fn derivative(self, mu: f64) -> f64 {
        match self {
            Self::NegativeBinomial { .. } => 1.0 / clean(mu),
            Self::GammaIdentity => 1.0,
        }
    }
    fn variance(self, mu: f64) -> f64 {
        match self {
            Self::NegativeBinomial { alpha } => {
                let mu = clean(mu);
                mu + alpha * mu.powi(2)
            }
            Self::GammaIdentity => mu.powi(2),
        }
    }
    fn deviance(self, y: &[f64], mu: &[f64], scale: f64) -> f64 {
        let values: Vec<_> = y
            .iter()
            .zip(mu)
            .map(|(&y, &mu)| {
                let ratio = clean(y / mu);
                let residual = match self {
                    Self::NegativeBinomial { alpha } => {
                        y * ratio.ln()
                            - (y + 1.0 / alpha) * ((y + 1.0 / alpha) / (mu + 1.0 / alpha)).ln()
                    }
                    Self::GammaIdentity => -ratio.ln() + (y - mu) / mu,
                };
                2.0 * residual / scale
            })
            .collect();
        numeric::sum(&values)
    }
    fn scale(self, y: &[f64], mu: &[f64], df_resid: usize) -> f64 {
        match self {
            Self::NegativeBinomial { .. } => 1.0,
            Self::GammaIdentity => {
                numeric::sum(
                    &y.iter()
                        .zip(mu)
                        .map(|(&y, &mu)| (y - mu).powi(2) / self.variance(mu))
                        .collect::<Vec<_>>(),
                ) / df_resid as f64
            }
        }
    }
}

#[derive(Debug, Serialize)]
pub struct Fit {
    pub params: Vec<f64>,
    pub mu: Vec<f64>,
    pub deviance: f64,
    pub scale: f64,
    pub converged: bool,
    pub iterations: usize,
}

pub fn fit(y: &[f64], design: &[Vec<f64>], offset: &[f64], family: Family) -> Result<Fit> {
    let n = y.len();
    ensure!(
        n > 0 && design.len() == n && offset.len() == n,
        "GLM observation shape mismatch"
    );
    let k = design[0].len();
    ensure!(
        k > 0 && design.iter().all(|r| r.len() == k),
        "GLM design shape mismatch"
    );
    ensure!(
        y.iter()
            .chain(offset)
            .chain(design.iter().flatten())
            .all(|v| v.is_finite()),
        "non-finite GLM data"
    );
    let x: Vec<f64> = design.iter().flatten().copied().collect();
    let rank = linalg::rank(&x, n, k)?;
    let mean = numeric::sum(y) / n as f64;
    let mut mu: Vec<_> = y.iter().map(|&v| (v + mean) / 2.0).collect();
    let mut eta: Vec<_> = mu.iter().map(|&m| family.link(m)).collect();
    let mut scale = family.scale(y, &mu, n - rank);
    let mut previous = family.deviance(y, &mu, scale);
    ensure!(!previous.is_nan(), "GLM initial deviance is NaN");
    let mut weighted_x = vec![0.0; n * k];
    let mut weighted_y = vec![0.0; n];
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 0..100 {
        for i in 0..n {
            let derivative = family.derivative(mu[i]);
            let weight = (1.0 / (derivative.powi(2) * family.variance(mu[i]))).sqrt();
            let working = eta[i] + derivative * (y[i] - mu[i]) - offset[i];
            ensure!(
                weight.is_finite() && working.is_finite(),
                "GLM weights or working response are non-finite"
            );
            weighted_y[i] = weight * working;
            for j in 0..k {
                weighted_x[i * k + j] = weight * x[i * k + j];
            }
        }
        let params = linalg::least_squares(&weighted_x, &weighted_y, n, k)?;
        let linear = linalg::matvec(&x, &params, n, k);
        for i in 0..n {
            eta[i] = linear[i] + offset[i];
            mu[i] = family.inverse(eta[i]);
        }
        let deviance = family.deviance(y, &mu, scale);
        scale = family.scale(y, &mu, n - rank);
        converged = (deviance - previous).abs() <= 1e-8;
        previous = deviance;
        iterations = iteration + 1;
        if converged {
            break;
        }
    }
    // statsmodels refits the final WLS using its pinv_extended cutoff 1e-15.
    let params = linalg::pinv_solve(&weighted_x, &weighted_y, n, k)?;
    let linear = linalg::matvec(&x, &params, n, k);
    for i in 0..n {
        mu[i] = family.inverse(linear[i] + offset[i]);
    }
    Ok(Fit {
        params,
        deviance: family.deviance(y, &mu, 1.0),
        mu,
        scale,
        converged,
        iterations,
    })
}
