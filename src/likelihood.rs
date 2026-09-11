// SplAdder v3.1.1 testing/likelihood.py and spladder_test.py (BSD-3-Clause).

unsafe extern "C" {
    fn ruspladder_lgamma(x: f64) -> f64;
    fn ruspladder_trigamma(x: f64) -> f64;
    fn ruspladder_chi2_cdf1(x: f64) -> f64;
    fn ruspladder_xlog1py(y: f64, x: f64) -> f64;
}

pub fn trigamma(x: f64) -> f64 {
    unsafe { ruspladder_trigamma(x) }
}

pub fn chi2_pvalue(statistic: f64) -> f64 {
    // Preserve 1-CDF, including cancellation at very small p-values.
    1.0 - unsafe { ruspladder_chi2_cdf1(statistic) }
}

pub fn nb_logpmf(y: f64, n: f64, p: f64) -> f64 {
    if n <= 0.0 || !(0.0..=1.0).contains(&p) || y.is_nan() {
        return f64::NAN;
    }
    if y < 0.0 || y.fract() != 0.0 {
        return f64::NEG_INFINITY;
    }
    let lgamma = |x| unsafe { ruspladder_lgamma(x) };
    let coefficient = lgamma(n + y) - lgamma(y + 1.0) - lgamma(n);
    coefficient + n * p.ln() + unsafe { ruspladder_xlog1py(y, -p) }
}

pub fn adjusted(dispersion: f64, design: &[Vec<f64>], response: &[f64], mu: &[f64]) -> f64 {
    let n = 1.0 / dispersion;
    let loglik = response
        .iter()
        .zip(mu)
        .map(|(&y, &mu)| nb_logpmf(y, n, n / (n + mu)))
        .sum::<f64>();
    let k = design[0].len();
    let information: Vec<f64> = (0..k)
        .flat_map(|j| {
            (0..k).map(move |i| {
                design
                    .iter()
                    .zip(mu)
                    .map(|(row, &mu)| (row[i] * (mu / (1.0 + mu * dispersion))) * row[j])
                    .sum::<f64>()
            })
        })
        .collect();
    // NumPy det uses LU and exponentiates the sum of log absolute diagonals.
    let determinant = crate::linalg::determinant(information, k);
    loglik - 0.5 * determinant.ln()
}

pub fn adjusted_shrink(
    dispersion: f64,
    design: &[Vec<f64>],
    response: &[f64],
    mu: &[f64],
    fitted: f64,
    prior: f64,
) -> f64 {
    adjusted(dispersion, design, response, mu)
        - (dispersion.ln() - fitted.ln()).powi(2) / (2.0 * prior.powi(2))
}
