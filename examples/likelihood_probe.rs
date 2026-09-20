use anyhow::Result;
use ruspladder::likelihood;
use serde::Deserialize;
use std::io::{self, BufRead};
#[derive(Deserialize)]
struct Request {
    design: Vec<Vec<f64>>,
    response: Vec<f64>,
    mu: Vec<f64>,
    dispersion: f64,
    fitted: f64,
    prior: f64,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let n = 1.0 / r.dispersion;
        let logpmf: Vec<_> = r
            .response
            .iter()
            .zip(&r.mu)
            .map(|(&y, &mu)| likelihood::nb_logpmf(y, n, n / (n + mu)))
            .collect();
        println!(
            "{}",
            serde_json::json!({"logpmf":logpmf,
            "adjusted":likelihood::adjusted(r.dispersion,&r.design,&r.response,&r.mu),
            "shrink":likelihood::adjusted_shrink(r.dispersion,&r.design,&r.response,&r.mu,r.fitted,r.prior),
            "trigamma":likelihood::trigamma(r.prior), "pvalue":likelihood::chi2_pvalue(r.prior)})
        );
    }
    Ok(())
}
