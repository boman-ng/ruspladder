use anyhow::Result;
use ruspladder::{
    likelihood,
    statistics::{self, TestOptions},
};
use serde::Deserialize;
use std::io::{self, BufRead};
#[derive(Deserialize)]
struct Request {
    counts: Vec<Vec<f64>>,
    null: Vec<Vec<f64>>,
    design: Vec<Vec<f64>>,
    sf: Vec<f64>,
    selected: Vec<bool>,
    options: TestOptions,
    parallel: usize,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(r.parallel)
            .build()?;
        let result = pool.install(|| -> Result<_> {
            let mut raw = statistics::estimate(&r.counts, &r.design, &r.sf, r.options.min_count)?;
            let estimated = raw.values.clone();
            let trend = statistics::fit_trend(&r.counts, &mut raw, &r.sf)?;
            let (adjusted, prior) = statistics::adjust(&r.counts, &r.design, &raw, &trend, &r.sf)?;
            let pvalues = statistics::test(&r.counts, &adjusted.values, &r.sf, &r.null, &r.design, &r.selected, r.options.max_zero_fraction)?;
            let final_result = statistics::run(&r.counts, &r.null, &r.design, &r.sf, &r.selected, r.options)?;
            Ok(serde_json::json!({"estimated":estimated,"raw":raw,"trend":trend,"adjusted":adjusted,"prior":prior,"pvalues":pvalues,"final":final_result,
                "trigamma":likelihood::trigamma((r.design.len() - r.design[0].len()) as f64 / 2.0)}))
        })?;
        println!("{result}");
    }
    Ok(())
}
