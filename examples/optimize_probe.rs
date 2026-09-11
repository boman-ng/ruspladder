use anyhow::Result;
use ruspladder::optimize;
use serde::Deserialize;
use std::io::{self, BufRead};
#[derive(Deserialize)]
struct Request {
    function: String,
    bounds: [f64; 2],
    center: f64,
    tolerance: f64,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let result = optimize::bounded(
            |x| {
                Ok(match r.function.as_str() {
                    "quadratic" => (x - r.center).powi(2),
                    "quartic" => (x - r.center).powi(4),
                    "linear" => x - r.center,
                    "flat" => 2.0,
                    "nan" => f64::NAN,
                    "absolute" => (x - r.center).abs(),
                    _ => anyhow::bail!("unknown objective"),
                })
            },
            r.bounds,
            r.tolerance,
        )?;
        println!("{}", serde_json::to_string(&result)?);
    }
    Ok(())
}
