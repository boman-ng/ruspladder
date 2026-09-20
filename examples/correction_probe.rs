use anyhow::Result;
use ruspladder::correction::{self, Correction};
use serde::Deserialize;
use std::io::{self, BufRead};
#[derive(Deserialize)]
struct Request {
    pvalues: Vec<Option<f64>>,
    method: Correction,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let value = match correction::adjust(
            &r.pvalues
                .iter()
                .map(|x| x.unwrap_or(f64::NAN))
                .collect::<Vec<_>>(),
            r.method,
        ) {
            Ok(values) => serde_json::json!({"values":values}),
            Err(error) => serde_json::json!({"error":error.to_string()}),
        };
        println!("{value}");
    }
    Ok(())
}
