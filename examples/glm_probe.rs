use anyhow::Result;
use ruspladder::glm::{self, Family};
use serde::Deserialize;
use std::io::{self, BufRead};
#[derive(Deserialize)]
struct Request {
    response: Vec<f64>,
    design: Vec<Vec<f64>>,
    offset: Vec<f64>,
    family: Family,
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let r: Request = serde_json::from_str(&line?)?;
        let output = match glm::fit(&r.response, &r.design, &r.offset, r.family) {
            Ok(result) => serde_json::json!({"fit":result}),
            Err(error) => serde_json::json!({"error":error.to_string()}),
        };
        println!("{output}");
    }
    Ok(())
}
