// Adapted from SciPy 1.13.1 optimize._optimize._minimize_scalar_bounded.
// BSD-3-Clause; see licenses/SciPy-BSD.txt. Keep its tolerances and update order.
use anyhow::{Result, ensure};
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct Minimum {
    pub x: f64,
    pub fun: f64,
    pub status: u8,
    pub evaluations: usize,
}

pub fn bounded(
    mut function: impl FnMut(f64) -> Result<f64>,
    bounds: [f64; 2],
    tolerance: f64,
) -> Result<Minimum> {
    ensure!(
        bounds.iter().all(|x| x.is_finite()) && bounds[0] <= bounds[1],
        "invalid scalar minimization bounds"
    );
    let sqrt_eps = 2.2e-16_f64.sqrt();
    let golden_mean = 0.5 * (3.0 - 5.0_f64.sqrt());
    let [mut a, mut b] = bounds;
    let mut fulc = a + golden_mean * (b - a);
    let mut nfc = fulc;
    let mut xf = fulc;
    let mut rat: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut fx = function(xf)?;
    let mut evaluations = 1;
    let mut fu = f64::INFINITY;
    let mut ffulc = fx;
    let mut fnfc = fx;
    let mut xm = 0.5 * (a + b);
    let mut tol1 = sqrt_eps * xf.abs() + tolerance / 3.0;
    let mut tol2 = 2.0 * tol1;
    let mut status = 0;
    let sign = |v: f64| {
        if v < 0.0 {
            -1.0
        } else if v >= 0.0 {
            1.0
        } else {
            f64::NAN
        }
    };
    while (xf - xm).abs() > tol2 - 0.5 * (b - a) {
        let mut golden = true;
        if e.abs() > tol1 {
            let r = (xf - nfc) * (fx - ffulc);
            let q = (xf - fulc) * (fx - fnfc);
            let mut p = (xf - fulc) * q - (xf - nfc) * r;
            let mut q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let r = e;
            e = rat;
            if p.abs() < (0.5 * q * r).abs() && p > q * (a - xf) && p < q * (b - xf) {
                golden = false;
                rat = (p + 0.0) / q;
                let x = xf + rat;
                if x - a < tol2 || b - x < tol2 {
                    rat = tol1 * sign(xm - xf);
                }
            }
        }
        if golden {
            e = if xf >= xm { a - xf } else { b - xf };
            rat = golden_mean * e;
        }
        let x = xf + sign(rat) * rat.abs().max(tol1);
        fu = function(x)?;
        evaluations += 1;
        if fu <= fx {
            if x >= xf {
                a = xf;
            } else {
                b = xf;
            }
            fulc = nfc;
            ffulc = fnfc;
            nfc = xf;
            fnfc = fx;
            xf = x;
            fx = fu;
        } else {
            if x < xf {
                a = x;
            } else {
                b = x;
            }
            if fu <= fnfc || nfc == xf {
                fulc = nfc;
                ffulc = fnfc;
                nfc = x;
                fnfc = fu;
            } else if fu <= ffulc || fulc == xf || fulc == nfc {
                fulc = x;
                ffulc = fu;
            }
        }
        xm = 0.5 * (a + b);
        tol1 = sqrt_eps * xf.abs() + tolerance / 3.0;
        tol2 = 2.0 * tol1;
        if evaluations >= 500 {
            status = 1;
            break;
        }
    }
    if xf.is_nan() || fx.is_nan() || fu.is_nan() {
        status = 2;
    }
    Ok(Minimum {
        x: xf,
        fun: fx,
        status,
        evaluations,
    })
}
