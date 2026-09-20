//! Native kernels used by the pinned NumPy reference (OpenBLAS ILP64).
use anyhow::{Result, ensure};
use std::sync::Once;

unsafe extern "C" {
    fn scipy_openblas_set_num_threads64_(threads: i32);
    fn scipy_dgetrf_64_(m: &i64, n: &i64, a: *mut f64, lda: &i64, pivots: *mut i64, info: &mut i64);
    fn scipy_LAPACKE_dgelsd64_(
        layout: i32,
        m: i64,
        n: i64,
        nrhs: i64,
        a: *mut f64,
        lda: i64,
        b: *mut f64,
        ldb: i64,
        singular: *mut f64,
        rcond: f64,
        rank: *mut i64,
    ) -> i64;
    fn scipy_LAPACKE_dgesdd64_(
        layout: i32,
        job: i8,
        m: i64,
        n: i64,
        a: *mut f64,
        lda: i64,
        singular: *mut f64,
        u: *mut f64,
        ldu: i64,
        vt: *mut f64,
        ldvt: i64,
    ) -> i64;
    fn scipy_cblas_dgemv64_(
        layout: i32,
        trans: i32,
        m: i64,
        n: i64,
        alpha: f64,
        a: *const f64,
        lda: i64,
        x: *const f64,
        incx: i64,
        beta: f64,
        y: *mut f64,
        incy: i64,
    );
    fn scipy_cblas_dgemm64_(
        layout: i32,
        transa: i32,
        transb: i32,
        m: i64,
        n: i64,
        k: i64,
        alpha: f64,
        a: *const f64,
        lda: i64,
        b: *const f64,
        ldb: i64,
        beta: f64,
        c: *mut f64,
        ldc: i64,
    );
}

pub fn least_squares(x: &[f64], y: &[f64], n: usize, k: usize) -> Result<Vec<f64>> {
    assert!(n > 0 && k > 0 && x.len() == n * k && y.len() == n);
    initialize();
    let mut a = x.to_vec();
    let mut b = vec![0.0; n.max(k)];
    b[..n].copy_from_slice(y);
    let mut singular = vec![0.0; n.min(k)];
    let mut rank = 0;
    let info = unsafe {
        scipy_LAPACKE_dgelsd64_(
            101,
            n as i64,
            k as i64,
            1,
            a.as_mut_ptr(),
            k as i64,
            b.as_mut_ptr(),
            1,
            singular.as_mut_ptr(),
            -1.0,
            &mut rank,
        )
    };
    ensure!(info == 0, "LAPACK least squares failed: {info}");
    b.truncate(k);
    Ok(b)
}

struct Svd {
    singular: Vec<f64>,
    u: Vec<f64>,
    vt: Vec<f64>,
}
fn svd(x: &[f64], n: usize, k: usize, vectors: bool) -> Result<Svd> {
    assert!(n > 0 && k > 0 && x.len() == n * k);
    initialize();
    let r = n.min(k);
    let mut a = x.to_vec();
    let mut result = Svd {
        singular: vec![0.0; r],
        u: vec![0.0; n * r],
        vt: vec![0.0; r * k],
    };
    let info = unsafe {
        scipy_LAPACKE_dgesdd64_(
            101,
            if vectors { b'S' } else { b'N' } as i8,
            n as i64,
            k as i64,
            a.as_mut_ptr(),
            k as i64,
            result.singular.as_mut_ptr(),
            result.u.as_mut_ptr(),
            r as i64,
            result.vt.as_mut_ptr(),
            k as i64,
        )
    };
    ensure!(info == 0, "LAPACK SVD failed: {info}");
    Ok(result)
}

pub fn rank(x: &[f64], n: usize, k: usize) -> Result<usize> {
    let s = svd(x, n, k, false)?.singular;
    Ok(s.iter()
        .filter(|&&x| x > s[0] * n.max(k) as f64 * f64::EPSILON)
        .count())
}

pub fn matvec(x: &[f64], y: &[f64], n: usize, k: usize) -> Vec<f64> {
    assert!(n > 0 && k > 0 && x.len() == n * k && y.len() == k);
    initialize();
    let mut result = vec![0.0; n];
    unsafe {
        scipy_cblas_dgemv64_(
            101,
            111,
            n as i64,
            k as i64,
            1.0,
            x.as_ptr(),
            k as i64,
            y.as_ptr(),
            1,
            0.0,
            result.as_mut_ptr(),
            1,
        );
    }
    result
}

pub fn pinv_solve(x: &[f64], y: &[f64], n: usize, k: usize) -> Result<Vec<f64>> {
    let mut s = svd(x, n, k, true)?;
    let r = n.min(k);
    let cutoff = 1e-15 * s.singular[0];
    for i in 0..r {
        let inv = if s.singular[i] > cutoff {
            1.0 / s.singular[i]
        } else {
            0.0
        };
        for j in 0..n {
            s.u[j * r + i] *= inv;
        }
    }
    let mut pinv = vec![0.0; k * n];
    unsafe {
        scipy_cblas_dgemm64_(
            101,
            112,
            112,
            k as i64,
            n as i64,
            r as i64,
            1.0,
            s.vt.as_ptr(),
            k as i64,
            s.u.as_ptr(),
            r as i64,
            0.0,
            pinv.as_mut_ptr(),
            n as i64,
        );
    }
    Ok(matvec(&pinv, y, k, n))
}

fn initialize() {
    static INIT: Once = Once::new();
    // Parallelism belongs to the independent event rows, not these tiny solves.
    INIT.call_once(|| unsafe { scipy_openblas_set_num_threads64_(1) });
}

pub fn determinant(mut column_major: Vec<f64>, n: usize) -> f64 {
    assert!(n > 0 && column_major.len() == n * n);
    initialize();
    let size = n as i64;
    let mut pivots = vec![0_i64; n];
    let mut info = 0;
    unsafe {
        scipy_dgetrf_64_(
            &size,
            &size,
            column_major.as_mut_ptr(),
            &size,
            pivots.as_mut_ptr(),
            &mut info,
        );
    }
    if info > 0 {
        return 0.0;
    }
    assert_eq!(info, 0);
    let mut sign = 1.0;
    let mut logdet = 0.0;
    for i in 0..n {
        if pivots[i] != i as i64 + 1 {
            sign = -sign;
        }
        let d = column_major[i * n + i];
        if d < 0.0 {
            sign = -sign;
        }
        logdet += d.abs().ln();
    }
    sign * logdet.exp()
}
