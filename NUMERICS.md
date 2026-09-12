# Statistical numerical compatibility

The baseline is SplAdder v3.1.1 with NumPy 2.2.6, SciPy 1.13.1 and
statsmodels 0.14.4 on Linux x86_64. Production code is Rust with native
OpenBLAS/LAPACK and small source-derived special-function kernels; it never
imports Python modules or starts a Python interpreter.

`glm.rs` follows statsmodels' NB/log and Gamma/identity IRLS, including the
initial means, working weights, previous-iteration scale, deviance stopping
rule, and final pseudoinverse refit. `linalg.rs` calls the same OpenBLAS 0.3.29
ILP64 kernels used by the pinned NumPy wheel: DGELSD, DGESDD, DGETRF and BLAS
matrix products. Independent rows use Rayon; BLAS computation uses one thread.
Set OPENBLAS_NUM_THREADS=1 before launching to prevent the native loader from
creating its default idle worker pool.

`optimize.rs` ports SciPy 1.13.1's bounded Brent implementation and default
termination/evaluation rules. `likelihood.rs` retains the reference's actual
Cox–Reid determinant computation. `statistics.rs` follows raw dispersion,
Gamma trend, empirical prior, shrunken dispersion, LRT and final isoform choice.

Numerical behavior depends on the reference reduction order: expression and
PSI rows use sequential addition, while contiguous count rows use pairwise
addition. Keep the LAPACK determinant followed by log in the Cox–Reid term;
a log-determinant changes dispersion estimates. SciPy's `xlog1py` uses the
Cephes `log1p` implementation shipped in `vendor/scipy-special`.
JSON probes require exact float round trips; production HDF5 is binary.

The native dependency bootstrap extracts only the three redistributable shared
libraries from a SHA-256-pinned NumPy wheel and validates each library's hash.
The 16.5 MB archive and 27 MB extracted libraries stay under the configured data
work root. See `licenses/NumPy-wheel-native-notices.txt`
for OpenBLAS, LAPACK, GCC runtime and quadmath notices. Development library paths are set through RUSPLADDER_BLAS_DIR. Release packaging
sets relative ELF library paths and includes the pinned native libraries beside
the executable, with their notices. The Linux amd64 archive requires glibc 2.36
or newer; Python is not a runtime dependency.

Primary references:

- [NumPy determinant/LAPACK contract](https://numpy.org/doc/2.2/reference/generated/numpy.linalg.det.html)
- [statsmodels GLM IRLS](https://github.com/statsmodels/statsmodels/blob/v0.14.4/statsmodels/genmod/generalized_linear_model.py)
- [statsmodels minimal WLS](https://github.com/statsmodels/statsmodels/blob/v0.14.4/statsmodels/regression/_tools.py)
- [SciPy bounded optimizer](https://github.com/scipy/scipy/blob/v1.13.1/scipy/optimize/_optimize.py)
- [SciPy xlogy/xlog1py](https://github.com/scipy/scipy/blob/v1.13.1/scipy/special/_xlogy.pxd)
- [SciPy Cephes unity.c](https://github.com/scipy/scipy/blob/v1.13.1/scipy/special/cephes/unity.c)
- [DLMF one-degree-of-freedom chi-square identity](https://dlmf.nist.gov/8.4.E1)
