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
creating its default idle worker pool. A launch probe observed one thread with
this setting and 64 without it, before any statistics ran.

`optimize.rs` ports SciPy 1.13.1's bounded Brent implementation and default
termination/evaluation rules. `likelihood.rs` retains the reference's actual
Cox–Reid determinant computation. `statistics.rs` follows raw dispersion,
Gamma trend, empirical prior, shrunken dispersion, LRT and final isoform choice.

Three observed compatibility problems determined the implementation:

1. A zero-heavy four-sample case produces an information matrix with entries
   around 1e-10 and 400. A nalgebra LU substitute changed the likelihood
   objective by about 2.2e-6 for identical nominal inputs. A full dispersion
   comparison failed at 2.185e-6 absolute error. NumPy documents that its
   determinant uses LAPACK DGETRF, so the replacement now reuses that kernel.
   The separate nalgebra least-squares implementation also required replacement
   to reproduce the subsequent iterations. The nalgebra dependency was removed.
2. JSON's default floating-point parser changed some fixture inputs by one ULP.
   Exact float round trips eliminated those changes. Replaying 152 fixed-input
   objectives from the failing dispersion trace then agreed exactly. This
   concerns the comparison transport; public production HDF5 is binary.
3. Remaining 1e-14 likelihood differences affected a later, ill-conditioned
   shrinkage optimization. Term-by-term comparison isolated `xlog1py`.
   Although `_xlogy.pxd` imports libc log1p, the SciPy wheel calls the renamed
   `cephes_log1p` from unity.c (confirmed by disassembly). Reusing this source
   approximation removed the differences, including all 113 objective
   evaluations from a second trace. No tolerances or optimizer settings changed.

The regression includes NB/Gamma families, zero-heavy and large counts,
convergence flags, intermediate dispersions, final p-values and isoform choices,
1/4-thread equality, six-decimal p-values and decisions at 0.01/0.05/0.1. These
are component comparisons, not acceptance of the complete CLI workflow.
Reference all-zero initial-deviance failures are counted separately.

The native dependency bootstrap extracts only the three redistributable shared
libraries from a SHA-256-pinned NumPy wheel and validates each library's hash.
The 16.5 MB archive and 27 MB extracted libraries stay under the configured data
work root. An independent installation was tested under
`../tmp/native-bootstrap-test`. See `licenses/NumPy-wheel-native-notices.txt`
for OpenBLAS, LAPACK, GCC runtime and quadmath notices. Native library paths are
set at build time through RUSPLADDER_BLAS_DIR; no portability beyond this tested
Linux x86_64 baseline is claimed.

Primary references:

- [NumPy determinant/LAPACK contract](https://numpy.org/doc/2.2/reference/generated/numpy.linalg.det.html)
- [statsmodels GLM IRLS](https://github.com/statsmodels/statsmodels/blob/v0.14.4/statsmodels/genmod/generalized_linear_model.py)
- [statsmodels minimal WLS](https://github.com/statsmodels/statsmodels/blob/v0.14.4/statsmodels/regression/_tools.py)
- [SciPy bounded optimizer](https://github.com/scipy/scipy/blob/v1.13.1/scipy/optimize/_optimize.py)
- [SciPy xlogy/xlog1py](https://github.com/scipy/scipy/blob/v1.13.1/scipy/special/_xlogy.pxd)
- [SciPy Cephes unity.c](https://github.com/scipy/scipy/blob/v1.13.1/scipy/special/cephes/unity.c)
- [DLMF one-degree-of-freedom chi-square identity](https://dlmf.nist.gov/8.4.E1)
