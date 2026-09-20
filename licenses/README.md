# Third-party redistribution materials

Ruspladder's BSD-3-Clause license does not replace component licenses. The
runtime archive contains this directory, the Rust dependency notices collected
from Cargo's locked sources, and the bundled system-library notices. NumPy,
SciPy, statsmodels and SplAdder attributions also cover source-derived work.

## Pinned numerical libraries

`native-runtime-provenance.json` records the NumPy 2.2.6 wheel URL, SHA-256
and three extracted library hashes consumed by `scripts/bootstrap-native.py`. `NumPy-wheel-native-notices.txt` retains that
wheel's notices: OpenBLAS 0.3.29/LAPACK, GPL-3.0 with GCC Runtime Library
Exception for libgfortran, and LGPL-2.1-or-later for libquadmath.
`LGPL-2.1.txt` is the complete license applicable to libquadmath, copied from
[GCC 4.8.5](https://github.com/gcc-mirror/gcc/blob/releases/gcc-4.8.5/COPYING.LIB).
The GPL text and GCC Runtime Library Exception are retained in the NumPy notices.

[NumPy's pinned build requirements](https://github.com/numpy/numpy/blob/v2.2.6/requirements/ci_requirements.txt)
select scipy-openblas64 0.3.29.0.0. Its
[build recipe](https://github.com/MacPython/openblas-libs/blob/v0.3.29.0.0/tools/build_steps.sh)
and [workflow](https://github.com/MacPython/openblas-libs/blob/v0.3.29.0.0/.github/workflows/posix.yml)
build OpenBLAS v0.3.29 with the ILP64/scipy symbol prefix in manylinux2014.
The workflow explicitly checks out OpenBLAS v0.3.29 before compiling, selecting
[`8795fc7985635de1ecf674b87e2008a15097ffab`](https://github.com/OpenMathLib/OpenBLAS/tree/8795fc7985635de1ecf674b87e2008a15097ffab)
for the numerical library build.

## GCC runtime source correspondence

`native-runtime-provenance.json` identifies the CentOS binary RPMs, their source
RPMs, hashes, original ELF build IDs and matching `.text`, `.rodata`, `.data` and
`.eh_frame` hashes for the two runtime libraries. The raw RPM library SHA-256
prefixes also match the first auditwheel hash in each wheel filename.
The binary RPM's `SOURCERPM` metadata selects:

| Bundled component | Distribution build | Corresponding source |
| --- | --- | --- |
| libquadmath | CentOS gcc 4.8.5-44.el7 | gcc-4.8.5-44.el7.src.rpm |
| libgfortran | CentOS gcc-libraries 8.3.1-2.1.1.el7 | gcc-libraries-8.3.1-2.1.1.el7.src.rpm |

The NumPy wheel changes ELF names/search paths through auditwheel. Ruspladder
further sets relative RPATHs with patchelf; it does not change runtime code or
numeric data. The full source RPMs include upstream source, distribution patches
and RPM build specifications. Source RPMs can be unpacked with `rpm2cpio` and
`cpio`; their specifications describe the distribution build. They are provided
as corresponding-source materials, not as prebuilt installation dependencies.

Packaging downloads and verifies these source RPMs and creates
`ruspladder-v<VERSION>-native-sources.tar.gz` beside the runtime archive, with its
own SHA-256 file. It includes this explanation and the provenance manifest.
Redistribute this source attachment alongside the runtime archive at the same
download location. Build caches and the large source archives are stored outside Git.

## Other native components

HTSlib, HDF5, libdeflate and xz/lzma sources are pinned through Cargo.lock and
retain their nested notices in the runtime package. Bundled zlib and bzip2 are
copied from the pinned build image with their distribution notices; their RPM
versions and source-package names are recorded in the package's build provenance.
The compiler's host runtime and glibc remain supplied by the target operating
system. x86-simd-sort's source revision and license are in `vendor/`; SciPy's
vendored special-function sources have their own source hash manifest.
