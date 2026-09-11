# ruspladder

An in-progress Rust migration of SplAdder v3.1.1. The intended scope is the
complete `prep`, `build`, and differential `test` workflow, excluding plots.
This is not yet a complete or validated replacement for SplAdder.

The production pipeline runs without Python. Python is a build/reference tool.
Public tabular and HDF5 results remain compatible;
Python pickle caches are outside the compatibility contract.

See [MIGRATION.md](MIGRATION.md) for requirements, progress, and verification.

On this machine, source worktrees and all dependencies, build caches, input data,
and run outputs live under `/home/wubw/data/ruspladder/`. Run Cargo through
`scripts/cargo.sh` to keep its cache and target directory on that disk.

The current numerical parity baseline targets Linux x86_64. Run
`python3 scripts/bootstrap-native.py` to install the pinned OpenBLAS 0.3.29
ILP64 libraries under the work directory (or use `scripts/bootstrap.sh` for
the complete reference/build setup). The bootstrap verifies the archive and
each native library; it does not install NumPy as a production dependency.
The binary links these native libraries and needs no Python interpreter.
Keep `OPENBLAS_NUM_THREADS=1` when launching it: event rows use Rayon, and
this also prevents OpenBLAS's loader from creating an idle 64-thread pool.
The library additionally fixes BLAS computation to one thread.

`test` now runs from native graph/event caches and public counted HDF5 files:

```sh
OPENBLAS_NUM_THREADS=1 /home/wubw/data/ruspladder/target/debug/ruspladder test \
  -o results -a sampleA1,sampleA2 -b sampleB1,sampleB2 --parallel 4
```

Annotation `prep` and differential `test` are available. The build command and
alignment preparation are still being wired; this example assumes the counted
inputs and native caches already exist. No raw-alignment-to-result benchmark is
claimed yet.

Upstream source: https://github.com/ratschlab/spladder/tree/v3.1.1
