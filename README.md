# ruspladder

A Rust migration of SplAdder v3.1.1: `prep`, `build`, and differential `test`,
including all six event types and nonvisual outputs. The production workflow
is compared against the pinned upstream source; see the validation boundaries
and observed upstream failures in [MIGRATION.md](MIGRATION.md) and
[COMPATIBILITY.md](COMPATIBILITY.md).

The production pipeline runs without Python. Python is a build/reference tool.
Public tabular and HDF5 results remain compatible;
Python pickle caches are outside the compatibility contract.

See [PERFORMANCE.md](PERFORMANCE.md) for the 4-thread / 8 GiB measurements.

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

Build the release binary:

```sh
bash scripts/cargo.sh build --release --bin ruspladder
```

`test` runs from native graph/event caches and public counted HDF5 files:

```sh
OPENBLAS_NUM_THREADS=1 /home/wubw/data/ruspladder/target/release/ruspladder test \
  -o results -a sampleA1,sampleA2 -b sampleB1,sampleB2 --parallel 4
```

`build` now connects annotation, BAM/CRAM graph generation, merging, graph counts,
event verification and all nonvisual outputs. For example:

```sh
OPENBLAS_NUM_THREADS=1 /home/wubw/data/ruspladder/target/release/ruspladder build \
  -a annotation.gtf -b sampleA.bam,sampleB.bam -o results --parallel 4
```

`prep`, direct or sparse-alignment `build`, and differential `test` are available.
Use `prep -a annotation.gtf -b sample.bam --sparse-bam --parallel 4` to create
bounded-memory public alignment summaries, then add `--sparse-bam` to `build`.
All commands accept up to 64 threads; validation uses 1/4 threads.
The constrained benchmark and output comparisons are reproducible with
`scripts/run_lifecycle_benchmarks.py`; run `bash scripts/check.sh` for the
complete source comparison suite in the prepared reference environment.

Upstream source: https://github.com/ratschlab/spladder/tree/v3.1.1
