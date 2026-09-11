# ruspladder

An in-progress Rust migration of SplAdder v3.1.1. The intended scope is the
complete `prep`, `build`, and differential `test` workflow, excluding plots.
This is not yet a complete or validated replacement for SplAdder.

The production pipeline must run without Python. Python is used only to generate
and compare reference results. Public tabular and HDF5 results remain compatible;
Python pickle caches are outside the compatibility contract.

See [MIGRATION.md](MIGRATION.md) for requirements, progress, and verification.

On this machine, source worktrees and all dependencies, build caches, input data,
and run outputs live under `/home/wubw/data/ruspladder/`. Run Cargo through
`scripts/cargo.sh` to keep its cache and target directory on that disk.

Upstream source: https://github.com/ratschlab/spladder/tree/v3.1.1
