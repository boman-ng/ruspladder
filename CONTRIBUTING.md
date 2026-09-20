# Contributing

Ruspladder is currently maintained by [boman-ng](https://github.com/boman-ng).
Report reproducible bugs or propose changes through
[Issues](https://github.com/boman-ng/ruspladder/issues) and pull requests.
For security issues use [SECURITY.md](SECURITY.md).

## Build and check

Use the pinned Docker build in [README.md](README.md#build-and-test-from-source)
with a fresh work directory. It runs formatting, Clippy, Rust unit tests, packaging,
the six public SplAdder comparison scenarios and release-evidence generation.
GitHub CI additionally runs the relocated package without Python or network access.

For local development, follow the same README prerequisites, then run:

```sh
export RUSPLADDER_WORK_ROOT="$HOME/data/ruspladder"
bash scripts/bootstrap.sh
python3 -m unittest discover -s scripts -p 'test_readiness.py'
bash scripts/cargo.sh fmt -- --check
bash scripts/cargo.sh clippy --locked --all-targets -- -D warnings
make test
```

The broader comparison suite is useful when changing analysis behavior:

```sh
python3 scripts/fetch_fixtures.py
bash scripts/check.sh
```

Public fixtures are downloaded outside the repository and checked against the
committed manifests. A checksum error is a failure to investigate, not a reason
to regenerate the expected manifest. Use a new work directory for independent
comparisons; stale result caches can otherwise hide computation changes.

## Scientific changes

[COMPATIBILITY.md](COMPATIBILITY.md) defines the supported scientific contracts,
known upstream behaviors and numerical tolerances. [NUMERICS.md](NUMERICS.md)
explains the pinned numerical kernels. A pull request should state the concrete
problem, resulting behavior, affected contracts, and checks actually performed.
Use the narrowest meaningful regression case; compare against the pinned SplAdder
reference when outputs or algorithms change. Do not weaken tolerances merely to
make a comparison pass. Intentional scientific differences need explicit rationale,
evidence and compatibility documentation.

Changes to native dependencies must update provenance, hashes and redistribution
materials, then pass package and scientific comparisons. See
[third-party materials](licenses/README.md). Preserve upstream copyright and
attribution. Contributions are distributed under the repository's BSD-3-Clause
license, with third-party components retaining their own terms.

## Reports and releases

Include the version, platform, exact command, error and a minimal public or
synthetic reproduction. Do not attach patient data, credentials or private
alignments. Mention whether annotations and result directories were fresh.

Release builds produce runtime and native-source archives, checksums and a
verification JSON containing the public comparison results and artifact hashes.
See [redistribution materials](licenses/README.md) for the bundled dependencies.
