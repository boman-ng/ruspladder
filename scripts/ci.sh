#!/usr/bin/env bash
set -euo pipefail
task_root="${RUSPLADDER_WORK_ROOT:-$HOME/data/ruspladder}"
python3 -c 'import sys; sys.path.insert(0, "scripts"); from resource_probe import snapshot; s = snapshot(); print(s); assert s["effective_memory_limit"] == 4 * 1024**3'
bash scripts/bootstrap.sh
bash scripts/cargo.sh fmt -- --check
bash scripts/cargo.sh clippy --locked --all-targets -- -D warnings
bash scripts/cargo.sh build --locked --release --bin ruspladder --example cache_export_probe
python3 scripts/package_release.py --out "$task_root/dist"
version=$(python3 -c 'import tomllib; print(tomllib.load(open("Cargo.toml", "rb"))["package"]["version"])')
package="$task_root/dist/ruspladder-v$version-linux-amd64"
make test BINARY="$package/ruspladder" TEST_WORK="$task_root/public-tests"
