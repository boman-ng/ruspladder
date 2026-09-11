#!/usr/bin/env bash
# Reuse the local image and locked reference environment in an enforced cgroup.
set -euo pipefail
task_root="${RUSPLADDER_WORK_ROOT:-/home/wubw/data/ruspladder}"
real_root=$(readlink -f "$task_root")
repo_root=$(cd "$(dirname "$0")/.." && pwd)
python_executable=$(readlink -f "$task_root/envs/reference/bin/python")
python_alias=$(readlink "$task_root/envs/reference/bin/python")
python_prefix=$(dirname "$(dirname "$python_executable")")
exec docker run --rm --network none --user "$(id -u):$(id -g)" \
    --cpuset-cpus "${RUSPLADDER_BENCH_CPUS:-0-3}" --cpus 4 \
    --memory 8g --memory-swap 8g \
    --mount "type=bind,source=$real_root,target=$real_root" \
    --mount "type=bind,source=$real_root,target=$task_root" \
    --mount "type=bind,source=$python_prefix,target=$python_prefix,readonly" \
    --mount "type=bind,source=$python_executable,target=$python_alias,readonly" \
    --workdir "$repo_root" \
    --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 \
    --env NUMBA_NUM_THREADS=4 --env RAYON_NUM_THREADS=4 \
    --env "NUMBA_CACHE_DIR=$task_root/cache/numba" --env "TMPDIR=$task_root/tmp" \
    --env "XDG_CACHE_HOME=$task_root/cache" --env "MPLCONFIGDIR=$task_root/cache/matplotlib" \
    --entrypoint "$1" \
    sha256:1c3316753323aa5edd4c2be762323794f4aa9f4335635dc4c7964ee99cb986d2 "${@:2}"
