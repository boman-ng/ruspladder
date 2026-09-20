#!/usr/bin/env bash
# Reuse the local image and locked reference environment in an enforced cgroup.
set -euo pipefail
task_root="${RUSPLADDER_WORK_ROOT:-$HOME/data/ruspladder}"
real_root=$(readlink -f "$task_root")
repo_root=$(cd "$(dirname "$0")/.." && pwd)
python_executable=$(readlink -f "$task_root/envs/reference/bin/python")
python_alias=$(readlink "$task_root/envs/reference/bin/python")
python_prefix=$(dirname "$(dirname "$python_executable")")
mounts=(--mount "type=bind,source=$real_root,target=$task_root")
if [[ "$real_root" != "$task_root" ]]; then
    mounts+=(--mount "type=bind,source=$real_root,target=$real_root")
fi
if [[ "$python_alias" = /* && "$python_alias" != "$python_executable" ]]; then
    mounts+=(--mount "type=bind,source=$python_executable,target=$python_alias,readonly")
fi
exec docker run --rm --network none --user "$(id -u):$(id -g)" \
    --cpuset-cpus "${RUSPLADDER_BENCH_CPUS:-0-3}" --cpus 4 \
    --memory 4g --memory-swap 4g \
    "${mounts[@]}" \
    --mount "type=bind,source=$python_prefix,target=$python_prefix,readonly" \
    --mount "type=bind,source=$repo_root,target=$repo_root,readonly" \
    --env "RUSPLADDER_WORK_ROOT=$task_root" \
    --workdir "$repo_root" \
    --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 \
    --env NUMBA_NUM_THREADS=4 --env RAYON_NUM_THREADS=4 \
    --env "NUMBA_CACHE_DIR=$task_root/cache/numba" --env "TMPDIR=$task_root/tmp" \
    --env "XDG_CACHE_HOME=$task_root/cache" --env "MPLCONFIGDIR=$task_root/cache/matplotlib" \
    --entrypoint "$1" \
    ruspladder-build:0.1.0 "${@:2}"
