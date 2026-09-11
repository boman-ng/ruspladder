#!/usr/bin/env bash
set -euo pipefail
task_root="${RUSPLADDER_WORK_ROOT:-/home/wubw/data/ruspladder}"
repo_root=$(cd "$(dirname "$0")/.." && pwd)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=4 RAYON_NUM_THREADS=4
export NUMBA_CACHE_DIR="$task_root/cache/numba"
export TMPDIR="$task_root/tmp"
run_root="$task_root/runs/p0-benchmark-$(date +%Y%m%dT%H%M%S)"
mkdir -p "$run_root"
printf '%s\n' "$run_root"
for run_label in reference-0 rust-0 rust-1 reference-1 reference-2 rust-2; do
    run_mode=${run_label%-*}
    bash "$repo_root/scripts/run_constrained.sh" \
        "$task_root/envs/reference/bin/python" "$repo_root/scripts/benchmark_sparse.py" \
        --mode "$run_mode" --probe "$task_root/target/release/examples/sparse_probe" \
        --bam "$task_root/fixtures/airway/SRR1039508_subset.bam" --chromosomes 1 \
        --work "$run_root/$run_label" > "$run_root/$run_label.log" 2>&1
done
bash "$repo_root/scripts/run_constrained.sh" \
    "$task_root/envs/reference/bin/python" "$repo_root/scripts/benchmark_sparse.py" \
    --mode rust --probe "$task_root/target/release/examples/sparse_probe" \
    --bam "$task_root/runs/p0-coverage-parallel/long-contigs.bam" --chromosomes contig0,contig1,contig2,contig3 \
    --work "$run_root/oom-reproducer-rust" > "$run_root/oom-reproducer-rust.log" 2>&1
printf 'Completed: %s\n' "$run_root"
