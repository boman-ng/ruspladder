#!/usr/bin/env bash
set -euo pipefail
task_root="${RUSPLADDER_WORK_ROOT:-/home/wubw/data/ruspladder}"
repo_root=$(cd "$(dirname "$0")/.." && pwd)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=4 RAYON_NUM_THREADS=4
export NUMBA_CACHE_DIR="$task_root/cache/numba"
export TMPDIR="$task_root/tmp"
cd "$repo_root"
scripts/cargo.sh fmt -- --check
scripts/cargo.sh clippy --all-targets -- -D warnings
scripts/cargo.sh build --examples
check_root="$task_root/runs/check-$(date +%Y%m%dT%H%M%S)"
mkdir -p "$check_root"
printf 'Evidence directory: %s\n' "$check_root"
reference_python="$task_root/envs/reference/bin/python"
probe_root="$task_root/target/debug/examples"
upstream="$task_root/upstream/spladder"
"$reference_python" scripts/compare_graphs.py "$probe_root/graph_probe" --failure "$check_root/graph-mismatch.json"
"$reference_python" scripts/compare_annotations.py "$probe_root/annotation_probe" --upstream "$upstream" --work "$check_root/annotation" --extra "$task_root/fixtures/airway/Homo_sapiens.GRCh37.75_subset.gtf"
"$reference_python" scripts/compare_reads.py "$probe_root/reads_probe" --upstream "$upstream" --work "$check_root/reads"
"$reference_python" scripts/compare_detectors.py "$probe_root/detect_probe" --upstream "$upstream" --work "$check_root/detectors"
"$reference_python" scripts/compare_collection.py "$probe_root/collect_probe" --upstream "$upstream" --work "$check_root/collection"
"$reference_python" scripts/compare_edits.py "$probe_root/edit_probe" --upstream "$upstream" --work "$check_root/edits"
"$reference_python" scripts/compare_augmentation.py "$probe_root/augment_probe" --upstream "$upstream" --work "$check_root/augmentation"
"$reference_python" scripts/compare_intron_edges.py "$probe_root/intron_probe" --upstream "$upstream" --work "$check_root/introns"
"$reference_python" scripts/compare_intron_filters.py "$probe_root/intron_filter_probe" --work "$check_root/intron-filters"
"$reference_python" scripts/compare_build_graphs.py "$probe_root/build_graph_probe" --upstream "$upstream" --work "$check_root/build-graphs" --airway "$task_root/fixtures/airway"
"$reference_python" scripts/compare_merge.py "$probe_root/merge_probe" --upstream "$upstream" --work "$check_root/merge"
"$reference_python" scripts/compare_counts.py "$probe_root/count_probe" --upstream "$upstream" --work "$check_root/counts"
"$reference_python" scripts/compare_verification.py "$probe_root/verify_probe" --upstream "$upstream" --work "$check_root/verification"
