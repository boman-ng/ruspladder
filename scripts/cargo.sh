#!/usr/bin/env bash
set -euo pipefail
task_root="${RUSPLADDER_WORK_ROOT:-/home/wubw/data/ruspladder}"
export CARGO_HOME="$task_root/cache/cargo"
export CARGO_TARGET_DIR="$task_root/target"
export TMPDIR="$task_root/tmp"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-16}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-4}"
export RUSPLADDER_BLAS_DIR="${RUSPLADDER_BLAS_DIR:-$task_root/native/openblas/lib}"
if [[ -d "$task_root/envs/build/bin" ]]; then
    export PATH="$task_root/envs/build/bin:$PATH"
    if [[ -z "${LIBCLANG_PATH:-}" ]]; then
        LIBCLANG_PATH=$("$task_root/envs/build/bin/python" -c 'from clang.cindex import Config; print(Config.library_path)')
        export LIBCLANG_PATH
        # The libclang wheel has no compiler resource headers. Reuse this
        # machine's C compiler headers for bindgen's stddef/stdarg includes.
        export BINDGEN_EXTRA_CLANG_ARGS="${BINDGEN_EXTRA_CLANG_ARGS:-} -isystem $(cc -print-file-name=include)"
    fi
fi
exec cargo "$@"
