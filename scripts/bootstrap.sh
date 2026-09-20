#!/usr/bin/env bash
set -euo pipefail
task_root="${RUSPLADDER_WORK_ROOT:-$HOME/data/ruspladder}"
repo_root=$(cd "$(dirname "$0")/.." && pwd)
export UV_CACHE_DIR="$task_root/cache/uv"
export TMPDIR="$task_root/tmp"
mkdir -p "$TMPDIR" "$task_root/envs" "$task_root/upstream" "$task_root/cache"
reference="$task_root/upstream/spladder"
if [[ ! -d "$reference/.git" ]]; then
    git clone --depth 1 --branch v3.1.1 https://github.com/ratschlab/spladder.git "$reference"
fi
test "$(git -C "$reference" rev-parse HEAD)" = 65ceec839b9ff0cf96703c1605ee43667662f410
if [[ ! -x "$task_root/envs/reference/bin/python" ]]; then
    uv venv "$task_root/envs/reference" --python 3.12
fi
uv pip install --python "$task_root/envs/reference/bin/python" -r "$repo_root/reference-requirements.lock"
uv pip install --python "$task_root/envs/reference/bin/python" --no-deps -e "$reference"
if [[ ! -x "$task_root/envs/build/bin/python" ]]; then
    uv venv "$task_root/envs/build" --python 3.13
fi
uv pip install --python "$task_root/envs/build/bin/python" cmake==4.4.3 libclang==18.1.1
"$task_root/envs/build/bin/python" "$repo_root/scripts/bootstrap-native.py"
