#!/usr/bin/env python3
"""Inspect effective cgroup-v1 memory limits and CPU affinity inside a Slurm job."""
import json
import os
from pathlib import Path


def snapshot():
    affinity = sorted(os.sched_getaffinity(0))
    entries = {}
    for line in Path("/proc/self/cgroup").read_text().splitlines():
        _, controllers, relative = line.split(":", 2)
        for controller in controllers.split(","):
            entries[controller] = relative
    if "memory" not in entries:
        raise RuntimeError("This machine's benchmark probe currently requires cgroup v1 memory control")
    root = Path("/sys/fs/cgroup/memory")
    current = root / entries["memory"].lstrip("/")
    limits = []
    while current != root:
        limits.append({"path": str(current),
                       "limit": int((current / "memory.limit_in_bytes").read_text()),
                       "peak": int((current / "memory.max_usage_in_bytes").read_text())})
        current = current.parent
    return {"cpu_affinity": affinity, "memory_limits": limits,
            "effective_memory_limit": min(x["limit"] for x in limits),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID")}


if __name__ == "__main__":
    result = snapshot()
    print(json.dumps(result, indent=2))
    assert len(result["cpu_affinity"]) == 4, "benchmark needs exactly 4 available logical CPUs"
    assert result["effective_memory_limit"] == 8 * 1024**3, "8 GiB hard memory limit was not enforced"
