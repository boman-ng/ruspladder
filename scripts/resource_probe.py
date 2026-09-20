#!/usr/bin/env python3
"""Inspect cgroup memory limits and CPU affinity in Slurm/Docker."""
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
    if "" in entries:
        for line in Path("/proc/self/mountinfo").read_text().splitlines():
            fields, filesystem = line.split(" - ", 1)
            before, after = fields.split(), filesystem.split()
            if after[0] == "cgroup2":
                root = Path(before[4])
                current = root / Path(entries[""]).relative_to(before[3])
                break
        else:
            raise RuntimeError("unified cgroup mount not found")
        limits = []
        while True:
            limit = (current / "memory.max").read_text().strip()
            if limit != "max":
                limits.append({"path": str(current), "limit": int(limit),
                               "peak": int((current / "memory.peak").read_text())})
            if current == root:
                break
            current = current.parent
        return {"cpu_affinity": affinity, "memory_limits": limits,
                "effective_memory_limit": min(x["limit"] for x in limits),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID")}
    if "memory" not in entries:
        raise RuntimeError("This machine's benchmark probe currently requires cgroup v1 memory control")
    # Docker mounts the container's cgroup as the filesystem root; Slurm
    # exposes the complete hierarchy. Resolve the actual mount root first.
    for line in Path("/proc/self/mountinfo").read_text().splitlines():
        fields, filesystem = line.split(" - ", 1)
        before, after = fields.split(), filesystem.split()
        if after[0] == "cgroup" and "memory" in after[2].split(","):
            root = Path(before[4])
            current = root / Path(entries["memory"]).relative_to(before[3])
            break
    else:
        raise RuntimeError("memory cgroup mount not found")
    limits = []
    while True:
        limits.append({"path": str(current),
                       "limit": int((current / "memory.limit_in_bytes").read_text()),
                       "peak": int((current / "memory.max_usage_in_bytes").read_text())})
        if current == root:
            break
        current = current.parent
    return {"cpu_affinity": affinity, "memory_limits": limits,
            "effective_memory_limit": min(x["limit"] for x in limits),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID")}


if __name__ == "__main__":
    result = snapshot()
    print(json.dumps(result, indent=2))
    assert len(result["cpu_affinity"]) == 4, "benchmark needs exactly 4 available logical CPUs"
    assert result["effective_memory_limit"] == 4 * 1024**3, "4 GiB hard memory limit was not enforced"
