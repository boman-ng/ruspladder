#!/usr/bin/env python3
"""Generate tiny valid BAMs with long contigs and run the unmodified reference."""
import argparse
import json
import multiprocessing
import resource
import time
from pathlib import Path
from types import SimpleNamespace

import pysam
from spladder.reads import summarize_chr


def summarize(task):
    bam, name = task
    options = SimpleNamespace(verbose=False, primary_only=True, var_aware=False, mm_tag="NM")
    start = time.monotonic()
    _, reads, minus, plus = summarize_chr(str(bam), name, options)
    return {"contig": name, "shape": list(reads.shape), "nonzeros": int(reads.nnz),
            "coverage_sum": int(reads.sum()), "introns": int(len(minus) + len(plus)),
            "wall_seconds": time.monotonic() - start,
            "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--length", type=int, default=250000000)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    bam = args.work / "long-contigs.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": f"contig{i}", "LN": args.length} for i in range(4)]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
        for i in range(4):
            for position in [1000, args.length - 1000]:
                record = pysam.AlignedSegment()
                record.query_name = f"read_{i}_{position}"
                record.query_sequence = "A" * 100
                record.query_qualities = pysam.qualitystring_to_array("I" * 100)
                record.reference_id = i
                record.reference_start = position
                record.mapping_quality = 60
                record.cigar = [(0, 50), (3, 100), (0, 50)]
                record.set_tag("NM", 0)
                output.write(record)
    pysam.index(str(bam))
    print(json.dumps({"input_bam_bytes": bam.stat().st_size, "contig_length": args.length,
                      "contigs": 4, "workers": args.workers}), flush=True)
    start = time.monotonic()
    with multiprocessing.Pool(args.workers) as pool:
        result = pool.map(summarize, [(bam, f"contig{i}") for i in range(4)])
    memory_group = next(line.split(":", 2)[2] for line in Path("/proc/self/cgroup").read_text().splitlines()
                        if "memory" in line.split(":", 2)[1].split(","))
    peak = int((Path("/sys/fs/cgroup/memory") / memory_group.lstrip("/") / "memory.max_usage_in_bytes").read_text())
    print(json.dumps({"wall_seconds": time.monotonic() - start, "workers": result,
                      "cgroup_peak_bytes": peak}, indent=2))


if __name__ == "__main__":
    main()
