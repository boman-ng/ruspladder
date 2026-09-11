#!/usr/bin/env python3
"""One isolated sparse-prep measurement inside a 4 CPU / 8 GiB Slurm step."""
import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path


def reference_worker(task):
    from spladder.reads import summarize_chr
    bam, chromosome, options = task
    return summarize_chr(bam, chromosome, options, usetmp=True)


def reference(args):
    import multiprocessing as mp
    from types import SimpleNamespace
    import h5py
    options = SimpleNamespace(verbose=False, primary_only=False, var_aware=False, mm_tag="NM", tmpdir=str(args.work), cram_ref=None)
    with h5py.File(args.work / "summary.hdf5", "w") as output, mp.Pool(4) as pool:
        # This is prep_sparse_bam_full's ordered temporary-file collection;
        # all scientific computation calls the unmodified summarize_chr.
        for _, path in pool.imap(reference_worker, [(str(args.bam), c, options) for c in args.chromosomes.split(",")]):
            with h5py.File(path, "r") as part:
                for key in part: output.create_dataset(key, data=part[key][:], compression="gzip")
            Path(path).unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["reference", "rust"], required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--chromosomes", required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    if args.worker:
        reference(args)
        return
    from resource_probe import snapshot
    before = snapshot()
    assert len(before["cpu_affinity"]) == 4
    assert before["effective_memory_limit"] == 8 * 1024**3
    request = dict(bam=str(args.bam), output=str(args.work / "summary.hdf5"),
                   chromosomes=args.chromosomes.split(","), reference=None,
                   options={}, parallel=4, window=1048576, unstranded=True)
    if args.mode == "reference":
        command = [sys.executable, __file__, *sys.argv[1:], "--worker"]
    else:
        command = [str(args.probe)]
    start = time.monotonic()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    process = subprocess.run(command, input=json.dumps(request) + "\n", text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed = time.monotonic() - start
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    after = snapshot()
    user = usage.ru_utime - usage_before.ru_utime
    system = usage.ru_stime - usage_before.ru_stime
    report = dict(mode=args.mode, wall_seconds=elapsed, user_seconds=user, system_seconds=system,
                  mean_cores=(user + system) / elapsed, utilization_4cpu_percent=100 * (user + system) / (4 * elapsed),
                  max_process_rss_kib=usage.ru_maxrss, input_blocks=usage.ru_inblock, output_blocks=usage.ru_oublock,
                  returncode=process.returncode, resources_before=before, resources_after=after)
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.work / "stderr.log").write_text(process.stderr)
    print(json.dumps(report, indent=2), flush=True)
    if process.returncode: raise SystemExit(process.returncode)


if __name__ == "__main__": main()
