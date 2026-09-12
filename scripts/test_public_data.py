#!/usr/bin/env python3
"""Replay SplAdder v3.1.1's six make-test scenarios without visualization.

Use its public BAM/CRAM and GTF files, run both complete CLIs from fresh
annotations, and reuse the existing scientific-output comparator.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from compare_lifecycle import compare


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--exporter", type=Path, required=True)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True)  # A new directory prevents cache-only passes.
    commit = subprocess.check_output(
        ["git", "-C", str(args.upstream), "rev-parse", "HEAD"], text=True).strip()
    assert commit == "65ceec839b9ff0cf96703c1605ee43667662f410"
    reference = Path(sys.executable).parent / "spladder"
    outputs = ["--output-conf-icgc", "--output-txt", "--output-txt-conf",
               "--output-gff3", "--output-struc", "--output-struc-conf",
               "--output-bed", "--output-conf-bed", "--output-conf-tcga"]
    event_types = "exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip"
    reports = []
    for name in ["basic-pos-merge", "basic-neg-merge", "basic-pos-single",
                 "basic-neg-single", "events-bam", "events-cram"]:
        root = args.work / name
        inputs = root / "inputs"
        inputs.mkdir(parents=True)
        basic = name.startswith("basic")
        single = name.endswith("single")
        source = args.upstream / "tests" / ("testcase_basic" if basic else "testcase_events") / "data"
        if basic:
            strand = name.split("-")[1]
            annotation = source / f"annotation_{strand}.gtf"
            bams = [source / f"align/{strand}_{i}.bam" for i in range(1, 2 if single else 6)]
            ref = None
            readlen = "15"
        else:
            ext = name.split("-")[1]
            annotation = source / "testcase_events_spladder.gtf"
            bams = [source / f"align/testcase_events_1_sample{i}.{ext}" for i in range(1, 21)]
            ref = source / "genome.fa" if ext == "cram" else None
            readlen = "50"
        local_bams = []
        for bam in bams:
            local = inputs / bam.name
            shutil.copy2(bam, local)
            index = Path(str(bam) + (".crai" if bam.suffix == ".cram" else ".bai"))
            shutil.copy2(index, inputs / index.name)
            local_bams.append(local)
        commands = []
        for mode, binary in [("reference", reference), ("rust", args.binary)]:
            local_annotation = inputs / f"{mode}.gtf"
            shutil.copy2(annotation, local_annotation)
            out = root / mode
            flags = ["build", "-a", str(local_annotation), "-b", ",".join(map(str, local_bams)),
                     "-o", str(out), "--parallel", "4", "--readlen", readlen,
                     "--merge-strat", "single" if single else "merge_graphs",
                     "--event-types", event_types, "--extract-ase"]
            flags += ["--output-conf-icgc"] if single else outputs
            if ref:
                flags += ["--reference", str(ref)]
            stages = [("build", flags)]
            if single:
                stages.append(("reuse-sparse", flags + ["--sparse-bam"]))
            if not basic:
                a, b = list(range(1, 11)), list(range(11, 21))
                if name.endswith("cram"):
                    a, b = [10, 8, 2, 1, 7, 6, 5, 3, 9, 4], [20, 13, 17, 11, 12, 19, 15, 14, 16, 18]
                group = lambda indices: ",".join(f"testcase_events_1_sample{i}" for i in indices)
                stages.append(("test", ["test", "-o", str(out), "-a", group(a), "-b", group(b),
                                        "--readlen", "50", "--event-types", "exon_skip",
                                        "--dpsi", "0", "--parallel", "4"]))
            for stage, flags in stages:
                command = [str(binary), *flags]
                with (root / f"{mode}-{stage}.log").open("w") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1200)
                commands.append(command)
        result = compare(root / "reference", root / "rust", args.exporter)
        reports.append(dict(case=name, commands=commands, **result))
        print("PASS:", name, result, flush=True)
    report = dict(reference_commit=commit, threads=4, cases=reports, result="pass")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
