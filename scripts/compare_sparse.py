#!/usr/bin/env python3
"""Compare sparse HDF5 arrays/dtypes to summarize_chr, across window boundaries."""
import argparse
import contextlib
import io
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pysam
from spladder.reads import summarize_chr
from compare_reads import create_bam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    bam = args.work / "cigar.bam"
    create_bam(bam)
    requests, expectations, reference_errors = [], [], []
    options_list = [{}, {"primary_only": True}, {"mapped": False}, {"spliced": False},
                    {"filter": dict(intron=350000, exon_len=12, mismatch=0, mincount=2)},
                    {"filter": dict(intron=350000, exon_len=12, mismatch=1, mincount=2), "var_aware": True}]
    source_files = [(bam, None, ["chr1", "MT", "missing"])]
    data = args.upstream / "tests/testcase_events/data"
    for ext in ["bam", "cram"]:
        path = data / f"align/testcase_events_1_sample1.{ext}"
        reference = data / "genome.fa" if ext == "cram" else None
        with pysam.AlignmentFile(str(path), reference_filename=str(reference) if reference else None) as source:
            source_files.append((path, reference, list(source.references)))
    for path, reference, chromosomes in source_files:
        for unstranded in [False, True]:
            for filters in options_list:
                opts = dict(primary_only=False, var_aware=False, mm_tag="NM", **{})
                opts.update(filters)
                source_options = SimpleNamespace(**opts, verbose=False, cram_ref=str(reference) if reference else None)
                expected = {}
                failed = False
                for chromosome in chromosomes:
                    with contextlib.redirect_stdout(io.StringIO()):
                        try:
                            _, coo, minus, plus = summarize_chr(str(path), chromosome, source_options,
                                filter=filters.get("filter"), mapped=filters.get("mapped", True),
                                spliced=filters.get("spliced", True), unstranded=unstranded)
                        except KeyError as error:
                            if str(error) != "'XM'" or not filters.get("var_aware"): raise
                            request = dict(bam=str(path), output=str(args.work / "invalid.hdf5"), chromosomes=[chromosome],
                                           reference=str(reference) if reference else None, options=opts, parallel=4, window=73, unstranded=unstranded)
                            failed_run = subprocess.run([str(args.probe)], input=json.dumps(request) + "\n", text=True, capture_output=True)
                            assert failed_run.returncode != 0 and "XM" in failed_run.stderr, failed_run
                            reference_errors.append(dict(bam=str(path), error=str(error), rust_error="missing XM", unstranded=unstranded))
                            failed = True
                            break
                    expected.update({chromosome + "_reads_row": coo.row.astype("uint8"), chromosome + "_reads_col": coo.col,
                                     chromosome + "_reads_dat": coo.data, chromosome + "_reads_shp": np.array(coo.shape),
                                     chromosome + "_introns_m": minus, chromosome + "_introns_p": plus})
                if failed: continue
                for parallel, window in [(1, 73), (4, 73), (4, 1024)]:
                    output = args.work / f"case-{len(requests)}.hdf5"
                    requests.append(dict(bam=str(path), output=str(output), chromosomes=chromosomes,
                                         reference=str(reference) if reference else None, options=opts,
                                         parallel=parallel, window=window, unstranded=unstranded))
                    expectations.append(expected)
    result = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests),
                            text=True, stdout=subprocess.PIPE, check=True)
    assert len(result.stdout.splitlines()) == len(requests)
    arrays = 0
    for i, (request, expected) in enumerate(zip(requests, expectations)):
        with h5py.File(request["output"], "r") as source:
            assert set(source) == set(expected), (i, set(source), set(expected))
            for key, want in expected.items():
                got = source[key][:]
                assert got.dtype == want.dtype, (i, key, got.dtype, want.dtype)
                assert got.shape == want.shape and np.array_equal(got, want), (i, key, got, want)
                arrays += 1
    report = dict(cases=len(requests), arrays=arrays, threads=[1, 4], windows=[73, 1024], reference_errors=reference_errors, result="pass")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(requests)} sparse BAM/CRAM summaries, {arrays} exact arrays and dtypes")


if __name__ == "__main__": main()
