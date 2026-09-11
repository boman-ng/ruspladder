#!/usr/bin/env python3
"""Differential CIGAR, strand, filtering and BAM/CRAM region evidence checks."""
import argparse
import json
import subprocess
from pathlib import Path

import pysam
from spladder.reads import get_reads, filter_read


def create_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "MT", "LN": 2000}]}
    cigars = ["50M", "20M50N30M", "10S20M2I10M50N20M", "20M5D50N30M",
              "20=5X50N30M", "13M50N37M", "12M50N38M", "20M50N30M",
              "20M50N30M", "20M50N30M", "20M50N30M"]
    with pysam.AlignmentFile(str(path), "wb", header=header) as output:
        for i, cigar in enumerate(cigars):
            read = pysam.AlignedSegment()
            read.query_name = f"read{i}"
            read.reference_id = 0
            read.reference_start = 100 + i * 20
            read.mapping_quality = 60
            read.cigarstring = cigar
            length = sum(n for op, n in read.cigartuples if op in [0, 1, 4, 7, 8])
            read.query_sequence = "A" * length
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            read.flag = 256 if i == 7 else 2048 if i == 8 else 0
            read.set_tag("NM", i % 3)
            read.set_tag("XM", i % 2)
            read.set_tag("XG", i % 3)
            read.set_tag("ZZ", i % 2)
            if i % 3: read.set_tag("XS", "+" if i % 3 == 1 else "-", value_type="A")
            output.write(read)
    pysam.index(str(path))


def reference(request):
    options = dict(filter=None, mapped=True, spliced=True, strand=None, primary_only=False,
                   var_aware=False, no_mm=False, mm_tag="NM")
    options.update(request.get("options", {}))
    coverage, plus, minus = get_reads(request["bam"], request["chromosome"], request["start"], request["stop"],
                                      collapse=True, cram_ref=request.get("reference"), **options)
    count = 0
    with pysam.AlignmentFile(request["bam"], reference_filename=request.get("reference")) as bam:
        if request["chromosome"] != "MT" and bam.get_tid(request["chromosome"]) >= 0:
            count = sum(not filter_read(read, **options) for read in
                        bam.fetch(request["chromosome"], request["start"], request["stop"], until_eof=True))
    return {"coverage": coverage[0].tolist(), "introns_plus": plus.tolist(),
            "introns_minus": minus.tolist(), "read_count": count}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    bam = args.work / "cigar.bam"
    create_bam(bam)
    requests = []
    for start, stop in [(0, 2000), (105, 180), (130, 320), (400, 500)]:
        for strand in [None, "+", "-"]:
            for options in [{}, {"primary_only": True}, {"mapped": False}, {"spliced": False},
                            {"filter": {"intron": 350000, "exon_len": 12, "mismatch": 0, "mincount": 2}},
                            {"filter": {"intron": 350000, "exon_len": 12, "mismatch": 1, "mincount": 2}, "var_aware": True},
                            {"filter": {"intron": 350000, "exon_len": 12, "mismatch": 0, "mincount": 2}, "no_mm": True},
                            {"filter": {"intron": 350000, "exon_len": 12, "mismatch": 0, "mincount": 2}, "mm_tag": "ZZ"}]:
                requests.append({"bam": str(bam), "chromosome": "chr1", "start": start, "stop": stop,
                                 "options": {**options, "strand": strand}})
    for contig in ["MT", "missing"]:
        requests.append({"bam": str(bam), "chromosome": contig, "start": 0, "stop": 200})
    data = args.upstream / "tests/testcase_events/data"
    for suffix in ["bam", "cram"]:
        path = data / f"align/testcase_events_1_sample1.{suffix}"
        reference_path = str(data / "genome.fa") if suffix == "cram" else None
        with pysam.AlignmentFile(str(path), reference_filename=reference_path) as source:
            contig = source.references[0]
            length = source.lengths[0]
        for start, stop in [(0, min(length, 5000)), (100, min(length, 1500))]:
            for filter_ in [None, {"intron": 350000, "exon_len": 13, "mismatch": 0, "mincount": 2}]:
                requests.append({"bam": str(path), "reference": reference_path, "chromosome": contig,
                                 "start": start, "stop": stop, "options": {"filter": filter_, "primary_only": True}})
    result = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests),
                            capture_output=True, text=True, check=True)
    actual = [json.loads(line) for line in result.stdout.splitlines()]
    assert len(actual) == len(requests)
    for i, (request, actual) in enumerate(zip(requests, actual)):
        expected = reference(request)
        if actual != expected:
            failure = args.work / "mismatch.json"
            failure.write_text(json.dumps({"case": i, "request": request, "expected": expected, "actual": actual}, indent=2))
            raise AssertionError(f"read evidence mismatch: {failure}")
    print(f"PASS: {len(requests)} BAM/CRAM region and filtering cases agree exactly")


if __name__ == "__main__":
    main()
