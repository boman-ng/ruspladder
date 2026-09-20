#!/usr/bin/env python3
"""Independently verify every original/normalized GTF row and ID placement."""
import argparse
import csv
import hashlib
import itertools
import json
import re
from collections import defaultdict
from pathlib import Path

ID = re.compile(r'\b(gene_id|transcript_id)\s+"([^"\n]*)"')


def audit(original, normalized):
    placements = {kind: defaultdict(dict) for kind in ("gene", "transcript")}
    inverse = {kind: {} for kind in placements}
    rows = changed = 0
    with original.open() as source, normalized.open() as target:
        for line_no, (a, b) in enumerate(itertools.zip_longest(source, target), 1):
            assert a is not None and b is not None, ("line count", line_no)
            if a.startswith("#") or not a.strip():
                assert a == b, ("comment", line_no)
                continue
            x, y = a.rstrip("\n").split("\t"), b.rstrip("\n").split("\t")
            assert len(x) == len(y) == 9 and x[:8] == y[:8], ("coordinates/fields", line_no)
            old, new = dict(ID.findall(x[8])), dict(ID.findall(y[8]))
            assert old.keys() == new.keys(), ("attributes", line_no)
            assert ID.sub(lambda m: m[1] + ' "ID"', x[8]) == ID.sub(lambda m: m[1] + ' "ID"', y[8]), ("other attributes", line_no)
            rows += 1
            changed += a != b
            for kind, key in [("gene", "gene_id"), ("transcript", "transcript_id")]:
                if key not in old:
                    continue
                scope = (x[0], x[6]) if kind == "gene" else (old["gene_id"], x[0], x[6])
                previous = placements[kind][old[key]].setdefault(scope, new[key])
                assert previous == new[key], ("one placement mapped twice", line_no)
                origin = (old[key], scope)
                assert inverse[kind].setdefault(new[key], origin) == origin, ("mixed coordinate frames", line_no)
    report = dict(rows=rows, changed_rows=changed)
    expected_mapping = set()
    for kind in placements:
        collisions = 0
        for original_id, scopes in placements[kind].items():
            if len(scopes) > 1:
                collisions += 1
                assert len(set(scopes.values())) == len(scopes), original_id
            else:
                assert next(iter(scopes.values())) == original_id, ("unnecessary rename", original_id)
            for scope, new_id in scopes.items():
                if new_id == original_id:
                    continue
                parent, chrom, strand = ("", *scope) if kind == "gene" else scope
                expected_mapping.add((kind, original_id, new_id, parent, chrom, strand))
        report[kind] = dict(original_ids=len(placements[kind]), placements=len(inverse[kind]), collision_ids=collisions)
    mapping = Path(str(original) + ".loci.tsv")
    with mapping.open() as f:
        reader = csv.reader(f, delimiter="\t")
        assert next(reader) == ["feature", "original_id", "locus_id", "gene_id", "seqname", "strand"]
        actual = [tuple(row) for row in reader]
    assert len(actual) == len(set(actual)) and set(actual) == expected_mapping, "mapping file"
    for label, path in [("original", original), ("normalized", normalized)]:
        with path.open("rb") as f:
            report[label + "_sha256"] = hashlib.file_digest(f, "sha256").hexdigest()
    return dict(result="pass", **report)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("original", type=Path)
    p.add_argument("normalized", type=Path)
    p.add_argument("--report", type=Path, required=True)
    a = p.parse_args()
    result = audit(a.original, a.normalized)
    a.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
