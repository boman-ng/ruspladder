#!/usr/bin/env python3
"""Compare upstream augmentation; inject identical tracks at its I/O boundary."""
import argparse
import contextlib
import copy
import io
import json
import random
import subprocess
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph
from spladder import editgraph
from compare_annotations import serialize
from compare_detectors import cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    retention = dict(min_retention_cov=5, min_retention_region=0.75,
                     min_retention_max_exon_fold_diff=4, min_retention_rel_cov=0.1,
                     max_retention_rel_cov=1.5)
    cassette = dict(min_cassette_cov=5, min_cassette_region=0.9, min_cassette_rel_diff=0.5)
    options = SimpleNamespace(verbose=False, sparse_bam=False, confidence=3, read_filter=None,
                              var_aware=False, primary_only=False, ignore_mismatches=False,
                              mm_tag="NM", ref_genome=None, intron_retention=retention, cassette_exon=cassette)
    rng = random.Random(461)
    inputs = []
    # Explicit consecutive retention paths plus a cassette between two introns.
    for n in [2, 3, 4, 8, 16]:
        gene = Gene(name="chain", start=10, stop=10 + n * 50, chr="chr1", strand="+")
        gene.exons = [np.array([[10 + i * 50, 30 + i * 50] for i in range(n)])]
        gene.splicegraph = Splicegraph(gene)
        inputs.append(gene)
    for item in cases(args.upstream):
        g = item["graph"]
        n = len(g["vertices"])
        if n == 0: continue
        gene = Gene(name="probe", start=min(v[0] for v in g["vertices"]),
                    stop=max(v[1] for v in g["vertices"]), chr="chr1", strand=item["strand"])
        gene.splicegraph.vertices = np.asarray(g["vertices"]).T
        gene.splicegraph.edges = np.zeros((n, n), dtype=int)
        for i, neighbors in enumerate(g["edges"]): gene.splicegraph.edges[i, neighbors] = 1
        gene.splicegraph.terminals = np.asarray(g["terminals"], dtype=int).T
        inputs.append(gene)
    requests, expected = [], []
    totals = dict(retention=0, cassette=0)
    for gene in inputs:
        vertices = gene.splicegraph.vertices.T.tolist()
        introns = []
        for a, b in zip(vertices, vertices[1:]):
            if b[0] - a[1] > 10:
                introns.extend([[a[1], a[1] + 3], [b[0] - 3, b[0]]])
        for mode in ["retention", "cassette"]:
            for profile in ["uniform", "islands", "random"]:
                track = np.full(gene.stop - gene.start, 10, dtype=np.uint64)
                if profile == "islands":
                    track[:] = 0
                    for a, b in zip(introns[::2], introns[1::2]): track[a[1] - gene.start:b[0] - gene.start] = 20
                elif profile == "random": track[:] = [rng.randrange(20) for _ in track]
                requests.append(dict(gene=serialize(gene), track=track.tolist(), introns=introns,
                                     retention=retention if mode == "retention" else None,
                                     cassette=cassette if mode == "cassette" else None))
                modified = copy.deepcopy(gene)
                modified.introns = [np.array([v + [10] for v in introns], dtype=int).reshape(-1, 3) for _ in range(2)]
                with patch.object(editgraph, "init_regions", return_value=([], options)), \
                     patch.object(editgraph, "add_reads_from_bam", return_value=track[np.newaxis, :]), \
                     contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    function = editgraph.insert_intron_retentions if mode == "retention" else editgraph.insert_cassette_exons
                    genes, count = function(np.array([modified], dtype=object), "injected.bam", options)
                expected.append(dict(gene=serialize(genes[0]), inserted=int(count)))
                totals[mode] += int(count)
    result = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests),
                            text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in result.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want != got:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2))
            raise AssertionError(f"augmentation mismatch: {path}")
    assert all(v > 0 for v in totals.values()), totals
    (args.work / "report.json").write_text(json.dumps(dict(cases=len(actual), inserted=totals, result="pass"), indent=2) + "\n")
    print(f"PASS: {len(actual)} coverage-dependent augmentation comparisons; {totals}")


if __name__ == "__main__": main()
