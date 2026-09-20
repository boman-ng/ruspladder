#!/usr/bin/env python3
"""Exercise intron-edge insertion, including boundary mutation by coverage queries."""
import argparse
import contextlib
import copy
import io
import json
import random
import subprocess
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
    rng = random.Random(1283)
    inputs = []
    for item in cases(args.upstream):
        graph = item["graph"]
        n = len(graph["vertices"])
        if n == 0: continue
        gene = Gene(name="probe", start=min(v[0] for v in graph["vertices"]), stop=max(v[1] for v in graph["vertices"]), chr="chr1", strand=item["strand"])
        gene.splicegraph.vertices = np.asarray(graph["vertices"], dtype=int).T
        gene.splicegraph.edges = np.zeros((n, n), dtype=int)
        for i, neighbors in enumerate(graph["edges"]): gene.splicegraph.edges[i, neighbors] = 1
        gene.splicegraph.terminals = np.asarray(graph["terminals"], dtype=int).T
        introns = [[a[1] + delta, b[0]] for a, b in zip(graph["vertices"], graph["vertices"][1:]) for delta in [-15, 0, 15]]
        introns += [[a[0] + 2, a[1] - 2] for a in graph["vertices"]]
        introns += [[a[1], b[0] + 15] for a, b in zip(graph["vertices"], graph["vertices"][1:])]
        introns = [v for v in introns if 0 <= v[0] < v[1]]
        rng.shuffle(introns)
        inputs.append(([gene], [introns[:12]]))
    # Adjacent gene connections must follow the original neighbor short-circuit.
    for strand in ["+", "-"]:
        genes = []
        for i, exons in enumerate([[[10, 60], [110, 160]], [[200, 260], [310, 360]]]):
            g = Gene(name=f"g{i}", start=exons[0][0], stop=exons[-1][1], chr="chr1", strand=strand)
            g.exons = [np.array(exons)]
            g.splicegraph = Splicegraph(g)
            genes.append(g)
        inputs.append((genes, [[[160, 200], [60, 110]], [[160, 200], [260, 310]]]))
    requests, expected = [], []
    totals = {}
    for genes, introns in inputs:
        for depth, append, retain in [(0, False, False), (10, True, True), (11, True, True), (50, False, True)]:
            config = dict(min_exon_len=10, vicinity_region=5, insert_intron_retention=retain,
                          gene_merges=False, append_new_terminal_exons=append, append_new_terminal_exons_len=30)
            request = dict(genes=[serialize(g) for g in genes], introns=introns, options=config, coverage_depth=depth)
            requests.append(request)
            modified = np.array(copy.deepcopy(genes), dtype=object)
            for g, pairs in zip(modified, introns):
                g.introns = np.array([np.array([v + [10] for v in pairs], dtype=int).reshape(-1, 3) for _ in range(2)])
            options = SimpleNamespace(verbose=False, sparse_bam=False, logfile="-", debug=False, intron_edges=config,
                                      read_filter=None, var_aware=False, primary_only=False, ignore_mismatches=False,
                                      mm_tag="NM", ref_genome=None)
            def coverage(blocks, *a, **kw): return np.full((1, blocks[0].stop - blocks[0].start), depth, dtype=np.uint64)
            with patch.object(editgraph, "add_reads_from_bam", side_effect=coverage), contextlib.redirect_stdout(io.StringIO()):
                result, counts = editgraph.insert_intron_edges(modified, "injected.bam", options)
            expected.append(dict(genes=[serialize(g) for g in result], inserted=counts))
            for k, v in counts.items(): totals[k] = totals.get(k, 0) + v
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want != got:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2))
            raise AssertionError(f"intron insertion mismatch: {path}")
    assert all(totals[k] > 0 for k in ["intron_in_exon", "alt_53_prime", "exon_skip", "new_terminal_exon"]), totals
    (args.work / "report.json").write_text(json.dumps(dict(cases=len(actual), inserted=totals, result="pass"), indent=2) + "\n")
    print(f"PASS: {len(actual)} intron-edge insertion cases; {totals}")


if __name__ == "__main__": main()
