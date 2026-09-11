#!/usr/bin/env python3
"""Exercise all six raw detectors against the actual upstream implementations."""
import argparse
import contextlib
import io
import json
import pickle
import random
import subprocess
from pathlib import Path

import numpy as np
from spladder.alt_splice import detect
from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph


def reference(request):
    gene = Gene(name="probe", start=0, stop=100000, chr="chr1", strand=request["strand"])
    graph = request["graph"]
    gene.splicegraph = Splicegraph()
    gene.splicegraph.vertices = np.asarray(graph["vertices"], dtype=np.int64).reshape(-1, 2).T
    n = len(graph["vertices"])
    gene.splicegraph.edges = np.zeros((n, n), dtype=np.int64)
    gene.splicegraph.terminals = np.asarray(graph["terminals"], dtype=np.int64).reshape(-1, 2).T
    for i, neighbors in enumerate(graph["edges"]):
        gene.splicegraph.edges[i, neighbors] = 1
    gene.to_sparse()
    genes = np.array([gene], dtype=object)
    indices = np.array([0], dtype=np.int64)
    result = {}
    with contextlib.redirect_stdout(io.StringIO()):
        for kind, function in [("exon_skip", detect.detect_exonskips), ("intron_retention", detect.detect_intronreten),
                               ("mutex_exons", detect.detect_xorexons), ("mult_exon_skip", detect.detect_multipleskips)]:
            ids, found = function(genes, indices, edge_limit=request["edge_limit"])
            assert len(ids) == len(found) and all(i == 0 for i in ids)
            if kind == "mult_exon_skip":
                result[kind] = [{"first": int(a), "skipped": skipped.tolist(), "last": int(b)} for a, skipped, b in found]
            else:
                result[kind] = [list(map(int, row)) for row in found]
        ids5, five, ids3, three = detect.detect_altprime(genes, indices, edge_limit=request["edge_limit"])
        assert all(i == 0 for i in list(ids5) + list(ids3))
        result["alt_5prime"] = [{"common": int(g["threeprimesite"]), "alternatives": g["fiveprimesites"].tolist()} for g in five]
        result["alt_3prime"] = [{"common": int(g["fiveprimesite"]), "alternatives": g["threeprimesites"].tolist()} for g in three]
    return result


def graph(vertices, pairs):
    neighbors = [set() for _ in vertices]
    for a, b in pairs:
        neighbors[a].add(b)
        neighbors[b].add(a)
    return {"vertices": vertices, "edges": [sorted(row) for row in neighbors],
            "terminals": [[not any(j <= i for j in row), not any(j >= i for j in row)] for i, row in enumerate(neighbors)]}


def cases(upstream):
    rng = random.Random(23)
    for n in [0, 1, 2, 3, 4, 8, 12]:
        for _ in range(12):
            vertices = [[10 * i, 10 * i + rng.randrange(1, 51)] for i in range(n)]
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.35]
            for strand in ["+", "-"]:
                yield {"graph": graph(vertices, pairs), "strand": strand, "edge_limit": 500}
    # Shortest-path and longest-path tie selection, and each limit's semantics.
    for n in [4, 6, 8]:
        vertices = [[20 * i, 20 * i + 10] for i in range(n)]
        for pairs in [[(i, i + 1) for i in range(n - 1)] + [(0, n - 1)],
                      [(i, j) for i in range(n) for j in range(i + 1, n)]]:
            for limit in [0, 1, n - 1, n, 500]:
                yield {"graph": graph(vertices, pairs), "strand": "+", "edge_limit": limit}
    for path in sorted((upstream / "tests").glob("*/results*/spladder/genes_graph_conf3.merge_graphs.pickle")):
        with path.open("rb") as source:
            genes = pickle.load(source, encoding="latin1")[0]
        for gene in genes:
            gene.from_sparse()
            sg = gene.splicegraph
            yield {"graph": {"vertices": sg.vertices.T.tolist(),
                             "edges": [np.where(row)[0].tolist() for row in sg.edges],
                             "terminals": sg.terminals.T.astype(bool).tolist()},
                   "strand": gene.strand, "edge_limit": 500}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    requests = list(cases(args.upstream))
    completed = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests),
                               capture_output=True, text=True, check=True)
    results = [json.loads(line) for line in completed.stdout.splitlines()]
    assert len(results) == len(requests)
    totals = {kind: 0 for kind in results[0]}
    reference_errors = []
    for index, (request, actual) in enumerate(zip(requests, results)):
        try:
            expected = reference(request)
        except ValueError as error:
            # Empty genes are removed by the annotation parser. Calling the
            # raw upstream multiple-skip wrapper on an empty graph is undefined
            # (hstack([])). Record that failure, never count it as parity.
            if request["graph"]["vertices"] or str(error) != "need at least one array to concatenate":
                raise
            reference_errors.append({"case": index, "error": str(error), "vertices": 0})
            assert not any(actual.values())
            continue
        if actual != expected:
            failure = args.work / "mismatch.json"
            failure.write_text(json.dumps({"case": index, "request": request, "expected": expected, "actual": actual}, indent=2))
            raise AssertionError(f"raw detector mismatch: {failure}")
        for kind, found in actual.items(): totals[kind] += len(found)
    assert all(count > 0 for count in totals.values()), totals
    compared = len(requests) - len(reference_errors)
    report = {"cases_compared": compared, "reference_errors": reference_errors,
              "event_totals": totals, "result": "pass for comparable cases"}
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {compared} raw detector cases; {len(reference_errors)} empty-graph reference failures recorded separately; event classes: {totals}")


if __name__ == "__main__":
    main()
