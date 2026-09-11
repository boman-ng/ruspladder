#!/usr/bin/env python3
"""Compare Rust graph operations with the actual pinned SplAdder classes.

Run with the reference Python after building `--example graph_probe`.
The JSON probe is a development-only interface, not a production file format.
"""
import argparse
import json
import random
import subprocess
from pathlib import Path

import numpy as np
from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph
from spladder.classes.segmentgraph import Segmentgraph


def reference(request):
    gene = Gene(name="test", start=0, stop=100000, chr="test", strand="+")
    gene.transcripts = [str(i) for i in range(len(request["transcripts"]))]
    gene.exons = [np.asarray(t, dtype=int).reshape(-1, 2) for t in request["transcripts"]]
    graph = Splicegraph(gene)
    for op in request.get("operations", []):
        kind = op["op"]
        if kind in ("sort", "uniquify", "update_terminals"):
            getattr(graph, kind)()
        elif kind == "subset":
            graph.subset(np.array(op["indices"], dtype=int))
        elif kind == "intron":
            graph.add_intron(np.array(op["left"], dtype=int), op["keep_end"],
                             np.array(op["right"], dtype=int), op["keep_start"])
        elif kind == "cassette":
            graph.add_cassette_exon(np.array(op["exon"], dtype=int),
                                   np.array(op["before"], dtype=int),
                                   np.array(op["after"], dtype=int))
        elif kind == "retention":
            graph.add_intron_retention(op["left"], op["right"])
        else:
            raise ValueError(kind)
    gene.splicegraph = graph
    gene.segmentgraph = Segmentgraph(gene)
    segments = gene.segmentgraph
    return {
        "graph": {
            "vertices": graph.vertices.T.tolist(),
            "edges": [np.where(row)[0].tolist() for row in graph.edges],
            "terminals": graph.terminals.T.astype(bool).tolist(),
        },
        "segments": {
            "segments": segments.segments.T.tolist(),
            "matches": [np.where(row)[0].tolist() for row in segments.seg_match],
            "edges": np.array(np.where(segments.seg_edges)).T.tolist(),
        },
        "non_alt": gene.get_non_alt_seg_ids().tolist(),
    }


def cases():
    yield {"transcripts": []}
    yield {"transcripts": [[[10, 20]], [[10, 20]], [[10, 30]]]}
    base = [[[0, 20], [40, 60], [80, 100]], [[0, 20], [80, 100]]]
    yield {"transcripts": base}
    for operations in [
        [{"op": "retention", "left": 0, "right": 1}],
        [{"op": "cassette", "exon": [65, 75], "before": [1], "after": [2]}, {"op": "sort"}],
        [{"op": "subset", "indices": [0, 2]}],
        [{"op": "intron", "left": [0], "right": [2], "keep_end": True, "keep_start": True}],
        [{"op": "uniquify"}, {"op": "update_terminals"}],
    ]:
        yield {"transcripts": base, "operations": operations}
    rng = random.Random(23)
    for _ in range(300):
        exons = [[i * 10, i * 10 + rng.randint(1, 25)] for i in range(rng.randint(2, 50))]
        txs = [sorted(rng.sample(exons, rng.randint(1, min(8, len(exons)))))
               for _ in range(rng.randint(1, 16))]
        yield {"transcripts": txs}
        yield {"transcripts": txs, "operations": [{"op": "uniquify"}, {"op": "update_terminals"}]}
    # Equal starts exercise the reference's argsort ordering, not just its sets.
    for count in [3, 15, 16, 17, 32, 64]:
        yield {"transcripts": [[[0, 100 + i], [200, 300]] for i in range(count)]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--failure", type=Path, required=True)
    args = parser.parse_args()
    inputs = list(cases())
    process = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in inputs),
                             text=True, capture_output=True, check=True)
    actual = [json.loads(line) for line in process.stdout.splitlines()]
    assert len(actual) == len(inputs), (len(actual), len(inputs))
    for index, (request, result) in enumerate(zip(inputs, actual)):
        expected = reference(request)
        if expected != result:
            args.failure.write_text(json.dumps({"index": index, "request": request,
                                                "expected": expected, "actual": result}, indent=2))
            raise AssertionError(f"graph mismatch in case {index}; see {args.failure}")
    print(f"PASS: {len(inputs)} graph/segment operation cases agree exactly with SplAdder")


if __name__ == "__main__":
    main()
