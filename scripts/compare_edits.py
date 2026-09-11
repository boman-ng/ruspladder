#!/usr/bin/env python3
"""Differential checks of labels, short exon removal, and duplicate merging."""
import argparse
import contextlib
import copy
import io
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph
from spladder.editgraph import remove_short_exons
from spladder.merge import merge_duplicate_exons
from compare_annotations import serialize
from compare_detectors import cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    remove = dict(terminal_short_extend=40, terminal_short_len=10, min_exon_len=50, min_exon_len_remove=10)
    options = SimpleNamespace(remove_exons=remove, verbose=False, logfile="-")
    inputs = list(cases(args.upstream))
    # Equal endpoints/duplicates exercise the source's positional-index quirk.
    for tx in [[[[10, 40]], [[20, 40]], [[20, 50]]],
               [[[10, 40], [60, 90]], [[20, 40], [60, 90]], [[60, 90]]],
               [[[10, 40]], [[10, 40]], [[10, 40]]]]:
        gene = Gene(name="test", start=10, stop=90, chr="chr1", strand="+")
        gene.exons = [np.array(x) for x in tx]
        gene.splicegraph = Splicegraph(gene)
        inputs.append({"graph": serialize(gene)["splicegraph"], "strand": "+"})
    requests, expected = [], []
    for item in inputs:
        graph = item["graph"]
        n = len(graph["vertices"])
        if n == 0: continue  # upstream remove_short_exons dereferences vertex 0
        gene = Gene(name="probe", start=min(e[0] for e in graph["vertices"]),
                    stop=max(e[1] for e in graph["vertices"]), chr="chr1", strand=item["strand"])
        gene.splicegraph.vertices = np.array(graph["vertices"], dtype=np.int64).T
        gene.splicegraph.edges = np.zeros((n, n), dtype=np.int64)
        for i, neighbors in enumerate(graph["edges"]): gene.splicegraph.edges[i, neighbors] = 1
        gene.splicegraph.terminals = np.array(graph["terminals"], dtype=np.int64).T
        for rm, merge, label in [(False, False, True), (True, False, False), (False, True, False), (True, True, True)]:
            requests.append(dict(gene=serialize(gene), remove=remove if rm else None, merge_duplicates=merge, label=label))
            modified = copy.deepcopy(gene)
            with contextlib.redirect_stdout(io.StringIO()):
                if rm: modified = remove_short_exons(np.array([modified], dtype=object), options)[0]
                if merge: modified = merge_duplicate_exons(np.array([modified], dtype=object), options)[0]
                if label: modified.label_alt()
            expected.append(serialize(modified))
    result = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests),
                            text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in result.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want != got:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2))
            raise AssertionError(f"graph edit mismatch: {path}")
    (args.work / "report.json").write_text(json.dumps(dict(cases=len(actual), result="pass"), indent=2) + "\n")
    print(f"PASS: {len(actual)} gene-label and graph-edit comparisons")


if __name__ == "__main__": main()
