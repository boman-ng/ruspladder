#!/usr/bin/env python3
"""Compare sample/chunk graph merging, edge support filtering and cache roundtrips."""
import argparse
import contextlib
import io
import json
import pickle
import subprocess
from pathlib import Path
from types import SimpleNamespace

from spladder.merge import merge_genes_by_splicegraph
from spladder.editgraph import filter_by_edgecount
from compare_annotations import serialize


def load_genes(path):
    with path.open("rb") as source: genes = pickle.load(source, encoding="latin1")[0]
    for gene in genes:
        gene.from_sparse()
        if not hasattr(gene, "introns_anno"): gene.populate_annotated_introns()
    return genes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    groups = []
    for directory in sorted((args.upstream / "tests").glob("*/results*/spladder")):
        files = sorted(p for p in directory.glob("genes_graph_conf3.*.pickle") if not ".merge_" in p.name)
        if files: groups.append(files)
    requests, expected = [], []
    options = SimpleNamespace(do_merge_all=False, confidence=3, outdir=str(args.work))
    for files in groups:
        for chunk, threshold in [(None, None), (None, 0), (None, 1), (None, 2), (None, len(files) + 1), (3, None), (3, 2)]:
            index = len(requests)
            requests.append(dict(samples=[[serialize(g) for g in load_genes(path)] for path in files],
                                 chunksize=chunk, min_count=threshold, cache=str(args.work / f"cache-{index}.hdf5")))
            with contextlib.redirect_stdout(io.StringIO()):
                inputs = files
                if chunk:
                    inputs = []
                    for i in range(0, len(files), chunk):
                        out = args.work / f"chunk-{index}-{i}.pickle"
                        merge_genes_by_splicegraph(options, merge_list=list(map(str, files[i:i + chunk])), fn_out=str(out))
                        inputs.append(out)
                output = args.work / f"merged-{index}.pickle"
                merge_genes_by_splicegraph(options, merge_list=list(map(str, inputs)), fn_out=str(output))
                genes = load_genes(output)
                if threshold is not None:
                    genes = filter_by_edgecount(genes, SimpleNamespace(sg_min_edge_count=threshold))
                    for gene in genes: gene.from_sparse()
            expected.append([serialize(g) for g in genes])
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want != got:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2))
            raise AssertionError(f"sample graph merge mismatch: {path}")
    (args.work / "report.json").write_text(json.dumps(dict(cases=len(actual), source_groups=len(groups), result="pass"), indent=2) + "\n")
    print(f"PASS: {len(actual)} sample/chunk merge, edge support, and graph cache comparisons")


if __name__ == "__main__": main()
