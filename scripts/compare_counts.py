#!/usr/bin/env python3
"""Compare segment means, covered positions, junction counts and gene mutations."""
import argparse
import contextlib
import copy
import io
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from spladder.count import count_graph_coverage
from compare_annotations import serialize
from compare_merge import load_genes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    cases = []
    basic = args.upstream / "tests/testcase_basic"
    for strand in ["pos", "neg"]:
        genes = load_genes(basic / f"results_merged_{strand}/spladder/genes_graph_conf3.merge_graphs.pickle")
        cases.append((genes, [basic / f"data/align/{strand}_{i}.bam" for i in [1, 2]], basic / f"data/genome_{strand}.fa"))
    events = args.upstream / "tests/testcase_events"
    for extension, directory in [("bam", "results_merged"), ("cram", "results_merged_cram")]:
        genes = load_genes(events / f"{directory}/spladder/genes_graph_conf3.merge_graphs.pickle")
        cases.append((genes, [events / f"data/align/testcase_events_1_sample{i}.{extension}" for i in [1, 2]], events / "data/genome.fa"))
    for chromosome in ["MT", "missing"]:
        genes, bams, reference = copy.deepcopy(cases[-2])
        for gene in genes: gene.chr = chromosome
        cases.append((genes, bams, reference))
    requests, expected = [], []
    for genes, bams, reference in cases:
        for primary in [False, True]:
            options = SimpleNamespace(sparse_bam=False, var_aware=False, primary_only=primary, mm_tag="NM", ref_genome=str(reference))
            modified = copy.deepcopy(genes)
            with contextlib.redirect_stdout(io.StringIO()):
                counts = count_graph_coverage(modified, list(map(str, bams)), options)
            found = [[dict(segments=c.segments.tolist(), seg_pos=c.seg_pos.tolist(), edges=c.edges.tolist(),
                           edges_float=c.edges.dtype.kind == "f") for c in sample] for sample in counts]
            for parallel in [1, 4]:
                requests.append(dict(genes=[serialize(g) for g in genes], bams=list(map(str, bams)), reference=str(reference),
                                     options=dict(primary_only=primary), parallel=parallel))
                expected.append(dict(genes=[serialize(g) for g in modified], counts=found))
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want != got:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2))
            raise AssertionError(f"graph counting mismatch: {path}")
    (args.work / "report.json").write_text(json.dumps(dict(cases=len(actual), threads=[1, 4], result="pass (exact)"), indent=2) + "\n")
    print(f"PASS: {len(actual)} exact graph counting cases at 1/4 threads")


if __name__ == "__main__": main()
