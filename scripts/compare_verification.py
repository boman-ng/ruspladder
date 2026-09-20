#!/usr/bin/env python3
"""Compare all event feature vectors, verification flags and PSI to SplAdder."""
import argparse
import copy
import json
import pickle
import subprocess
import warnings
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
from spladder import settings
from spladder.alt_splice import verify
from spladder.classes.event import Event
from spladder.classes.gene import Gene
from spladder.classes.splicegraph import Splicegraph
from spladder.classes.segmentgraph import Segmentgraph
from spladder.helpers import compute_psi
from compare_annotations import serialize
from compare_collection import event_value, KINDS
from compare_merge import load_genes


def sources(upstream):
    for directory in sorted((upstream / "tests").glob("*/results_merged*")):
        graph = directory / "spladder/genes_graph_conf3.merge_graphs.pickle"
        if not graph.exists(): continue
        genes = load_genes(graph)
        with h5py.File(directory / "spladder/genes_graph_conf3.merge_graphs.count.hdf5", "r") as counts:
            for kind in KINDS:
                path = directory / f"merge_graphs_{kind}_C3.pickle"
                if not path.exists(): continue
                with path.open("rb") as source: events = pickle.load(source, encoding="latin1")
                for event in events:
                    gene = genes[event.gene_idx]
                    seg = np.where(counts["gene_ids_segs"][:, 0] == event.gene_idx)[0]
                    edge = np.where(counts["gene_ids_edges"][:, 0] == event.gene_idx)[0]
                    for sample in [0, counts["segments"].shape[1] - 1]:
                        yield gene, event, counts["segments"][seg, sample], np.c_[counts["edge_idx"][edge], counts["edges"][edge, sample]], counts["seg_pos"][seg, sample]
    gene = Gene(name="multi", chr="chr1", strand="+", start=10, stop=200)
    gene.exons = [np.array([[10, 30], [60, 80], [110, 130], [180, 200]]), np.array([[10, 30], [180, 200]])]
    gene.transcripts = ["long", "short"]
    gene.splicegraph = Splicegraph(gene)
    gene.segmentgraph = Segmentgraph(gene)
    gene.populate_annotated_introns()
    event = Event("mult_exon_skip", chr="chr1", strand="+")
    event.exons1, event.exons2 = gene.exons[1], gene.exons[0]
    event.gene_name, event.gene_idx, event.id = np.array(["multi"]), 0, 1
    event.set_annotation_flag(gene.introns_anno)
    n = gene.segmentgraph.segments.shape[1]
    a, b = np.where(gene.segmentgraph.seg_edges)
    yield gene, event, np.array([20., 3., 7., 30.]), np.c_[a * n + b, np.array([3, 2, 4, 5])], np.array([20., 20., 20., 20.])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    requests, expected = [], []
    totals = {kind: 0 for kind in KINDS}
    for gene, event, segments, edges, positions in sources(args.upstream):
        for profile in ["actual", "zero", "threshold", "empty", "invalid"]:
            cs, ce, cp = segments.copy(), edges.copy(), positions.copy()
            ev = copy.deepcopy(event)
            if profile == "zero": cs[:] = 0; ce[:, 1] = 0; cp[:] = 0
            elif profile == "threshold": cs[:] = 10; ce[:, 1] = 3; cp[:] = (gene.segmentgraph.segments[1] - gene.segmentgraph.segments[0]) * 0.9
            elif profile == "empty": ce = np.zeros((0, 2))
            elif profile == "invalid": ev.exons1[0, 0] = -1
            for annotation_support in [False, True]:
                options = SimpleNamespace(confidence=3, readlen=50, use_anno_support=annotation_support, psi_min_reads=10)
                settings.default_settings(options)
                settings.set_confidence_level(options)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    if not len(ce): flags, info = verify.verify_empty(ev.event_type)
                    elif ev.event_type == "intron_retention": flags, info = verify.verify_intron_retention(ev, gene, cs, ce, cp, options)
                    else:
                        name = "alt_prime" if ev.event_type.startswith("alt_") else ev.event_type
                        flags, info = getattr(verify, "verify_" + name)(ev, gene, cs, ce, options)
                    psi, a, b = compute_psi(info[np.newaxis, :], ev.event_type, options)
                assert not np.isinf(info).any() and not np.isinf(psi).any()
                config = dict(use_anno_support=annotation_support, exon_skip=options.exon_skip, mult_exon_skip=options.mult_exon_skip,
                              intron_retention={k: v for k, v in options.intron_retention.items() if k != "read_filter"},
                              alt_prime=options.alt_prime, mutex_exons=options.mutex_exons)
                requests.append(dict(gene=serialize(gene), event=event_value(ev),
                    counts=dict(segments=cs.tolist(), seg_pos=cp.tolist(), edges=[[int(x), int(y)] for x, y in ce], edges_float=True),
                    options=config, min_reads=10))
                expected.append(dict(verified=list(map(bool, flags)), info=info, psi=np.array([psi[0], a[0], b[0]])))
                totals[ev.event_type] += 1
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    max_difference = 0.0
    for i, (want, got) in enumerate(zip(expected, actual)):
        try:
            assert want["verified"] == got["verified"]
            for key in ["info", "psi"]:
                values = np.asarray(got[key], dtype=float)
                np.testing.assert_allclose(values, want[key], atol=1e-10, rtol=1e-8, equal_nan=True)
                finite = np.isfinite(values) & np.isfinite(want[key])
                if finite.any(): max_difference = max(max_difference, float(np.abs(values[finite] - want[key][finite]).max()))
        except AssertionError:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected={k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in want.items()}, actual=got), indent=2))
            raise AssertionError(f"event verification mismatch: {path}")
    assert all(totals.values()), totals
    report = dict(cases=len(actual), kinds=totals, max_absolute_difference=max_difference, atol=1e-10, rtol=1e-8, result="pass")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(actual)} event verification/PSI cases; max absolute difference {max_difference}; {totals}")


if __name__ == "__main__": main()
