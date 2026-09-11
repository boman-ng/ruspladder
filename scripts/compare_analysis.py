#!/usr/bin/env python3
"""Compare complete analyze_events public HDF5, with bounded batches and threads."""
import argparse
import contextlib
import copy
import json
import os
import pickle
import subprocess
import warnings
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
from spladder import settings
from spladder.alt_splice.analyze import analyze_events
from compare_annotations import serialize
from compare_collection import event_value, KINDS
from compare_count_io import compare_hdf5
from compare_merge import load_genes
from compare_verification import sources


def multi_source(upstream, root):
    """Use the already checked explicit four-exon skip, with two count samples."""
    gene, event, segments, edges, positions = list(sources(upstream))[-1]
    (root / "spladder").mkdir(parents=True, exist_ok=True)
    with (root / "spladder/genes_graph_conf3.merge_graphs.pickle").open("wb") as f:
        pickle.dump((np.array([gene], dtype=object), {}), f, -1)
    with (root / "merge_graphs_mult_exon_skip_C3.pickle").open("wb") as f:
        pickle.dump(np.array([event], dtype=object), f, -1)
    path = root / "spladder/genes_graph_conf3.merge_graphs.count.hdf5"
    with h5py.File(path, "w") as f:
        f["samples"] = np.array(["multi_1", "multi_2"], dtype="S")
        f["strains"] = h5py.SoftLink("/samples")
        f["segments"] = np.c_[segments, segments * 2]
        f["seg_pos"] = np.c_[positions, positions]
        f["edges"] = np.c_[edges[:, 1], edges[:, 1] * 2]
        f["edge_idx"] = edges[:, 0]
        f["gene_ids_segs"] = np.zeros((len(segments), 1), dtype=int)
        f["gene_ids_edges"] = np.zeros((len(edges), 1), dtype=int)
        f["gene_names"] = np.array([[gene.name]], dtype="S")
        f["seg_len"] = (gene.segmentgraph.segments[1] - gene.segmentgraph.segments[0])[:, None]
    return root


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "report.json").unlink(missing_ok=True)
    requests, expected, paths = [], [], []
    kinds = dict.fromkeys(KINDS, 0)
    directories = sorted((args.upstream / "tests").glob("*/results_merged*"))
    directories.append(multi_source(args.upstream, args.work / "multi-input"))
    for directory in directories:
        graph = directory / "spladder/genes_graph_conf3.merge_graphs.pickle"
        if not graph.exists(): continue
        genes = load_genes(graph)
        count_path = directory / "spladder/genes_graph_conf3.merge_graphs.count.hdf5"
        with h5py.File(count_path) as f: samples = f["samples"][:].astype(str)
        for kind in KINDS:
            event_path = directory / f"merge_graphs_{kind}_C3.pickle"
            if not event_path.exists(): continue
            with event_path.open("rb") as f: source = pickle.load(f, encoding="latin1")
            for profile in ["actual", "empty", "batches"]:
                if profile == "batches" and not len(source): continue
                events = copy.deepcopy(source)
                if profile == "empty": events = np.array([], dtype=object)
                elif profile == "batches":
                    events = np.array([copy.deepcopy(source[i % len(source)]) for i in range(261)], dtype=object)[::-1]
                    for i, event in enumerate(events): event.id = i + 1
                for parallel in [1, 4]:
                    root = args.work / str(len(requests))
                    (root / "spladder").mkdir(parents=True, exist_ok=True)
                    with (root / "spladder/genes_graph_conf3.merge_graphs.pickle").open("wb") as f: pickle.dump((genes, {}), f, -1)
                    linked = root / "spladder" / count_path.name
                    if not linked.exists(): os.link(count_path, linked)
                    with (root / event_path.name).open("wb") as f: pickle.dump(events, f, -1)
                    output = root / f"merge_graphs_{kind}_C3.counts.hdf5"
                    output.unlink(missing_ok=True)
                    options = SimpleNamespace(outdir=str(root), confidence=3, readlen=50, merge="merge_graphs", validate_sg=False,
                        samples=samples, use_anno_support=False, psi_min_reads=10, parallel=parallel, verbose=False, compress_text=False)
                    settings.default_settings(options)
                    settings.set_confidence_level(options)
                    for flag in ["txt", "struc", "bed", "gff3", "confirmed_txt", "confirmed_struc", "confirmed_bed", "confirmed_gff3", "confirmed_tcga", "confirmed_icgc"]:
                        setattr(options, "output_" + flag, False)
                    with (root / "reference.log").open("w") as log, contextlib.redirect_stdout(log), warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        analyze_events(kind, [], options)
                    with (root / f"merge_graphs_{kind}_C3.confirmed.pickle").open("rb") as f: expected.append(pickle.load(f).tolist())
                    config = dict(use_anno_support=False, exon_skip=options.exon_skip, mult_exon_skip=options.mult_exon_skip,
                        intron_retention=options.intron_retention, alt_prime=options.alt_prime, mutex_exons=options.mutex_exons)
                    requests.append(dict(genes=[serialize(g) for g in genes], events=[event_value(e) for e in events], kind=kind, counts=str(count_path),
                        output=str(root / "rust.hdf5"), samples=samples.tolist(), sample_idx=list(range(len(samples))), options=config, min_reads=10, parallel=parallel))
                    paths.append((output, root / "rust.hdf5"))
                    kinds[kind] += len(events)
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line)["confirmed"] for line in run.stdout.splitlines()]
    assert actual == expected
    arrays = sum(compare_hdf5(a, b) for a, b in paths)
    assert all(kinds.values()), kinds
    report = dict(cases=len(actual), arrays=arrays, events=kinds, threads=[1, 4], result="pass")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {report}")


if __name__ == "__main__": main()
