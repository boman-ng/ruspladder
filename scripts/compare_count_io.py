#!/usr/bin/env python3
"""Compare public graph-count/expression HDF5, including shape, dtype and links."""
import argparse
import contextlib
import copy
import json
import pickle
import subprocess
import warnings
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
from spladder.count import count_graph_coverage_wrapper, compute_gene_expression, get_size_factors
from compare_annotations import serialize
from compare_merge import load_genes


def compare_hdf5(want, got):
    compared = 0
    with h5py.File(want) as a, h5py.File(got) as b:
        assert sorted(a) == sorted(b), (want, sorted(a), sorted(b))
        for key in a:
            x, y = a[key], b[key]
            assert x.shape == y.shape and x.dtype == y.dtype, (key, x.shape, y.shape, x.dtype, y.dtype)
            link_a, link_b = a.get(key, getlink=True), b.get(key, getlink=True)
            assert type(link_a) == type(link_b), key
            if isinstance(link_a, h5py.SoftLink): assert link_a.path == link_b.path, key
            if x.dtype.kind == "f": np.testing.assert_allclose(x[:], y[:], atol=1e-10, rtol=1e-8, equal_nan=True, err_msg=key)
            else: np.testing.assert_array_equal(x[:], y[:], err_msg=key)
            assert x.compression == y.compression, (key, x.compression, y.compression)
            compared += 1
    return compared


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("normalization_probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "report.json").unlink(missing_ok=True)
    cases = []
    basic = args.upstream / "tests/testcase_basic"
    for strand in ["pos", "neg"]:
        cases.append((load_genes(basic / f"results_merged_{strand}/spladder/genes_graph_conf3.merge_graphs.pickle"),
                      [basic / f"data/align/{strand}_{i}.bam" for i in [1, 2]], basic / f"data/genome_{strand}.fa"))
    events = args.upstream / "tests/testcase_events"
    for ext, directory in [("bam", "results_merged"), ("cram", "results_merged_cram")]:
        cases.append((load_genes(events / f"{directory}/spladder/genes_graph_conf3.merge_graphs.pickle"),
                      [events / f"data/align/testcase_events_1_sample{i}.{ext}" for i in [1, 2]], events / "data/genome.fa"))
    genes, bams, reference = copy.deepcopy(cases[2])
    for gene in genes: gene.chr = "missing"
    cases.append((genes, bams, reference))
    requests, comparisons = [], []
    for source, bams, reference in cases:
        for mode in ["all", "single", "merge_single", "reordered"]:
            for parallel in [1, 4]:
                root = args.work / str(len(requests))
                root.mkdir(exist_ok=True)
                genes = copy.deepcopy(source)
                graph_path = root / "genes.pickle"
                with graph_path.open("wb") as out: pickle.dump((genes, {}), out, -1)
                labels = np.array(["sample_1", "sample_long_2"])
                selected = [1, 0, 1] if mode == "reordered" else None
                options = SimpleNamespace(sparse_bam=False, var_aware=False, primary_only=False, mm_tag="NM",
                    ref_genome=str(reference), samples=labels, merge="single" if mode == "single" else "merge_graphs", readlen=50, verbose=False)
                qmode = "single" if mode == "merge_single" else "all"
                with (root / "reference.log").open("w") as log, contextlib.redirect_stdout(log), warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    count_graph_coverage_wrapper(str(graph_path), str(root / "reference.counts.h5"), list(map(str, bams)), options, sample_idx=1, qmode=qmode)
                    if mode in ["all", "reordered"]:
                        compute_gene_expression(options, str(graph_path), str(root / "reference.counts.h5"), str(root / "reference.expression.h5"), sample_idx=selected)
                requested_bams = [bams[1]] if mode == "single" else [bams[0]] if mode == "merge_single" else bams
                requests.append(dict(genes=[serialize(g) for g in source], bams=list(map(str, requested_bams)), reference=str(reference),
                    samples=labels.tolist(), counts=str(root / "rust.counts.h5"), expression=str(root / "rust.expression.h5") if mode in ["all", "reordered"] else None,
                    sample_idx=selected, readlen=50, parallel=parallel))
                comparisons.append((root / "reference.counts.h5", root / "rust.counts.h5"))
                if mode in ["all", "reordered"]: comparisons.append((root / "reference.expression.h5", root / "rust.expression.h5"))
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    assert len(run.stdout.splitlines()) == len(requests)
    arrays = sum(compare_hdf5(a, b) for a, b in comparisons)
    rng = np.random.default_rng(412)
    inputs, expected = [], []
    reference_errors = 0
    for rows in [1, 2, 7, 129, 1001]:
        for samples in [1, 2, 4, 8]:
            for profile in ["zero", "integer", "fractional", "large"]:
                counts = np.zeros((rows, samples)) if profile == "zero" else rng.integers(0, 1000, (rows, samples)).astype(float)
                if profile == "fractional": counts *= rng.random(counts.shape)
                if profile == "large": counts *= 1e8
                for kind in ["geomean", "tc", "uq"]:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        try:
                            expected.append(get_size_factors(counts, SimpleNamespace(verbose=False), kind))
                        except NameError as error:
                            assert kind == "uq" and "scoreatpercentile" in str(error)
                            expected.append(None)
                            reference_errors += 1
                    inputs.append(dict(counts=counts.tolist(), samples=samples, kind=kind))
    run = subprocess.run([str(args.normalization_probe)], input="".join(json.dumps(r) + "\n" for r in inputs), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want is None:
            assert "scoreatpercentile" in got["error"]
        else:
            np.testing.assert_allclose(np.asarray(got["values"], dtype=float), want, atol=1e-10, rtol=1e-8, equal_nan=True, err_msg=f"normalization case {i}")
    report = dict(count_cases=len(requests), files=len(comparisons), arrays=arrays, normalizations=len(actual) - reference_errors, upstream_uq_failures=reference_errors, threads=[1, 4], result="pass")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {report}")


if __name__ == "__main__": main()
