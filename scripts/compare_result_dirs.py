#!/usr/bin/env python3
"""Compare full native result directories, or Python public files and pickles.

Native comparisons require exact float bits. Python comparisons retain the
project's numerical tolerance (atol=1e-10, rtol=1e-8), recording largest errors.
"""
import argparse
import gzip
import json
import pickle
import subprocess
from pathlib import Path

import h5py
import numpy as np


def compare(reference, native, exporter=None):
    report = dict(hdf5_files=0, datasets=0, text_files=0, graphs=0, events=0,
                  max_abs_float_error=0.0, python_reference=exporter is not None)
    reference_files = {p.relative_to(reference) for p in reference.rglob("*") if p.is_file()}
    native_files = {p.relative_to(native) for p in native.rglob("*") if p.is_file()}
    def public(p):
        return p.suffix in (".txt", ".gz", ".bed", ".gff3") or p.name.endswith((".count.hdf5", ".counts.hdf5", ".gene_exp.hdf5"))
    assert {p for p in reference_files if public(p)} == {p for p in native_files if public(p)}, "public file inventory"
    if exporter is None:
        assert {p for p in reference_files if p.suffix == ".hdf5"} == {p for p in native_files if p.suffix == ".hdf5"}, "HDF5 file inventory"
    for relative in sorted(reference_files):
        source, target = reference / relative, native / relative
        if source.suffix == ".hdf5":
            with h5py.File(source) as a, h5py.File(target) as b:
                names_a, names_b = [], []
                a.visit(names_a.append)
                b.visit(names_b.append)
                assert names_a == names_b, (relative, "object names")
                for name in ["", *names_a]:
                    x, y = a[name] if name else a, b[name] if name else b
                    assert sorted(x.attrs) == sorted(y.attrs), (relative, name, "attribute names")
                    for key in x.attrs:
                        np.testing.assert_array_equal(x.attrs[key], y.attrs[key])
                    if isinstance(x, h5py.Group):
                        assert sorted(x) == sorted(y), (relative, name, "group members")
                        for key in x:
                            u, v = x.get(key, getlink=True), y.get(key, getlink=True)
                            assert type(u) is type(v), (relative, name, key, "link type")
                            if isinstance(u, h5py.SoftLink):
                                assert u.path == v.path, (relative, name, key, "link target")
                for name in names_a:
                    x, y = a[name], b[name]
                    assert type(x) is type(y), (relative, name, "object type")
                    if not isinstance(x, h5py.Dataset):
                        continue
                    assert x.shape == y.shape and x.dtype == y.dtype, (relative, name, "shape/dtype")
                    assert x.compression == y.compression, (relative, name, "compression")
                    u, v = x[()], y[()]
                    if x.dtype.kind == "f":
                        if exporter is None:
                            assert np.asarray(u).tobytes() == np.asarray(v).tobytes(), (relative, name, "float bits")
                        else:
                            np.testing.assert_allclose(u, v, rtol=1e-8, atol=1e-10, equal_nan=True, err_msg=f"{relative}:{name}")
                            finite = np.isfinite(u) & np.isfinite(v)
                            delta = np.abs(np.asarray(u)[finite] - np.asarray(v)[finite])
                            report["max_abs_float_error"] = max(report["max_abs_float_error"], float(delta.max(initial=0)))
                    else:
                        np.testing.assert_array_equal(u, v, err_msg=f"{relative}:{name}")
                    report["datasets"] += 1
            report["hdf5_files"] += 1
        elif source.suffix in (".txt", ".gz", ".bed", ".gff3"):
            read = lambda p: gzip.open(p, "rb").read() if p.suffix == ".gz" else p.read_bytes()
            assert read(source) == read(target), (relative, "text")
            report["text_files"] += 1
        elif exporter and source.suffix == ".pickle" and not source.name.endswith(".confirmed.pickle"):
            from compare_annotations import serialize
            from compare_collection import event_value
            from compare_merge import load_genes
            if relative.parts[0] == "spladder":
                request = dict(graph=str(target.with_suffix(".hdf5")))
                expected = [serialize(g) for g in load_genes(source)]
                report["graphs"] += 1
            else:
                request = dict(events=str(target.with_suffix(".events.hdf5")))
                with source.open("rb") as f:
                    expected = [event_value(e) for e in pickle.load(f)]
                report["events"] += 1
            actual = subprocess.check_output([str(exporter)], input=json.dumps(request) + "\n", text=True)
            assert json.loads(actual) == expected, (relative, "serialized graph/events")
    assert report["datasets"] > 0 and report["text_files"] > 0, "incomplete output"
    return dict(result="pass", **report)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("reference", type=Path)
    p.add_argument("native", type=Path)
    p.add_argument("--exporter", type=Path)
    p.add_argument("--report", type=Path, required=True)
    a = p.parse_args()
    report = compare(a.reference, a.native, a.exporter)
    a.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)
