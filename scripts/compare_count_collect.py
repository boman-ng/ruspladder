#!/usr/bin/env python3
"""Check single-quantification collection, first-file types, links and maxshape."""
import argparse
import contextlib
import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
from spladder.count import collect_single_quantification_results
from compare_count_io import compare_hdf5


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "report.json").unlink(missing_ok=True)
    requests, pairs = [], []
    for selected in [[0, 2], [8, 10], [16, 18], [24, 26], [32, 16], [16, 32]]:
        for order in [[0, 1], [1, 0, 1]]:
            root = args.work / str(len(requests))
            (root / "spladder").mkdir(parents=True, exist_ok=True)
            paths = []
            for i, source in enumerate(selected):
                path = root / f"spladder/genes_graph_conf3.merge_graphs.part_{i}.count.hdf5"
                shutil.copyfile(args.counts / str(source) / "reference.counts.h5", path)
                # The collection contract explicitly truncates names to S255.
                if i == 1:
                    with h5py.File(path, "r+") as f:
                        values = f["samples"][:].astype(str)
                        values = np.array([v + "x" * 260 for v in values], dtype="S")
                        del f["samples"]
                        f["samples"] = values
                paths.append(path)
            options = SimpleNamespace(outdir=str(root), confidence=3, merge="merge_graphs", verbose=False, samples=np.array(["part_0", "part_1"]), validate_sg=False)
            expected = root / "reference.h5"
            with (root / "reference.log").open("w") as log, contextlib.redirect_stdout(log):
                collect_single_quantification_results(str(expected), order, options)
            actual = root / "rust.h5"
            pairs.append((expected, actual))
            requests.append(dict(paths=[str(paths[i]) for i in order], output=str(actual)))
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    assert len(run.stdout.splitlines()) == len(requests)
    arrays = sum(compare_hdf5(a, b) for a, b in pairs)
    for a, b in pairs:
        with h5py.File(a) as x, h5py.File(b) as y:
            for key in x: assert x[key].maxshape == y[key].maxshape, (key, x[key].maxshape, y[key].maxshape)
    report = dict(cases=len(requests), arrays=arrays, result="pass")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {report}")


if __name__ == "__main__": main()
