#!/usr/bin/env python3
"""Compare the isoform counts and sample/PSI ordering actually used by test CLI."""
import argparse
import contextlib
import io
import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
from spladder.alt_splice.quantify import quantify_from_counted_events
from compare_collection import KINDS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "report.json").unlink(missing_ok=True)
    paths = sorted((args.upstream / "tests").glob("*/results_merged*/merge_graphs_*_C3.counts.hdf5"))
    paths.extend(sorted(args.analysis.glob("[0-9]*/merge_graphs_mult_exon_skip_C3.counts.hdf5")))
    requests, expected = [], []
    kinds = dict.fromkeys(KINDS, 0)
    for path in paths:
        kind = path.name.removeprefix("merge_graphs_").removesuffix("_C3.counts.hdf5")
        with h5py.File(path) as f:
            if "conf_idx" not in f: continue
            samples = len(f["samples"])
        for profile in ["actual", "fractional", "filter_idx"]:
            local = args.work / f"counts-{len(requests)}.h5"
            shutil.copyfile(path, local)
            if profile == "fractional":
                with h5py.File(local, "r+") as f:
                    f["event_counts"][:, 1:, :] = f["event_counts"][:, 1:, :] * 0.75
            if profile == "filter_idx":
                with h5py.File(local, "r+") as f:
                    f["filter_idx"] = f["conf_idx"][:] + 100
                    labels = np.array([str(i) + ".npz" for i in range(samples)], dtype="S")
                    del f["samples"]
                    f["samples"] = labels
            for a, b in [([0], [samples - 1]), (list(reversed(range(samples))), list(range(samples))), ([], list(reversed(range(samples))))]:
                for high_memory in [False, True]:
                    with contextlib.redirect_stdout(io.StringIO()):
                        cov, psi, gene, event, ids, labels = quantify_from_counted_events(str(local), np.array(a, dtype=int), np.array(b, dtype=int), kind, SimpleNamespace(use_exon_counts=False, verbose=False), gen_event_ids=False, high_mem=high_memory)
                    assert ids is None
                    expected.append(dict(coverage=cov, psi=psi, gene_idx=gene, event_idx=event, samples=labels.tolist()))
                    requests.append(dict(path=str(local), group_a=a, group_b=b, kind=kind))
                    kinds[kind] += 1
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        assert got["samples"] == want["samples"], i
        for key in ["coverage", "psi", "gene_idx", "event_idx"]:
            np.testing.assert_array_equal(np.asarray(got[key], dtype=float), np.asarray(want[key], dtype=float), err_msg=f"case {i} {key}")
    assert all(kinds.values()), kinds
    report = dict(cases=len(actual), kinds=kinds, result="pass (exact)")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {report}")


if __name__ == "__main__": main()
