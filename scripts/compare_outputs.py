#!/usr/bin/env python3
"""Compare every non-visual event text format byte-for-byte after decompression."""
import argparse
import contextlib
import gzip
import json
import pickle
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
from spladder.alt_splice import write
from compare_collection import event_value, KINDS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "report.json").unlink(missing_ok=True)
    directories = sorted((args.upstream / "tests").glob("*/results_merged*"))
    # The last six analysis cases are the explicit multi-exon fixture.
    directories.extend(sorted(args.analysis.glob("[0-9]*"), key=lambda x: int(x.name))[-6:-5])
    requests, pairs = [], []
    totals = dict.fromkeys(KINDS, 0)
    for directory in directories:
        for kind in KINDS:
            event_path = directory / f"merge_graphs_{kind}_C3.pickle"
            count_path = directory / f"merge_graphs_{kind}_C3.counts.hdf5"
            if not event_path.exists(): continue
            with event_path.open("rb") as f: events = pickle.load(f, encoding="latin1")
            if not len(events): continue
            with h5py.File(count_path) as f:
                samples = f["samples"][:].astype(str)
                confirmed = f["conf_idx"][:] if "conf_idx" in f else np.array([], dtype=int)
            for profile in ["actual", "invalid"]:
                counts = args.work / f"counts-{len(requests)}.hdf5"
                shutil.copyfile(count_path, counts)
                if profile == "invalid":
                    with h5py.File(counts, "r+") as f: f["event_counts"][:, 0, :] = 0; f["psi"][:] = np.nan
                for select in [np.arange(len(events)), confirmed, np.array([], dtype=int)]:
                    for format in ["txt", "structured", "bed", "gff3", "gtf", "tcga", "icgc"]:
                        for compressed in ([False, True] if format in ["txt", "structured", "tcga", "icgc"] else [False]):
                            root = args.work / str(len(requests))
                            root.mkdir(exist_ok=True)
                            suffix = ".gz" if compressed else ""
                            expected, actual = root / ("reference.txt" + suffix), root / ("rust.txt" + suffix)
                            with (root / "reference.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                                if format == "txt": write.write_events_txt(str(expected), samples, events, str(counts), event_idx=select)
                                elif format == "structured": write.write_events_structured(str(expected), events, str(counts), idx=select)
                                elif format == "bed": write.write_events_bed(str(expected), events, idx=select)
                                elif format in ["gff3", "gtf"]: write.write_events_gff3(str(expected), events, idx=select, as_gtf=format == "gtf")
                                else: getattr(write, "write_events_" + format)(str(expected), samples, events, str(counts), event_idx=select)
                            requests.append(dict(events=[event_value(e) for e in events], format=format, counts=str(counts), output=str(actual), samples=samples.tolist(), indices=select.tolist()))
                            pairs.append((expected, actual, compressed))
                            totals[kind] += 1
    with (args.work / "rust.log").open("w") as log:
        run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, stderr=log, check=True)
    assert len(run.stdout.splitlines()) == len(requests)
    for expected, actual, compressed in pairs:
        read = lambda p: gzip.open(p, "rb").read() if compressed else p.read_bytes()
        a, b = read(expected), read(actual)
        if a != b:
            raise AssertionError(f"text mismatch: {expected}, {actual}\nexpected: {a[:2000]!r}\nactual: {b[:2000]!r}")
    assert all(totals.values()), totals
    report = dict(files=len(pairs), kinds=totals, comparison="exact decompressed bytes", result="pass")
    (args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {report}")


if __name__ == "__main__": main()
