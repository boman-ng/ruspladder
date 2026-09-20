#!/usr/bin/env python3
"""Compare every SplAdder correction choice and its observed failure paths."""
import argparse
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from spladder.spladder_test import adj_pval

parser = argparse.ArgumentParser()
parser.add_argument("probe", type=Path)
parser.add_argument("--work", type=Path, required=True)
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
(args.work / "report.json").unlink(missing_ok=True)
rng = np.random.default_rng(429)
requests, expected = [], []
errors = {}
for n in [0, 1, 2, 7, 128, 1001]:
    for profile in ["uniform", "small", "ties", "zero", "one", "missing", "all_missing"]:
        p = rng.random(n)
        if profile == "small": p *= 0.005
        if profile == "ties": p = np.round(p, 1)
        if profile == "zero": p[:] = 0
        if profile == "one": p[:] = 1
        if profile == "missing": p[::3] = np.nan
        if profile == "all_missing": p[:] = np.nan
        for method in ["BH", "Bonferroni", "Holm", "Hochberg", "Hommel", "BY", "TSBH"]:
            try: expected.append(adj_pval(p, SimpleNamespace(correction=method)))
            except (ZeroDivisionError, ValueError) as error:
                assert not np.any(~np.isnan(p)) or method == "TSBH", (method, error)
                expected.append(None)
                key = type(error).__name__ + ": " + str(error)
                errors[key] = errors.get(key, 0) + 1
            requests.append(dict(pvalues=[None if np.isnan(x) else x for x in p], method=method))
run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
actual = [json.loads(line) for line in run.stdout.splitlines()]
assert len(actual) == len(expected)
passed = 0
for i, (want, got) in enumerate(zip(expected, actual)):
    if want is None: assert "error" in got, i
    else:
        values = np.asarray(got["values"], dtype=float)
        np.testing.assert_allclose(values, want, atol=1e-8, rtol=1e-6, equal_nan=True, err_msg=f"case {i}")
        np.testing.assert_array_equal(np.round(values, 6), np.round(want, 6))
        for alpha in [0.01, 0.05, 0.1]: np.testing.assert_array_equal(values <= alpha, want <= alpha)
        passed += 1
report = dict(cases=passed, upstream_failures=errors, result="pass")
(args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"PASS: {report}")
