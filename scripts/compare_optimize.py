#!/usr/bin/env python3
"""Check SciPy bounded-minimizer values, statuses and evaluation counts."""
import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
parser = argparse.ArgumentParser()
parser.add_argument("probe", type=Path)
parser.add_argument("--work", type=Path, required=True)
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
(args.work / "report.json").unlink(missing_ok=True)
requests, expected = [], []
for bounds in [[0, 10], [-100, 100], [1, 1], [1e-7, 2e-7]]:
    for center in [-100, 0, 1e-8, 0.1, 5, 9.99999, 100]:
        for kind in ["quadratic", "quartic", "linear", "flat", "nan", "absolute"]:
            function = lambda x: {"quadratic": lambda: (x-center)**2, "quartic": lambda: (x-center)**4, "linear": lambda: x-center,
                "flat": lambda: 2., "nan": lambda: np.nan, "absolute": lambda: abs(x-center)}[kind]()
            result = minimize_scalar(function, bounds=bounds, method="bounded", options=dict(xatol=1e-5))
            requests.append(dict(function=kind, bounds=bounds, center=center, tolerance=1e-5))
            expected.append(dict(x=result.x, fun=result.fun, status=result.status, evaluations=result.nfev))
run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
actual = [json.loads(line) for line in run.stdout.splitlines()]
assert len(actual) == len(expected)
for i, (want, got) in enumerate(zip(expected, actual)):
    assert want["status"] == got["status"] and want["evaluations"] == got["evaluations"], (i, want, got)
    np.testing.assert_allclose(np.array([got["x"], got["fun"]], dtype=float), [want["x"], want["fun"]], atol=1e-12, rtol=1e-12, equal_nan=True, err_msg=f"case {i}")
report = dict(cases=len(actual), result="pass")
(args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"PASS: {report}")
