#!/usr/bin/env python3
"""Compare the complete raw/trend/shrunk dispersion and LRT numerical sequence."""
import argparse
import contextlib
import io
import json
import subprocess
import warnings
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.special import polygamma
from spladder.spladder_test import estimate_dispersion, fit_dispersion, adjust_dispersion, test_count, calculate_varPrior, run_testing
parser = argparse.ArgumentParser()
parser.add_argument("probe", type=Path)
parser.add_argument("--work", type=Path, required=True)
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
(args.work / "report.json").unlink(missing_ok=True)
rng = np.random.default_rng(381)
requests, expected = [], []
for samples in [4, 8, 16]:
    n = 2 * samples
    design = np.zeros((n, 4))
    design[:, 0] = 1
    design[:samples//2, 1:3] = 1
    design[samples:samples+samples//2, 2] = 1
    design[samples:, 3] = 1
    null = design[:, [0, 2, 3]]
    sf = rng.uniform(0.5, 2, n)
    for profile in ["regular", "filtered"]:
        means = np.exp(rng.uniform(2, 6, 60))
        dispersions = rng.uniform(0.01, 1, 60)
        counts = np.array([rng.negative_binomial(1/a, (1/a)/(1/a+mu*sf)) for a, mu in zip(dispersions, means)]).astype(float)
        if profile == "filtered": counts[0, :n//2] = 0; counts[1, :3*n//4] = 0; counts[2] = 1
        selected = np.ones(len(counts), dtype=bool); selected[::7] = False
        options = SimpleNamespace(min_count=10, max_0_frac=0.5, verbose=False, parallel=1, diagnose_plots=False)
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            raw, raw_conv = estimate_dispersion(counts, design, sf, options, "exon_skip")
            estimated = raw.copy()
            fitted, coefficients, indices = fit_dispersion(counts, raw, raw_conv, sf, options, design, "exon_skip")
            adjusted, adjusted_conv = adjust_dispersion(counts, design, raw, fitted, indices, sf, options, "exon_skip")
            pvalues = test_count(counts, adjusted, sf, null, design, options, selected)
            trigamma = float(polygamma(1, (n-4)/2))
            prior = float(calculate_varPrior(raw, fitted, indices, trigamma))
            final_p, final_cov, final_raw, final_adj = run_testing(counts, null, design, sf, options, "exon_skip", selected)
        want = dict(estimated=estimated.ravel(), raw=dict(values=raw.ravel(), converged=raw_conv.ravel()), trend=dict(values=fitted.ravel(), coefficients=coefficients, indices=indices),
            adjusted=dict(values=adjusted.ravel(), converged=adjusted_conv.ravel()), pvalues=pvalues.ravel(), trigamma=trigamma, prior=prior,
            final=dict(pvalues=final_p, coverage=final_cov, dispersion_raw=final_raw.ravel(), dispersion_adjusted=final_adj.ravel()))
        for parallel in [1, 4]:
            requests.append(dict(counts=counts.tolist(), null=null.tolist(), design=design.tolist(), sf=sf.tolist(), selected=selected.tolist(), options=dict(min_count=10, max_zero_fraction=0.5), parallel=parallel))
            expected.append(want)
run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
actual = [json.loads(line) for line in run.stdout.splitlines()]
assert len(actual) == len(expected)
def compare(want, got, key=""):
    if isinstance(want, dict):
        for k in want: compare(want[k], got[k], key+"."+k)
    else:
        values = np.asarray(got, dtype=float)
        if key.endswith(("converged", "indices")): np.testing.assert_array_equal(values, want, err_msg=key)
        else: np.testing.assert_allclose(values, want, atol=1e-8, rtol=1e-6, equal_nan=True, err_msg=key)
        if key.endswith(".pvalues"):
            np.testing.assert_array_equal(np.round(values, 6), np.round(want, 6))
            for alpha in [0.01, 0.05, 0.1]: np.testing.assert_array_equal(values <= alpha, want <= alpha)
for i, (want, got) in enumerate(zip(expected, actual)):
    try: compare(want, got)
    except AssertionError:
        default = lambda x: x.tolist() if isinstance(x, np.ndarray) else x
        (args.work / "mismatch.json").write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), default=default, indent=2) + "\n")
        raise
for i in range(0, len(actual), 2): assert actual[i] == actual[i+1], f"thread nondeterminism: {i}"
report = dict(cases=len(actual), rows_per_case=60, threads=[1, 4], result="pass")
(args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"PASS: {report}")
