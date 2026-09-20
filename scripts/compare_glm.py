#!/usr/bin/env python3
"""Compare the IRLS paths SplAdder uses: NB/log and Gamma/identity."""
import argparse
import json
import subprocess
import warnings
from pathlib import Path
import numpy as np
import statsmodels.api as sm

parser = argparse.ArgumentParser()
parser.add_argument("probe", type=Path)
parser.add_argument("--work", type=Path, required=True)
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
(args.work / "report.json").unlink(missing_ok=True)
rng = np.random.default_rng(1245)
requests, expected = [], []
for n in [4, 8, 16, 40, 128]:
    design = np.zeros((n, 4))
    design[:, 0] = 1
    design[:n//4, 1:3] = 1
    design[n//2:n//2+n//4, 2] = 1
    design[n//2:, 3] = 1
    offset = np.log(rng.uniform(0.25, 4, n))
    for profile in ["poisson", "dispersed", "zeros", "all_zero", "constant", "large"]:
        y = rng.poisson(10, n).astype(float)
        if profile == "dispersed": y = rng.negative_binomial(2, 0.1, n).astype(float)
        if profile == "zeros": y[:n//3] = 0
        if profile == "all_zero": y[:] = 0
        if profile == "constant": y[:] = 3
        if profile == "large": y *= 100000
        for alpha in [0.00001, 0.01, 0.1, 2, 9.9]:
            for full in [False, True]:
                x = design if full else design[:, [0, 2, 3]]
                requests.append(dict(response=y.tolist(), design=x.tolist(), offset=offset.tolist(), family=dict(kind="negative_binomial", alpha=alpha)))
for n in [8, 30, 128, 1001]:
    for profile in ["regular", "outlier", "small", "constant", "zero"]:
        means = np.exp(rng.uniform(0, 10, n))
        x = np.c_[1 / means, np.ones(n)]
        y = (1.5 / means + 0.05) * rng.gamma(2, 0.5, n)
        if profile == "outlier": y[:2] *= 100
        if profile == "small": y *= 1e-5
        if profile == "constant": y[:] = 0.1
        if profile == "zero": y[:] = 0
        requests.append(dict(response=y.tolist(), design=x.tolist(), offset=[0.] * n, family=dict(kind="gamma_identity")))
failures = {}
for request in requests:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        family = request["family"]
        model_family = sm.families.NegativeBinomial(alpha=family["alpha"]) if family["kind"] == "negative_binomial" else sm.families.Gamma(sm.families.links.Identity())
        try:
            fit = sm.GLM(np.array(request["response"]), np.array(request["design"]), offset=np.array(request["offset"]), family=model_family).fit()
            expected.append(dict(params=fit.params.tolist(), mu=fit.mu.tolist(), deviance=float(fit.deviance), scale=float(fit.scale), converged=bool(fit.converged), iterations=fit.fit_history["iteration"]))
        except (ValueError, np.linalg.LinAlgError) as error:
            expected.append(None)
            failures[type(error).__name__ + ": " + str(error)] = failures.get(type(error).__name__ + ": " + str(error), 0) + 1
run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
actual = [json.loads(line) for line in run.stdout.splitlines()]
assert len(actual) == len(expected)
passed, iteration_differences = 0, []
for i, (want, got) in enumerate(zip(expected, actual)):
    try:
        if want is None: assert "error" in got
        else:
            got = got["fit"]
            assert got["converged"] == want["converged"]
            for key in ["params", "mu", "deviance", "scale"]:
                np.testing.assert_allclose(got[key], want[key], atol=1e-8, rtol=1e-6, equal_nan=True, err_msg=key)
            if got["iterations"] != want["iterations"]: iteration_differences.append(dict(case=i, expected=want["iterations"], actual=got["iterations"]))
            passed += 1
    except (AssertionError, KeyError):
        mismatch = args.work / "mismatch.json"
        mismatch.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2) + "\n")
        raise AssertionError(f"GLM mismatch: {mismatch}")
report = dict(cases=passed, upstream_failures=failures, iteration_differences=iteration_differences, result="pass")
(args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"PASS: {report}")
