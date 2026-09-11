#!/usr/bin/env python3
"""Check negative-binomial, Cox-Reid, prior, trigamma and chi-square kernels."""
import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.stats import nbinom, chi2
from scipy.special import polygamma
from spladder.spladder_test import adj_loglikelihood_scalar, adj_loglikelihood_shrink_scalar_onedisper
parser = argparse.ArgumentParser()
parser.add_argument("probe", type=Path)
parser.add_argument("--work", type=Path, required=True)
args = parser.parse_args()
args.work.mkdir(parents=True, exist_ok=True)
(args.work / "report.json").unlink(missing_ok=True)
rng = np.random.default_rng(255)
requests, expected = [], []
for n in [8, 20, 100]:
    design = np.c_[np.ones(n), np.arange(n) % 2, np.arange(n) % 3 == 0]
    for dispersion in [1e-5, 0.01, 0.1, 2., 9.9]:
        for profile in ["regular", "zero", "large"]:
            response = rng.poisson(30, n).astype(float)
            mu = rng.uniform(0.1, 100, n)
            if profile == "zero": response[:] = 0
            if profile == "large": response *= 1e5; mu *= 1e5
            prior = rng.uniform(0.01, 100)
            fitted = rng.uniform(0.01, 10)
            requests.append(dict(design=design.tolist(), response=response.tolist(), mu=mu.tolist(), dispersion=dispersion, fitted=fitted, prior=prior))
            expected.append(dict(logpmf=nbinom.logpmf(response, 1/dispersion, (1/dispersion)/(1/dispersion+mu)).tolist(),
                adjusted=float(adj_loglikelihood_scalar(dispersion, design, response, mu, 1.)),
                shrink=float(adj_loglikelihood_shrink_scalar_onedisper(dispersion, design, response, mu, fitted, prior, 1.)),
                trigamma=float(polygamma(1, prior)), pvalue=float(1-chi2.cdf(prior, 1))))
run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
actual = [json.loads(line) for line in run.stdout.splitlines()]
assert len(actual) == len(expected)
for i, (want, got) in enumerate(zip(expected, actual)):
    for key in want:
        np.testing.assert_allclose(got[key], want[key], atol=1e-8, rtol=1e-6, equal_nan=True, err_msg=f"case {i} {key}")
report = dict(cases=len(actual), result="pass")
(args.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"PASS: {report}")
