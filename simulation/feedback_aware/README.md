# Feedback-Aware Tuning of Recursive Q-Learning: simulation module

This module reproduces the six simulation tables in the current main paper and Supplementary Material. The numerical core is the uploaded v17 implementation, with internal version **1.1.0**; it has not been altered by this packaging audit. The packaging revision is recorded separately in `PACKAGING_REVISION`.

**This is not the complete manuscript repository.** The primary, response-gated ADHD diagnostic requires its own analysis scripts and data instructions and is not supplied by this two-stage simulation module. Older Pretest-Stein R sensitivity experiments are not the method evaluated in the current FA tables. The numerical core was audited locally and is distributed here as the current feedback-aware simulation module.

## Fresh environment and reproducible commands

The audit reran every production configuration with **Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0**. Use a fresh environment created by Python 3.11 or later. Do not copy a virtual environment from another computer. The original upload included a Python-3.7.4 environment whose workers stopped because NumPy was unavailable; that environment and those failed outputs are deliberately excluded here.

From this directory, using an installed Python 3.13 interpreter:

```sh
python3.13 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-tested.txt
python -m unittest discover -s tests -v
python run_feedback_tuning.py smoke --workers 1 --output results/smoke
python run_feedback_tuning.py paper --workers 3 --output results/paper
python review/gaussian_audit.py --output results/supplement_gaussian
```

A compatible `python3.11` or `python3.12` can replace `python3.13` when creating the environment, but this audit's exact replay was performed on 3.13.5. Check `python --version` inside the environment. `requirements-tested.txt` pins the NumPy and SciPy versions used in the replay.

The paper run uses four local signals, sample sizes 250/1000/4000, and 5,000 independent datasets per cell. The log-RSS and quadratic variants share the same 60,000 datasets; this gives 120,000 FA analyses, not 120,000 independent datasets. Seed 20260909 and per-cell streams are recorded in `run_manifest.json`. The library contains the nine penalty pairs in `{0.5,1,2}^2` and an all-wide reference. The tuning parameter is fixed at T=2 and the reference weights are uniform.

Completed jobs can be resumed. Changes to source files or run settings require a new output directory, rather than mixing incompatible chunks. Set `--workers 1` for serial execution or `--criterion logrss` for only the main finite-sample variant.

## Which output reproduces which table?

| Manuscript table | Command/output |
|---|---|
| Main Table 1: Gaussian risks | `results/paper/gaussian/gaussian_results.csv` |
| Main Table 2: paired risk differences | `results/paper/logrss/main_table_gaps.csv` |
| Supplement Table S1: full risks and risk estimates | `results/paper/logrss/smart_results.csv` |
| Supplement Table S2: criterion variants | `results/paper/logrss/` and `results/paper/wald/` |
| Supplement Table S3: adaptation terms | `results/supplement_gaussian/table_s3_adaptation.csv` |
| Supplement Table S4: fixed-T sensitivity | `results/supplement_gaussian/table_s4_tuning.csv` |

`review/gaussian_audit.py` was newly supplied during the audit to fill the reproduction entry point named in the supplement. It is not presented as a recovered historical script. It uses the unchanged numerical core and recomputes the supplement's T=1,2,4 calculations at two quadrature orders.

`reference_results/` contains selected original production summaries used as replay baselines. A fresh full run in the audit reproduced all 1,674 compared numerical entries exactly. All 224 printed numerical results and Monte Carlo standard errors in the six tables match at their displayed precision. Selected audit evidence is in `validation/`; checksums for the distributed numerical core and retained reference summaries are in `SHA256SUMS.txt`.

## Model, target and interpretation

`fit_smart` takes observed `(X1,A1,X2,A2,Y1,Y2)` only. Unknown local signals and true Q-functions are used by the simulator and population-risk evaluator, not by the fitter. Each specification performs its own backward Q-learning fit. The final output is a weighted combination of completed stage-1 fitted Q-functions, with no new regression after combination.

The implementation is for the fully observed, two-stage E3-prime model. A general fixed-K transfer theorem in the paper does not make this a generic K-stage software package. Gaussian risk identities are exact under their smooth-model assumptions; the empirical SMART risk estimate is first-order, not finite-sample unbiased. The risk is Q-function prediction risk, not treatment-regime value, and there is no uniform improvement claim over Akaike weighting or hard AIC.

The source parameter `temperature` is the fixed tuning parameter **T** in the paper. Internal `candidate_maps` and `score` names are retained for backward compatibility: the maps represent complete Q-learning specifications, and the scores are the risk estimates S_j. `wald` is the manuscript's quadratic implementation. These code names do not change the mathematical definitions.

Important CSV naming conventions:

- `fa_minus_akaike`: **FA minus Akaike**; positive means FA has greater risk. The publication tables use the opposite sign.
- `main_table_gaps.csv`: **Akaike minus FA**; positive favours FA.
- `risk_akaike_limit`: legacy name for the **finite-sample empirical risk** of the `(1,1)` specification; this is not the asymptotic limit.
- `limit_akaike`: the separate asymptotic Akaike benchmark.
- `reported_scaled_sure`: reported first-order estimate of full n-scaled stage-1 risk, including the common component.

The `(1,1)` specification is exactly ordinary recursive Akaike weighting in the log-RSS implementation when safeguards are inactive; the quadratic version agrees at first order. Hard AIC is an external Gaussian comparator, not a differentiated library member.

A prespecified Gram failure returns a common zero prediction. Its prediction risk must remain in unconditional risk summaries, and the availability of risk estimates must be reported separately. No such failure occurred in the reproduced paper runs. Independent stress checks of action changes and coefficient caps validate implementation mechanics, not new theoretical coverage.

## Core source hashes

See `validation/core_hashes.json`. The estimator, simulator, Gaussian verifier and runner have the same bytes as in the uploaded v17 source. Only packaging, documentation, validation evidence and the supplemental-table entry point were added or updated.
