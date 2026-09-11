# Feedback-Aware Tuning of Recursive Q-Learning

This repository contains reproducibility materials for the manuscript

**Feedback-Aware Tuning of Recursive Q-Learning**.

The current feedback-aware (FA) implementation is distributed in
`simulation/feedback_aware/`. The repository also retains the earlier
v1.0.0 simulation material under `simulation/main/` and
`simulation/extensions/` so that the archived v1.0.0 release remains
traceable. Those legacy simulations are not used to generate the current
FA tables.

## Current FA simulation

The audited FA module reproduces the six numerical tables in the current
main paper and Supplementary Material. Its numerical core is unchanged from
the audited implementation; the packaging adds documentation, validation
evidence, and the supplemental Gaussian-audit entry point.

From `simulation/feedback_aware/`, use Python 3.11 or later:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python run_feedback_tuning.py smoke --workers 1 --output results/smoke
python run_feedback_tuning.py paper --workers 3 --output results/paper
python review/gaussian_audit.py --output results/supplement_gaussian
```

The paper run uses four local signals, sample sizes 250, 1000 and 4000,
and 5,000 independent datasets per cell. The log-RSS and quadratic
implementations use the same 60,000 datasets. The library contains the nine
penalty pairs in `{0.5,1,2}^2` plus the all-wide reference, with fixed
`T = 2` and uniform prior weights.

`CODE_TO_MANUSCRIPT.md` gives the exact script-to-table mapping.

## Simulated ADHD diagnostic

The diagnostic uses the simulated `adhd` data distributed with
`DTRlearn2` 1.1. The primary analysis is:

```bash
cd application/adhd
Rscript run_adhd_sel_ma.R main
```

The application is a propagation diagnostic, not a clinical validation of
the FA tuning algorithm.

## Validation

The FA module contains unit tests, saved reference outputs, independent
implementation checks, and a numerical audit of the Gaussian risk identity.
The production simulation was independently replayed during the packaging
audit and the manuscript table values were checked against the reproduced
outputs.

## Versions

- **v1.0.0**: initial repository release for the earlier
  prediction-risk-propagation / pretest--Stein manuscript.
- **v2.0.0**: current repository release for
  *Feedback-Aware Tuning of Recursive Q-Learning*.

The FA simulation module keeps its own internal numerical-core version
(`simulation/feedback_aware/VERSION`) so that repository-release versioning
is distinct from the implementation version used in the audit.

## Data and citation

No original participant-level clinical data are included. The ADHD example
uses a publicly distributed simulated dataset. Please use `CITATION.cff` for
software citation. Zenodo metadata are provided in `.zenodo.json` for the
v2.0.0 release.
