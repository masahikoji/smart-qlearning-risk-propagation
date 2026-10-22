# Reproducibility materials

This directory contains journal-neutral numerical validation programs, archived outputs,
metadata, and replicate-level arrays supporting the recursive Q-learning risk-propagation analyses.

## Recursive Gaussian local-limit calculation

- `compute_recursive_limit.py`
- `recursive_limit_validation.csv`
- `recursive_limit_validation_raw.txt`

These files contain the independent finite-dimensional Gaussian recursion used to compare
finite-sample risk gaps with their recursive local limits.

## Smooth-boundary and exact-tie validation

The current release-level boundary-validation program is:

- `run_boundary_validation_exact.py`

Archived full-precision summary:

- `boundary_validation.csv`

Replicate-level arrays:

- `smooth_boundary_raw_n250.npz`
- `smooth_boundary_raw_n500.npz`
- `smooth_boundary_raw_n1000.npz`
- `exact_tie_raw_n250.npz`
- `exact_tie_raw_n500.npz`
- `exact_tie_raw_n1000.npz`

The three smooth-boundary arrays are the original completed-run arrays and retain the
recorded final-run SHA-256 identities. The exact-tie arrays were deterministically regenerated
from the final boundary-validation program using the same seed schedule and numerical settings;
their means and Monte Carlo standard errors agree with the archived full-precision summary.
See `boundary_validation_metadata.json` for provenance and hashes.

## Smooth-boundary limiting constants

- `compute_boundary_theory_constants.py`
- `boundary_theory_constants.json`

The script implements the displayed one-dimensional expectation for the smooth-boundary
limit and reproduces the rounded constants 5.920 for AIC selection and 4.527 for Akaike
model averaging.

## Other numerical code

The principal two- and three-stage simulations and decomposition code are maintained under
`../simulation/main/`. The simulated ADHD analysis is maintained under `../application/adhd/`.
