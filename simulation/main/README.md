# Main SMART simulations

This directory contains the two- and three-stage simulations used to study recursive AIC selection and Akaike weighting, generated-response propagation, treatment switching, and the exact risk decomposition.

The fitted procedures are ordinary backward least-squares Q-learning regressions. Numerical integration is used only in a known-truth scoring calculation for one three-stage target; it is not part of the fitted Q-learning procedure.

## Run order

```bash
Rscript smoke_test.R
Rscript calibrate_scenarios.R
Rscript run_population_checks.R

Rscript run_two_stage.R paper
Rscript run_three_stage.R paper

Rscript run_secondary.R 2 paper TRUE
Rscript run_secondary.R 3 paper TRUE

Rscript run_decomposition.R 2 paper "" TRUE
Rscript run_decomposition.R 3 paper "" TRUE

Rscript summarize.R
```

The paper profile uses 20,000 repeated SMART samples for the principal simulations. A shorter pilot can be run by setting `SMART_R`, for example

```bash
SMART_R=1000 Rscript run_two_stage.R paper
```

Compatible completed chunks are reused when a run is extended. Do not change the chunk size, evaluation-sample size, base seed, or scenario definition when extending an existing result set.

## Main outputs

Production results are written to `results_v2/`. The summary directory contains the files used most often in the manuscript:

- `risk_summary.csv`
- `risk_gap_summary.csv`
- `model_weight_summary.csv`
- `diagnostic_summary.csv`
- `secondary_summary.csv`
- `decomposition_summary.csv`
- `scenario_calibration.csv`
- `population_orthogonality_checks.csv`

The decomposition calculation uses paired conditionally independent descendant trajectories after exact one-step Rao--Blackwellisation of the local Bellman source. It is a diagnostic layer and does not change the fitted Q-learning estimator.

## Scenario families

The simulations include separated treatment gaps, smooth treatment boundaries, exact ties, weak and strong state transport, generated-response criterion perturbation, local terminal-signal grids, and fixed-misspecification robustness settings. Scenario definitions are in `R/00_config.R`.

## Validation

`smoke_test.R` checks the analytic moment formulas, information-criterion identities, representative recursive fits, the paired-descendant decomposition, and the regime-regret identity. Independent Python checks are in `tests/`.
