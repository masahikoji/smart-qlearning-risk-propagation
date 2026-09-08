# Prediction-risk propagation in SMART Q-learning

This repository contains the code used for the numerical work in the manuscript

**Prediction-Risk Propagation and Pretest--Stein Improvement in Q-Learning for Sequential Multiple Assignment Randomised Trials**.

The code is split into the main simulation study, the additional simulations for the finite-library and pretest--Stein results, and the simulated ADHD SMART illustration. The analysis code uses the same recursive pseudo-outcome construction for model selection and model averaging as described in the manuscript.

## Software

The simulation code uses base R. The ADHD analysis additionally requires the `DTRlearn2` package; the reported analysis was run with `DTRlearn2` 1.1 under R 4.3.0. The reference session information is saved in `application/adhd/sessionInfo_reference.txt`.

Python is not required for the analysis. It is used only by the independent reference checks in `simulation/main/tests/`.

## Repository layout

| Directory | Contents |
| --- | --- |
| `simulation/main/` | Two- and three-stage SMART simulations, generated-response diagnostics, treatment-boundary experiments, and exact error-decomposition checks |
| `simulation/extensions/` | Canonical pretest--Stein calculations, transported effective-rank simulations, criterion-feedback control, fixed finite-library AIC experiment, and the focused E2'' experiment |
| `application/adhd/` | Primary ADHD analysis, candidate-library sensitivity analysis, and exploratory broad-block pretest--Stein diagnostics |

Each directory has its own README with the run commands and output files.

## Quick validation

After installing R and, for the ADHD part, `DTRlearn2`, run

```bash
bash run_smoke_tests.sh
```

The Python reference checks can also be run independently:

```bash
cd simulation/main/tests
python3 reference_v2_1_validation.py
python3 reference_v2_2_decomposition_validation.py
python3 reference_v2_3_paired_decomposition.py
python3 reference_v2_4_regret_validation.py
python3 reference_primary_truth_quadrature.py
```

## Reproducing the paper runs

The production simulations are deliberately separated from the short smoke profiles so that pilot chunks cannot be mixed with the paper runs.

Main SMART simulations:

```bash
cd simulation/main
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

Additional simulations:

```bash
cd simulation/extensions
Rscript smoke_test.R
Rscript smoke_test_E2pp.R
Rscript run_additional_simulation.R all paper FALSE
Rscript run_E2pp_stein_simulation.R all paper FALSE
```

ADHD illustration:

```bash
cd application/adhd
Rscript run_adhd_sel_ma.R all 2000 8
Rscript run_adhd_candidate_library_sensitivity.R all
Rscript smoke_test_stein.R
Rscript run_adhd_stein_additional.R
```

The bootstrap core count in the first ADHD command may be changed for the local machine. The production simulation engines are chunked and can be restarted with the same command.

## Data

The illustrative data are the simulated `adhd` data set distributed with `DTRlearn2`. No original participant-level clinical data are included in this repository.

## Citation and archival metadata

`CITATION.cff` and `.zenodo.json` are included for a GitHub release archived through Zenodo. The software DOI and the final article citation can be added after the first Zenodo release is created.
