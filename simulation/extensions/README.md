# Additional simulations

This directory contains the simulations added for the pretest--Stein and fixed finite-library results. The code uses base R and is independent of the main simulation engine.

## General extension study

`run_additional_simulation.R` has four modes:

- `normal`: canonical Gaussian risk calculations for hard AIC selection, Akaike weighting, and pretest--Stein rules;
- `transport`: two-stage SMART simulations with a population transported Gram matrix of prespecified effective rank;
- `feedback`: a criterion-feedback control outside the simple effective-rank preservation theorem;
- `finite`: a four-model local AIC experiment for the fixed finite-library result.

Run

```bash
Rscript smoke_test.R
Rscript run_additional_simulation.R all paper FALSE
```

Production outputs are written under `results/`. The main summary files are

```text
results/normal_mean/normal_mean_risk_grid.csv
results/transport/transport_risk_summary.csv
results/transport/transport_limit_theory.csv
results/transport/transport_pairwise_gaps.csv
results/feedback/feedback_summary.csv
results/finite_library/finite_library_summary.csv
results/finite_library/finite_library_limit_theory.csv
```

The long simulations are chunked and restartable. The runner stores source-code checksums with the chunks and stops if incompatible chunks are found.

## Focused E2'' experiment

`run_E2pp_stein_simulation.R` implements the Gaussian two-stage construction used to verify the terminal pretest--Stein result and its survival under multidirectional versus collinear backward transport.

The paper run uses `n=500`, 20,000 paired finite-sample replicates at each signal value, and 2,000,000 Gaussian-limit draws. Run

```bash
Rscript smoke_test_E2pp.R
Rscript run_E2pp_stein_simulation.R all paper FALSE
```

The main outputs are written to `results/E2pp_stein/`, including

```text
E2pp_main_table.csv
E2pp_pairwise_gaps.csv
E2pp_risk_summary.csv
E2pp_limit_gaps.csv
E2pp_limit_risks.csv
E2pp_transport_validation.csv
```

Parallelism can be changed without editing the code, for example

```bash
E2PP_CORES=8 Rscript run_E2pp_stein_simulation.R all paper FALSE
```
