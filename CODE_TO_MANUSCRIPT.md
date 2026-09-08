# Code-to-manuscript map

This file lists the scripts that produce the numerical results discussed in the manuscript and Supplementary Material.

| Manuscript component | Script | Main output |
| --- | --- | --- |
| Recursive AIC selection versus Akaike weighting | `simulation/main/run_two_stage.R`, `simulation/main/run_three_stage.R` | `simulation/main/results_v2/summary/risk_gap_summary.csv` |
| Generated-response coefficient and criterion perturbation | `simulation/main/run_two_stage.R`, `simulation/main/run_three_stage.R` | `simulation/main/results_v2/summary/diagnostic_summary.csv` |
| Treatment switching and regime regret | `simulation/main/run_secondary.R` | `simulation/main/results_v2/summary/secondary_summary.csv` |
| Exact error-decomposition diagnostics | `simulation/main/run_decomposition.R` | `simulation/main/results_v2/summary/decomposition_summary.csv` |
| Population orthogonality benchmark | `simulation/main/run_population_checks.R` | `simulation/main/results_v2/summary/population_orthogonality_checks.csv` |
| Canonical pretest--Stein risk calculations | `simulation/extensions/run_additional_simulation.R normal paper` | `simulation/extensions/results/normal_mean/normal_mean_risk_grid.csv` |
| Transported effective-rank experiment | `simulation/extensions/run_additional_simulation.R transport paper` | `simulation/extensions/results/transport/transport_risk_summary.csv` |
| Criterion-feedback control | `simulation/extensions/run_additional_simulation.R feedback paper` | `simulation/extensions/results/feedback/feedback_summary.csv` |
| Fixed finite-library AIC experiment | `simulation/extensions/run_additional_simulation.R finite paper` | `simulation/extensions/results/finite_library/finite_library_summary.csv` |
| Focused E2'' simulation and Gaussian-limit benchmark | `simulation/extensions/run_E2pp_stein_simulation.R all paper FALSE` | `simulation/extensions/results/E2pp_stein/E2pp_main_table.csv` |
| Primary simulated ADHD analysis | `application/adhd/run_adhd_sel_ma.R main` | `application/adhd/results_adhd_sel_ma/` |
| ADHD candidate-library sensitivity | `application/adhd/run_adhd_candidate_library_sensitivity.R all` | `application/adhd/results_candidate_library_sensitivity/` |
| ADHD broad-block pretest--Stein sensitivity | `application/adhd/run_adhd_stein_additional.R` | `application/adhd/results_adhd_stein/` |

The production simulations use fixed seeds and write their session information or run metadata to the result directories. Short smoke or debug profiles are stored separately from the paper profiles.
