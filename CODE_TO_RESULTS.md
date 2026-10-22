# Code-to-results map

This file maps scientific result types to the code and archived outputs that generate or validate them.
It deliberately avoids journal-specific table and section numbers.

| Scientific result | Primary code / source | Main archived output |
| --- | --- | --- |
| Two-stage finite-sample recursive risk comparisons | `simulation/main/run_two_stage.R` | simulation summary outputs under `simulation/main/` results directories |
| Three-stage finite-sample recursive risk comparisons | `simulation/main/run_three_stage.R` | simulation summary outputs under `simulation/main/` results directories |
| Cross-stage linear/nonlinear decomposition | `simulation/main/run_decomposition.R` | decomposition summaries and validation outputs |
| Consolidated simulation summaries | `simulation/main/summarize.R` | consolidated CSV summaries |
| Recursive Gaussian local-limit calculation | `reproducibility/compute_recursive_limit.py` | `reproducibility/recursive_limit_validation.csv` and raw text output |
| Smooth-boundary and exact-tie conditional-integration validation | `reproducibility/run_boundary_validation_exact.py` | `reproducibility/boundary_validation.csv` and replicate-level NPZ arrays |
| Smooth-boundary theoretical constants | `reproducibility/compute_boundary_theory_constants.py` | `reproducibility/boundary_theory_constants.json` |
| Criterion-feedback / coefficient-perturbation diagnostics | `simulation/main/` diagnostic routines | diagnostic summary CSV |
| Regime-value regret | current `simulation/main/` implementation | secondary-outcome summary CSV |
| Simulated ADHD primary recursive selection/averaging analysis | `application/adhd/run_adhd_sel_ma.R main` | `application/adhd/results_adhd_sel_ma/` |
| ADHD candidate-library sensitivity | `application/adhd/run_adhd_candidate_library_sensitivity.R all` | `application/adhd/results_candidate_library_sensitivity/` |
| Risk-ordering / backward-reversal figure | `manuscript_assets/make_risk_reversal_figure.py` | vector figure PDF |

## Release rule

Before a release is archived, every result explicitly described as being available in the
reproducibility archive should have a corresponding file in the release.
