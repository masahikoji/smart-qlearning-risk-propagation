# Code-to-results map

This file maps scientific result types to the code and archived outputs that generate or validate them. It deliberately avoids journal-specific table and section numbers.

| Scientific result | Primary code / source | Main archived output |
| --- | --- | --- |
| Two-stage finite-sample recursive risk comparisons | `simulation/main/run_two_stage.R` | simulation summary outputs under `simulation/main/` results directories |
| Three-stage finite-sample recursive risk comparisons | `simulation/main/run_three_stage.R` | simulation summary outputs under `simulation/main/` results directories |
| Cross-stage linear/nonlinear decomposition | `simulation/main/run_decomposition.R` | decomposition summaries and validation outputs |
| Consolidated simulation summaries | `simulation/main/summarize.R` | consolidated CSV summaries |
| Recursive Gaussian local-limit calculation | `reproducibility/compute_recursive_limit.py` | archived recursive-limit CSV and raw output |
| Smooth-boundary exact conditional-integration validation | `reproducibility/run_table_s7_exact_aligned.py` or its release-equivalent path | full-precision CSV, metadata JSON, and archived replicate-level arrays |
| Smooth-boundary theoretical constants | theoretical-constant script in `reproducibility/` | archived constant summary |
| Exact-tie nonlinear contribution validation | exact conditional-integration code in `reproducibility/` | archived exact-tie summary output |
| Criterion-feedback / coefficient-perturbation diagnostics | `simulation/main/` diagnostic routines | diagnostic summary CSV |
| Regime-value regret | current `simulation/main/` implementation | secondary-outcome summary CSV |
| Simulated ADHD primary recursive selection/averaging analysis | `application/adhd/run_adhd_sel_ma.R main` | `application/adhd/results_adhd_sel_ma/` |
| ADHD candidate-library sensitivity | `application/adhd/run_adhd_candidate_library_sensitivity.R all` | `application/adhd/results_candidate_library_sensitivity/` |
| Risk-ordering / backward-reversal figure | `manuscript_assets/make_risk_reversal_figure.py` | vector figure PDF |

## Release rule

Before a release is archived, every result explicitly described as being available in the reproducibility archive should have a corresponding file in the release. In particular, if the research report states that full-precision estimates, metadata, or replicate-level arrays are archived, those files must be included in the Zenodo release even if they are not practical to keep in the Git repository itself.
