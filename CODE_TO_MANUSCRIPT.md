# Code-to-manuscript map

This file maps the current **Feedback-Aware Tuning of Recursive Q-Learning** manuscript and Supplementary Material to the code used for each numerical result.

| Manuscript component | Script / source | Main output |
| --- | --- | --- |
| Main Table 1: Gaussian risks and expected risk estimates | `simulation/feedback_aware/run_feedback_tuning.py paper` | `results/paper/gaussian/gaussian_results.csv` |
| Main Table 2: paired Akaike-minus-FA risk differences | `simulation/feedback_aware/run_feedback_tuning.py paper` | `results/paper/logrss/main_table_gaps.csv` |
| Supplement Table S1: finite-sample full risks and reported risk estimates | `simulation/feedback_aware/run_feedback_tuning.py paper` | `results/paper/logrss/smart_results.csv` |
| Supplement Table S2: log-RSS and quadratic criterion variants | `simulation/feedback_aware/run_feedback_tuning.py paper` | `results/paper/logrss/` and `results/paper/wald/` |
| Supplement Table S3: adaptation-term audit | `simulation/feedback_aware/review/gaussian_audit.py` | `results/supplement_gaussian/table_s3_adaptation.csv` |
| Supplement Table S4: fixed-`T` sensitivity | `simulation/feedback_aware/review/gaussian_audit.py` | `results/supplement_gaussian/table_s4_tuning.csv` |
| Primary simulated ADHD propagation diagnostic | `application/adhd/run_adhd_sel_ma.R main` | `application/adhd/results_adhd_sel_ma/` |
| ADHD candidate-library sensitivity | `application/adhd/run_adhd_candidate_library_sensitivity.R all` | `application/adhd/results_candidate_library_sensitivity/` |
| ADHD broad-block exploratory sensitivity | `application/adhd/run_adhd_stein_additional.R` | `application/adhd/results_adhd_stein/` |

The current FA tables are generated only by `simulation/feedback_aware/`. The folders `simulation/main/` and `simulation/extensions/` are retained as legacy v1.0.0 materials and are not used to validate the current FA method.
