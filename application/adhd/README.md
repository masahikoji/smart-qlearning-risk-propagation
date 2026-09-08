# Simulated ADHD SMART illustration

This directory reproduces the application based on the simulated `DTRlearn2::adhd` data. The source design rerandomises nonresponders only, so the stage-2 analysis uses records with `r == 0`; responders retain their observed final outcome in the stage-1 pseudo-outcome.

The reported analysis was run with R 4.3.0 and `DTRlearn2` 1.1. The corresponding session information is included in `sessionInfo_reference.txt`.

Install the data package if needed:

```r
install.packages("DTRlearn2")
```

## Primary analysis

```bash
Rscript run_adhd_sel_ma.R main
```

The primary analysis compares recursive AIC selection and Akaike weighting, with AICc and BIC sensitivity calculations. Results are written to `results_adhd_sel_ma/`.

The nonparametric bootstrap stability diagnostic is run separately:

```bash
Rscript run_adhd_sel_ma.R bootstrap 2000 8
```

or together with the primary analysis:

```bash
Rscript run_adhd_sel_ma.R all 2000 8
```

## Candidate-library sensitivity

```bash
Rscript run_adhd_candidate_library_sensitivity.R all
```

This analysis keeps the source-aligned Q-learning structure while varying the candidate effect modifiers and functional form. Results are written to `results_candidate_library_sensitivity/`.

## Broad-block pretest--Stein sensitivity

The primary M22-versus-M23 comparison adds only one parameter and therefore lies outside the `d >= 3` pretest--Stein result. The additional analysis uses only broad nested pairs already present in the candidate-library sensitivity study.

Run

```bash
Rscript smoke_test_stein.R
Rscript run_adhd_stein_additional.R
```

Results are written to `results_adhd_stein/`. The reported effective rank is a same-data diagnostic and is not used to choose the pretest--Stein tuning constant.
