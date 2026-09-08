#!/usr/bin/env Rscript
ca <- commandArgs(trailingOnly = FALSE)
hit <- grep("^--file=", ca, value = TRUE)
root <- if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)) else normalizePath(getwd())
source(file.path(root, "R", "adhd_sel_ma_functions.R"))
source(file.path(root, "R", "adhd_candidate_sensitivity_functions.R"))
source(file.path(root, "R", "adhd_stein_additional_functions.R"))
dat <- load_adhd_data()
centers <- compute_centers(dat)
dat_c <- apply_centers(dat, centers)
pairs <- make_block_pairs()
expected_d <- c(expanded_E22_E29 = 6L, expanded_E20_E29 = 7L,
                spline_N22_N26 = 7L, spline_N20_N26 = 8L)
# These likelihood-ratio values are implied by the AIC values in the validated
# candidate-library sensitivity output supplied with the existing analysis.
expected_lambda <- c(
  expanded_E22_E29 = 11.917915781814,
  expanded_E20_E29 = 47.178047038176,
  spline_N22_N26 = 12.663717621009,
  spline_N20_N26 = 51.320981460347
)
for (nm in names(pairs)) {
  nr <- stage2_eligible(dat_c)
  bf <- fit_nested_block(dat_c[nr, , drop = FALSE], pairs[[nm]]$narrow, pairs[[nm]]$wide)
  stopifnot(bf$d == expected_d[[nm]])
  stopifnot(abs(bf$aic_identity_error) < 1e-8)
  stopifnot(abs(bf$lambda - expected_lambda[[nm]]) < 1e-5)
  z <- run_one_adhd_block(dat_c, nm, pairs[[nm]], a_values = c(1, 2))
  stopifnot(all(is.finite(z$summary$lambda)), all(is.finite(z$summary$wide_weight)))
  stopifnot(all(z$summary$wide_weight >= 0 & z$summary$wide_weight <= 1))
}
cat("All ADHD Pretest--Stein additional-analysis smoke checks passed.\n")
