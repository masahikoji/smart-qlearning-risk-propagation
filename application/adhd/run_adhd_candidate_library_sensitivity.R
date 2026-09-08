#!/usr/bin/env Rscript


args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  script_dir <- dirname(script_path)
} else {
  script_dir <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(script_dir)

args <- commandArgs(trailingOnly = TRUE)
criterion_arg <- if (length(args) >= 1L) args[[1]] else "all"
allowed <- c("all", "AIC", "AICc", "BIC")
if (!criterion_arg %in% allowed) {
  stop("First argument must be one of: all, AIC, AICc, BIC")
}
criteria <- if (criterion_arg == "all") c("AIC", "AICc", "BIC") else criterion_arg

source(file.path("R", "adhd_sel_ma_functions.R"))
source(file.path("R", "adhd_candidate_sensitivity_functions.R"))

outdir <- file.path(script_dir, "results_candidate_library_sensitivity")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat <- load_adhd_data()
elig <- stage2_eligible(dat)
resp <- dat$r == 1

cat("ADHD SMART candidate-library sensitivity analysis\n")
cat("================================================\n")
cat("Program directory:", script_dir, "\n")
cat("Output directory :", outdir, "\n")
cat("Criteria         :", paste(criteria, collapse = ", "), "\n")
cat("N total          :", nrow(dat), "\n")
cat("N stage-2 eligible nonresponders:", sum(elig), "\n")
cat("N responders     :", sum(resp), "\n\n")

# Source-aligned coding/reference check inherited from the primary program.
centers <- compute_centers(dat)
write.csv(
  data.frame(variable = names(centers), center = as.numeric(centers)),
  file.path(outdir, "centering_constants.csv"), row.names = FALSE
)

ref <- reference_reproduction(dat)
write.csv(ref, file.path(outdir, "d3c_reference_reproduction.csv"), row.names = FALSE)
max_ref_diff <- attr(ref, "max_abs_difference")
if (!is.finite(max_ref_diff)) stop("Reference-reproduction check returned a non-finite discrepancy.")
if (max_ref_diff > 0.03) {
  warning(
    "Source-aligned full-model coefficient discrepancy exceeds 0.03 (max = ",
    signif(max_ref_diff, 6), "). Inspect data/package version before interpreting sensitivity results."
  )
}

# Save exact prespecified library definitions before fitting anything.
def_table <- library_definition_table()
write.csv(def_table, file.path(outdir, "candidate_library_definitions.csv"), row.names = FALSE)

libs <- adhd_sensitivity_formula_libraries(response_stage2 = "y", response_stage1 = "pseudo")

all_results <- lapply(names(libs), function(variant) {
  vv <- libs[[variant]]
  rr <- setNames(lapply(criteria, function(cr) {
    cat("Running", variant, "-", cr, "...\n")
    run_recursive_pair_custom(
      train = dat,
      eval = dat,
      criterion = cr,
      stage2_formulas = vv$stage2,
      stage1_formulas = vv$stage1,
      analysis_variant = variant
    )
  }), criteria)
  rr
})
names(all_results) <- names(libs)

# Collect candidate-level IC weights, recursive coefficients, propagation metrics,
# decision summaries, margin diagnostics, and subject-level outputs.
ic_rows <- list()
coef_rows <- list()
prop_rows <- list()
decision_rows <- list()
margin_rows <- list()
subject_rows <- list()
fixed_coef_rows <- list()
fixed_diag_rows <- list()
small_rows <- list()
idx <- 1L

for (variant in names(all_results)) {
  for (cr in names(all_results[[variant]])) {
    res <- all_results[[variant]][[cr]]

    ic_rows[[idx]] <- sensitivity_ic_table(res)

    cc <- model_coefficients_table(res)
    cc$analysis_variant <- variant
    coef_rows[[idx]] <- cc[, c("analysis_variant", setdiff(names(cc), "analysis_variant"))]

    prop_rows[[idx]] <- sensitivity_propagation_summary(res, dat)

    d1 <- decision_summary(
      res$stage1_sel$actions, res$stage1_ma$actions,
      ids = dat$id,
      strata = paste0("o13=", dat$o13),
      criterion = cr, stage = 1L
    )
    d1$analysis_variant <- variant

    d2 <- decision_summary(
      res$stage2_sel$stage2_eval[elig, , drop = FALSE],
      res$stage2_ma$stage2_eval[elig, , drop = FALSE],
      ids = dat$id[elig],
      strata = paste0("a1=", dat$a1[elig], ",o22=", dat$o22[elig]),
      criterion = cr, stage = 2L
    )
    d2$analysis_variant <- variant
    decision_rows[[idx]] <- rbind(d1, d2)

    m1 <- decision_margin_diagnostics(
      res$stage1_sel$actions, res$stage1_ma$actions,
      criterion = cr, stage = 1L, analysis_variant = variant
    )
    m1$comparison <- "recursive_SEL_vs_MA"

    fd <- fixed_stage1_diagnostic(res)
    mf <- decision_margin_diagnostics(
      fd$action_sel, fd$action_ma,
      criterion = cr, stage = 1L, analysis_variant = variant
    )
    mf$comparison <- "fixed_primary_M11_downstream_only"
    margin_rows[[idx]] <- rbind(m1, mf)

    subject_rows[[idx]] <- sensitivity_subject_table(res, dat)

    fc <- fd$coef_table
    fc$analysis_variant <- variant
    fc$criterion <- cr
    fixed_coef_rows[[idx]] <- fc[, c("analysis_variant", "criterion", setdiff(names(fc), c("analysis_variant", "criterion")))]

    fdiag <- fd$diagnostics
    fdiag$analysis_variant <- variant
    fixed_diag_rows[[idx]] <- fdiag[, c("analysis_variant", setdiff(names(fdiag), "analysis_variant"))]

    # Candidate-dimensional diagnostics: one row per stage/candidate/target fit.
    it <- sensitivity_ic_table(res)
    ss <- unique(it[, c(
      "analysis_variant", "stage", "procedure_target", "model", "n",
      "p_beta", "k_ic", "residual_df", "aicc_correction"
    )])
    ss$p_beta_over_n <- ss$p_beta / ss$n
    ss$residual_df_fraction <- ss$residual_df / ss$n
    small_rows[[idx]] <- ss

    idx <- idx + 1L
  }
}

ic_table <- do.call(rbind, ic_rows)
coef_table <- do.call(rbind, coef_rows)
prop_table <- do.call(rbind, prop_rows)
decision_table <- do.call(rbind, decision_rows)
margin_table <- do.call(rbind, margin_rows)
subject_table <- do.call(rbind, subject_rows)
fixed_coef_table <- do.call(rbind, fixed_coef_rows)
fixed_diag_table <- do.call(rbind, fixed_diag_rows)
small_table <- do.call(rbind, small_rows)

# Stable ordering for easy manuscript checking.
prop_table <- prop_table[order(prop_table$criterion, prop_table$analysis_variant), ]
ic_table <- ic_table[order(ic_table$criterion, ic_table$analysis_variant, -ic_table$stage, ic_table$procedure_target, ic_table$model), ]
decision_table <- decision_table[order(decision_table$criterion, decision_table$analysis_variant, decision_table$stage, decision_table$stratum), ]
margin_table <- margin_table[order(margin_table$criterion, margin_table$analysis_variant, margin_table$comparison), ]

write.csv(ic_table, file.path(outdir, "sensitivity_candidate_model_ic_weights.csv"), row.names = FALSE)
write.csv(coef_table, file.path(outdir, "sensitivity_recursive_SEL_MA_coefficients.csv"), row.names = FALSE)
write.csv(prop_table, file.path(outdir, "sensitivity_propagation_summary.csv"), row.names = FALSE)
write.csv(decision_table, file.path(outdir, "sensitivity_decision_summary.csv"), row.names = FALSE)
write.csv(margin_table, file.path(outdir, "sensitivity_decision_margin_diagnostics.csv"), row.names = FALSE)
write.csv(subject_table, file.path(outdir, "sensitivity_subject_level_stage1.csv"), row.names = FALSE)
write.csv(fixed_coef_table, file.path(outdir, "sensitivity_fixed_M11_coefficients.csv"), row.names = FALSE)
write.csv(fixed_diag_table, file.path(outdir, "sensitivity_fixed_M11_generated_response_diagnostics.csv"), row.names = FALSE)
write.csv(small_table, file.path(outdir, "sensitivity_small_sample_diagnostics.csv"), row.names = FALSE)

# Validate that the sensitivity program reproduces the current primary AIC
# analysis before changing the candidate library.
if ("AIC" %in% criteria) {
  val <- primary_reference_validation(all_results$primary_reference$AIC)
  write.csv(val, file.path(outdir, "primary_AIC_reproduction_validation.csv"), row.names = FALSE)
  if (!all(val$pass)) {
    warning(
      "One or more primary-AIC reproduction checks fell outside the prespecified tolerance. ",
      "Inspect primary_AIC_reproduction_validation.csv before interpreting sensitivity results."
    )
  }
}

# Save complete fit objects for later inspection without rerunning.
saveRDS(
  list(
    results = all_results,
    library_definitions = def_table,
    criteria = criteria,
    source_reference_max_abs_difference = max_ref_diff
  ),
  file.path(outdir, "sensitivity_full_results.rds")
)

write_sensitivity_plots(all_results, outdir)

# Concise human-readable report.
sink(file.path(outdir, "sensitivity_analysis_summary.txt"))
cat("ADHD SMART candidate-library sensitivity analysis\n")
cat("================================================\n\n")
cat("Data: DTRlearn2::adhd; source-aligned stage-2 eligibility r == 0.\n")
cat("N total:", nrow(dat), "\n")
cat("N stage-2 eligible:", sum(elig), "\n")
cat("Source-reference reproduction max absolute coefficient discrepancy:", max_ref_diff, "\n\n")
cat("Sensitivity libraries\n")
cat("---------------------\n")
for (variant in names(libs)) {
  cat("*", variant, "\n  ", libs[[variant]]$description, "\n")
  cat("  Stage-2 candidates:", length(libs[[variant]]$stage2), "\n")
  cat("  Stage-1 candidates:", length(libs[[variant]]$stage1), "\n")
}
cat("\nCore interpretation rule\n")
cat("------------------------\n")
cat("The sensitivity analysis is descriptive. Robustness is assessed by whether the\n")
cat("generated-response pathway persists across prespecified candidate libraries:\n")
cat("stage-2 model uncertainty -> stage-1 pseudo-outcome difference -> fixed-M11\n")
cat("coefficient/LR perturbation -> full recursive stage-1 decision comparison.\n")
cat("A disappearance of the pseudo-outcome difference when one model receives nearly\n")
cat("all information-criterion weight is also theoretically coherent and should not be\n")
cat("treated as a failed sensitivity analysis.\n\n")

for (cr in criteria) {
  cat("================================================\n")
  cat("Criterion:", cr, "\n")
  cat("================================================\n")
  ss <- prop_table[prop_table$criterion == cr, , drop = FALSE]
  for (i in seq_len(nrow(ss))) {
    x <- ss[i, ]
    cat("\nVariant:", x$analysis_variant, "\n")
    cat("Stage 2 selected:", x$stage2_selected_model,
        "; top MA model:", x$stage2_top_MA_model,
        "; top weight:", signif(x$stage2_top_weight, 5),
        "; effective models:", signif(x$stage2_effective_models, 5), "\n")
    cat("Pseudo-outcome RMS SEL-MA:", signif(x$pseudo_diff_rms, 6),
        "; max |difference|:", signif(x$pseudo_diff_max_abs, 6), "\n")
    cat("Fixed primary M11 sqrt(n) coefficient shift:", signif(x$fixed_M11_sqrt_n_gamma_shift, 6),
        "; LR shift:", signif(x$fixed_M11_LR_shift, 6), "\n")
    cat("Stage 1 SEL selected:", x$stage1_SEL_selected_model,
        "; top MA model:", x$stage1_MA_top_model,
        "; recursive discordance:", x$stage1_recursive_discordance_n, "/", nrow(dat), "\n")
    cat("Stage 1 min |contrast| SEL:", signif(x$stage1_min_abs_contrast_SEL, 6),
        "; max |SEL-MA contrast difference|:",
        signif(x$stage1_max_abs_SEL_MA_contrast_difference, 6), "\n")
  }
  cat("\n")
}
sink()

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nSensitivity analysis complete.\n")
cat("Key summary:", file.path(outdir, "sensitivity_propagation_summary.csv"), "\n")
cat("Text report:", file.path(outdir, "sensitivity_analysis_summary.txt"), "\n")
