#!/usr/bin/env Rscript

ca <- commandArgs(trailingOnly = FALSE)
hit <- grep("^--file=", ca, value = TRUE)
root <- if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)) else normalizePath(getwd())
setwd(root)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1]] else "main"
boot_B <- if (length(args) >= 2L) as.integer(args[[2]]) else 2000L
nphys <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.na(nphys) || nphys < 1L) nphys <- 1L
boot_cores <- if (length(args) >= 3L) as.integer(args[[3]]) else max(1L, min(8L, nphys))

if (!mode %in% c("main", "bootstrap", "all")) {
  stop("mode must be one of: main, bootstrap, all")
}
if (is.na(boot_B) || boot_B < 1L) stop("bootstrap B must be a positive integer.")
if (is.na(boot_cores) || boot_cores < 1L) stop("bootstrap cores must be a positive integer.")

source(file.path(root, "R", "adhd_sel_ma_functions.R"))

outdir <- file.path(root, "results_adhd_sel_ma")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat <- load_adhd_data()
elig <- stage2_eligible(dat)
resp <- dat$r == 1

cat("Loaded DTRlearn2::adhd\n")
cat("N total:", nrow(dat), "\n")
cat("N responders:", sum(resp), "\n")
cat("N source-aligned stage-2 eligible nonresponders:", sum(elig), "\n")
cat("Note: DTRlearn2::adhd contains coded a2 values for responders; source-aligned analysis ignores them.\n")

# Basic data summary and the DTRlearn2 coding/source-design distinction.
centers <- compute_centers(dat)
write.csv(data.frame(variable = names(centers), center = as.numeric(centers)),
          file.path(outdir, "centering_constants.csv"), row.names = FALSE)

data_summary <- data.frame(
  item = c(
    "N_total", "N_stage2_eligible_nonresponders", "N_responders",
    "A1_minus1", "A1_plus1",
    "A2_minus1_all_rows", "A2_plus1_all_rows",
    "A2_minus1_stage2_eligible", "A2_plus1_stage2_eligible",
    "A2_minus1_responders_ignored", "A2_plus1_responders_ignored",
    "mean_y", "sd_y"
  ),
  value = c(
    nrow(dat), sum(elig), sum(resp),
    sum(dat$a1 == -1, na.rm = TRUE), sum(dat$a1 == 1, na.rm = TRUE),
    sum(dat$a2 == -1, na.rm = TRUE), sum(dat$a2 == 1, na.rm = TRUE),
    sum(dat$a2[elig] == -1), sum(dat$a2[elig] == 1),
    sum(dat$a2[resp] == -1, na.rm = TRUE), sum(dat$a2[resp] == 1, na.rm = TRUE),
    mean(dat$y, na.rm = TRUE), sd(dat$y, na.rm = TRUE)
  )
)
write.csv(data_summary, file.path(outdir, "data_summary.csv"), row.names = FALSE)

sink(file.path(outdir, "data_structure_note.txt"))
cat("Source-aligned ADHD SMART analysis\n")
cat("==================================\n\n")
cat("DTRlearn2::adhd contains a2 coded as -1/+1 for all 150 rows, including responders.\n")
cat("The source d3center ADHD SMART / PROC QLEARN analysis treats only nonresponders as stage-2 randomized/eligible.\n")
cat("This program therefore uses stage-2 eligibility r == 0 and observed a2, and ignores responder a2 values.\n")
cat("For responders, the stage-1 pseudo-outcome is the observed final y.\n\n")
cat("Counts in this data set:\n")
print(data_summary)
sink()

# Reproduce the source-aligned full Q-learning model before doing SEL/MA.
ref <- reference_reproduction(dat)
write.csv(ref, file.path(outdir, "d3c_reference_reproduction.csv"), row.names = FALSE)
max_ref_diff <- attr(ref, "max_abs_difference")
cat("Max absolute difference from published d3center coefficient references:",
    format(max_ref_diff, digits = 6), "\n")
if (is.finite(max_ref_diff) && max_ref_diff > 0.03) {
  warning("The source-aligned full-model coefficients differ by >0.03 from the d3center reference. Inspect before interpretation.")
}

criteria <- c("AIC", "BIC", "AICc")
variant_defs <- list(
  primary_four_model = list(
    library = "four_model",
    description = "Primary prespecified four-model stage-2 library M20/M21/M22/M23"
  ),
  nested_M22_M23 = list(
    library = "nested_M22_M23",
    description = "Theory-oriented one-df nested reduction: M22 versus M23"
  ),
  strong_M20_M23 = list(
    library = "strong_M20_M23",
    description = "Strong-separation sensitivity: M20 versus M23"
  )
)

all_results <- lapply(variant_defs, function(v) {
  setNames(lapply(criteria, function(cr) {
    run_recursive_pair(dat, dat, criterion = cr, stage2_library = v$library)
  }), criteria)
})

primary <- all_results$primary_four_model

# Combined outputs across the primary and both prespecified reductions.
all_ic <- list()
all_coef <- list()
all_decision <- list()
all_subject <- list()
all_fixed_coef <- list()
all_fixed_diag <- list()
all_margin <- list()
all_nested <- list()

for (variant in names(all_results)) {
  for (cr in criteria) {
    z <- all_results[[variant]][[cr]]

    t2 <- z$stage2_lib$table
    t2$analysis_variant <- variant
    t2$stage <- 2L
    t2$procedure_target <- "common_terminal_fit"

    t1s <- z$stage1_sel$lib$table
    t1s$analysis_variant <- variant
    t1s$stage <- 1L
    t1s$procedure_target <- paste0(cr, "_SEL_pseudooutcome")

    t1m <- z$stage1_ma$lib$table
    t1m$analysis_variant <- variant
    t1m$stage <- 1L
    t1m$procedure_target <- paste0(cr, "_MA_pseudooutcome")

    all_ic[[paste(variant, cr, sep = "_")]] <- rbind(t2, t1s, t1m)

    ct <- model_coefficients_table(z)
    ct$analysis_variant <- variant
    all_coef[[paste(variant, cr, sep = "_")]] <- ct[, c("analysis_variant", setdiff(names(ct), "analysis_variant"))]

    ds1 <- decision_summary(
      z$stage1_sel$actions, z$stage1_ma$actions,
      ids = dat$id, strata = paste0("o13=", dat$o13), criterion = cr, stage = 1L
    )
    ds1$analysis_variant <- variant

    ds2 <- decision_summary(
      z$stage2_sel$stage2_eval[elig, ], z$stage2_ma$stage2_eval[elig, ],
      ids = dat$id[elig], strata = paste0("a1=", dat$a1[elig], ",o22=", dat$o22[elig]),
      criterion = cr, stage = 2L
    )
    ds2$analysis_variant <- variant
    all_decision[[paste(variant, cr, sep = "_")]] <- rbind(ds1, ds2)

    fz <- fixed_stage1_diagnostic(z)
    ft <- fz$coef_table
    ft$analysis_variant <- variant
    ft$criterion <- cr
    all_fixed_coef[[paste(variant, cr, sep = "_")]] <-
      ft[, c("analysis_variant", "criterion", setdiff(names(ft), c("analysis_variant", "criterion")))]

    fd <- fz$diagnostics
    fd$analysis_variant <- variant
    all_fixed_diag[[paste(variant, cr, sep = "_")]] <-
      fd[, c("analysis_variant", setdiff(names(fd), "analysis_variant"))]

    all_subject[[paste(variant, cr, sep = "_")]] <- data.frame(
      id = dat$id,
      analysis_variant = variant,
      criterion = cr,
      o11 = dat$o11,
      o12 = dat$o12,
      o13 = dat$o13,
      o14 = dat$o14,
      a1_observed = dat$a1,
      response = dat$r,
      stage1_Q_minus1_SEL = z$stage1_sel$actions$q_minus,
      stage1_Q_plus1_SEL = z$stage1_sel$actions$q_plus,
      stage1_contrast_SEL = z$stage1_sel$actions$contrast,
      stage1_recommendation_SEL = ifelse(z$stage1_sel$actions$action == 1, "BMOD", "MED"),
      stage1_Q_minus1_MA = z$stage1_ma$actions$q_minus,
      stage1_Q_plus1_MA = z$stage1_ma$actions$q_plus,
      stage1_contrast_MA = z$stage1_ma$actions$contrast,
      stage1_recommendation_MA = ifelse(z$stage1_ma$actions$action == 1, "BMOD", "MED"),
      stage1_SEL_MA_discordant = z$stage1_sel$actions$action != z$stage1_ma$actions$action,
      pseudooutcome_stage2_SEL = z$stage2_sel$pseudo_eval,
      pseudooutcome_stage2_MA = z$stage2_ma$pseudo_eval,
      pseudooutcome_SEL_minus_MA = z$stage2_sel$pseudo_eval - z$stage2_ma$pseudo_eval,
      fixed_stage1_contrast_from_stage2_SEL = fz$action_sel$contrast,
      fixed_stage1_contrast_from_stage2_MA = fz$action_ma$contrast,
      fixed_stage1_recommendation_from_stage2_SEL = ifelse(fz$action_sel$action == 1, "BMOD", "MED"),
      fixed_stage1_recommendation_from_stage2_MA = ifelse(fz$action_ma$action == 1, "BMOD", "MED"),
      fixed_stage1_downstream_only_discordant = fz$action_sel$action != fz$action_ma$action,
      stringsAsFactors = FALSE
    )

    m1 <- decision_margin_diagnostics(
      z$stage1_sel$actions, z$stage1_ma$actions,
      criterion = cr, stage = 1L, analysis_variant = variant
    )
    m1$comparison <- "recursive_SEL_vs_MA"
    mfix <- decision_margin_diagnostics(
      fz$action_sel, fz$action_ma,
      criterion = cr, stage = 1L, analysis_variant = variant
    )
    mfix$comparison <- "fixed_stage1_downstream_only"
    all_margin[[paste(variant, cr, sep = "_")]] <- rbind(m1, mfix)

    if (variant == "nested_M22_M23") {
      all_nested[[paste(variant, cr, sep = "_")]] <-
        nested_stage2_diagnostic(z, "M22", "M23", analysis_variant = variant)
    }
    if (variant == "strong_M20_M23") {
      all_nested[[paste(variant, cr, sep = "_")]] <-
        nested_stage2_diagnostic(z, "M20", "M23", analysis_variant = variant)
    }
  }
}

ic_table <- do.call(rbind, all_ic)
coef_table <- do.call(rbind, all_coef)
decision_table <- do.call(rbind, all_decision)
subject_table <- do.call(rbind, all_subject)
fixed_coef_table <- do.call(rbind, all_fixed_coef)
fixed_diag_table <- do.call(rbind, all_fixed_diag)
margin_table <- do.call(rbind, all_margin)
nested_table <- do.call(rbind, all_nested)

write.csv(ic_table, file.path(outdir, "all_variants_candidate_model_ic_weights.csv"), row.names = FALSE)
write.csv(coef_table, file.path(outdir, "all_variants_recursive_SEL_MA_coefficients.csv"), row.names = FALSE)
write.csv(decision_table, file.path(outdir, "all_variants_decision_summary.csv"), row.names = FALSE)
write.csv(subject_table, file.path(outdir, "all_variants_subject_level_stage1_decisions.csv"), row.names = FALSE)
write.csv(fixed_coef_table, file.path(outdir, "all_variants_downstream_only_fixed_stage1_coefficients.csv"), row.names = FALSE)
write.csv(fixed_diag_table, file.path(outdir, "all_variants_downstream_only_generated_response_diagnostics.csv"), row.names = FALSE)
write.csv(margin_table, file.path(outdir, "all_variants_decision_margin_diagnostics.csv"), row.names = FALSE)
write.csv(nested_table, file.path(outdir, "nested_stage2_comparison_diagnostics.csv"), row.names = FALSE)

# Primary four-model files retain simple names for manuscript use.
primary_ic <- subset(ic_table, analysis_variant == "primary_four_model")
primary_coef <- subset(coef_table, analysis_variant == "primary_four_model")
primary_decision <- subset(decision_table, analysis_variant == "primary_four_model")
primary_subject <- subset(subject_table, analysis_variant == "primary_four_model")
primary_fixed_coef <- subset(fixed_coef_table, analysis_variant == "primary_four_model")
primary_fixed_diag <- subset(fixed_diag_table, analysis_variant == "primary_four_model")
primary_margin <- subset(margin_table, analysis_variant == "primary_four_model")

write.csv(primary_ic, file.path(outdir, "candidate_model_ic_weights.csv"), row.names = FALSE)
write.csv(primary_coef, file.path(outdir, "recursive_SEL_MA_coefficients.csv"), row.names = FALSE)
write.csv(primary_decision, file.path(outdir, "decision_summary.csv"), row.names = FALSE)
write.csv(primary_subject, file.path(outdir, "subject_level_stage1_decisions.csv"), row.names = FALSE)
write.csv(primary_fixed_coef, file.path(outdir, "downstream_only_fixed_stage1_coefficients.csv"), row.names = FALSE)
write.csv(primary_fixed_diag, file.path(outdir, "downstream_only_generated_response_diagnostics.csv"), row.names = FALSE)
write.csv(primary_margin, file.path(outdir, "decision_margin_diagnostics.csv"), row.names = FALSE)

# Explicit small-sample diagnostics from the primary library.
small_sample <- unique(primary_ic[, c("stage", "model", "n", "p_beta", "k_ic", "residual_df", "aicc_correction")])
small_sample$p_beta_over_n <- small_sample$p_beta / small_sample$n
small_sample$residual_df_fraction <- small_sample$residual_df / small_sample$n
small_sample <- small_sample[order(small_sample$stage, small_sample$model), ]
write.csv(small_sample, file.path(outdir, "small_sample_diagnostics.csv"), row.names = FALSE)

# Main descriptive plots use the prespecified four-model primary analysis.
write_base_plots(primary, dat, outdir)

# Compact text summary.
sink(file.path(outdir, "analysis_summary.txt"))
cat("ADHD SMART recursive model selection vs model averaging\n")
cat("=====================================================\n\n")
cat("Publicly available simulated ADHD SMART illustrative analysis.\n")
cat("N =", nrow(dat), "; source-aligned stage-2 nonresponders =", sum(elig), "\n")
cat("Responder a2 values present in DTRlearn2 are ignored in the source-aligned analysis.\n\n")
cat("Published d3center reference reproduction max absolute coefficient difference:", max_ref_diff, "\n\n")

for (variant in names(all_results)) {
  cat("=====================================================\n")
  cat("Analysis variant:", variant, "\n")
  cat(variant_defs[[variant]]$description, "\n")
  cat("=====================================================\n\n")
  for (cr in criteria) {
    z <- all_results[[variant]][[cr]]
    cat("---", cr, "---\n")
    cat("Stage 2 selected model:", z$stage2_lib$selected_model, "\n")
    cat("Stage 2 weights:\n")
    print(round(z$stage2_lib$weights, 5))
    cat("Stage 1 SEL selected model:", z$stage1_sel$lib$selected_model, "\n")
    cat("Stage 1 MA weights:\n")
    print(round(z$stage1_ma$lib$weights, 5))

    d1 <- decision_summary(z$stage1_sel$actions, z$stage1_ma$actions, criterion = cr, stage = 1L)
    cat("Stage 1 recursive SEL-MA decision discordance:", d1$discordant_n, "/", d1$n,
        " (", round(100 * d1$discordant_rate, 2), "%)\n", sep = "")

    fd <- fixed_stage1_diagnostic(z)$diagnostics
    cat("Downstream-only fixed-stage-1 decision discordance:", fd$fixed_stage1_decision_discordance_n,
        "/", fd$n, " (", round(100 * fd$fixed_stage1_decision_discordance_rate, 2), "%)\n", sep = "")
    cat("Generated-response RMS SEL-MA shift:", signif(fd$pseudo_diff_rms, 6), "\n")
    cat("Stage-1 wide-only coefficient SEL-MA shift:", signif(fd$gamma_shift_direct, 6), "\n")
    cat("sqrt(n) x coefficient shift:", signif(fd$sqrt_n_gamma_shift_direct, 6), "\n")
    cat("Stage-1 LR SEL-MA shift:", signif(fd$lr_shift, 6), "\n")
    cat("FWL identity error:", signif(fd$gamma_identity_error, 4), "\n")

    md <- decision_margin_diagnostics(z$stage1_sel$actions, z$stage1_ma$actions,
                                      criterion = cr, stage = 1L, analysis_variant = variant)
    cat("Minimum |stage-1 contrast|, SEL:", signif(md$min_abs_contrast_SEL, 6), "\n")
    cat("Minimum |stage-1 contrast|, MA :", signif(md$min_abs_contrast_MA, 6), "\n")
    cat("Maximum |SEL-MA contrast difference|:", signif(md$max_abs_SEL_MA_contrast_difference, 6), "\n\n")
  }
}
cat("Primary application library: four prespecified stage-2 models M20/M21/M22/M23.\n")
cat("M22-vs-M23 is a theory-oriented one-df nested reduction, not a data-selected primary library.\n")
cat("M20-vs-M23 is retained as a strong-separation sensitivity analysis.\n")
cat("AIC is theory-aligned; AICc is a small-sample sensitivity; BIC is a criterion sensitivity.\n")
cat("Do not interpret the observed LR statistic as an estimate of the local parameter delta.\n")
cat("This application is descriptive and is not used to verify asymptotic n*Delta limits.\n")
sink()

cat("Main analysis complete. Results are in:", normalizePath(outdir), "\n")

if (mode %in% c("bootstrap", "all")) {
  cat("Running AIC bootstrap stability diagnostic for the primary four-model library with B =",
      boot_B, "and cores =", boot_cores, "\n")
  b <- bootstrap_stability(dat, B = boot_B, criterion = "AIC", cores = boot_cores,
                           stage2_library = "four_model")
  write.csv(b$subject, file.path(outdir, "bootstrap_AIC_primary_subject_stability.csv"), row.names = FALSE)
  write.csv(b$replicate, file.path(outdir, "bootstrap_AIC_primary_replicate_summary.csv"), row.names = FALSE)
  write.csv(b$failures, file.path(outdir, "bootstrap_AIC_primary_failures.csv"), row.names = FALSE)

  sink(file.path(outdir, "bootstrap_AIC_primary_summary.txt"))
  cat("Algorithmic stability diagnostic only; not a formal confidence interval.\n")
  cat("Primary library: four-model M20/M21/M22/M23.\n")
  cat("Requested B:", boot_B, "\n")
  cat("Successful B:", nrow(b$replicate), "\n")
  cat("Failed B:", nrow(b$failures), "\n")
  cat("Median recursive SEL-MA discordance rate:", median(b$replicate$SEL_MA_discordance_rate), "\n")
  cat("IQR recursive SEL-MA discordance rate:",
      paste(quantile(b$replicate$SEL_MA_discordance_rate, c(0.25, 0.75)), collapse = " to "), "\n")
  cat("Stage-2 selected-model frequencies:\n")
  print(prop.table(table(b$replicate$stage2_selected_model)))
  cat("Stage-1 SEL selected-model frequencies:\n")
  print(prop.table(table(b$replicate$stage1_SEL_selected_model)))
  sink()
  cat("Bootstrap stability diagnostic complete.\n")
}

# Reproducibility metadata.
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()
