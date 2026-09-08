#!/usr/bin/env Rscript
ca <- commandArgs(trailingOnly = FALSE)
hit <- grep("^--file=", ca, value = TRUE)
root <- if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)) else normalizePath(getwd())
source(file.path(root, "R", "adhd_sel_ma_functions.R"))
source(file.path(root, "R", "adhd_candidate_sensitivity_functions.R"))
source(file.path(root, "R", "adhd_stein_additional_functions.R"))
setwd(root)
outdir <- file.path(root, "results_adhd_stein")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- load_adhd_data()
centers <- compute_centers(dat)
dat_c <- apply_centers(dat, centers)
pairs <- make_block_pairs()

all_summary <- list(); all_subject <- list(); all_eigen <- list(); all_diff <- list()
k <- 1L
for (nm in names(pairs)) {
  cat("Running block:", nm, "\n")
  z <- run_one_adhd_block(dat_c, nm, pairs[[nm]], a_values = c(1, 2))
  all_summary[[k]] <- z$summary
  all_subject[[k]] <- z$subject
  all_eigen[[k]] <- z$eigen
  all_diff[[k]] <- pairwise_rule_differences(z$subject)
  k <- k + 1L
}
summary_df <- do.call(rbind, all_summary)
subject_df <- do.call(rbind, all_subject)
eigen_df <- do.call(rbind, all_eigen)
diff_df <- do.call(rbind, all_diff)

utils::write.csv(summary_df, file.path(outdir, "adhd_stein_block_summary.csv"), row.names = FALSE)
threshold_df <- unique(summary_df[, c(
  "pair", "role", "stage2_n", "d", "lambda", "AIC_threshold",
  "AIC_diff_wide_minus_narrow", "AIC_identity_error",
  "stage2_wide_selected_by_AIC", "transported_reff_hat"
)])
utils::write.csv(threshold_df, file.path(outdir, "adhd_block_threshold_diagnostics.csv"), row.names = FALSE)
utils::write.csv(eigen_df, file.path(outdir, "adhd_transport_effective_rank.csv"), row.names = FALSE)
utils::write.csv(diff_df, file.path(outdir, "adhd_stein_vs_selection_propagation.csv"), row.names = FALSE)
utils::write.csv(subject_df, file.path(outdir, "adhd_stein_subject_level.csv"), row.names = FALSE)
utils::write.csv(data.frame(variable = names(centers), center = as.numeric(centers)),
                 file.path(outdir, "centering_constants.csv"), row.names = FALSE)

sink(file.path(outdir, "analysis_summary.txt"))
cat("ADHD Pretest--Stein additional sensitivity\n")
cat("==========================================\n\n")
cat("This is an exploratory sensitivity based only on nested pairs formed from candidates\n")
cat("that already existed in the prespecified candidate-library sensitivity code. No new\n")
cat("d>=3 effect-modifier block was created solely to make the Stein theorem applicable.\n\n")
cat("Threshold diagnostics:\n")
print(threshold_df)
cat("\nRule-specific summaries:\n")
print(summary_df)
cat("\nTransported effective-rank diagnostics:\n")
print(eigen_df)
cat("\nSelection-versus-alternative propagation differences:\n")
print(diff_df)
cat("\nImportant: transported_reff_hat is a same-data empirical projection diagnostic.\n")
cat("It must not be used as if it were the population M in the theorem, and a is not\n")
cat("adaptively chosen from this estimate. The primary M22-vs-M23 analysis remains d=1\n")
cat("and therefore remains outside the d>=3 Pretest--Stein dominance result.\n")
sink()
base::writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("Done. Results written to", outdir, "\n")
