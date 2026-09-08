#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1]] else "all"
profile <- if (length(args) >= 2L) args[[2]] else "paper"
overwrite <- if (length(args) >= 3L) toupper(args[[3]]) %in% c("TRUE", "T", "1", "YES") else FALSE
valid_modes <- c("all", "normal", "transport", "feedback", "finite", "summarize", "smoke")
if (!mode %in% valid_modes) stop("mode must be one of: ", paste(valid_modes, collapse = ", "))
if (!profile %in% c("paper", "smoke", "debug")) stop("profile must be paper, smoke, or debug")

ca <- commandArgs(trailingOnly = FALSE)
hit <- grep("^--file=", ca, value = TRUE)
root <- if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)) else normalizePath(getwd())
source(file.path(root, "R", "additional_sim_functions.R"))
setwd(root)
control <- add_control(profile)
# Keep short pilot/debug chunks physically separate from production chunks.
result_dir_name <- if (profile == "paper") "results" else paste0("results_", profile)
outroot <- safe_dir(file.path(root, result_dir_name))

log_msg("mode=", mode, "; profile=", profile, "; cores=", control$cores, "; results=", outroot)
run_smoke_checks()

if (mode %in% c("smoke")) {
  log_msg("Algebraic smoke checks passed.")
  quit(save = "no", status = 0)
}
if (mode %in% c("normal", "all")) {
  log_msg("Running canonical Gaussian-limit checks")
  run_normal_mean_grid(control, safe_dir(file.path(outroot, "normal_mean")))
}
if (mode %in% c("transport", "all")) {
  log_msg("Running two-stage transported effective-rank SMART simulations")
  run_transport_grid(control, safe_dir(file.path(outroot, "transport")), overwrite = overwrite)
}
if (mode %in% c("feedback", "all")) {
  log_msg("Running criterion-feedback negative control")
  run_feedback_grid(control, safe_dir(file.path(outroot, "feedback")), overwrite = overwrite)
}
if (mode %in% c("finite", "all")) {
  log_msg("Running fixed finite-library AIC simulations")
  run_finite_library_grid(control, safe_dir(file.path(outroot, "finite_library")), overwrite = overwrite)
}
if (mode == "summarize") {
  tr <- file.path(outroot, "transport")
  if (file.exists(file.path(tr, "transport_manifest.csv"))) summarize_transport_grid(tr)
}
base::writeLines(utils::capture.output(utils::sessionInfo()), file.path(outroot, "sessionInfo.txt"))
log_msg("Done.")
