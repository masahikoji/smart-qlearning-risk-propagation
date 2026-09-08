#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1]] else "all"
profile <- if (length(args) >= 2L) args[[2]] else "paper"
overwrite <- if (length(args) >= 3L) toupper(args[[3]]) %in% c("TRUE", "T", "1", "YES") else FALSE

valid_modes <- c("all", "simulation", "theory", "summarize", "smoke")
if (!mode %in% valid_modes) stop("mode must be one of: ", paste(valid_modes, collapse = ", "))
if (!profile %in% c("paper", "smoke", "debug")) stop("profile must be paper, smoke, or debug")

ca <- commandArgs(trailingOnly = FALSE)
hit <- grep("^--file=", ca, value = TRUE)
root <- if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)) else normalizePath(getwd())
source(file.path(root, "R", "e2pp_stein_functions.R"))
setwd(root)

control <- e2pp_control(profile)
result_base <- if (profile == "paper") "results" else paste0("results_", profile)
outdir <- e2pp_safe_dir(file.path(root, result_base, "E2pp_stein"))

version_file <- file.path(root, "E2PP_VERSION")
version <- if (file.exists(version_file)) trimws(readLines(version_file, warn = FALSE)[1]) else "unknown"
e2pp_log("E2'' version=", version, "; mode=", mode, "; profile=", profile,
          "; n=", control$n, "; R=", control$R, "; cores=", control$cores,
          "; results=", outdir)

e2pp_smoke_checks(control)
e2pp_log("Algebraic and finite-sample smoke checks passed.")

if (mode == "smoke") quit(save = "no", status = 0)

if (mode %in% c("theory", "all")) {
  e2pp_log("Running Gaussian-limit benchmarks with R_theory=", control$R_theory)
  e2pp_run_limit_theory(control, outdir)
}

if (mode %in% c("simulation", "all")) {
  e2pp_log("Running focused E2'' finite-sample simulation")
  e2pp_run_simulation_grid(control, outdir, root, overwrite = overwrite)
}

if (mode == "summarize") {
  e2pp_summarise(control, outdir)
}

base::writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, "sessionInfo.txt"))
e2pp_log("Done.")
