args0 <- commandArgs(trailingOnly = FALSE)
filearg <- sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])
root <- dirname(normalizePath(filearg))
Sys.setenv(SMART_PROGRAM_DIR = root); setwd(root)
source(file.path(root, "R", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1L) args[1] else "paper"
overwrite <- if (length(args) >= 2L) isTRUE(as.logical(args[2])) else FALSE
run_theory_checks()
run_scenario_set(2L, profile = profile, overwrite = overwrite)
run_scenario_set(3L, profile = profile, overwrite = overwrite)
s <- summarize_all(); plot_scaled_gaps(s)
