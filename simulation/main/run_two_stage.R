args0 <- commandArgs(trailingOnly = FALSE)
filearg <- sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])
root <- dirname(normalizePath(filearg))
Sys.setenv(SMART_PROGRAM_DIR = root); setwd(root)
source(file.path(root, "R", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1L) args[1] else "paper"
ids <- if (length(args) >= 2L && nzchar(args[2])) strsplit(args[2], ",", fixed = TRUE)[[1]] else NULL
overwrite <- if (length(args) >= 3L) isTRUE(as.logical(args[3])) else FALSE
run_theory_checks()
run_scenario_set(2L, profile = profile, ids = ids, overwrite = overwrite)
