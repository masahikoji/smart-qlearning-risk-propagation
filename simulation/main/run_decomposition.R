args0 <- commandArgs(trailingOnly = FALSE)
filearg <- sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])
root <- dirname(normalizePath(filearg))
Sys.setenv(SMART_PROGRAM_DIR = root); setwd(root)
source(file.path(root, "R", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1L) as.integer(args[1]) else 2L
profile <- if (length(args) >= 2L) args[2] else "paper"
ids <- if (length(args) >= 3L && nzchar(args[3])) strsplit(args[3], ",", fixed = TRUE)[[1]] else NULL
overwrite <- if (length(args) >= 4L) isTRUE(as.logical(args[4])) else FALSE
scens <- smart_scenarios(K)
if (is.null(ids)) ids <- principal_decomposition_ids(K)
scens <- scens[vapply(scens, function(s) s$id %in% ids, logical(1))]
ctl <- smart_run_control(profile)
for (s in scens) for (n in ctl$n_values) run_decomposition_scenario_n(s, n, profile, target_stages = 1L, overwrite = overwrite)
