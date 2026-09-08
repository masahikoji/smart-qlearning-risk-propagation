args0 <- commandArgs(trailingOnly = FALSE)
filearg <- sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])
root <- dirname(normalizePath(filearg))
Sys.setenv(SMART_PROGRAM_DIR = root); setwd(root)
source(file.path(root, "R", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1L) as.integer(args[1]) else 2L
profile <- if (length(args) >= 2L) args[2] else "paper"
do_regret <- if (length(args) >= 3L) as.logical(args[3]) else TRUE
ids <- if (length(args) >= 4L && nzchar(args[4])) strsplit(args[4], ",", fixed = TRUE)[[1]] else NULL
overwrite <- if (length(args) >= 5L) isTRUE(as.logical(args[5])) else FALSE
scens <- smart_scenarios(K)
if (!is.null(ids)) scens <- scens[vapply(scens, function(s) s$id %in% ids, logical(1))]
ctl <- smart_run_control(profile)
for (s in scens) for (n in ctl$n_values) run_secondary_scenario_n(s, n, profile, do_regret, overwrite)
