# Source all simulation modules in dependency order.
.root <- Sys.getenv("SMART_PROGRAM_DIR", unset = getwd())
mods <- c(
  "00_config.R", "01_utils.R", "02_features.R", "03_fit.R",
  "04_truth_dgp.R", "05_validation.R", "06_eval_cache.R",
  "07_recursive_fit.R", "08_runner.R", "09_summary.R",
  "10_secondary_metrics.R", "11_decomposition.R", "12_theory_checks.R",
  "13_calibration.R", "14_population_checks.R"
)
for (m in mods) source(file.path(.root, "R", m), local = .GlobalEnv)
rm(.root, mods)
