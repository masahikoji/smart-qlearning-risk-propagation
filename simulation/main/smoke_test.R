args0 <- commandArgs(trailingOnly = FALSE)
filearg <- sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1])
root <- dirname(normalizePath(filearg))
Sys.setenv(SMART_PROGRAM_DIR = root); setwd(root)
source(file.path(root, "R", "load_all.R"))
run_theory_checks()
ctl <- smart_run_control("debug")
ids <- c(
  "2s_sep_weak_terminal_d1", "2s_smooth_terminal_local",
  "2s_exact_tie_terminal", "2s_criterion_aligned", "2s_four_model_sensitivity",
  "3s_sep_weak_terminal_d1", "3s_smooth_terminal_local",
  "3s_exact_tie_terminal", "3s_criterion_aligned", "3s_four_model_sensitivity"
)
for (sid in ids) {
  s <- scenario_by_id(sid); n <- ctl$n_values[1]
  validate_scenario(s, n, ctl)
  ec <- build_evaluation_cache(s, n, ctl)
  z <- simulate_one_replicate(s, n, replicate_seed(ctl$seed, s$id, n, 1L), ec)
  risks <- z[grep("\\.risk$", names(z))]
  if (!all(is.finite(risks)) || any(risks < -1e-10)) stop("Non-finite/negative risk in ", sid)
  ie <- z[grep("identity_error$", names(z))]
  if (length(ie) && max(abs(ie), na.rm = TRUE) > 1e-7) stop("IC identity failure in ", sid)
  log_msg("Smoke replicate passed: ", sid)
}
# Direct paired-descendant decomposition checks.  Local sources are
# Rao-Blackwellized exactly at one step; paired descendants are used only for
# multi-stage transport products.
for (sid in c("2s_sep_weak_terminal_d1", "2s_smooth_terminal_local",
              "3s_sep_strong_terminal_d1", "3s_smooth_terminal_local")) {
  s <- scenario_by_id(sid); n <- 250L
  ec <- build_evaluation_cache(s, n, ctl)
  z <- simulate_one_replicate(s, n, replicate_seed(ctl$seed, s$id, n, 3L), ec)
  row <- matrix(z, nrow = 1L, dimnames = list(NULL, names(z)))
  cache <- build_decomposition_cache(s, n, ec$data, target_stage = 1L,
                                     N_decomp = 80L, base_seed = ctl$seed,
                                     n_pairs = 2L)
  d <- paired_metrics_batch(s, n, row, "AIC_MA", cache, target_stage = 1L)
  if (!all(is.finite(d))) stop("Non-finite paired decomposition metric in ", sid)
  if (d[1L, "diag_s1"] < -1e-12) {
    stop("Target-stage local diagonal must be nonnegative in ", sid)
  }
  if (abs(d[1L, "diag_s1"] - d[1L, "A_to_s1"]) > 1e-10) {
    stop("First partial linear risk must equal the target-stage diagonal in ", sid)
  }
  if (!is.finite(d[1L, "full_risk_reconstructed"]) ||
      !is.finite(d[1L, "reconstruction_error"])) {
    stop("Invalid reconstructed risk in ", sid)
  }
  log_msg("Paired-descendant exact-local decomposition check passed: ", sid)
}

# Regime-regret check. The estimator uses the exact performance-
# difference identity and must be nonnegative; the true optimal policy has zero
# regret pathwise on the same exogenous trajectories.
for (sid in c("2s_smooth_terminal_local", "3s_smooth_terminal_local")) {
  s <- scenario_by_id(sid); n <- 250L
  ec <- build_evaluation_cache(s, n, ctl)
  z <- simulate_one_replicate(s, n, replicate_seed(ctl$seed, s$id, n, 7L), ec)
  exog <- simulate_value_exogenous(s$K, 120L,
                                   secondary_value_seed(ctl$seed, s$id, n, 7L))
  oracle_err <- true_policy_regret_check(s, n, exog, ec$quad)
  if (!is.finite(oracle_err) || oracle_err > 1e-10) {
    stop("Oracle regret check failed in ", sid)
  }
  rr <- vapply(names(procedure_definitions()), function(proc) {
    regime_regret_gap_one(s, n, z, proc, exog, ec$quad)
  }, numeric(1))
  if (!all(is.finite(rr)) || any(rr < -1e-12)) {
    stop("Invalid optimality-gap regime regret in ", sid)
  }
  log_msg("Optimality-gap regime-regret check passed: ", sid)
}

log_msg("All smoke tests passed.")
