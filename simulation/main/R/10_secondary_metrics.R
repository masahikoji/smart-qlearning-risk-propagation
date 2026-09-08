# Treatment switching and full-regime value regret.
# Regret uses the finite-horizon performance-difference identity with an
# independent evaluation sample shared across procedures within each replicate.

extract_beta_from_row <- function(row, proc, stage, union) {
  nm <- paste0(proc, ".s", stage, ".coef_", union)
  z <- as.numeric(row[nm])
  names(z) <- union
  z
}

compute_switching_metrics <- function(scenario, n, mat, eval_cache, batch_size = 25L) {
  procs <- names(procedure_definitions())
  R <- nrow(mat)
  out <- matrix(NA_real_, R, 0)
  for (s in seq_len(scenario$K)) {
    cache <- eval_cache$risk[[paste0("stage", s)]]
    C <- cache$contrast_matrix
    for (r in procs) {
      cn <- paste0(r, ".s", s, ".coef_", cache$union)
      B <- mat[, cn, drop = FALSE]
      colnames(B) <- cache$union
      metric <- numeric(R)
      freq_plus <- numeric(R)
      starts <- seq.int(1L, R, by = batch_size)
      for (st in starts) {
        en <- min(R, st + batch_size - 1L)
        Dhat <- C %*% t(B[st:en, , drop = FALSE])
        plus <- Dhat >= 0
        freq_plus[st:en] <- colMeans(plus)
        if (scenario$boundary == "exact_tie" && s == scenario$K) {
          metric[st:en] <- NA_real_
        } else {
          truth_plus <- cache$true_opt == 1
          metric[st:en] <- colMeans(plus != truth_plus)
        }
      }
      out <- cbind(out, metric, freq_plus)
      colnames(out)[(ncol(out)-1L):ncol(out)] <- c(
        paste0(r, ".s", s, ".switch_prob"),
        paste0(r, ".s", s, ".plus_action_freq"))
    }
  }
  out
}

# Separate deterministic seed stream for regime-value evaluation.
secondary_value_seed <- function(base_seed, scenario_id, n, rep_id) {
  z <- as.double(base_seed) + 1300000033 +
    700001 * scenario_hash(scenario_id) +
    3001 * as.integer(n) + 7919 * as.integer(rep_id)
  as.integer(z %% 2147483000)
}

simulate_value_exogenous <- function(K, N_value, seed) {
  set.seed(seed)
  out <- data.frame(X1 = runif(N_value, -1, 1))
  if (K >= 2L) {
    for (s in 2:K) out[[paste0("U", s)]] <- runif(N_value, -1, 1)
  }
  out
}

# Exact finite-horizon performance-difference representation:
#
#   V^* - V^d = E_d sum_s { V_s^*(H_s) - Q_s^*(H_s,d_s(H_s)) }.
#
# Each summand is a nonnegative optimality gap.  We therefore estimate regret by
# simulating only the fitted policy forward and accumulating true optimality gaps
# along that trajectory.  This is preferable to subtracting an analytic true
# value from a separate Monte Carlo policy value: the latter can leave a common
# finite-evaluation offset even when the fitted and optimal policies coincide.
regime_regret_gap_one <- function(scenario, n, row, proc, eval_exog, quad,
                                  stage1_truth = NULL, tol = 1e-10) {
  K <- scenario$K
  N <- nrow(eval_exog)
  d <- data.frame(X1 = eval_exog$X1)
  regret_path <- numeric(N)

  for (s in seq_len(K)) {
    spec <- stage_model_spec(scenario, s)
    beta <- extract_beta_from_row(row, proc, s, spec$union)

    qhat_plus <- predict_embedded_q(beta, d, s, 1)
    qhat_minus <- predict_embedded_q(beta, d, s, -1)
    a <- ifelse(qhat_plus >= qhat_minus, 1, -1)

    if (s == 1L && !is.null(stage1_truth)) {
      qstar_plus <- stage1_truth$plus
      qstar_minus <- stage1_truth$minus
    } else {
      qstar_plus <- true_q(scenario, n, s, d, 1, quad)
      qstar_minus <- true_q(scenario, n, s, d, -1, quad)
    }
    vstar <- pmax(qstar_plus, qstar_minus)
    qstar_taken <- ifelse(a == 1, qstar_plus, qstar_minus)
    gap <- vstar - qstar_taken
    if (any(gap < -tol, na.rm = TRUE)) {
      stop("Negative true optimality gap beyond numerical tolerance in ",
           scenario$id, ", n=", n, ", stage=", s, ", proc=", proc)
    }
    regret_path <- regret_path + pmax(gap, 0)

    d[[paste0("A", s)]] <- a
    if (s < K) {
      u <- eval_exog[[paste0("U", s + 1L)]]
      d[[paste0("X", s + 1L)]] <- scenario$rho_X * d[[paste0("X", s)]] +
        (1 - scenario$rho_X) * u + scenario$xi_X * a
    }
  }

  mean(regret_path)
}

# Diagnostic identity check: the true optimal policy must have zero regret
# pathwise under the same performance-difference representation (up to roundoff).
true_policy_regret_check <- function(scenario, n, eval_exog, quad, tol = 1e-10) {
  K <- scenario$K
  N <- nrow(eval_exog)
  d <- data.frame(X1 = eval_exog$X1)
  total <- numeric(N)
  for (s in seq_len(K)) {
    qp <- true_q(scenario, n, s, d, 1, quad)
    qm <- true_q(scenario, n, s, d, -1, quad)
    a <- ifelse(qp >= qm, 1, -1)
    gap <- pmax(qp, qm) - ifelse(a == 1, qp, qm)
    if (max(abs(gap), na.rm = TRUE) > tol) {
      stop("Oracle policy regret identity failure in ", scenario$id,
           ", n=", n, ", stage=", s)
    }
    total <- total + pmax(gap, 0)
    d[[paste0("A", s)]] <- a
    if (s < K) {
      u <- eval_exog[[paste0("U", s + 1L)]]
      d[[paste0("X", s + 1L)]] <- scenario$rho_X * d[[paste0("X", s)]] +
        (1 - scenario$rho_X) * u + scenario$xi_X * a
    }
  }
  max(abs(total), na.rm = TRUE)
}

compute_regime_regret <- function(scenario, n, mat, eval_cache, control,
                                  cores = control$cores) {
  procs <- names(procedure_definitions())
  R <- nrow(mat)

  # Independent evaluation trajectories are generated for every fitted SMART
  # replicate.  The same trajectories are shared across procedures within that
  # replicate (common random numbers), so method contrasts are efficient while
  # replicate-level MCSE remains valid.  With R=20000, N_value=500 gives
  # 10 million evaluation trajectories per procedure/scenario-n combination.
  nv <- Sys.getenv("SMART_N_VALUE", unset = "")
  N_value <- if (nzchar(nv)) as.integer(nv) else min(500L, as.integer(control$N_eval))
  if (!is.finite(N_value) || N_value < 1L) stop("SMART_N_VALUE must be a positive integer")

  eval_one <- function(i) {
    exog <- simulate_value_exogenous(
      scenario$K, N_value,
      secondary_value_seed(control$seed, scenario$id, n, mat[i, "rep_id"])
    )
    d1 <- data.frame(X1 = exog$X1)
    stage1_truth <- list(
      plus = true_q(scenario, n, 1L, d1, 1, eval_cache$quad),
      minus = true_q(scenario, n, 1L, d1, -1, eval_cache$quad)
    )
    vapply(procs, function(proc) {
      regime_regret_gap_one(scenario, n, mat[i, ], proc, exog, eval_cache$quad,
                            stage1_truth = stage1_truth)
    }, numeric(1))
  }

  if (cores > 1L && .Platform$OS.type != "windows") {
    vals <- parallel::mclapply(seq_len(R), eval_one, mc.cores = cores,
                               mc.preschedule = TRUE, mc.set.seed = FALSE)
  } else {
    vals <- lapply(seq_len(R), eval_one)
  }
  out <- do.call(rbind, vals)
  colnames(out) <- paste0(procs, ".regime_regret")

  attr(out, "N_value") <- N_value
  attr(out, "regret_method") <- "optimality_gap_forward_mc_v1"
  attr(out, "evaluation_design") <- "independent_by_fit_shared_across_procedures"
  out
}

run_secondary_scenario_n <- function(scenario, n, profile = "paper",
                                     do_regret = TRUE, overwrite = FALSE) {
  outdir <- scenario_result_dir(scenario, n)
  metafile <- file.path(outdir, "metadata.rds")
  if (file.exists(metafile)) {
    meta <- readRDS(metafile)
    scenario <- meta$scenario
    control <- meta$control
  } else {
    control <- smart_run_control(profile)
  }
  mat <- read_scenario_matrix(scenario$id, n)
  eval_cache <- build_evaluation_cache(scenario, n, control)
  outfile <- file.path(outdir, "secondary_metrics.rds")
  if (file.exists(outfile) && !overwrite) {
    log_msg("Skipping existing secondary metrics: ", outfile)
    return(invisible(readRDS(outfile)))
  }
  log_msg("Computing treatment-switching metrics for ", scenario$id, ", n=", n)
  sw <- compute_switching_metrics(scenario, n, mat, eval_cache)
  rr <- NULL
  rr_meta <- list(N_value = NA_integer_, regret_method = NA_character_,
                  evaluation_design = NA_character_)
  if (do_regret) {
    log_msg("Computing regime regret by optimality-gap forward MC for ",
            scenario$id, ", n=", n)
    rr <- compute_regime_regret(scenario, n, mat, eval_cache, control)
    rr_meta <- list(
      N_value = attr(rr, "N_value"),
      regret_method = attr(rr, "regret_method"),
      evaluation_design = attr(rr, "evaluation_design")
    )
  }
  ans <- list(rep_id = mat[, "rep_id"], switching = sw, regret = rr,
              regret_meta = rr_meta, scenario = scenario$id, n = n)
  saveRDS(ans, outfile, compress = "xz")
  invisible(ans)
}
