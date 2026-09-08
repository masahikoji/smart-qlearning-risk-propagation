# Monte Carlo verification of the exact Q-error decomposition.
# One-step Bellman sources are evaluated from closed-form conditional moments.
# Paired conditionally independent descendants are used for multi-stage
# transport products and switching terms.

# Decomposition seeds

decomposition_seed <- function(base_seed, scenario_id, n, target_stage,
                               action, pair_id, branch_id) {
  akey <- if (action > 0) 1L else 2L
  z <- as.double(base_seed) + 1100000003 +
    700001 * scenario_hash(scenario_id) +
    1009 * as.integer(n) +
    101 * as.integer(target_stage) +
    37 * akey + 17 * as.integer(pair_id) + 7 * as.integer(branch_id)
  as.integer(z %% 2147483000)
}

# Exact one-step local sources

# Exact conditional expectation of max_a Qhat_{s+1}(H_{s+1},a) given H_s,A_s.
# Qhat is a cubic prognostic polynomial plus a quadratic treatment contrast.
expected_fitted_value_next <- function(scenario, beta_next, data, stage, action) {
  N <- nrow(data)
  if (length(action) == 1L) action <- rep(action, N)
  if (length(action) != N) stop("action length mismatch")
  next_stage <- stage + 1L
  center <- scenario$rho_X * data[[paste0("X", stage)]] + scenario$xi_X * action
  h <- 1 - scenario$rho_X

  p0 <- rep(coef_or_zero(beta_next, "int"), N)
  if (next_stage > 1L) {
    for (k in seq_len(next_stage - 1L)) {
      ak <- if (k == stage) action else data[[paste0("A", k)]]
      p0 <- p0 + coef_or_zero(beta_next, paste0("pastA", k)) * ak
    }
  }
  pcoef <- c(0,
             coef_or_zero(beta_next, "X"),
             coef_or_zero(beta_next, "X_sq"),
             coef_or_zero(beta_next, "X_cu"))
  ep <- p0 + mean_poly_uniform(center, h, pcoef)
  c0 <- coef_or_zero(beta_next, "A")
  c1 <- coef_or_zero(beta_next, "AX")
  c2 <- coef_or_zero(beta_next, "AX2")
  ep + mean_abs_quadratic_uniform(center, h, c0, c1, c2)
}

local_zeta_one <- function(scenario, n, beta_s, beta_next,
                           data, stage, action) {
  qhat <- predict_embedded_q(beta_s, data, stage, action)
  past <- if (stage == 1L) numeric(0) else
    as.matrix(data[paste0("A", seq_len(stage - 1L))])
  mu <- stage_reward_mean(scenario, n, stage,
                          data[[paste0("X", stage)]], past, action)
  if (stage == scenario$K) return(qhat - mu)
  qhat - mu - expected_fitted_value_next(
    scenario, beta_next, data, stage, action)
}

switching_remainder_one <- function(beta_s, data, stage, true_action) {
  qp <- predict_embedded_q(beta_s, data, stage, 1)
  qm <- predict_embedded_q(beta_s, data, stage, -1)
  qd <- predict_embedded_q(beta_s, data, stage, true_action)
  pmax(pmax(qp, qm) - qd, 0)
}

# Conditionally independent descendants

simulate_optimal_descendant <- function(scenario, n, base_data, target_stage,
                                        target_action, seed) {
  K <- scenario$K
  N <- nrow(base_data)
  if (target_stage < 1L || target_stage > K) stop("Invalid target_stage")
  if (length(target_action) == 1L) target_action <- rep(target_action, N)
  if (length(target_action) != N) stop("target_action length mismatch")

  d <- base_data
  d[[paste0("A", target_stage)]] <- target_action
  if (target_stage < K) {
    set.seed(seed)
    for (s in target_stage:(K - 1L)) {
      u <- runif(N, -1, 1)
      d[[paste0("U", s + 1L)]] <- u
      d[[paste0("X", s + 1L)]] <-
        scenario$rho_X * d[[paste0("X", s)]] +
        (1 - scenario$rho_X) * u +
        scenario$xi_X * d[[paste0("A", s)]]
      # All future actions are the designated true optimal actions.  In the
      # two-/three-stage DGPs these comparisons are analytic; no quadrature is
      # required here.
      d[[paste0("A", s + 1L)]] <- true_opt_action(
        scenario, n, s + 1L, d, quad = NULL, tie = 1)
    }
  }
  d
}

build_decomposition_cache <- function(scenario, n, eval_data, target_stage,
                                      N_decomp, base_seed, n_pairs = 1L) {
  if (N_decomp > nrow(eval_data)) stop("N_decomp exceeds evaluation cache size")
  base <- eval_data[seq_len(N_decomp), , drop = FALSE]
  actions <- c(-1, 1)
  pairs <- vector("list", n_pairs)
  for (p in seq_len(n_pairs)) {
    pa <- vector("list", length(actions)); names(pa) <- as.character(actions)
    for (a in actions) {
      branches <- vector("list", 2L)
      for (b in 1:2) {
        sd <- decomposition_seed(base_seed, scenario$id, n, target_stage,
                                 a, p, b)
        branches[[b]] <- simulate_optimal_descendant(
          scenario, n, base, target_stage, a, sd)
      }
      pa[[as.character(a)]] <- branches
    }
    pairs[[p]] <- pa
  }
  list(base = base, actions = actions, pairs = pairs,
       N = N_decomp, n_pairs = n_pairs, target_stage = target_stage)
}

# Coefficient extraction and source matrices

coef_matrix_from_rows <- function(scenario, rows, proc, stage) {
  spec <- stage_model_spec(scenario, stage)
  nm <- paste0(proc, ".s", stage, ".coef_", spec$union)
  B <- as.matrix(rows[, nm, drop = FALSE])
  storage.mode(B) <- "double"
  colnames(B) <- spec$union
  B
}

beta_from_matrix_row <- function(B, i) {
  z <- as.numeric(B[i, , drop = TRUE])
  names(z) <- colnames(B)
  z
}

# Exact local-source and switching-remainder matrices for all fits in a batch on
# one descendant path. Rows are evaluation histories; columns are fitted SMART
# replicates.
source_matrices_for_path <- function(scenario, n, rows, proc,
                                     path_data, target_stage) {
  K <- scenario$K
  Rb <- nrow(rows); N <- nrow(path_data)
  B <- lapply(target_stage:K, function(s)
    coef_matrix_from_rows(scenario, rows, proc, s))
  names(B) <- paste0("s", target_stage:K)

  Z <- vector("list", K - target_stage + 1L)
  names(Z) <- paste0("s", target_stage:K)
  Rho <- if (target_stage < K) {
    z <- vector("list", K - target_stage)
    names(z) <- paste0("s", (target_stage + 1L):K)
    z
  } else list()

  for (s in target_stage:K) {
    spec <- stage_model_spec(scenario, s)
    action <- path_data[[paste0("A", s)]]
    bs <- B[[paste0("s", s)]]
    Fact <- feature_matrix(path_data, s, action, spec$union)
    qhat <- Fact %*% t(bs)

    past <- if (s == 1L) numeric(0) else
      as.matrix(path_data[paste0("A", seq_len(s - 1L))])
    mu <- stage_reward_mean(scenario, n, s,
                            path_data[[paste0("X", s)]], past, action)
    M <- qhat - mu

    if (s < K) {
      bn <- B[[paste0("s", s + 1L)]]
      EV <- matrix(NA_real_, nrow = N, ncol = Rb)
      for (j in seq_len(Rb)) {
        EV[, j] <- expected_fitted_value_next(
          scenario, beta_from_matrix_row(bn, j), path_data, s, action)
      }
      M <- M - EV
    }
    Z[[paste0("s", s)]] <- M

    if (s > target_stage) {
      Fp <- feature_matrix(path_data, s, 1, spec$union)
      Fm <- feature_matrix(path_data, s, -1, spec$union)
      qp <- Fp %*% t(bs)
      qm <- Fm %*% t(bs)
      Rho[[paste0("s", s)]] <- pmax(pmax(qp, qm) - qhat, 0)
    }
  }
  list(Z = Z, rho = Rho)
}

col_mean_prod <- function(A, B) colMeans(A * B)

# One row per fitted replicate. Paired branches are conditionally independent
# given the common target history/action. Off-diagonal source products are
# symmetrized to reduce numerical Monte Carlo noise.
paired_metrics_batch <- function(scenario, n, rows, proc, cache, target_stage) {
  K <- scenario$K
  Rb <- nrow(rows)
  ss <- target_stage:K

  metric_names <- c(
    paste0("diag_s", ss),
    paste0("A_to_s", ss),
    if (length(ss) > 1L) unlist(lapply(seq_along(ss), function(i) {
      if (i == length(ss)) return(character(0))
      paste0("cov_s", ss[i], "_s", ss[(i + 1L):length(ss)])
    })) else character(0),
    if (target_stage < K) paste0("nonlinear_source_sq_s", (target_stage + 1L):K) else character(0),
    "linear_risk", "nonlinear_sq", "linear_nonlinear_cross2",
    "full_risk_direct", "full_risk_reconstructed", "reconstruction_error"
  )
  out <- matrix(0, nrow = Rb, ncol = length(metric_names),
                dimnames = list(NULL, metric_names))

  weight <- 1 / (length(cache$actions) * cache$n_pairs)
  for (p in seq_len(cache$n_pairs)) {
    for (a in cache$actions) {
      paths <- cache$pairs[[p]][[as.character(a)]]
      r1 <- source_matrices_for_path(
        scenario, n, rows, proc, paths[[1L]], target_stage)
      r2 <- source_matrices_for_path(
        scenario, n, rows, proc, paths[[2L]], target_stage)

      L1 <- matrix(0, nrow = cache$N, ncol = Rb)
      L2 <- matrix(0, nrow = cache$N, ncol = Rb)
      for (s in ss) {
        z1 <- r1$Z[[paste0("s", s)]]
        z2 <- r2$Z[[paste0("s", s)]]
        out[, paste0("diag_s", s)] <- out[, paste0("diag_s", s)] +
          weight * col_mean_prod(z1, z2)
        L1 <- L1 + z1; L2 <- L2 + z2
        out[, paste0("A_to_s", s)] <- out[, paste0("A_to_s", s)] +
          weight * col_mean_prod(L1, L2)
      }

      if (length(ss) > 1L) {
        for (ii in seq_len(length(ss) - 1L)) {
          for (jj in (ii + 1L):length(ss)) {
            s <- ss[ii]; u <- ss[jj]
            x <- 0.5 * (
              col_mean_prod(r1$Z[[paste0("s", s)]], r2$Z[[paste0("s", u)]]) +
              col_mean_prod(r1$Z[[paste0("s", u)]], r2$Z[[paste0("s", s)]])
            )
            nm <- paste0("cov_s", s, "_s", u)
            out[, nm] <- out[, nm] + weight * x
          }
        }
      }

      N1 <- matrix(0, nrow = cache$N, ncol = Rb)
      N2 <- matrix(0, nrow = cache$N, ncol = Rb)
      if (target_stage < K) {
        for (s in (target_stage + 1L):K) {
          q1 <- r1$rho[[paste0("s", s)]]
          q2 <- r2$rho[[paste0("s", s)]]
          out[, paste0("nonlinear_source_sq_s", s)] <-
            out[, paste0("nonlinear_source_sq_s", s)] +
            weight * col_mean_prod(q1, q2)
          N1 <- N1 + q1; N2 <- N2 + q2
        }
      }

      out[, "linear_risk"] <- out[, "linear_risk"] +
        weight * col_mean_prod(L1, L2)
      out[, "nonlinear_sq"] <- out[, "nonlinear_sq"] +
        weight * col_mean_prod(N1, N2)
      out[, "linear_nonlinear_cross2"] <-
        out[, "linear_nonlinear_cross2"] +
        weight * colMeans(L1 * N2 + N1 * L2)
      out[, "full_risk_reconstructed"] <-
        out[, "full_risk_reconstructed"] +
        weight * col_mean_prod(L1 + N1, L2 + N2)
    }
  }

  # Reference full risk comes from the already-completed primary run, evaluated
  # on its common N_eval=100000 scoring sample.  This avoids re-evaluating Q^*
  # inside the decomposition diagnostic.
  risk_name <- paste0(proc, ".s", target_stage, ".risk")
  if (!(risk_name %in% colnames(rows))) stop("Missing primary risk column: ", risk_name)
  out[, "full_risk_direct"] <- as.numeric(rows[, risk_name])
  out[, "reconstruction_error"] <-
    out[, "full_risk_direct"] - out[, "full_risk_reconstructed"]
  out
}

# Public driver

run_decomposition_scenario_n <- function(scenario, n, profile = "paper",
                                         target_stages = 1L,
                                         overwrite = FALSE) {
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

  nd <- Sys.getenv("SMART_N_DECOMP", unset = "")
  N_decomp <- if (nzchar(nd)) as.integer(nd) else min(2000L, control$N_eval)
  rd <- Sys.getenv("SMART_R_DECOMP", unset = "")
  R_decomp <- if (nzchar(rd)) as.integer(rd) else min(2000L, nrow(mat))
  pd <- Sys.getenv("SMART_PAIRS_DECOMP", unset = "")
  N_pairs <- if (nzchar(pd)) as.integer(pd) else 1L
  bd <- Sys.getenv("SMART_BATCH_DECOMP", unset = "")
  batch_size <- if (nzchar(bd)) as.integer(bd) else 25L
  cd <- Sys.getenv("SMART_CORES_DECOMP", unset = "")
  cores_decomp <- if (nzchar(cd)) as.integer(cd) else min(8L, control$cores)
  if (!is.finite(N_decomp) || N_decomp < 50L) stop("SMART_N_DECOMP must be >= 50")
  if (!is.finite(R_decomp) || R_decomp < 1L) stop("SMART_R_DECOMP must be >= 1")
  if (!is.finite(N_pairs) || N_pairs < 1L) stop("SMART_PAIRS_DECOMP must be >= 1")
  if (!is.finite(batch_size) || batch_size < 1L) stop("SMART_BATCH_DECOMP must be >= 1")
  if (!is.finite(cores_decomp) || cores_decomp < 1L) stop("SMART_CORES_DECOMP must be >= 1")

  mat <- mat[seq_len(min(R_decomp, nrow(mat))), , drop = FALSE]
  R_decomp <- nrow(mat)
  outfile <- file.path(outdir, "decomposition_metrics.rds")
  if (file.exists(outfile) && !overwrite) {
    old <- readRDS(outfile)
    if (identical(old$method, "paired_descendant_exact_local_v1")) return(invisible(old))
    stop("Existing decomposition file was created by an older method. Re-run with overwrite=TRUE.")
  }

  procs_env <- Sys.getenv("SMART_PROCS_DECOMP", unset = "")
  procs <- if (nzchar(procs_env)) {
    z <- strsplit(procs_env, ",", fixed = TRUE)[[1L]]
    trimws(z)
  } else names(procedure_definitions())
  bad <- setdiff(procs, names(procedure_definitions()))
  if (length(bad)) stop("Unknown SMART_PROCS_DECOMP procedure(s): ", paste(bad, collapse=", "))

  # A decomposition-specific evaluation sample is generated directly; unlike
  # build_evaluation_cache(), this does not evaluate Q^* and therefore does not
  # invoke any quadrature.  The seed is the same deterministic evaluation seed
  # used elsewhere in the project.
  eval_data <- simulate_smart_data(
    scenario, N_decomp, seed = evaluation_seed(control$seed, scenario$id, n),
    include_rewards = FALSE)

  ans <- list(); kk <- 0L
  for (t in target_stages) {
    if (t < 1L || t > scenario$K) stop("Invalid target stage")
    cache <- build_decomposition_cache(
      scenario, n, eval_data, t, N_decomp,
      base_seed = control$seed, n_pairs = N_pairs)

    starts <- seq.int(1L, R_decomp, by = batch_size)
    for (r in procs) {
      log_msg("Decomposition ", scenario$id, " n=", n, " target=", t,
              " proc=", r, " R=", R_decomp, " N=", N_decomp,
              " pairs=", N_pairs, " method=paired-descendant-exact-local")
      one_block <- function(bb) {
        st <- starts[bb]; en <- min(R_decomp, st + batch_size - 1L)
        paired_metrics_batch(
          scenario, n, mat[st:en, , drop = FALSE], r, cache, t)
      }
      if (cores_decomp > 1L && .Platform$OS.type != "windows") {
        blocks <- parallel::mclapply(
          seq_along(starts), one_block, mc.cores = cores_decomp,
          mc.preschedule = TRUE, mc.set.seed = FALSE)
      } else {
        blocks <- lapply(seq_along(starts), one_block)
      }
      M <- do.call(rbind, blocks)
      kk <- kk + 1L
      ans[[kk]] <- list(target_stage = t, procedure = r, metrics = M)
    }
  }

  result <- list(
    scenario = scenario$id, n = n, R_decomp = R_decomp,
    N_decomp = N_decomp, N_pairs = N_pairs, cores_decomp = cores_decomp,
    method = "paired_descendant_exact_local_v1",
    quadrature_order = NA_integer_,
    integration = paste0(
      "closed-form one-step local sources; paired conditionally independent ",
      "optimal-policy descendants for multi-stage transport; no quadrature"),
    evaluation_note = paste0(
      "Full-risk reference is the primary N_eval risk. Transported source products use ",
      "a common fixed descendant grid across fitted procedures/replicates."),
    entries = ans)
  saveRDS(result, outfile, compress = "xz")
  invisible(result)
}
