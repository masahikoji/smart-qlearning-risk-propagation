# Configuration for the SMART Q-learning model-selection/model-averaging simulations.
# Every numerical choice not fixed by theory is centralized here.

smart_run_control <- function(profile = c("paper", "smoke", "debug")) {
  profile <- match.arg(profile)
  detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  cores_default <- if (length(detected_cores) != 1L || is.na(detected_cores)) 1L else max(1L, as.integer(detected_cores) - 2L)
  if (profile == "paper") {
    ctl <- list(
      n_values = c(250L, 500L, 1000L),
      R = 20000L,
      N_eval = 100000L,
      chunk_size = 250L,
      cores = cores_default,
      quadrature_order = 48L,
      save_coefficients = TRUE,
      stop_on_gap_failure = TRUE,
      gap_grid = 801L,
      gap_validation_floor = 0.20,
      seed = 260906L
    )
  } else if (profile == "smoke") {
    ctl <- list(
      n_values = 250L,
      R = 20L,
      N_eval = 5000L,
      chunk_size = 10L,
      cores = 1L,
      quadrature_order = 24L,
      save_coefficients = TRUE,
      stop_on_gap_failure = TRUE,
      gap_grid = 201L,
      gap_validation_floor = 0.10,
      seed = 260906L
    )
  } else {
    ctl <- list(
      n_values = 250L,
      R = 2L,
      N_eval = 500L,
      chunk_size = 2L,
      cores = 1L,
      quadrature_order = 16L,
      save_coefficients = TRUE,
      stop_on_gap_failure = TRUE,
      gap_grid = 101L,
      gap_validation_floor = 0.05,
      seed = 260906L
    )
  }

  int_override <- function(name, old) {
    z <- Sys.getenv(name, unset = "")
    if (nzchar(z)) as.integer(z) else old
  }
  ctl$R <- int_override("SMART_R", ctl$R)
  ctl$N_eval <- int_override("SMART_N_EVAL", ctl$N_eval)
  ctl$chunk_size <- int_override("SMART_CHUNK", ctl$chunk_size)
  ctl$cores <- int_override("SMART_CORES", ctl$cores)
  ctl$quadrature_order <- int_override("SMART_QUAD", ctl$quadrature_order)
  ctl
}

base_reward_parameters <- function(K) {
  stopifnot(K %in% c(2L, 3L))
  beta0 <- rep(0.25, K)
  beta1 <- rep(0.50, K)
  beta_past <- vector("list", K)
  beta_past[[1L]] <- numeric(0)
  if (K >= 2L) beta_past[[2L]] <- 0.20
  if (K >= 3L) beta_past[[3L]] <- c(0.15, 0.15)
  list(beta0 = beta0, beta1 = beta1, beta_past = beta_past,
       sigma = rep(1.0, K))
}

# Even state moments for the symmetric transition
# X_{s+1}=rho X_s+(1-rho)U_{s+1}+xi A_s.
state_even_moments <- function(rho_X, xi_X, K) {
  m2 <- numeric(K); m4 <- numeric(K)
  m2[1L] <- 1/3; m4[1L] <- 1/5
  if (K >= 2L) {
    h <- 1 - rho_X
    for (s in 2:K) {
      m2[s] <- rho_X^2 * m2[s-1L] + h^2/3 + xi_X^2
      m4[s] <- rho_X^4 * m4[s-1L] + h^4/5 + xi_X^4 +
        6 * rho_X^2 * h^2 * m2[s-1L] / 3 +
        6 * rho_X^2 * xi_X^2 * m2[s-1L] +
        6 * h^2 * xi_X^2 / 3
    }
  }
  list(m2 = m2, m4 = m4)
}

make_scenario <- function(id, K, family = "standard", candidate_count = 2L,
                          rho_X = 0.20, xi_X = 0.10,
                          boundary = "separated",
                          departure = "local",
                          b = NULL, c_fixed = NULL,
                          tau0 = NULL, tau1 = NULL,
                          criterion_variant = NA_character_,
                          accumulation_role = "general",
                          target_delta_terminal = NA_real_,
                          notes = "") {
  stopifnot(K %in% c(2L, 3L))
  if (is.null(b)) b <- rep(0, K)
  if (is.null(c_fixed)) c_fixed <- rep(0, K)
  if (is.null(tau0)) tau0 <- rep(1.35, K)
  if (is.null(tau1)) tau1 <- rep(0.20, K)
  stopifnot(length(b) == K, length(c_fixed) == K,
            length(tau0) == K, length(tau1) == K)
  out <- c(
    list(id = id, K = as.integer(K), family = family,
         candidate_count = as.integer(candidate_count),
         rho_X = rho_X, xi_X = xi_X, boundary = boundary,
         departure = departure, b = as.numeric(b),
         c_fixed = as.numeric(c_fixed), tau0 = as.numeric(tau0),
         tau1 = as.numeric(tau1), criterion_variant = criterion_variant,
         accumulation_role = accumulation_role,
         target_delta_terminal = as.numeric(target_delta_terminal), notes = notes),
    base_reward_parameters(as.integer(K))
  )
  class(out) <- c("smart_scenario", "list")
  if (is.finite(out$target_delta_terminal) && out$departure == "local") {
    mm <- state_even_moments(out$rho_X, out$xi_X, out$K)
    gK <- mm$m4[out$K] - mm$m2[out$K]^2
    if (!(gK > 0)) stop("Nonpositive terminal residual variance in ", id)
    out$b[out$K] <- out$target_delta_terminal * out$sigma[out$K] / sqrt(gK)
  }
  out
}

# Main scenarios. The delta-grid scenarios are included explicitly because the
# one-stage AIC gap changes sign around delta=2.3648; using b=2 everywhere would
# leave all terminal comparisons in the same positive-gap region.
smart_scenarios <- function(K = c(2L, 3L)) {
  K <- as.integer(match.arg(as.character(K), c("2", "3")))

  if (K == 2L) {
    list(
      make_scenario(
        "2s_sep_weak_terminal_d1", 2L,
        rho_X = 0.20, xi_X = 0.10, boundary = "separated",
        target_delta_terminal = 1.0, accumulation_role = "orthogonal",
        notes = "Primary weak-transport accumulation benchmark; terminal local departure only."
      ),
      make_scenario(
        "2s_sep_strong_terminal_d1", 2L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        target_delta_terminal = 1.0, accumulation_role = "orthogonal",
        notes = "Primary strong-transport accumulation benchmark; terminal local departure only."
      ),
      make_scenario(
        "2s_sep_stage1_stage2_local", 2L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        b = c(2.0, 0), target_delta_terminal = 1.0,
        accumulation_role = "multiple_source",
        notes = "Local departures at both stages; finite-sample cross-stage covariance is estimated without a sign claim."
      ),
      make_scenario(
        "2s_smooth_terminal_local", 2L,
        rho_X = 0.30, xi_X = 0.00, boundary = "smooth",
        target_delta_terminal = 1.0,
        tau0 = c(1.35, 0.00), tau1 = c(0.20, 1.00),
        accumulation_role = "boundary",
        notes = "Smooth boundary only at stage 2; X^3 is included only for this smooth-boundary basis."
      ),
      make_scenario(
        "2s_exact_tie_terminal", 2L,
        rho_X = 0.30, xi_X = 0.00, boundary = "exact_tie",
        b = c(0, 0), tau0 = c(1.35, 0.00), tau1 = c(0.20, 0.00),
        accumulation_role = "boundary",
        notes = "Exact treatment tie at stage 2."
      ),
      make_scenario(
        "2s_criterion_aligned", 2L, family = "criterion", candidate_count = 2L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        target_delta_terminal = 1.0, criterion_variant = "aligned",
        accumulation_role = "criterion",
        notes = "Stage-1 wide-only direction is X1^2; separated gaps only."
      ),
      make_scenario(
        "2s_criterion_control", 2L, family = "criterion", candidate_count = 2L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        target_delta_terminal = 1.0, criterion_variant = "control",
        accumulation_role = "criterion",
        notes = "Stage-1 wide-only direction is A1*X1^2; direct first-order population projection is zero."
      ),
      make_scenario(
        "2s_sep_fixed_sensitivity", 2L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        departure = "fixed", c_fixed = c(0.15, 0.15),
        accumulation_role = "sensitivity",
        notes = "Fixed-effect sensitivity away from the root-n local regime."
      ),
      make_scenario(
        "2s_four_model_sensitivity", 2L, candidate_count = 4L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        b = c(2.0, 2.0), accumulation_role = "sensitivity",
        notes = "Four-candidate robustness analysis; M2 omits AX despite nonzero tau1 and is intentionally outside the local-correct-specification theory."
      ),
      make_scenario(
        "2s_delta_d05", 2L, rho_X = 0.50, xi_X = 0.30,
        boundary = "separated", target_delta_terminal = 0.5,
        accumulation_role = "delta_grid", notes = "Terminal local-signal grid: delta=0.5."
      ),
      make_scenario(
        "2s_delta_d24", 2L, rho_X = 0.50, xi_X = 0.30,
        boundary = "separated", target_delta_terminal = 2.4,
        accumulation_role = "delta_grid", notes = "Terminal local-signal grid: delta=2.4, near the sign-change point."
      ),
      make_scenario(
        "2s_delta_d30", 2L, rho_X = 0.50, xi_X = 0.30,
        boundary = "separated", target_delta_terminal = 3.0,
        accumulation_role = "delta_grid", notes = "Terminal local-signal grid: delta=3.0."
      )
    )
  } else {
    list(
      make_scenario(
        "3s_sep_weak_terminal_d1", 3L,
        rho_X = 0.20, xi_X = 0.10, boundary = "separated",
        target_delta_terminal = 1.0, accumulation_role = "orthogonal",
        notes = "Primary weak-transport three-stage accumulation benchmark."
      ),
      make_scenario(
        "3s_sep_strong_terminal_d1", 3L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        target_delta_terminal = 1.0, accumulation_role = "orthogonal",
        notes = "Primary strong-transport three-stage accumulation benchmark."
      ),
      make_scenario(
        "3s_sep_stage2_stage3_local", 3L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        b = c(0, 2.0, 0), target_delta_terminal = 1.0,
        accumulation_role = "multiple_source",
        notes = "Local departures at stages 2 and 3; covariance is estimated explicitly without an orthogonality claim."
      ),
      make_scenario(
        "3s_smooth_terminal_local", 3L,
        rho_X = 0.30, xi_X = 0.00, boundary = "smooth",
        target_delta_terminal = 1.0,
        tau0 = c(1.50, 1.50, 0.00), tau1 = c(0.15, 0.15, 1.00),
        accumulation_role = "boundary",
        notes = "Smooth boundary confined to stage 3."
      ),
      make_scenario(
        "3s_exact_tie_terminal", 3L,
        rho_X = 0.30, xi_X = 0.00, boundary = "exact_tie",
        b = c(0, 0, 0), tau0 = c(1.50, 1.50, 0.00),
        tau1 = c(0.15, 0.15, 0.00), accumulation_role = "boundary",
        notes = "Exact tie confined to stage 3."
      ),
      make_scenario(
        "3s_criterion_aligned", 3L, family = "criterion", candidate_count = 2L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        target_delta_terminal = 1.0, criterion_variant = "aligned",
        accumulation_role = "criterion",
        notes = "Stage 3 -> stage 2 -> stage 1 criterion-perturbation cascade."
      ),
      make_scenario(
        "3s_criterion_control", 3L, family = "criterion", candidate_count = 2L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        target_delta_terminal = 1.0, criterion_variant = "control",
        accumulation_role = "criterion",
        notes = "Criterion-score control; direct first-order population projection is zero."
      ),
      make_scenario(
        "3s_sep_fixed_sensitivity", 3L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        departure = "fixed", c_fixed = c(0, 0.15, 0.15),
        accumulation_role = "sensitivity", notes = "Fixed-effect sensitivity."
      ),
      make_scenario(
        "3s_four_model_sensitivity", 3L, candidate_count = 4L,
        rho_X = 0.50, xi_X = 0.30, boundary = "separated",
        b = c(0, 2.0, 2.0), accumulation_role = "sensitivity",
        notes = "Four-candidate robustness analysis; M2 is intentionally outside the local-correct-specification theory when tau1 is nonzero."
      ),
      make_scenario(
        "3s_delta_d05", 3L, rho_X = 0.50, xi_X = 0.30,
        boundary = "separated", target_delta_terminal = 0.5,
        accumulation_role = "delta_grid", notes = "Terminal local-signal grid: delta=0.5."
      ),
      make_scenario(
        "3s_delta_d24", 3L, rho_X = 0.50, xi_X = 0.30,
        boundary = "separated", target_delta_terminal = 2.4,
        accumulation_role = "delta_grid", notes = "Terminal local-signal grid: delta=2.4."
      ),
      make_scenario(
        "3s_delta_d30", 3L, rho_X = 0.50, xi_X = 0.30,
        boundary = "separated", target_delta_terminal = 3.0,
        accumulation_role = "delta_grid", notes = "Terminal local-signal grid: delta=3.0."
      )
    )
  }
}

resolve_c <- function(scenario, n) {
  if (scenario$departure == "fixed") out <- scenario$c_fixed else out <- scenario$b / sqrt(n)
  if (scenario$boundary == "exact_tie") out[scenario$K] <- 0
  out
}

scenario_by_id <- function(id) {
  all <- c(smart_scenarios(2L), smart_scenarios(3L))
  idx <- which(vapply(all, function(z) identical(z$id, id), logical(1)))
  if (length(idx) != 1L) stop("Unknown or duplicated scenario id: ", id)
  all[[idx]]
}

# Decomposition is intentionally run only for the scenarios that diagnose the
# mechanisms promised in the manuscript. Other ids can still be supplied manually.
principal_decomposition_ids <- function(K) {
  if (K == 2L) {
    c("2s_sep_strong_terminal_d1", "2s_sep_stage1_stage2_local",
      "2s_smooth_terminal_local", "2s_exact_tie_terminal")
  } else {
    c("3s_sep_strong_terminal_d1", "3s_sep_stage2_stage3_local",
      "3s_smooth_terminal_local", "3s_exact_tie_terminal")
  }
}
