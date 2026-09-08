# Scenario-level validation against the theoretical restrictions in Section 3.

validate_scenario <- function(scenario, n, control, quad = gauss_legendre(control$quadrature_order),
                              verbose = TRUE) {
  K <- scenario$K
  if (scenario$boundary == "smooth") {
    if (abs(scenario$xi_X) > 1e-14) stop("Smooth-boundary scenarios require xi_X=0 in this program.")
    if (!(scenario$rho_X > 0 && scenario$rho_X < 0.5)) stop("Smooth-boundary scenarios require 0<rho_X<1/2.")
  }
  if (scenario$family == "criterion" && scenario$boundary != "separated") {
    stop("Criterion-perturbation scenarios are restricted to separated treatment gaps.")
  }

  # Validate all stages that are supposed to be separated.  For smooth/exact tie,
  # only the terminal stage is exempt; upstream stages remain separated.
  separated_stages <- if (scenario$boundary %in% c("smooth", "exact_tie")) {
    if (K > 1L) seq_len(K - 1L) else integer(0)
  } else seq_len(K)

  mins <- rep(NA_real_, K)
  for (s in separated_stages) {
    pa <- past_action_grid(s)
    min_gap <- Inf
    for (r in seq_len(nrow(pa))) {
      past <- if (s == 1L) numeric(0) else pa[r, ]
      su <- state_support_given_past(scenario, s, past)
      x <- seq(su[1], su[2], length.out = control$gap_grid)
      dat <- data.frame(X1 = rep(0, length(x)))
      dat[[paste0("X", s)]] <- x
      if (s > 1L) {
        for (j in seq_len(s - 1L)) dat[[paste0("A", j)]] <- rep(past[j], length(x))
      }
      qp <- true_q(scenario, n, s, dat, 1, quad)
      qm <- true_q(scenario, n, s, dat, -1, quad)
      min_gap <- min(min_gap, min(abs(qp - qm)))
    }
    mins[s] <- min_gap
    if (min_gap < control$gap_validation_floor) {
      msg <- sprintf("Scenario %s, n=%d: stage-%d minimum true gap %.6f < validation floor %.6f",
                     scenario$id, n, s, min_gap, control$gap_validation_floor)
      if (isTRUE(control$stop_on_gap_failure)) stop(msg) else warning(msg)
    }
  }

  # Terminal boundary checks.
  if (scenario$boundary == "smooth") {
    cvec <- resolve_c(scenario, n)
    if (abs(scenario$tau0[K]) > 1e-12 || scenario$tau1[K] <= 0) {
      stop("Smooth boundary is coded at x=0 and requires terminal tau0=0, tau1>0.")
    }
    # Ensure the local quadratic does not introduce an extra root over the reachable support.
    all_su <- apply(past_action_grid(K), 1L, function(pa) state_support_given_past(scenario, K, pa))
    max_abs_x <- max(abs(all_su))
    if (abs(cvec[K]) * max_abs_x >= scenario$tau1[K]) {
      stop("Local quadratic is too large: it can create an unintended extra terminal treatment boundary.")
    }
  }
  if (scenario$boundary == "exact_tie") {
    cvec <- resolve_c(scenario, n)
    if (max(abs(c(scenario$tau0[K], scenario$tau1[K], cvec[K]))) > 1e-12) {
      stop("Exact-tie scenario does not have identically zero terminal treatment contrast.")
    }
  }

  if (verbose) {
    sg <- paste(sprintf("s%d=%s", seq_len(K),
                        ifelse(is.na(mins), "boundary", format(round(mins, 4), nsmall = 4))),
                collapse = ", ")
    log_msg("Validated ", scenario$id, " at n=", n, " [", sg, "]")
  }
  invisible(list(min_gap = mins))
}

validate_nested_ic_identity <- function(fits, n, tol = 1e-8) {
  d <- nested_ic_diagnostics(fits, n)
  if (is.null(d)) return(TRUE)
  if (!is.finite(d$aic_identity_error) || !is.finite(d$bic_identity_error)) return(TRUE)
  ok <- abs(d$aic_identity_error) < tol && abs(d$bic_identity_error) < tol
  if (!ok) stop("AIC/BIC nested-model identity failed.")
  TRUE
}

validate_conditional_second_moment <- function(scenario, N = 200000L, seed = 123) {
  set.seed(seed)
  x1 <- runif(N, -1, 1)
  a1 <- sample(c(-1, 1), N, TRUE)
  u2 <- runif(N, -1, 1)
  x2 <- scenario$rho_X * x1 + (1 - scenario$rho_X) * u2 + scenario$xi_X * a1
  target <- scenario$rho_X^2 * x1^2 + 2 * scenario$rho_X * scenario$xi_X * x1 * a1 +
    scenario$xi_X^2 + (1 - scenario$rho_X)^2 / 3
  # Moment condition checked through residual correlation/mean, not pointwise MC noise.
  res <- x2^2 - target
  c(mean_residual = mean(res), cor_residual_target = cor(res, target))
}
