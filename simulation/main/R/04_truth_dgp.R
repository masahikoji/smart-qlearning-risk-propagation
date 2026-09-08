# Data-generating mechanism and true Q/V functions.

quad_value <- function(x, a0, a1, a2) a0 + a1 * x + a2 * x^2
quad_antiderivative <- function(x, a0, a1, a2) a0 * x + 0.5 * a1 * x^2 + (a2 / 3) * x^3

real_quad_roots <- function(a0, a1, a2, tol = 1e-12) {
  if (abs(a2) < tol) {
    if (abs(a1) < tol) return(numeric(0))
    return(-a0 / a1)
  }
  disc <- a1^2 - 4 * a2 * a0
  if (disc < -tol) return(numeric(0))
  if (abs(disc) <= tol) return(-a1 / (2 * a2))
  s <- sqrt(max(0, disc))
  sort(c((-a1 - s) / (2 * a2), (-a1 + s) / (2 * a2)))
}

# E|a0 + a1 X + a2 X^2| for X ~ Uniform(center-halfwidth, center+halfwidth).
# Roots are global for a given polynomial, so the calculation vectorizes over center.
mean_abs_quadratic_uniform <- function(center, halfwidth, a0, a1, a2) {
  if (halfwidth <= 0) stop("halfwidth must be positive")
  l <- center - halfwidth
  u <- center + halfwidth
  roots <- real_quad_roots(a0, a1, a2)
  cuts <- c(-Inf, roots, Inf)
  integ <- rep(0, length(center))
  for (j in seq_len(length(cuts) - 1L)) {
    lo <- pmax(l, cuts[j])
    hi <- pmin(u, cuts[j + 1L])
    ok <- hi > lo
    if (!any(ok)) next
    mid <- (lo[ok] + hi[ok]) / 2
    sg <- sign(quad_value(mid, a0, a1, a2))
    sg[sg == 0] <- 1
    piece <- sg * (quad_antiderivative(hi[ok], a0, a1, a2) -
                     quad_antiderivative(lo[ok], a0, a1, a2))
    # Numerical roundoff can make a tiny negative value at a tangential root.
    integ[ok] <- integ[ok] + pmax(piece, 0)
  }
  integ / (2 * halfwidth)
}

stage_reward_mean <- function(scenario, n, stage, x, past_actions, action) {
  cvec <- resolve_c(scenario, n)
  out <- scenario$beta0[stage] + scenario$beta1[stage] * x
  if (stage > 1L) {
    bp <- scenario$beta_past[[stage]]
    if (is.null(dim(past_actions))) {
      pa <- matrix(past_actions, ncol = stage - 1L, byrow = TRUE)
    } else pa <- past_actions
    if (nrow(pa) == 1L && length(x) > 1L) pa <- pa[rep(1L, length(x)), , drop = FALSE]
    out <- out + as.numeric(pa %*% bp)
  }
  out + action * (scenario$tau0[stage] + scenario$tau1[stage] * x + cvec[stage] * x^2)
}

true_q2_terminal <- function(scenario, n, x2, a1, a2) {
  stopifnot(scenario$K == 2L)
  past <- cbind(a1)
  stage_reward_mean(scenario, n, 2L, x2, past, a2)
}

true_v2_terminal <- function(scenario, n, x2, a1) {
  qplus <- true_q2_terminal(scenario, n, x2, a1, 1)
  qminus <- true_q2_terminal(scenario, n, x2, a1, -1)
  pmax(qplus, qminus)
}

true_q1_two <- function(scenario, n, x1, a1) {
  cvec <- resolve_c(scenario, n)
  y1 <- stage_reward_mean(scenario, n, 1L, x1, numeric(0), a1)
  center <- scenario$rho_X * x1 + scenario$xi_X * a1
  h <- 1 - scenario$rho_X
  # V2 = beta20 + beta21 X2 + beta22 A1 + |tau20+tau21 X2+c2 X2^2|.
  ev2 <- scenario$beta0[2L] + scenario$beta1[2L] * center +
    scenario$beta_past[[2L]][1L] * a1 +
    mean_abs_quadratic_uniform(center, h, scenario$tau0[2L],
                               scenario$tau1[2L], cvec[2L])
  y1 + ev2
}

true_q3_terminal <- function(scenario, n, x3, a1, a2, a3) {
  stopifnot(scenario$K == 3L)
  stage_reward_mean(scenario, n, 3L, x3, cbind(a1, a2), a3)
}

true_v3_terminal <- function(scenario, n, x3, a1, a2) {
  pmax(true_q3_terminal(scenario, n, x3, a1, a2, 1),
       true_q3_terminal(scenario, n, x3, a1, a2, -1))
}

true_q2_three <- function(scenario, n, x2, a1, a2) {
  cvec <- resolve_c(scenario, n)
  y2 <- stage_reward_mean(scenario, n, 2L, x2, cbind(a1), a2)
  center3 <- scenario$rho_X * x2 + scenario$xi_X * a2
  h <- 1 - scenario$rho_X
  # E[V3 | H2,A2] is analytic because V3 is linear prognostic part + abs(quadratic).
  ev3 <- scenario$beta0[3L] + scenario$beta1[3L] * center3 +
    scenario$beta_past[[3L]][1L] * a1 + scenario$beta_past[[3L]][2L] * a2 +
    mean_abs_quadratic_uniform(center3, h, scenario$tau0[3L],
                               scenario$tau1[3L], cvec[3L])
  y2 + ev3
}

true_v2_three <- function(scenario, n, x2, a1) {
  pmax(true_q2_three(scenario, n, x2, a1, 1),
       true_q2_three(scenario, n, x2, a1, -1))
}

true_q1_three <- function(scenario, n, x1, a1, quad = gauss_legendre(32L),
                          chunk = 10000L) {
  if (length(a1) == 1L) a1 <- rep(a1, length(x1))
  if (length(a1) != length(x1)) stop("a1 length mismatch in true_q1_three")
  y1 <- stage_reward_mean(scenario, n, 1L, x1, numeric(0), a1)
  center2 <- scenario$rho_X * x1 + scenario$xi_X * a1
  h <- 1 - scenario$rho_X
  ans <- numeric(length(x1))
  ids <- split(seq_along(x1), ceiling(seq_along(x1) / chunk))
  for (ii in ids) {
    cen <- center2[ii]
    # Matrix rows = histories, columns = quadrature nodes.
    x2mat <- outer(cen, rep(1, length(quad$nodes))) +
      h * outer(rep(1, length(cen)), quad$nodes)
    a1mat <- matrix(a1[ii], nrow = length(ii), ncol = length(quad$nodes))
    vp <- true_q2_three(scenario, n, as.vector(x2mat), as.vector(a1mat), 1)
    vm <- true_q2_three(scenario, n, as.vector(x2mat), as.vector(a1mat), -1)
    vmat <- matrix(pmax(vp, vm), nrow = length(ii), ncol = length(quad$nodes))
    ans[ii] <- as.numeric(vmat %*% quad$prob_weights)
  }
  y1 + ans
}

true_q <- function(scenario, n, stage, data, action, quad = gauss_legendre(32L)) {
  if (scenario$K == 2L) {
    if (stage == 2L) return(true_q2_terminal(scenario, n, data$X2, data$A1, action))
    if (stage == 1L) return(true_q1_two(scenario, n, data$X1, action))
  } else {
    if (stage == 3L) return(true_q3_terminal(scenario, n, data$X3, data$A1, data$A2, action))
    if (stage == 2L) return(true_q2_three(scenario, n, data$X2, data$A1, action))
    if (stage == 1L) return(true_q1_three(scenario, n, data$X1, action, quad = quad))
  }
  stop("Invalid stage/K combination")
}

true_v <- function(scenario, n, stage, data, quad = gauss_legendre(32L)) {
  pmax(true_q(scenario, n, stage, data, 1, quad),
       true_q(scenario, n, stage, data, -1, quad))
}

true_opt_action <- function(scenario, n, stage, data, quad = gauss_legendre(32L), tie = 1) {
  qp <- true_q(scenario, n, stage, data, 1, quad)
  qm <- true_q(scenario, n, stage, data, -1, quad)
  d <- qp - qm
  out <- ifelse(d > 0, 1, ifelse(d < 0, -1, tie))
  as.numeric(out)
}

simulate_smart_data <- function(scenario, n, seed = NULL, include_rewards = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  K <- scenario$K
  dat <- data.frame(X1 = runif(n, -1, 1))
  for (s in seq_len(K)) {
    dat[[paste0("A", s)]] <- sample(c(-1, 1), n, replace = TRUE)
    if (include_rewards) {
      past <- if (s == 1L) numeric(0) else as.matrix(dat[paste0("A", seq_len(s - 1L))])
      mu <- stage_reward_mean(scenario, n, s, dat[[paste0("X", s)]], past,
                              dat[[paste0("A", s)]])
      dat[[paste0("Y", s)]] <- mu + rnorm(n, 0, scenario$sigma[s])
    }
    if (s < K) {
      u <- runif(n, -1, 1)
      dat[[paste0("U", s + 1L)]] <- u
      dat[[paste0("X", s + 1L)]] <- scenario$rho_X * dat[[paste0("X", s)]] +
        (1 - scenario$rho_X) * u + scenario$xi_X * dat[[paste0("A", s)]]
    }
  }
  dat
}

# Exact conditional support of X_s given a fixed past-action sequence.
state_support_given_past <- function(scenario, stage, past_actions) {
  l <- -1; u <- 1
  if (stage == 1L) return(c(l, u))
  stopifnot(length(past_actions) == stage - 1L)
  for (q in seq_len(stage - 1L)) {
    lnew <- scenario$rho_X * l - (1 - scenario$rho_X) + scenario$xi_X * past_actions[q]
    unew <- scenario$rho_X * u + (1 - scenario$rho_X) + scenario$xi_X * past_actions[q]
    l <- lnew; u <- unew
  }
  c(l, u)
}

# ---- Exact conditional polynomial integrations used by decomposition diagnostics ----

coef_or_zero <- function(beta, nm) {
  if (nm %in% names(beta) && is.finite(beta[[nm]])) as.numeric(beta[[nm]]) else 0
}

mean_poly_uniform <- function(center, halfwidth, coef) {
  # coef[k+1] multiplies x^k, k=0,...,3 (higher entries are ignored here).
  c0 <- if (length(coef) >= 1L) coef[1L] else 0
  c1 <- if (length(coef) >= 2L) coef[2L] else 0
  c2 <- if (length(coef) >= 3L) coef[3L] else 0
  c3 <- if (length(coef) >= 4L) coef[4L] else 0
  ex1 <- center
  ex2 <- center^2 + halfwidth^2 / 3
  ex3 <- center^3 + center * halfwidth^2
  c0 + c1 * ex1 + c2 * ex2 + c3 * ex3
}

mean_quadratic_uniform <- function(center, halfwidth, a0, a1, a2) {
  a0 + a1 * center + a2 * (center^2 + halfwidth^2 / 3)
}

# E[sign(C_*(X)) C_hat(X)] for X uniform on a moving interval.
# Both C_* and C_hat are quadratic. The integral is exact after splitting at
# the real roots of C_*. If C_* is identically zero, tie_sign is used.
mean_signed_quadratic_uniform <- function(center, halfwidth,
                                          chat0, chat1, chat2,
                                          cstar0, cstar1, cstar2,
                                          tie_sign = 1, tol = 1e-12) {
  if (halfwidth <= 0) stop("halfwidth must be positive")
  center <- as.numeric(center)
  true_is_tie <- max(abs(c(cstar0, cstar1, cstar2)), na.rm = TRUE) < tol
  if (true_is_tie) {
    return(tie_sign * mean_quadratic_uniform(center, halfwidth,
                                             chat0, chat1, chat2))
  }

  l <- center - halfwidth
  u <- center + halfwidth
  roots <- real_quad_roots(cstar0, cstar1, cstar2, tol = tol)
  cuts <- c(-Inf, roots[is.finite(roots)], Inf)
  integ <- numeric(length(center))
  for (j in seq_len(length(cuts) - 1L)) {
    lo0 <- cuts[j]
    hi0 <- cuts[j + 1L]
    lo <- pmax(l, lo0)
    hi <- pmin(u, hi0)
    ok <- hi > lo
    if (!any(ok)) next
    mid <- if (is.finite(lo0) && is.finite(hi0)) {
      (lo0 + hi0) / 2
    } else if (!is.finite(lo0) && is.finite(hi0)) {
      hi0 - 1
    } else if (is.finite(lo0) && !is.finite(hi0)) {
      lo0 + 1
    } else 0
    sg <- sign(quad_value(mid, cstar0, cstar1, cstar2))
    if (sg == 0) sg <- tie_sign
    integ[ok] <- integ[ok] + sg *
      (quad_antiderivative(hi[ok], chat0, chat1, chat2) -
         quad_antiderivative(lo[ok], chat0, chat1, chat2))
  }
  integ / (2 * halfwidth)
}

poly4_antiderivative <- function(x, b0, b1, b2, b3, b4) {
  b0*x + b1*x^2/2 + b2*x^3/3 + b3*x^4/4 + b4*x^5/5
}

# Exact conditional mean of the binary-treatment switching remainder when both
# the fitted half-contrast C_hat(x) and true half-contrast C_*(x) are quadratic.
# The integral is split at all real roots, so the O(n^{-1/2}) boundary bump is
# resolved exactly rather than by a fixed quadrature grid.
mean_switching_remainder_uniform <- function(center, halfwidth,
                                             chat0, chat1, chat2,
                                             cstar0 = NULL, cstar1 = NULL, cstar2 = NULL,
                                             true_sign = NULL, tie_sign = 1,
                                             return_square = FALSE,
                                             tol = 1e-12) {
  if (halfwidth <= 0) stop("halfwidth must be positive")
  center <- as.numeric(center)
  if (!is.null(true_sign)) {
    if (length(true_sign) == 1L) true_sign <- rep(true_sign, length(center))
    if (length(true_sign) != length(center)) stop("true_sign length mismatch")
    eabs <- mean_abs_quadratic_uniform(center, halfwidth, chat0, chat1, chat2)
    emean <- mean_quadratic_uniform(center, halfwidth, chat0, chat1, chat2)
    ans <- pmax(eabs - true_sign * emean, 0)
    if (!return_square) return(ans)
    # For squared rho with constant true sign we still split at C_hat roots.
    cstar0 <- cstar1 <- cstar2 <- NULL
  }

  true_is_tie <- is.null(cstar0) ||
    max(abs(c(cstar0, cstar1, cstar2)), na.rm = TRUE) < tol
  roots_hat <- real_quad_roots(chat0, chat1, chat2, tol = tol)
  roots_true <- if (true_is_tie) numeric(0) else real_quad_roots(cstar0, cstar1, cstar2, tol = tol)
  roots <- sort(unique(c(roots_hat[is.finite(roots_hat)], roots_true[is.finite(roots_true)])))
  cuts <- c(-Inf, roots, Inf)
  l <- center - halfwidth; u <- center + halfwidth
  integ <- numeric(length(center))

  b0 <- chat0^2
  b1 <- 2*chat0*chat1
  b2 <- chat1^2 + 2*chat0*chat2
  b3 <- 2*chat1*chat2
  b4 <- chat2^2

  for (j in seq_len(length(cuts)-1L)) {
    lo0 <- cuts[j]; hi0 <- cuts[j+1L]
    lo <- pmax(l, lo0); hi <- pmin(u, hi0)
    ok <- hi > lo
    if (!any(ok)) next
    mid <- if (is.finite(lo0) && is.finite(hi0)) {
      (lo0 + hi0)/2
    } else if (!is.finite(lo0) && is.finite(hi0)) {
      hi0 - 1
    } else if (is.finite(lo0) && !is.finite(hi0)) {
      lo0 + 1
    } else 0
    sh <- sign(quad_value(mid, chat0, chat1, chat2)); if (sh == 0) sh <- 1
    st <- if (true_is_tie) tie_sign else {
      q <- sign(quad_value(mid, cstar0, cstar1, cstar2)); if (q == 0) tie_sign else q
    }
    k <- sh - st
    if (abs(k) < tol) next
    if (!return_square) {
      piece <- k * (quad_antiderivative(hi[ok], chat0, chat1, chat2) -
                      quad_antiderivative(lo[ok], chat0, chat1, chat2))
    } else {
      piece <- k^2 * (poly4_antiderivative(hi[ok], b0,b1,b2,b3,b4) -
                        poly4_antiderivative(lo[ok], b0,b1,b2,b3,b4))
    }
    integ[ok] <- integ[ok] + piece
  }
  pmax(integ / (2*halfwidth), 0)
}
