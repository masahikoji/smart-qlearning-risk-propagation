# Independent algebraic/numerical checks used before launching a long run.

check_abs_quadratic_integral <- function(tol = 1e-8) {
  centers <- c(-0.7, -0.2, 0, 0.4, 0.9)
  h <- 0.6; a0 <- -0.15; a1 <- 0.7; a2 <- 0.4
  fast <- mean_abs_quadratic_uniform(centers, h, a0, a1, a2)
  slow <- vapply(centers, function(m) {
    integrate(function(x) abs(a0 + a1*x + a2*x^2), m-h, m+h,
              subdivisions = 500L, rel.tol = 1e-12)$value / (2*h)
  }, numeric(1))
  err <- max(abs(fast - slow))
  if (err > tol) stop("Absolute-quadratic integral check failed: ", err)
  err
}

check_E3_constants <- function(tol = 5e-10) {
  rho <- 0.5; omega <- 0.3; su <- 1-rho
  # X2=rho X1+su U2+omega A1, with independent centered uniforms and A1^2=1.
  EX2_2 <- (rho^2 + su^2)/3 + omega^2
  # Fourth moment of a sum of independent centered variables.
  EX2_4 <- rho^4/5 + su^4/5 + omega^4 +
    6*(rho^2/3)*(su^2/3) + 6*(rho^2/3)*omega^2 + 6*(su^2/3)*omega^2
  g2 <- EX2_4 - EX2_2^2
  lambda2 <- ((4/3)*rho^2*omega^2)/g2
  target_g2 <- 89/900
  target_lambda <- 27/89
  if (abs(g2-target_g2)>tol || abs(lambda2-target_lambda)>tol) {
    stop(sprintf("E3 constants failed: g2=%.12f lambda=%.12f", g2, lambda2))
  }
  c(g2 = g2, lambda_pi_sq = lambda2)
}


check_signed_quadratic_integral <- function(tol = 1e-9) {
  centers <- c(-0.65, -0.1, 0.35, 0.8)
  h <- 0.55
  ch <- c(0.08, -0.45, 0.22)
  cs <- c(-0.03, 0.9, 0.12)
  fast <- mean_signed_quadratic_uniform(centers, h,
                                        ch[1], ch[2], ch[3],
                                        cs[1], cs[2], cs[3])
  roots <- real_quad_roots(cs[1], cs[2], cs[3])
  numeric_ref <- function(m) {
    lo <- m - h; hi <- m + h
    rr <- roots[is.finite(roots) & roots > lo & roots < hi]
    edges <- c(lo, rr, hi)
    val <- 0
    for (j in seq_len(length(edges) - 1L)) {
      a <- edges[j]; b <- edges[j + 1L]
      if (!(b > a)) next
      val <- val + integrate(function(x) {
        chat <- ch[1] + ch[2]*x + ch[3]*x^2
        cstar <- cs[1] + cs[2]*x + cs[3]*x^2
        ifelse(cstar >= 0, chat, -chat)
      }, a, b, rel.tol = 1e-12, abs.tol = 1e-14,
      subdivisions = 1000L, stop.on.error = TRUE)$value
    }
    val / (2*h)
  }
  slow <- vapply(centers, numeric_ref, numeric(1))
  err <- max(abs(fast - slow))
  if (err > tol) stop("Signed-quadratic integral check failed: ", err)
  err
}

run_theory_checks <- function() {
  a <- check_abs_quadratic_integral()
  b <- check_E3_constants()
  c <- check_switching_integral()
  d <- check_aic_gap_constants()
  e <- check_signed_quadratic_integral()
  log_msg("Theory checks passed. abs-quadratic max error=", format(a, scientific = TRUE),
          "; switching max error=", format(c, scientific = TRUE),
          "; signed-quadratic max error=", format(e, scientific = TRUE),
          "; E3 g2=", format(b["g2"], digits = 10),
          "; lambda_pi^2=", format(b["lambda_pi_sq"], digits = 10),
          "; barDelta(1)=", format(d["delta1"], digits=8),
          "; barDelta(3)=", format(d["delta3"], digits=8))
  invisible(list(abs_quad_error = a, switching_error=c, signed_quad_error=e,
                 E3 = b, AIC_gap=d))
}

check_switching_integral <- function(tol = 1e-9) {
  centers <- c(-0.2, 0, 0.25)
  h <- 0.7
  ch <- c(0.03, 0.9, 0.15)
  cs <- c(0, 1.0, 0.1)
  fast <- mean_switching_remainder_uniform(centers, h, ch[1], ch[2], ch[3],
                                           cs[1], cs[2], cs[3])
  fast_sq <- mean_switching_remainder_uniform(centers, h, ch[1], ch[2], ch[3],
                                              cs[1], cs[2], cs[3],
                                              return_square = TRUE)

  # A single adaptive integrate() call can completely miss the narrow interval
  # between a fitted-contrast root and a true-contrast root.  That is exactly the
  # boundary phenomenon this routine is meant to test.  Build an independent
  # numerical reference by splitting only at the polynomial roots and integrating
  # the original positive-part integrand on each resulting interval.  The actual
  # production routine above uses polynomial antiderivatives, not integrate().
  all_roots <- sort(unique(c(real_quad_roots(ch[1], ch[2], ch[3]),
                             real_quad_roots(cs[1], cs[2], cs[3]))))
  numeric_ref <- function(m, square = FALSE) {
    lo <- m - h; hi <- m + h
    rr <- all_roots[is.finite(all_roots) & all_roots > lo & all_roots < hi]
    edges <- c(lo, rr, hi)
    val <- 0
    for (j in seq_len(length(edges) - 1L)) {
      a <- edges[j]; b <- edges[j + 1L]
      if (!(b > a)) next
      val <- val + integrate(function(x) {
        chat <- ch[1] + ch[2]*x + ch[3]*x^2
        cstar <- cs[1] + cs[2]*x + cs[3]*x^2
        st <- ifelse(cstar >= 0, 1, -1)
        rr <- pmax(abs(chat) - st*chat, 0)
        if (square) rr^2 else rr
      }, a, b, rel.tol = 1e-12, abs.tol = 1e-14,
      subdivisions = 1000L, stop.on.error = TRUE)$value
    }
    val / (2*h)
  }
  slow <- vapply(centers, numeric_ref, numeric(1), square = FALSE)
  slow_sq <- vapply(centers, numeric_ref, numeric(1), square = TRUE)
  err <- max(abs(fast-slow), abs(fast_sq-slow_sq))
  if (err > tol) {
    stop(sprintf(paste0("Exact switching integral check failed: %.16g; ",
                        "mean exact=%s; mean split-numeric=%s; ",
                        "square exact=%s; square split-numeric=%s"),
                 err, paste(format(fast, digits=16), collapse=","),
                 paste(format(slow, digits=16), collapse=","),
                 paste(format(fast_sq, digits=16), collapse=","),
                 paste(format(slow_sq, digits=16), collapse=",")))
  }
  err
}


check_aic_gap_constants <- function(tol = 2e-6) {
  d1 <- bar_delta_aic(1)
  d3 <- bar_delta_aic(3)
  if (abs(d1-0.3829919)>tol || abs(d3+0.2118179)>tol) {
    stop(sprintf("AIC gap constants failed: delta1=%.8f delta3=%.8f", d1, d3))
  }
  c(delta1=d1, delta3=d3)
}
