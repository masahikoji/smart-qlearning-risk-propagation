# Population direction checks for the accumulation-orthogonal benchmarks.
# These are algebraic checks of the function directions, not finite-sample
# estimates of the random cross-covariance.

population_orthogonality_row <- function(scenario, target_stage, source_stage) {
  if (scenario$boundary != "separated" || scenario$candidate_count != 2L ||
      scenario$family != "standard") return(NULL)
  if (source_stage <= target_stage) return(NULL)
  mm <- state_even_moments(scenario$rho_X, scenario$xi_X, scenario$K)

  if (source_stage == target_stage + 1L) {
    # z_t=A_t(X_t^2-E X_t^2); transported residual direction from t+1 is
    # rho^2(X_t^2-E X_t^2)+2 rho xi X_t A_t. By symmetry and balanced
    # randomization their inner product is exactly zero.
    gt <- mm$m4[target_stage] - mm$m2[target_stage]^2
    transport_norm <- scenario$rho_X^4 * gt +
      4 * scenario$rho_X^2 * scenario$xi_X^2 * mm$m2[target_stage]
    return(data.frame(
      scenario=scenario$id, target_stage=target_stage, source_stage=source_stage,
      analytic_inner_product=0, target_direction_norm_sq=gt,
      transported_direction_norm_sq=transport_norm,
      rationale="balanced A and symmetric X: E[A q(X)^2]=E[X q(X)]=0",
      stringsAsFactors=FALSE))
  }

  if (scenario$K == 3L && target_stage == 1L && source_stage == 3L) {
    # After two transports, the stage-3 quadratic residual direction is a linear
    # combination of q1=X1^2-E X1^2, X1*A1, X1 and A1. Its inner product with
    # z1=A1*q1 is zero term by term under symmetry/balanced randomization.
    g1 <- mm$m4[1L] - mm$m2[1L]^2
    return(data.frame(
      scenario=scenario$id, target_stage=1L, source_stage=3L,
      analytic_inner_product=0, target_direction_norm_sq=g1,
      transported_direction_norm_sq=NA_real_,
      rationale="two-step transported polynomial is orthogonal to A1*(X1^2-E X1^2) term by term",
      stringsAsFactors=FALSE))
  }
  NULL
}

population_orthogonality_table <- function(K = NULL) {
  scens <- if (is.null(K)) c(smart_scenarios(2L), smart_scenarios(3L)) else smart_scenarios(K)
  rows <- list(); k <- 0L
  for (sc in scens) {
    if (sc$accumulation_role != "orthogonal") next
    for (t in seq_len(sc$K - 1L)) {
      for (s in (t + 1L):sc$K) {
        z <- population_orthogonality_row(sc, t, s)
        if (!is.null(z)) { k <- k + 1L; rows[[k]] <- z }
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

write_population_orthogonality_checks <- function() {
  out <- population_orthogonality_table()
  f <- file.path(smart_dir("results", "summary"), "population_orthogonality_checks.csv")
  write.csv(out, f, row.names=FALSE)
  invisible(out)
}
