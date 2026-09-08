# Pre-run calibration summaries and theoretical reference values.

terminal_added_g_exact <- function(scenario) {
  if (scenario$candidate_count != 2L) return(NA_real_)
  mm <- state_even_moments(scenario$rho_X, scenario$xi_X, scenario$K)
  mm$m4[scenario$K] - mm$m2[scenario$K]^2
}

terminal_delta_exact <- function(scenario) {
  if (scenario$departure != "local" || scenario$candidate_count != 2L) return(NA_real_)
  g <- terminal_added_g_exact(scenario)
  scenario$b[scenario$K] * sqrt(g) / scenario$sigma[scenario$K]
}

one_step_transport_factor_exact <- function(scenario) {
  if (scenario$K < 2L || scenario$candidate_count != 2L) return(NA_real_)
  mm <- state_even_moments(scenario$rho_X, scenario$xi_X, scenario$K)
  s <- scenario$K - 1L
  gprev <- mm$m4[s] - mm$m2[s]^2
  gK <- mm$m4[scenario$K] - mm$m2[scenario$K]^2
  num <- scenario$rho_X^4 * gprev +
    4 * scenario$rho_X^2 * scenario$xi_X^2 * mm$m2[s]
  num / gK
}

# One-dimensional AIC local-risk gap from the supplement.
bar_delta_aic <- function(delta, rel.tol = 1e-10) {
  if (!is.finite(delta)) return(NA_real_)
  f <- function(z, mode) {
    xi <- delta + z
    w <- if (mode == "SEL") as.numeric(xi^2 > 2) else 1/(1 + exp(1 - xi^2/2))
    (w * xi - delta)^2 * dnorm(z)
  }
  # Split hard-selection integration at its two discontinuity points.
  cut <- sort(c(-sqrt(2)-delta, sqrt(2)-delta))
  es <- integrate(function(z) f(z, "SEL"), -Inf, cut[1], rel.tol=rel.tol,
                  subdivisions=1000L)$value +
    integrate(function(z) f(z, "SEL"), cut[1], cut[2], rel.tol=rel.tol,
              subdivisions=1000L)$value +
    integrate(function(z) f(z, "SEL"), cut[2], Inf, rel.tol=rel.tol,
              subdivisions=1000L)$value
  em <- integrate(function(z) f(z, "MA"), -Inf, Inf, rel.tol = rel.tol,
                  subdivisions = 1000L)$value
  es - em
}

terminal_theory_n_gap <- function(scenario, criterion) {
  if (scenario$departure != "local" || scenario$candidate_count != 2L) return(NA_real_)
  delta <- terminal_delta_exact(scenario)
  if (!is.finite(delta)) return(NA_real_)
  if (criterion == "AIC") return(scenario$sigma[scenario$K]^2 * bar_delta_aic(delta))
  if (criterion == "BIC") return(0)
  NA_real_
}

terminal_added_direction_calibration <- function(scenario, n = NULL, N = NULL, seed = NULL) {
  g <- terminal_added_g_exact(scenario)
  delta <- terminal_delta_exact(scenario)
  data.frame(g_added = g, sqrt_g = sqrt(g), delta_terminal = delta,
             one_step_transport_factor = one_step_transport_factor_exact(scenario))
}

scenario_calibration_table <- function(K = NULL, n_values = c(250L, 500L, 1000L),
                                       N_direction = NULL, profile = "paper") {
  scens <- if (is.null(K)) c(smart_scenarios(2L), smart_scenarios(3L)) else smart_scenarios(K)
  ctl <- smart_run_control(profile)
  rows <- list(); k <- 0L
  for (sc in scens) {
    dircal <- terminal_added_direction_calibration(sc)
    for (n in n_values) {
      v <- validate_scenario(sc, as.integer(n), ctl, verbose = FALSE)
      k <- k + 1L
      rows[[k]] <- data.frame(
        scenario = sc$id, K = sc$K, n = n, boundary = sc$boundary,
        rho_X = sc$rho_X, xi_X = sc$xi_X,
        min_gap_s1 = v$min_gap[1L],
        min_gap_s2 = if (sc$K >= 2L) v$min_gap[2L] else NA_real_,
        min_gap_s3 = if (sc$K >= 3L) v$min_gap[3L] else NA_real_,
        g_added_terminal = dircal$g_added,
        delta_terminal = dircal$delta_terminal,
        one_step_transport_factor = dircal$one_step_transport_factor,
        theory_n_delta_AIC_terminal = terminal_theory_n_gap(sc, "AIC"),
        theory_n_delta_BIC_terminal = terminal_theory_n_gap(sc, "BIC"),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

scenario_manifest_table <- function(K = NULL) {
  scens <- if (is.null(K)) c(smart_scenarios(2L), smart_scenarios(3L)) else smart_scenarios(K)
  rows <- lapply(scens, function(sc) {
    data.frame(
      scenario = sc$id, K = sc$K, family = sc$family,
      candidate_count = sc$candidate_count, boundary = sc$boundary,
      departure = sc$departure, rho_X = sc$rho_X, xi_X = sc$xi_X,
      b = paste(format(sc$b, digits=10), collapse = ";"),
      target_delta_terminal = sc$target_delta_terminal,
      delta_terminal_exact = terminal_delta_exact(sc),
      c_fixed = paste(sc$c_fixed, collapse = ";"),
      tau0 = paste(sc$tau0, collapse = ";"),
      tau1 = paste(sc$tau1, collapse = ";"),
      sigma = paste(sc$sigma, collapse = ";"),
      criterion_variant = ifelse(is.na(sc$criterion_variant), "", sc$criterion_variant),
      accumulation_role = sc$accumulation_role,
      notes = sc$notes,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
