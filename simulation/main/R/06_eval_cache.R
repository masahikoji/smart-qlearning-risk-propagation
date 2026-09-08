# Evaluation-sample caches and quadratic-form risk evaluation.

build_risk_cache <- function(scenario, n, stage, eval_data, quad) {
  spec <- stage_model_spec(scenario, stage)
  xp <- feature_matrix(eval_data, stage, 1, spec$union)
  xm <- feature_matrix(eval_data, stage, -1, spec$union)
  qp <- true_q(scenario, n, stage, eval_data, 1, quad)
  qm <- true_q(scenario, n, stage, eval_data, -1, quad)
  Phi <- rbind(xp, xm)
  qstar <- c(qp, qm)
  # q_t(a|h)=1/2 and N histories -> mean over the 2N stacked rows.
  G <- crossprod(Phi) / nrow(Phi)
  h <- as.numeric(crossprod(Phi, qstar) / nrow(Phi))
  names(h) <- spec$union
  cc <- mean(qstar^2)
  Cmat <- xp - xm
  true_diff <- qp - qm
  list(stage = stage, union = spec$union, G = G, h = h, c = cc,
       contrast_matrix = Cmat, true_diff = true_diff,
       true_opt = ifelse(true_diff >= 0, 1, -1))
}

risk_from_cache <- function(beta, cache) {
  b <- beta[cache$union]
  as.numeric(crossprod(b, cache$G %*% b) - 2 * crossprod(b, cache$h) + cache$c)
}

qdiff_rms_from_cache <- function(beta1, beta0, cache) {
  d <- beta1[cache$union] - beta0[cache$union]
  sqrt(max(0, as.numeric(crossprod(d, cache$G %*% d))))
}

build_evaluation_cache <- function(scenario, n, control) {
  seed <- evaluation_seed(control$seed, scenario$id, n)
  eval_data <- simulate_smart_data(scenario, control$N_eval, seed = seed,
                                   include_rewards = FALSE)
  quad <- gauss_legendre(control$quadrature_order)
  caches <- lapply(seq_len(scenario$K), function(s) build_risk_cache(scenario, n, s, eval_data, quad))
  names(caches) <- paste0("stage", seq_len(scenario$K))
  list(data = eval_data, risk = caches, quad = quad, seed = seed)
}
