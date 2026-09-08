# Additional simulations for the SMART model-uncertainty manuscript.

script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

env_int <- function(name, default) {
  z <- Sys.getenv(name, unset = "")
  if (!nzchar(z)) return(as.integer(default))
  out <- suppressWarnings(as.integer(z))
  if (is.na(out)) stop(name, " must be an integer")
  out
}

env_num <- function(name, default) {
  z <- Sys.getenv(name, unset = "")
  if (!nzchar(z)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(z))
  if (!is.finite(out)) stop(name, " must be numeric")
  out
}

env_chr <- function(name, default) {
  z <- Sys.getenv(name, unset = "")
  if (nzchar(z)) z else default
}

safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

add_control <- function(profile = c("paper", "smoke", "debug")) {
  profile <- match.arg(profile)
  nphys <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (length(nphys) != 1L || is.na(nphys) || nphys < 1L) nphys <- 1L
  if (profile == "paper") {
    out <- list(
      n_values = c(250L, 500L, 1000L),
      R_smart = 20000L,
      R_feedback = 10000L,
      R_finite = 20000L,
      R_normal = 3000000L,
      R_theory = 2000000L,
      chunk_size = 250L,
      normal_chunk = 100000L,
      cores = max(1L, min(8L, nphys)),
      seed = 260908L
    )
  } else if (profile == "smoke") {
    out <- list(
      n_values = c(250L, 1000L),
      R_smart = 100L,
      R_feedback = 100L,
      R_finite = 100L,
      R_normal = 20000L,
      R_theory = 30000L,
      chunk_size = 25L,
      normal_chunk = 10000L,
      cores = 1L,
      seed = 260908L
    )
  } else {
    out <- list(
      n_values = 250L,
      R_smart = 10L,
      R_feedback = 10L,
      R_finite = 10L,
      R_normal = 5000L,
      R_theory = 5000L,
      chunk_size = 5L,
      normal_chunk = 5000L,
      cores = 1L,
      seed = 260908L
    )
  }
  out$R_smart <- env_int("SMART_ADD_R", out$R_smart)
  out$R_feedback <- env_int("SMART_ADD_R_FEEDBACK", out$R_feedback)
  out$R_finite <- env_int("SMART_ADD_R_FINITE", out$R_finite)
  out$R_normal <- env_int("SMART_ADD_R_NORMAL", out$R_normal)
  out$R_theory <- env_int("SMART_ADD_R_THEORY", out$R_theory)
  out$chunk_size <- env_int("SMART_ADD_CHUNK", out$chunk_size)
  out$normal_chunk <- env_int("SMART_ADD_NORMAL_CHUNK", out$normal_chunk)
  out$cores <- env_int("SMART_ADD_CORES", out$cores)
  out
}

stable_logistic <- function(x) {
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  ex <- exp(x[!pos])
  out[!pos] <- ex / (1 + ex)
  out
}

stable_softmax <- function(logw) {
  m <- max(logw)
  z <- exp(logw - m)
  z / sum(z)
}

mc_mean_se <- function(x) {
  ok <- is.finite(x)
  x <- x[ok]
  if (!length(x)) return(c(mean = NA_real_, se = NA_real_, sd = NA_real_, n = 0))
  c(mean = mean(x), se = if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_,
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_, n = length(x))
}

# Pretest--Stein weights and exact normal-mean identity

stein_weight <- function(F, d, c = 2 * d, a = d - 2) {
  if (d < 3L) stop("Pretest--Stein dominance requires d >= 3")
  if (!(a > 0 && a <= min(c, 2 * (d - 2)))) {
    stop("a must satisfy 0 < a <= min(c, 2(d-2))")
  }
  out <- numeric(length(F))
  idx <- is.finite(F) & F > c
  out[idx] <- 1 - a / F[idx]
  out
}

selection_weight <- function(F, d, c = 2 * d) as.numeric(F > c)
akaike_weight <- function(F, d) stable_logistic(F / 2 - d)

inv_tail_moment_chisq <- function(c, df, ncp) {
  stats::integrate(
    function(x) stats::dchisq(x, df = df, ncp = ncp) / x,
    lower = c, upper = Inf, rel.tol = 1e-10, subdivisions = 2000L,
    stop.on.error = TRUE
  )$value
}

stein_exact_risk_difference <- function(d, delta_norm, a, c = 2 * d) {
  if (d < 3L) return(NA_real_)
  if (!(a > 0 && a <= min(c, 2 * (d - 2)))) return(NA_real_)
  e_inv <- inv_tail_moment_chisq(c, d, delta_norm^2)
  # By the coarea formula, (2a/sqrt(c)) int_{||x||^2=c} phi_delta dS
  # equals 4a times the noncentral chi-square density at c.
  a * (a - 2 * (d - 2)) * e_inv -
    4 * a * stats::dchisq(c, df = d, ncp = delta_norm^2)
}

make_delta <- function(delta_norm, d, direction = c("maxeig", "balanced", "mineig")) {
  direction <- match.arg(direction)
  if (delta_norm == 0) return(rep(0, d))
  if (direction == "balanced") return(rep(delta_norm / sqrt(d), d))
  out <- rep(0, d)
  out[if (direction == "maxeig") 1L else d] <- delta_norm
  out
}

normal_mean_chunk <- function(d, delta, nmc, seed, c, a_values) {
  set.seed(seed)
  x <- matrix(stats::rnorm(nmc * d), nrow = nmc, ncol = d)
  x <- sweep(x, 2L, delta, "+")
  F <- rowSums(x^2)
  delta_mat <- matrix(delta, nrow = nmc, ncol = d, byrow = TRUE)
  rule_weights <- list(SEL = selection_weight(F, d, c), AIC_MA = akaike_weight(F, d))
  if (d >= 3L) {
    for (a in a_values) {
      if (a <= min(c, 2 * (d - 2))) {
        rule_weights[[paste0("ST_a", format(a, trim = TRUE, scientific = FALSE))]] <-
          stein_weight(F, d, c, a)
      }
    }
  }
  losses <- lapply(rule_weights, function(w) {
    err <- x * w - delta_mat
    rowSums(err^2)
  })
  sel_loss <- losses$SEL
  out <- lapply(names(rule_weights), function(nm) {
    w <- rule_weights[[nm]]
    loss <- losses[[nm]]
    gap <- sel_loss - loss
    data.frame(rule = nm, sum_loss = sum(loss), sum_loss2 = sum(loss^2),
               sum_gap = sum(gap), sum_gap2 = sum(gap^2), n = nmc,
               mean_weight = mean(w), wide_prob = mean(F > c), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

run_normal_mean_grid <- function(control, outdir) {
  safe_dir(outdir)
  dims <- c(1L, 3L, 5L)
  delta_grid <- c(0, 0.5, 1, 2, 3, 4, 5, 8)
  rows <- list(); idx <- 1L
  for (d in dims) {
    avec <- if (d == 3L) c(1, 2) else if (d == 5L) c(1, 3, 6) else numeric(0)
    cthr <- 2 * d
    for (dn in delta_grid) {
      delta <- make_delta(dn, d, "maxeig")
      nleft <- control$R_normal
      part <- list(); p <- 1L
      while (nleft > 0L) {
        m <- min(control$normal_chunk, nleft)
        part[[p]] <- normal_mean_chunk(
          d, delta, m,
          control$seed + 100000L * d + 1000L * round(10 * dn) + p,
          cthr, avec
        )
        nleft <- nleft - m; p <- p + 1L
      }
      z <- do.call(rbind, part)
      for (rule in unique(z$rule)) {
        q <- z[z$rule == rule, , drop = FALSE]
        N <- sum(q$n)
        sm <- sum(q$sum_loss); sm2 <- sum(q$sum_loss2)
        mn <- sm / N
        vr <- max(0, (sm2 - N * mn^2) / max(1, N - 1))
        gmn <- sum(q$sum_gap) / N
        gvr <- max(0, (sum(q$sum_gap2) - N * gmn^2) / max(1, N - 1))
        a_num <- if (grepl("^ST_a", rule)) as.numeric(sub("^ST_a", "", rule)) else NA_real_
        exact_st_minus_sel <- if (is.finite(a_num)) stein_exact_risk_difference(d, dn, a_num, cthr) else NA_real_
        rows[[idx]] <- data.frame(
          d = d, delta_norm = dn, c = cthr, rule = rule, a = a_num,
          risk = mn, mcse = sqrt(vr / N),
          SEL_minus_rule = gmn, mcse_SEL_minus_rule = sqrt(gvr / N),
          mean_weight = weighted.mean(q$mean_weight, q$n),
          wide_selection_probability = weighted.mean(q$wide_prob, q$n),
          exact_ST_minus_SEL = exact_st_minus_sel,
          exact_SEL_minus_ST = if (is.finite(exact_st_minus_sel)) -exact_st_minus_sel else NA_real_,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  ans <- do.call(rbind, rows)
  utils::write.csv(ans, file.path(outdir, "normal_mean_risk_grid.csv"), row.names = FALSE)
  ans
}


# Gaussian limit under a transported quadratic loss

transport_limit_chunk <- function(delta, M, nmc, seed, a_values, c = 2 * length(delta)) {
  d <- length(delta)
  set.seed(seed)
  xi <- matrix(stats::rnorm(nmc * d), nrow = nmc, ncol = d)
  xi <- sweep(xi, 2L, delta, "+")
  F <- rowSums(xi^2)
  dm <- matrix(delta, nrow = nmc, ncol = d, byrow = TRUE)
  weights <- list(SEL = selection_weight(F, d, c), AIC_MA = akaike_weight(F, d))
  if (d >= 3L) {
    for (a in unique(a_values)) {
      if (a > 0 && a <= min(c, 2 * (d - 2))) {
        weights[[paste0("ST_a", format(a, trim = TRUE, scientific = FALSE))]] <-
          stein_weight(F, d, c, a)
      }
    }
  }
  losses <- lapply(weights, function(w) {
    er <- xi * w - dm
    l2 <- rowSums(er^2)
    l1 <- rowSums((er %*% M) * er)
    list(stage2 = l2, stage1 = l1)
  })
  sel2 <- losses$SEL$stage2
  sel1 <- losses$SEL$stage1
  do.call(rbind, lapply(names(losses), function(nm) {
    z2 <- losses[[nm]]$stage2
    z1 <- losses[[nm]]$stage1
    g2 <- sel2 - z2
    g1 <- sel1 - z1
    data.frame(
      rule = nm, n = nmc,
      sum_l2 = sum(z2), sum2_l2 = sum(z2^2),
      sum_l1 = sum(z1), sum2_l1 = sum(z1^2),
      sum_g2 = sum(g2), sum2_g2 = sum(g2^2),
      sum_g1 = sum(g1), sum2_g1 = sum(g1^2),
      stringsAsFactors = FALSE
    )
  }))
}

transport_limit_theory <- function(delta, M, control, seed, a_values = c(1, 2)) {
  left <- control$R_theory; parts <- list(); k <- 1L
  while (left > 0L) {
    m <- min(control$normal_chunk, left)
    parts[[k]] <- transport_limit_chunk(delta, M, m, seed + k, a_values)
    left <- left - m; k <- k + 1L
  }
  z <- do.call(rbind, parts)
  moment_from_sums <- function(sm, sm2, nn) {
    N <- sum(nn); S <- sum(sm); S2 <- sum(sm2); mu <- S / N
    vr <- max(0, (S2 - N * mu^2) / max(1, N - 1))
    c(mean = mu, se = sqrt(vr / N))
  }
  do.call(rbind, lapply(unique(z$rule), function(rule) {
    q <- z[z$rule == rule, , drop = FALSE]
    l2 <- moment_from_sums(q$sum_l2, q$sum2_l2, q$n)
    l1 <- moment_from_sums(q$sum_l1, q$sum2_l1, q$n)
    g2 <- moment_from_sums(q$sum_g2, q$sum2_g2, q$n)
    g1 <- moment_from_sums(q$sum_g1, q$sum2_g1, q$n)
    data.frame(rule = rule,
      theory_stage2_block_risk = l2["mean"], theory_stage2_block_mcse = l2["se"],
      theory_stage1_transported_risk = l1["mean"], theory_stage1_transported_mcse = l1["se"],
      theory_stage2_SEL_minus_rule = g2["mean"], theory_stage2_gap_mcse = g2["se"],
      theory_stage1_SEL_minus_rule = g1["mean"], theory_stage1_gap_mcse = g1["se"],
      stringsAsFactors = FALSE)
  }))
}

# Two-stage SMART with a d-dimensional local added block

rho_from_transport_eigenvalue <- function(m) {
  if (m <= 0) return(0)
  if (m >= 1) return(1)
  k <- sqrt(m / (1 - m))
  k / (1 + k)
}

transport_profile <- function(d = 3L, target_reff = d, mmax = 0.8) {
  if (!(target_reff >= 1 && target_reff <= d)) stop("target_reff must be in [1,d]")
  if (!(mmax > 0 && mmax < 1)) stop("mmax must be in (0,1)")
  m <- rep(if (d > 1L) mmax * (target_reff - 1) / (d - 1) else mmax, d)
  m[1L] <- mmax
  if (target_reff == 1 && d > 1L) m[-1L] <- 0
  rho <- vapply(m, rho_from_transport_eigenvalue, numeric(1))
  reff <- sum(m) / max(m)
  list(m = m, rho = rho, reff = reff, M = diag(m, nrow = d))
}

stage2_sd <- function(rho) sqrt((rho^2 + (1 - rho)^2) / 3)

make_union_names <- function(d) {
  c("(Intercept)", paste0("Z", seq_len(d)), "A1", "A2", paste0("A2Z", seq_len(d)))
}

make_stage1_names <- function(d) c("(Intercept)", paste0("X1", seq_len(d)), "A1")

ols_fit <- function(X, y, tol = 1e-10) {
  q <- qr(X, tol = tol)
  if (q$rank < ncol(X)) stop("Rank-deficient design in additional simulation")
  G <- crossprod(X) / nrow(X)
  evmin <- tryCatch(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values),
                    error = function(e) -Inf)
  if (!is.finite(evmin) || evmin < tol) {
    stop("Empirical Gram matrix is numerically singular in additional simulation")
  }
  beta <- as.numeric(qr.coef(q, y))
  names(beta) <- colnames(X)
  res <- as.numeric(y - X %*% beta)
  list(beta = beta, rss = max(sum(res^2), .Machine$double.xmin),
       residuals = res, rank = q$rank, min_gram_eigenvalue = evmin)
}

embed_beta <- function(beta, union_names) {
  nm <- names(beta)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
    stop("All coefficient vectors must have nonempty names before embedding")
  }
  if (anyDuplicated(nm)) stop("Duplicate coefficient names before embedding")
  bad <- setdiff(nm, union_names)
  if (length(bad)) {
    stop("Coefficient name(s) outside the declared union: ", paste(bad, collapse = ", "))
  }
  out <- setNames(rep(0, length(union_names)), union_names)
  out[nm] <- beta
  out
}

squared_named_error <- function(beta, truth, label = "coefficient risk") {
  nb <- names(beta); nt <- names(truth)
  if (is.null(nb) || is.null(nt)) stop(label, ": both vectors must be named")
  if (anyDuplicated(nb) || anyDuplicated(nt)) stop(label, ": duplicate coefficient names")
  if (!setequal(nb, nt)) {
    miss <- setdiff(nt, nb); extra <- setdiff(nb, nt)
    stop(label, ": coefficient-name mismatch; missing={", paste(miss, collapse = ","),
         "}; extra={", paste(extra, collapse = ","), "}")
  }
  beta <- beta[nt]
  sum((beta - truth)^2)
}

simulate_transport_data <- function(n, d, delta, rho, sigma1 = 1, sigma2 = 1,
                                    tau1 = 1.25, tau2 = 2.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- matrix(stats::runif(n * d, -1, 1), nrow = n, ncol = d)
  x1s <- sqrt(3) * x1
  a1 <- sample(c(-1, 1), n, replace = TRUE)
  u2 <- matrix(stats::runif(n * d, -1, 1), nrow = n, ncol = d)
  x2 <- sweep(x1, 2L, rho, "*") + sweep(u2, 2L, 1 - rho, "*")
  s2 <- stage2_sd(rho)
  z2 <- sweep(x2, 2L, s2, "/")
  a2 <- sample(c(-1, 1), n, replace = TRUE)
  local_coef <- sigma2 * delta / sqrt(n)
  tr2 <- tau2 + as.numeric(z2 %*% local_coef)
  y2 <- a2 * tr2 + stats::rnorm(n, 0, sigma2)
  y1 <- tau1 * a1 + stats::rnorm(n, 0, sigma1)
  list(n = n, d = d, X1 = x1s, X2 = z2, A1 = a1, A2 = a2, Y1 = y1, Y2 = y2,
       delta = delta, rho = rho, m = rho^2 / (rho^2 + (1 - rho)^2),
       sigma1 = sigma1, sigma2 = sigma2, tau1 = tau1, tau2 = tau2,
       local_coef = local_coef, true_contrast_half = tr2)
}

stage2_designs <- function(dat, subset_k = NULL) {
  d <- dat$d
  Z <- dat$X2
  colnames(Z) <- paste0("Z", seq_len(d))
  AZ <- Z * dat$A2
  colnames(AZ) <- paste0("A2Z", seq_len(d))
  base <- cbind(`(Intercept)` = 1, Z, A1 = dat$A1, A2 = dat$A2)
  if (is.null(subset_k)) {
    wide <- cbind(base, AZ)
  } else if (subset_k <= 0L) {
    wide <- base
  } else {
    wide <- cbind(base, AZ[, seq_len(subset_k), drop = FALSE])
  }
  list(base = base, wide = wide, AZ = AZ)
}

fit_nested_stage2 <- function(dat) {
  ds <- stage2_designs(dat)
  f0 <- ols_fit(ds$base, dat$Y2)
  f1 <- ols_fit(ds$wide, dat$Y2)
  dadd <- dat$d
  lambda <- dat$n * log(f0$rss / f1$rss)
  lambda <- max(0, lambda)
  union <- make_union_names(dat$d)
  list(
    narrow = embed_beta(f0$beta, union),
    wide = embed_beta(f1$beta, union),
    rss0 = f0$rss, rss1 = f1$rss, lambda = lambda, dadd = dadd
  )
}

aggregate_stage2_beta <- function(fit, rule, a = NA_real_) {
  d <- fit$dadd; F <- fit$lambda; cthr <- 2 * d
  w <- switch(rule,
    SEL = selection_weight(F, d, cthr),
    AIC_MA = akaike_weight(F, d),
    ST = stein_weight(F, d, cthr, a),
    stop("Unknown stage-2 rule: ", rule)
  )
  beta <- fit$narrow + w * (fit$wide - fit$narrow)
  list(beta = beta, weight = as.numeric(w), lambda = F)
}

stage2_predict <- function(beta, Z, A1, A2) {
  d <- ncol(Z)
  out <- beta["(Intercept)"] + as.numeric(Z %*% beta[paste0("Z", seq_len(d))]) +
    beta["A1"] * A1 + beta["A2"] * A2
  out + as.numeric((Z * A2) %*% beta[paste0("A2Z", seq_len(d))])
}

stage2_true_beta <- function(dat) {
  out <- setNames(rep(0, length(make_union_names(dat$d))), make_union_names(dat$d))
  out["A2"] <- dat$tau2
  out[paste0("A2Z", seq_len(dat$d))] <- dat$local_coef
  out
}

stage2_integrated_risk <- function(beta, dat) {
  b0 <- stage2_true_beta(dat)
  # Population evaluation Gram is identity for the chosen standardized basis.
  squared_named_error(beta, b0, "stage-2 integrated risk")
}

make_stage1_pseudo <- function(dat, beta2) {
  qp <- stage2_predict(beta2, dat$X2, dat$A1, 1)
  qm <- stage2_predict(beta2, dat$X2, dat$A1, -1)
  dat$Y1 + pmax(qp, qm)
}

fit_stage1_common <- function(dat, pseudo) {
  X <- cbind(`(Intercept)` = 1, dat$X1, A1 = dat$A1)
  colnames(X)[seq.int(2L, 1L + dat$d)] <- paste0("X1", seq_len(dat$d))
  ols_fit(X, pseudo)$beta
}

stage1_true_beta <- function(dat) {
  out <- setNames(rep(0, length(make_stage1_names(dat$d))), make_stage1_names(dat$d))
  out["(Intercept)"] <- dat$tau2
  out["A1"] <- dat$tau1
  out[paste0("X1", seq_len(dat$d))] <-
    dat$sigma2 * dat$delta * sqrt(dat$m) / sqrt(dat$n)
  out
}

stage1_integrated_risk <- function(beta, dat) {
  b0 <- stage1_true_beta(dat)
  squared_named_error(beta, b0, "stage-1 integrated risk")
}

true_gap_validation <- function(n, d, delta, rho, tau2 = 2, sigma2 = 1) {
  sd2 <- stage2_sd(rho)
  # X2_j is supported on [-1,1], so |Z_j| <= 1/sd2_j.
  worst_local <- sigma2 / sqrt(n) * sum(abs(delta) / sd2)
  c(tau2 = tau2, worst_local = worst_local, minimum_half_contrast = tau2 - worst_local)
}

simulate_transport_replicate <- function(n, d, delta_norm, target_reff,
                                         direction = "maxeig", seed = 1L,
                                         a_values = c(1), tau2 = 2.0) {
  tp <- transport_profile(d, target_reff)
  delta <- make_delta(delta_norm, d, direction)
  gv <- true_gap_validation(n, d, delta, tp$rho, tau2 = tau2)
  if (gv["minimum_half_contrast"] <= 0) {
    stop("Scenario does not have a uniformly separated terminal treatment gap")
  }
  dat <- simulate_transport_data(n, d, delta, tp$rho, tau2 = tau2, seed = seed)
  fit <- fit_nested_stage2(dat)
  rules <- list(SEL = list(rule = "SEL", a = NA_real_),
                AIC_MA = list(rule = "AIC_MA", a = NA_real_))
  if (d >= 3L) {
    for (a in unique(a_values)) {
      if (a > 0 && a <= min(2 * d, 2 * (d - 2))) {
        rules[[paste0("ST_a", format(a, trim = TRUE, scientific = FALSE))]] <- list(rule = "ST", a = a)
      }
    }
  }
  ans <- list(); k <- 1L
  for (nm in names(rules)) {
    ag <- aggregate_stage2_beta(fit, rules[[nm]]$rule, rules[[nm]]$a)
    pseudo <- make_stage1_pseudo(dat, ag$beta)
    b1 <- fit_stage1_common(dat, pseudo)
    contrast_hat <- ag$beta["A2"] + as.numeric(dat$X2 %*% ag$beta[paste0("A2Z", seq_len(d))])
    ans[[k]] <- data.frame(
      rule = nm,
      risk_stage2 = stage2_integrated_risk(ag$beta, dat),
      risk_stage1 = stage1_integrated_risk(b1, dat),
      lambda2 = fit$lambda,
      wide_weight = ag$weight,
      fitted_switch_rate = mean(contrast_hat < 0),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  out <- do.call(rbind, ans)
  out$n <- n; out$d <- d; out$delta_norm <- delta_norm; out$direction <- direction
  out$target_reff <- target_reff; out$reff_exact <- tp$reff
  out$m_max <- max(tp$m); out$m_min <- min(tp$m)
  out$m_eigenvalues <- paste(format(tp$m, digits = 8), collapse = ";")
  out$rho <- paste(format(tp$rho, digits = 8), collapse = ";")
  out$minimum_half_contrast <- gv["minimum_half_contrast"]
  out
}

# Criterion-feedback negative control

fit_stage1_feedback <- function(dat, pseudo) {
  d <- dat$d
  # The narrow model deliberately excludes X1_1; the wide model adds it.
  keep <- if (d > 1L) 2:d else integer(0)
  X0 <- cbind(`(Intercept)` = 1,
              if (length(keep)) dat$X1[, keep, drop = FALSE] else NULL,
              A1 = dat$A1)
  if (length(keep)) colnames(X0)[seq.int(2L, 1L + length(keep))] <- paste0("X1", keep)
  X1 <- cbind(X0, dat$X1[, 1L])
  colnames(X1)[ncol(X1)] <- paste0("X1", 1L)
  f0 <- ols_fit(X0, pseudo); f1 <- ols_fit(X1, pseudo)
  lambda <- max(0, dat$n * log(f0$rss / f1$rss))
  use_wide <- as.numeric(lambda > 2)
  union <- make_stage1_names(d)
  b0 <- embed_beta(f0$beta, union)
  b1 <- embed_beta(f1$beta, union)
  beta <- b0 + use_wide * (b1 - b0)
  list(beta = beta, lambda = lambda, selected_wide = use_wide)
}

simulate_feedback_replicate <- function(n, delta_norm, target_reff, seed,
                                        a = 1, direction = "maxeig") {
  d <- 3L
  tp <- transport_profile(d, target_reff)
  delta <- make_delta(delta_norm, d, direction)
  dat <- simulate_transport_data(n, d, delta, tp$rho, seed = seed)
  fit2 <- fit_nested_stage2(dat)
  proc <- list(
    SEL = aggregate_stage2_beta(fit2, "SEL"),
    ST = aggregate_stage2_beta(fit2, "ST", a)
  )
  rows <- lapply(names(proc), function(nm) {
    ps <- make_stage1_pseudo(dat, proc[[nm]]$beta)
    f1 <- fit_stage1_feedback(dat, ps)
    data.frame(
      rule = nm,
      risk_stage1 = stage1_integrated_risk(f1$beta, dat),
      lambda1 = f1$lambda,
      stage1_selected_wide = f1$selected_wide,
      lambda2 = fit2$lambda,
      stage2_weight = proc[[nm]]$weight,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$n <- n; out$d <- d; out$delta_norm <- delta_norm; out$target_reff <- target_reff
  out$reff_exact <- tp$reff; out$a <- a; out
}

# Fixed finite-library AIC experiment

fit_finite_library_stage2 <- function(dat) {
  ds <- stage2_designs(dat)
  d <- dat$d
  union <- make_union_names(d)
  fits <- vector("list", d + 1L)
  names(fits) <- paste0("M", 0:d)
  for (k in 0:d) {
    X <- if (k == 0L) ds$base else cbind(ds$base, ds$AZ[, seq_len(k), drop = FALSE])
    f <- ols_fit(X, dat$Y2)
    fits[[k + 1L]] <- list(beta = embed_beta(f$beta, union), rss = f$rss, kadd = k)
  }
  rss0 <- fits[[1L]]$rss
  ic_rel <- vapply(fits, function(f) 2 * f$kadd - dat$n * log(rss0 / f$rss), numeric(1))
  w <- stable_softmax(-0.5 * ic_rel)
  j <- which.min(ic_rel)[1L]
  bsel <- fits[[j]]$beta
  bma <- Reduce(`+`, Map(function(f, ww) ww * f$beta, fits, w))
  list(fits = fits, ic_rel = ic_rel, weights = w, selected = names(fits)[j],
       beta_sel = bsel, beta_ma = bma)
}

simulate_finite_library_replicate <- function(n, delta, seed = 1L) {
  d <- length(delta)
  tp <- transport_profile(d, d)
  dat <- simulate_transport_data(n, d, delta, tp$rho, seed = seed)
  ff <- fit_finite_library_stage2(dat)
  data.frame(
    risk_SEL = stage2_integrated_risk(ff$beta_sel, dat),
    risk_MA = stage2_integrated_risk(ff$beta_ma, dat),
    selected = ff$selected,
    w_M0 = ff$weights[1L], w_M1 = if (d >= 1) ff$weights[2L] else NA_real_,
    w_M2 = if (d >= 2) ff$weights[3L] else NA_real_,
    w_M3 = if (d >= 3) ff$weights[4L] else NA_real_,
    stringsAsFactors = FALSE
  )
}

finite_library_limit_chunk <- function(delta, nmc, seed) {
  d <- length(delta)
  set.seed(seed)
  xi <- matrix(stats::rnorm(nmc * d), nrow = nmc, ncol = d)
  xi <- sweep(xi, 2L, delta, "+")
  # Nested chain M_k contains coordinates 1,...,k. Common AIC constants cancel.
  partial <- cbind(0, t(apply(xi^2, 1L, cumsum)))
  pen <- matrix(rep(2 * (0:d), each = nmc), nrow = nmc)
  ic <- pen - partial
  sel <- max.col(-ic, ties.method = "first") - 1L
  min_ic <- apply(ic, 1L, min)
  ww <- exp(-0.5 * (ic - min_ic))
  ww <- ww / rowSums(ww)
  est_sel <- matrix(0, nrow = nmc, ncol = d)
  est_ma <- matrix(0, nrow = nmc, ncol = d)
  for (j in seq_len(d)) {
    est_sel[, j] <- xi[, j] * as.numeric(sel >= j)
    est_ma[, j] <- xi[, j] * rowSums(ww[, (j + 1L):(d + 1L), drop = FALSE])
  }
  dm <- matrix(delta, nrow = nmc, ncol = d, byrow = TRUE)
  ls <- rowSums((est_sel - dm)^2)
  lm <- rowSums((est_ma - dm)^2)
  data.frame(sum_sel = sum(ls), sum2_sel = sum(ls^2), sum_ma = sum(lm), sum2_ma = sum(lm^2), n = nmc)
}

finite_library_limit <- function(delta, control, seed = 1L) {
  left <- control$R_theory; parts <- list(); k <- 1L
  while (left > 0L) {
    m <- min(control$normal_chunk, left)
    parts[[k]] <- finite_library_limit_chunk(delta, m, seed + k)
    left <- left - m; k <- k + 1L
  }
  z <- do.call(rbind, parts); N <- sum(z$n)
  summarize_sum <- function(sm, sm2) {
    mn <- sum(sm) / N
    vr <- max(0, (sum(sm2) - N * mn^2) / max(1, N - 1))
    c(mean = mn, se = sqrt(vr / N))
  }
  sel <- summarize_sum(z$sum_sel, z$sum2_sel)
  ma <- summarize_sum(z$sum_ma, z$sum2_ma)
  c(SEL = unname(sel["mean"]),
    SEL_se = unname(sel["se"]),
    MA = unname(ma["mean"]),
    MA_se = unname(ma["se"]))
}

# Chunked restartable runners

scenario_slug <- function(...) {
  x <- paste(..., sep = "__")
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

rep_seed <- function(base, key, i) {
  raw <- utf8ToInt(key)
  h <- if (length(raw)) sum((seq_along(raw) * raw) %% 1000003L) else 1L
  as.integer((base + h * 1009 + i * 104729) %% 2147483000L + 1L)
}

additional_engine_md5 <- function() {
  root <- script_dir()
  fs <- c(
    file.path(root, "R", "additional_sim_functions.R"),
    file.path(root, "run_additional_simulation.R"),
    file.path(root, "VERSION")
  )
  fs <- fs[file.exists(fs)]
  unname(as.character(tools::md5sum(fs)))
}

run_chunked <- function(key, R, chunk_size, cores, outdir, fun, overwrite = FALSE) {
  d <- safe_dir(file.path(outdir, "raw", key))
  metafile <- file.path(d, "metadata.rds")
  existing <- sort(list.files(d, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = TRUE))
  md5_now <- additional_engine_md5()

  if (overwrite) {
    if (length(existing)) unlink(existing, force = TRUE)
    if (file.exists(metafile)) unlink(metafile, force = TRUE)
    existing <- character(0)
  } else if (length(existing)) {
    if (!file.exists(metafile)) {
      stop("Existing chunks have no metadata for ", key, ". Rerun with overwrite=TRUE.")
    }
    old <- readRDS(metafile)
    if (!identical(as.integer(old$chunk_size), as.integer(chunk_size))) {
      stop("Chunk size changed for ", key, ". Keep SMART_ADD_CHUNK unchanged or rerun with overwrite=TRUE.")
    }
    if (!identical(old$engine_md5, md5_now)) {
      stop("Simulation source changed since chunks were created for ", key,
           ". Rerun with overwrite=TRUE to avoid mixing code versions.")
    }
  }

  saveRDS(list(key = key, chunk_size = as.integer(chunk_size), requested_R = as.integer(R),
               engine_md5 = md5_now, timestamp = as.character(Sys.time()),
               R_version = R.version.string), metafile)

  starts <- seq.int(1L, R, by = chunk_size)
  for (st in starts) {
    en <- min(R, st + chunk_size - 1L)
    f <- file.path(d, sprintf("chunk_%06d_%06d.rds", st, en))
    if (file.exists(f)) next
    existing_names <- list.files(d, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = FALSE)
    if (length(existing_names)) {
      rg <- strcapture("^chunk_([0-9]+)_([0-9]+)\\.rds$", existing_names,
                       proto = list(lo = integer(), hi = integer()))
      if (nrow(rg) && any(rg$lo <= en & rg$hi >= st)) {
        stop("Existing chunk boundaries overlap requested block ", st, "-", en,
             " for ", key, ". Rerun with overwrite=TRUE.")
      }
    }
    ids <- st:en
    one <- function(i) fun(i)
    z <- if (cores > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
    } else lapply(ids, one)
    z <- Map(function(obj, id) {
      obj$rep_id <- id
      obj
    }, z, ids)
    ans <- do.call(rbind, z)
    tmp <- paste0(f, ".tmp")
    saveRDS(ans, tmp, compress = "gzip")
    if (!file.rename(tmp, f)) stop("Could not finalize chunk: ", f)
  }
  invisible(TRUE)
}

read_chunks <- function(outdir, key) {
  d <- file.path(outdir, "raw", key)
  fs <- sort(list.files(d, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = TRUE))
  if (!length(fs)) stop("No chunks for ", key)
  z <- do.call(rbind, lapply(fs, readRDS))
  if ("rule" %in% names(z)) {
    if (anyDuplicated(paste(z$rep_id, z$rule, sep = "::"))) {
      stop("Duplicate replicate/rule rows for ", key)
    }
  } else if (anyDuplicated(z$rep_id)) {
    stop("Duplicate replicate ids for ", key)
  }
  metafile <- file.path(d, "metadata.rds")
  if (file.exists(metafile)) {
    target_R <- readRDS(metafile)$requested_R
    if (is.finite(target_R)) z <- z[z$rep_id <= target_R, , drop = FALSE]
  }
  z[order(z$rep_id), , drop = FALSE]
}

run_transport_grid <- function(control, outdir, overwrite = FALSE) {
  safe_dir(outdir)
  d <- 3L
  reffs <- c(3, 2, 1.4, 1)
  dns <- c(0, 1, 2, 3, 4, 6)
  direction <- env_chr("SMART_ADD_DIRECTION", "maxeig")
  a_values <- c(1, 2)
  # Independent Gaussian-limit calculation of the terminal and transported
  # selection-minus-rule gaps. This is the direct numerical target of the
  # first-order theorem and is shared across finite sample sizes.
  theory_rows <- list(); kt <- 1L
  for (rr in reffs) for (dn in dns) {
    tp <- transport_profile(d, rr)
    delta <- make_delta(dn, d, direction)
    th <- transport_limit_theory(delta, tp$M, control,
                                 seed = control$seed + 700000L + kt * 1000L,
                                 a_values = a_values)
    th$target_reff <- rr; th$delta_norm <- dn; th$direction <- direction
    theory_rows[[kt]] <- th; kt <- kt + 1L
  }
  theory_df <- do.call(rbind, theory_rows)
  utils::write.csv(theory_df, file.path(outdir, "transport_limit_theory.csv"), row.names = FALSE)

  manifest <- list(); q <- 1L
  for (n in control$n_values) for (rr in reffs) for (dn in dns) {
    key <- scenario_slug("transport", paste0("n", n), paste0("r", rr), paste0("dn", dn), direction)
    manifest[[q]] <- data.frame(key = key, n = n, d = d, target_reff = rr,
                                delta_norm = dn, direction = direction, stringsAsFactors = FALSE)
    q <- q + 1L
    run_chunked(key, control$R_smart, control$chunk_size, control$cores, outdir,
      fun = function(i) simulate_transport_replicate(
        n, d, dn, rr, direction,
        seed = rep_seed(control$seed, key, i), a_values = a_values
      ), overwrite = overwrite)
  }
  man <- do.call(rbind, manifest)
  utils::write.csv(man, file.path(outdir, "transport_manifest.csv"), row.names = FALSE)
  summarize_transport_grid(outdir, man)
}

summarize_transport_grid <- function(outdir, manifest = NULL) {
  if (is.null(manifest)) manifest <- utils::read.csv(file.path(outdir, "transport_manifest.csv"), stringsAsFactors = FALSE)
  rows <- list(); k <- 1L
  for (ii in seq_len(nrow(manifest))) {
    mm <- manifest[ii, ]; z <- read_chunks(outdir, mm$key)
    for (rule in unique(z$rule)) {
      zz <- z[z$rule == rule, , drop = FALSE]
      r2 <- mc_mean_se(zz$risk_stage2); r1 <- mc_mean_se(zz$risk_stage1)
      rows[[k]] <- data.frame(
        key = mm$key, n = mm$n, d = mm$d, target_reff = mm$target_reff,
        delta_norm = mm$delta_norm, direction = mm$direction, rule = rule,
        risk_stage2 = r2["mean"], mcse_stage2 = r2["se"], n_risk_stage2 = mm$n * r2["mean"],
        risk_stage1 = r1["mean"], mcse_stage1 = r1["se"], n_risk_stage1 = mm$n * r1["mean"],
        mean_weight = mean(zz$wide_weight), mean_lambda2 = mean(zz$lambda2),
        mean_switch_rate = mean(zz$fitted_switch_rate), reff_exact = zz$reff_exact[1],
        m_max = zz$m_max[1], m_min = zz$m_min[1], m_eigenvalues = zz$m_eigenvalues[1],
        rho = zz$rho[1],
        stringsAsFactors = FALSE
      ); k <- k + 1L
    }
  }
  out <- do.call(rbind, rows)
  utils::write.csv(out, file.path(outdir, "transport_risk_summary.csv"), row.names = FALSE)
  # Paired selection-minus-rule gaps make the theorem visible directly and use
  # the common-random-number pairing within each fitted SMART replicate.
  pairs <- list(); kp <- 1L
  for (ii in seq_len(nrow(manifest))) {
    mm <- manifest[ii, ]; z <- read_chunks(outdir, mm$key)
    zs <- z[z$rule == "SEL", c("rep_id", "risk_stage2", "risk_stage1", "reff_exact")]
    for (arule in setdiff(unique(z$rule), "SEL")) {
      za <- z[z$rule == arule, c("rep_id", "risk_stage2", "risk_stage1")]
      names(zs)[2:3] <- c("r2_sel", "r1_sel")
      names(za)[2:3] <- c("r2_alt", "r1_alt")
      pp <- merge(zs, za, by = "rep_id", sort = TRUE)
      d2 <- pp$r2_sel - pp$r2_alt; d1 <- pp$r1_sel - pp$r1_alt
      q2 <- mc_mean_se(d2); q1 <- mc_mean_se(d1)
      aa <- suppressWarnings(as.numeric(sub("^ST_a", "", arule)))
      reff <- pp$reff_exact[1]
      pairs[[kp]] <- data.frame(
        key = mm$key, n = mm$n, d = mm$d, target_reff = mm$target_reff,
        delta_norm = mm$delta_norm, direction = mm$direction, alt_rule = arule, a = aa,
        gap_stage2_SEL_minus_ALT = q2["mean"], mcse_gap_stage2 = q2["se"], n_gap_stage2 = mm$n * q2["mean"],
        gap_stage1_SEL_minus_ALT = q1["mean"], mcse_gap_stage1 = q1["se"], n_gap_stage1 = mm$n * q1["mean"],
        reff_exact = reff,
        upstream_uniform_guarantee = grepl("^ST_a", arule) && is.finite(aa) && reff > 2 && aa <= 2 * reff - 4,
        stringsAsFactors = FALSE
      )
      kp <- kp + 1L
      names(zs)[2:3] <- c("risk_stage2", "risk_stage1")
    }
  }
  wide <- do.call(rbind, pairs)
  thfile <- file.path(outdir, "transport_limit_theory.csv")
  if (file.exists(thfile)) {
    th <- utils::read.csv(thfile, stringsAsFactors = FALSE)
    names(th)[names(th) == "rule"] <- "alt_rule"
    keep <- c("target_reff", "delta_norm", "direction", "alt_rule",
              "theory_stage2_SEL_minus_rule", "theory_stage2_gap_mcse",
              "theory_stage1_SEL_minus_rule", "theory_stage1_gap_mcse")
    wide <- merge(wide, th[, keep],
                  by = c("target_reff", "delta_norm", "direction", "alt_rule"), all.x = TRUE)
  }
  utils::write.csv(wide, file.path(outdir, "transport_pairwise_gaps.csv"), row.names = FALSE)
  out
}

run_feedback_grid <- function(control, outdir, overwrite = FALSE) {
  safe_dir(outdir)
  manifest <- list(); k <- 1L
  for (n in control$n_values) for (rr in c(3, 2, 1.4)) for (dn in c(1, 3, 6)) {
    key <- scenario_slug("feedback", paste0("n", n), paste0("r", rr), paste0("dn", dn))
    manifest[[k]] <- data.frame(key = key, n = n, target_reff = rr, delta_norm = dn, stringsAsFactors = FALSE); k <- k + 1L
    run_chunked(key, control$R_feedback, control$chunk_size, control$cores, outdir,
      fun = function(i) simulate_feedback_replicate(n, dn, rr, rep_seed(control$seed + 77L, key, i), a = 1),
      overwrite = overwrite)
  }
  man <- do.call(rbind, manifest)
  utils::write.csv(man, file.path(outdir, "feedback_manifest.csv"), row.names = FALSE)
  rows <- list(); k <- 1L
  for (ii in seq_len(nrow(man))) {
    z <- read_chunks(outdir, man$key[ii])
    for (rule in unique(z$rule)) {
      q <- z[z$rule == rule, , drop = FALSE]
      rr <- mc_mean_se(q$risk_stage1)
      rows[[k]] <- data.frame(key = man$key[ii], n = man$n[ii], target_reff = man$target_reff[ii],
        delta_norm = man$delta_norm[ii], rule = rule, risk_stage1 = rr["mean"], mcse = rr["se"],
        n_risk_stage1 = man$n[ii] * rr["mean"], mean_lambda1 = mean(q$lambda1),
        stage1_wide_selection_probability = mean(q$stage1_selected_wide), mean_stage2_weight = mean(q$stage2_weight),
        stringsAsFactors = FALSE); k <- k + 1L
    }
  }
  out <- do.call(rbind, rows)
  utils::write.csv(out, file.path(outdir, "feedback_summary.csv"), row.names = FALSE)
  pairs <- list(); kp <- 1L
  for (ii in seq_len(nrow(man))) {
    z <- read_chunks(outdir, man$key[ii])
    zs <- z[z$rule == "SEL", c("rep_id", "risk_stage1", "lambda1", "stage1_selected_wide")]
    za <- z[z$rule == "ST", c("rep_id", "risk_stage1", "lambda1", "stage1_selected_wide")]
    names(zs)[2:4] <- c("risk_sel", "lambda_sel", "wide_sel")
    names(za)[2:4] <- c("risk_st", "lambda_st", "wide_st")
    pp <- merge(zs, za, by = "rep_id", sort = TRUE)
    dg <- mc_mean_se(pp$risk_sel - pp$risk_st)
    pairs[[kp]] <- data.frame(
      key = man$key[ii], n = man$n[ii], target_reff = man$target_reff[ii], delta_norm = man$delta_norm[ii],
      n_gap_stage1_SEL_minus_ST = man$n[ii] * dg["mean"], n_mcse_gap = man$n[ii] * dg["se"],
      mean_lambda1_SEL_minus_ST = mean(pp$lambda_sel - pp$lambda_st),
      upstream_selection_disagreement_probability = mean(pp$wide_sel != pp$wide_st),
      stringsAsFactors = FALSE
    ); kp <- kp + 1L
  }
  utils::write.csv(do.call(rbind, pairs), file.path(outdir, "feedback_pairwise_SEL_vs_ST.csv"), row.names = FALSE)
  out
}

run_finite_library_grid <- function(control, outdir, overwrite = FALSE) {
  safe_dir(outdir)
  d <- 3L
  delta_list <- list(
    zero = c(0, 0, 0),
    equal_1 = rep(1 / sqrt(3), 3),
    equal_2 = rep(2 / sqrt(3), 3),
    heterogeneous = c(1.5, 1.0, 0.5),
    maxeig_3 = c(3, 0, 0)
  )
  manifest <- list(); k <- 1L; theory <- list(); kt <- 1L
  for (nm in names(delta_list)) {
    dl <- delta_list[[nm]]
    th <- finite_library_limit(dl, control, seed = control$seed + 9000L + kt)
    theory[[kt]] <- data.frame(delta_case = nm, delta = paste(format(dl, digits = 8), collapse = ";"),
      theory_block_SEL = unname(th["SEL"]), theory_block_SEL_mcse = unname(th["SEL_se"]),
      theory_block_MA = unname(th["MA"]), theory_block_MA_mcse = unname(th["MA_se"]),
      theory_common_nuisance = d + 3,
      theory_total_SEL = d + 3 + unname(th["SEL"]),
      theory_total_MA = d + 3 + unname(th["MA"]),
      stringsAsFactors = FALSE)
    kt <- kt + 1L
    for (n in control$n_values) {
      key <- scenario_slug("finite", nm, paste0("n", n))
      manifest[[k]] <- data.frame(key = key, n = n, delta_case = nm,
                                  delta = paste(format(dl, digits = 8), collapse = ";"), stringsAsFactors = FALSE); k <- k + 1L
      run_chunked(key, control$R_finite, control$chunk_size, control$cores, outdir,
        fun = function(i) simulate_finite_library_replicate(n, dl, rep_seed(control$seed + 123L, key, i)),
        overwrite = overwrite)
    }
  }
  man <- do.call(rbind, manifest); thdf <- do.call(rbind, theory)
  utils::write.csv(man, file.path(outdir, "finite_library_manifest.csv"), row.names = FALSE)
  utils::write.csv(thdf, file.path(outdir, "finite_library_limit_theory.csv"), row.names = FALSE)
  rows <- list(); k <- 1L
  for (ii in seq_len(nrow(man))) {
    z <- read_chunks(outdir, man$key[ii])
    rs <- mc_mean_se(z$risk_SEL); rm <- mc_mean_se(z$risk_MA)
    rg <- mc_mean_se(z$risk_SEL - z$risk_MA)
    selprob <- prop.table(table(factor(z$selected, levels = paste0("M", 0:3))))
    rows[[k]] <- data.frame(
      key = man$key[ii], n = man$n[ii], delta_case = man$delta_case[ii], delta = man$delta[ii],
      n_risk_SEL = man$n[ii] * rs["mean"], n_mcse_SEL = man$n[ii] * rs["se"],
      n_risk_MA = man$n[ii] * rm["mean"], n_mcse_MA = man$n[ii] * rm["se"],
      selection_minus_averaging_n_gap = man$n[ii] * rg["mean"],
      n_mcse_gap = man$n[ii] * rg["se"],
      p_M0 = selprob[1], p_M1 = selprob[2], p_M2 = selprob[3], p_M3 = selprob[4],
      mean_w_M0 = mean(z$w_M0), mean_w_M1 = mean(z$w_M1), mean_w_M2 = mean(z$w_M2), mean_w_M3 = mean(z$w_M3),
      stringsAsFactors = FALSE); k <- k + 1L
  }
  out <- merge(do.call(rbind, rows), thdf, by = c("delta_case", "delta"), all.x = TRUE)
  utils::write.csv(out, file.path(outdir, "finite_library_summary.csv"), row.names = FALSE)
  out
}

run_smoke_checks <- function() {
  # 0. Naming contract used by the recursive feedback experiment.
  u1 <- make_stage1_names(3L)
  if (!identical(u1, c("(Intercept)", "X11", "X12", "X13", "A1"))) {
    stop("Stage-1 union-name contract failed")
  }
  tb <- setNames(seq_along(u1), u1)
  if (!identical(names(embed_beta(tb, u1)), u1)) stop("Coefficient embedding contract failed")

  # 1. Effective-rank construction is exact by design.
  targets <- c(3, 2, 1.4, 1)
  got <- vapply(targets, function(r) transport_profile(3, r)$reff, numeric(1))
  if (max(abs(got - targets)) > 1e-12) stop("Effective-rank construction failed")
  # 2. AIC convexity condition for the two standard d=3 Stein choices.
  F <- c(6.1, 10, 100)
  for (a in c(1, 2)) {
    w <- stein_weight(F, 3, 6, a)
    if (any(w < 0 | w > 1)) stop("Stein weights left the simplex")
  }
  # 3. The exact Stein risk difference must be strictly negative.
  for (dn in c(0, 1, 3, 8)) {
    if (!(stein_exact_risk_difference(3, dn, 1, 6) < 0)) stop("Exact Stein risk identity check failed")
  }
  # 4. Separated-gap validation for all main d=3 transport scenarios.
  for (rr in targets) for (dn in c(0, 1, 2, 3, 4, 6)) {
    tp <- transport_profile(3, rr)
    g <- true_gap_validation(250, 3, make_delta(dn, 3, "maxeig"), tp$rho)
    if (!(g["minimum_half_contrast"] > 0)) stop("Gap validation failed")
  }
  invisible(TRUE)
}
