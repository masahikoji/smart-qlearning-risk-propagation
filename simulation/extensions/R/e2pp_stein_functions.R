# Focused E2'' pretest--Stein simulation.

.e2pp_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit)) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

e2pp_log <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

e2pp_safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

e2pp_env_int <- function(name, default) {
  z <- Sys.getenv(name, unset = "")
  if (!nzchar(z)) return(as.integer(default))
  out <- suppressWarnings(as.integer(z))
  if (length(out) != 1L || is.na(out) || out < 1L) stop(name, " must be a positive integer")
  out
}

e2pp_control <- function(profile = c("paper", "smoke", "debug")) {
  profile <- match.arg(profile)
  nphys <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (length(nphys) != 1L || is.na(nphys) || nphys < 1L) nphys <- 1L

  if (profile == "paper") {
    out <- list(
      n = 500L,
      R = 20000L,
      R_theory = 2000000L,
      chunk_size = 250L,
      theory_chunk = 100000L,
      cores = max(1L, min(8L, nphys)),
      seed = 260909L,
      s_grid = c(0, 1, 2, 3, 4, 6),
      rho = 0.6,
      tau0 = 2,
      beta = c(0.5, -0.3, 0.2),
      a_values = c(1, 2)
    )
  } else if (profile == "smoke") {
    out <- list(
      n = 120L,
      R = 80L,
      R_theory = 30000L,
      chunk_size = 20L,
      theory_chunk = 10000L,
      cores = 1L,
      seed = 260909L,
      s_grid = c(0, 4),
      rho = 0.6,
      tau0 = 2,
      beta = c(0.5, -0.3, 0.2),
      a_values = c(1, 2)
    )
  } else {
    out <- list(
      n = 80L,
      R = 8L,
      R_theory = 5000L,
      chunk_size = 4L,
      theory_chunk = 5000L,
      cores = 1L,
      seed = 260909L,
      s_grid = c(0, 4),
      rho = 0.6,
      tau0 = 2,
      beta = c(0.5, -0.3, 0.2),
      a_values = c(1, 2)
    )
  }

  out$R <- e2pp_env_int("E2PP_R", out$R)
  out$R_theory <- e2pp_env_int("E2PP_R_THEORY", out$R_theory)
  out$chunk_size <- e2pp_env_int("E2PP_CHUNK", out$chunk_size)
  out$theory_chunk <- e2pp_env_int("E2PP_THEORY_CHUNK", out$theory_chunk)
  out$cores <- e2pp_env_int("E2PP_CORES", out$cores)
  out
}

e2pp_stable_logistic <- function(x) {
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  ex <- exp(x[!pos])
  out[!pos] <- ex / (1 + ex)
  out
}

e2pp_selection_weight <- function(F, d = 3L, c = 2 * d) as.numeric(F > c)

e2pp_akaike_weight <- function(F, d = 3L) e2pp_stable_logistic(F / 2 - d)

e2pp_stein_weight <- function(F, a, d = 3L, c = 2 * d) {
  if (d < 3L) stop("Pretest--Stein requires d >= 3")
  if (!(a > 0 && a <= min(c, 2 * (d - 2)))) {
    stop("a must satisfy 0 < a <= min(c, 2(d-2))")
  }
  out <- numeric(length(F))
  idx <- is.finite(F) & F > c
  out[idx] <- 1 - a / F[idx]
  out
}

e2pp_sym_inv_sqrt <- function(S, tol = 1e-12) {
  ee <- eigen(S, symmetric = TRUE)
  if (any(!is.finite(ee$values)) || min(ee$values) <= tol) {
    stop("Sigma_Z is not positive definite")
  }
  ee$vectors %*% diag(1 / sqrt(ee$values), nrow = length(ee$values)) %*% t(ee$vectors)
}

e2pp_transport_spec <- function(type = c("multidirectional", "collinear"), rho = 0.6) {
  type <- match.arg(type)
  d <- 3L
  I3 <- diag(d)
  one <- rep(1, d)
  if (type == "multidirectional") {
    Sigma <- I3
    B <- rho * I3
  } else {
    Sigma <- (1 - rho^2) * I3 + rho^2 * tcrossprod(one)
    # E(Z | X1) = B X1 = rho * X11 * 1.
    B <- rho * outer(one, c(1, 0, 0))
  }
  Sinv2 <- e2pp_sym_inv_sqrt(Sigma)
  A <- Sinv2 %*% B
  M <- A %*% t(A)
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  lmax <- max(ev)
  reff <- if (lmax > 0) sum(ev) / lmax else 0
  list(type = type, d = d, rho = rho, Sigma = Sigma, Sinv2 = Sinv2,
       B = B, transported_map = A, M = M, eigenvalues_M = ev, reff = reff)
}

e2pp_mc_mean_se <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, se = NA_real_, sd = NA_real_, n = 0))
  c(mean = mean(x),
    se = if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_,
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
    n = length(x))
}

e2pp_ols <- function(X, y, tol = 1e-10) {
  fit <- stats::lm.fit(x = X, y = y, tol = tol)
  if (fit$rank < ncol(X)) stop("Rank-deficient OLS design in E2'' simulation")
  beta <- as.numeric(fit$coefficients)
  if (any(!is.finite(beta))) stop("Non-finite OLS coefficient in E2'' simulation")
  names(beta) <- colnames(X)
  rss <- sum(fit$residuals^2)
  if (!is.finite(rss) || rss <= 0) stop("Invalid RSS in E2'' simulation")
  list(beta = beta, rss = rss, residuals = as.numeric(fit$residuals))
}

e2pp_embed_narrow <- function(beta0) {
  out <- setNames(rep(0, 9L), c("(Intercept)", paste0("Z", 1:3), "A1", "A2", paste0("A2Z", 1:3)))
  out[names(beta0)] <- beta0
  out
}

e2pp_generate_latent <- function(n, seed) {
  set.seed(seed)
  list(
    X1 = matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L),
    U = matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L),
    A1 = sample(c(-1, 1), n, replace = TRUE),
    A2 = sample(c(-1, 1), n, replace = TRUE),
    eps = stats::rnorm(n)
  )
}

e2pp_make_dataset <- function(latent, s, control, spec) {
  n <- control$n
  delta <- rep(s / sqrt(3), 3L)
  cn <- as.numeric(spec$Sinv2 %*% delta) / sqrt(n)
  if (spec$type == "multidirectional") {
    Z <- spec$rho * latent$X1 + sqrt(1 - spec$rho^2) * latent$U
  } else {
    # Add rho * X11 to every column explicitly.
    Z <- sweep(sqrt(1 - spec$rho^2) * latent$U, 1L,
               spec$rho * latent$X1[, 1L], "+")
  }
  colnames(Z) <- paste0("Z", 1:3)
  beta <- control$beta
  half_contrast_true <- control$tau0 + as.numeric(Z %*% cn)
  Y2 <- 1 + as.numeric(Z %*% beta) + 0.3 * latent$A1 +
    latent$A2 * half_contrast_true + latent$eps
  list(X1 = latent$X1, Z = Z, A1 = latent$A1, A2 = latent$A2, Y2 = Y2,
       delta = delta, cn = cn, half_contrast_true = half_contrast_true, spec = spec)
}

e2pp_fit_stage2 <- function(dat) {
  n <- nrow(dat$Z)
  AZ <- dat$Z * dat$A2
  colnames(AZ) <- paste0("A2Z", 1:3)
  X0 <- cbind(`(Intercept)` = 1, dat$Z, A1 = dat$A1, A2 = dat$A2)
  X1 <- cbind(X0, AZ)
  f0 <- e2pp_ols(X0, dat$Y2)
  f1 <- e2pp_ols(X1, dat$Y2)
  lambda <- max(0, n * log(f0$rss / f1$rss))
  list(narrow = e2pp_embed_narrow(f0$beta), wide = f1$beta,
       rss0 = f0$rss, rss1 = f1$rss, lambda = lambda)
}

e2pp_rule_weights <- function(lambda, a_values = c(1, 2)) {
  out <- c(SEL = e2pp_selection_weight(lambda, 3L, 6),
           AIC_MA = e2pp_akaike_weight(lambda, 3L))
  for (a in a_values) {
    nm <- paste0("ST_a", format(a, trim = TRUE, scientific = FALSE))
    out[nm] <- e2pp_stein_weight(lambda, a = a, d = 3L, c = 6)
  }
  out
}

e2pp_stage2_true_beta <- function(dat, control) {
  setNames(c(1, control$beta, 0.3, control$tau0, dat$cn),
           c("(Intercept)", paste0("Z", 1:3), "A1", "A2", paste0("A2Z", 1:3)))
}

e2pp_stage2_risk <- function(beta_hat, dat, control) {
  truth <- e2pp_stage2_true_beta(dat, control)
  e <- beta_hat[names(truth)] - truth
  S <- dat$spec$Sigma
  as.numeric(e[1]^2 + crossprod(e[2:4], S %*% e[2:4]) +
               e[5]^2 + e[6]^2 + crossprod(e[7:9], S %*% e[7:9]))
}

e2pp_stage2_value <- function(beta_hat, dat) {
  base <- beta_hat["(Intercept)"] +
    as.numeric(dat$Z %*% beta_hat[paste0("Z", 1:3)]) +
    beta_hat["A1"] * dat$A1
  half_contrast <- beta_hat["A2"] +
    as.numeric(dat$Z %*% beta_hat[paste0("A2Z", 1:3)])
  list(value = base + abs(half_contrast), half_contrast = half_contrast)
}

e2pp_stage1_true_beta <- function(dat, control) {
  # Ignoring the exponentially rare true stage-2 sign reversal, which is
  # o(n^{-m}) for every fixed m in the present Gaussian local sequence.
  bx <- as.numeric(t(dat$spec$B) %*% (control$beta + dat$cn))
  setNames(c(1 + control$tau0, bx, 0.3),
           c("(Intercept)", paste0("X1", 1:3), "A1"))
}

e2pp_stage1_risk <- function(beta_hat, dat, control) {
  truth <- e2pp_stage1_true_beta(dat, control)
  e <- beta_hat[names(truth)] - truth
  # E[(1,X1,A1)(1,X1,A1)^T] = I_5 under the evaluation law.
  sum(e^2)
}

e2pp_fit_stage1_many <- function(dat, pseudo_matrix) {
  X <- cbind(`(Intercept)` = 1, dat$X1, A1 = dat$A1)
  colnames(X)[2:4] <- paste0("X1", 1:3)
  q <- qr(X, tol = 1e-10)
  if (q$rank < ncol(X)) stop("Rank-deficient stage-1 design in E2'' simulation")
  bh <- qr.coef(q, pseudo_matrix)
  if (is.null(dim(bh))) bh <- matrix(bh, ncol = 1L)
  rownames(bh) <- colnames(X)
  bh
}

e2pp_simulate_transport <- function(latent, s, control, spec) {
  dat <- e2pp_make_dataset(latent, s, control, spec)
  fit2 <- e2pp_fit_stage2(dat)
  weights <- e2pp_rule_weights(fit2$lambda, control$a_values)
  nr <- length(weights)
  rule_names <- names(weights)

  beta2 <- matrix(NA_real_, nrow = 9L, ncol = nr,
                  dimnames = list(names(fit2$wide), rule_names))
  pseudo <- matrix(NA_real_, nrow = control$n, ncol = nr,
                   dimnames = list(NULL, rule_names))
  risk2 <- switch_fit <- numeric(nr)

  for (j in seq_len(nr)) {
    b <- fit2$narrow + weights[j] * (fit2$wide - fit2$narrow)
    beta2[, j] <- b
    vv <- e2pp_stage2_value(b, dat)
    pseudo[, j] <- vv$value
    risk2[j] <- e2pp_stage2_risk(b, dat, control)
    switch_fit[j] <- mean(vv$half_contrast < 0)
  }

  b1mat <- e2pp_fit_stage1_many(dat, pseudo)
  risk1 <- vapply(seq_len(nr), function(j) {
    bh <- b1mat[, j]
    names(bh) <- rownames(b1mat)
    e2pp_stage1_risk(bh, dat, control)
  }, numeric(1))

  true_switch <- mean(dat$half_contrast_true < 0)
  data.frame(
    transport = spec$type,
    s = s,
    rule = rule_names,
    risk_stage2 = risk2,
    risk_stage1 = risk1,
    lambda2 = fit2$lambda,
    wide_weight = as.numeric(weights),
    fitted_switch_rate = switch_fit,
    true_switch_rate = true_switch,
    reff_exact = spec$reff,
    stringsAsFactors = FALSE
  )
}

e2pp_simulate_replicate <- function(rep_id, s, control, seed) {
  latent <- e2pp_generate_latent(control$n, seed)
  specs <- list(
    e2pp_transport_spec("multidirectional", control$rho),
    e2pp_transport_spec("collinear", control$rho)
  )
  z <- do.call(rbind, lapply(specs, function(sp) e2pp_simulate_transport(latent, s, control, sp)))
  z$rep_id <- rep_id
  z$n <- control$n
  z
}

# Gaussian-limit benchmarks

e2pp_inv_tail_moment_chisq <- function(c, df, ncp) {
  stats::integrate(function(x) stats::dchisq(x, df = df, ncp = ncp) / x,
                   lower = c, upper = Inf, rel.tol = 1e-10,
                   subdivisions = 2000L, stop.on.error = TRUE)$value
}

e2pp_exact_terminal_ST_minus_SEL <- function(s, a, d = 3L, c = 6) {
  e_inv <- e2pp_inv_tail_moment_chisq(c, d, s^2)
  a * (a - 2 * (d - 2)) * e_inv -
    4 * a * stats::dchisq(c, df = d, ncp = s^2)
}

e2pp_limit_chunk <- function(s, control, nmc, seed) {
  set.seed(seed)
  d <- 3L
  delta <- rep(s / sqrt(3), d)
  xi <- matrix(stats::rnorm(nmc * d), nrow = nmc, ncol = d)
  xi <- sweep(xi, 2L, delta, "+")
  F <- rowSums(xi^2)
  dm <- matrix(delta, nrow = nmc, ncol = d, byrow = TRUE)
  weights <- cbind(
    SEL = e2pp_selection_weight(F, d, 6),
    AIC_MA = e2pp_akaike_weight(F, d),
    ST_a1 = e2pp_stein_weight(F, 1, d, 6),
    ST_a2 = e2pp_stein_weight(F, 2, d, 6)
  )
  specs <- list(
    e2pp_transport_spec("multidirectional", control$rho),
    e2pp_transport_spec("collinear", control$rho)
  )

  out <- list(); k <- 1L
  for (j in seq_len(ncol(weights))) {
    w <- weights[, j]
    er <- xi * w - dm
    l2 <- rowSums(er^2)
    for (sp in specs) {
      l1 <- rowSums((er %*% sp$M) * er)
      out[[k]] <- data.frame(
        transport = sp$type,
        rule = colnames(weights)[j],
        nmc = nmc,
        sum_l2 = sum(l2), sum2_l2 = sum(l2^2),
        sum_l1 = sum(l1), sum2_l1 = sum(l1^2),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

e2pp_moment_from_sums <- function(sumx, sumx2, n) {
  N <- sum(n); S <- sum(sumx); S2 <- sum(sumx2)
  mu <- S / N
  vr <- max(0, (S2 - N * mu^2) / max(1, N - 1))
  c(mean = mu, se = sqrt(vr / N))
}

e2pp_run_limit_theory <- function(control, outdir) {
  e2pp_safe_dir(outdir)
  rows <- list(); krow <- 1L
  for (s in control$s_grid) {
    left <- control$R_theory
    parts <- list(); kp <- 1L
    while (left > 0L) {
      m <- min(control$theory_chunk, left)
      parts[[kp]] <- e2pp_limit_chunk(s, control, m,
                                      control$seed + 800000L + 10000L * round(s) + kp)
      left <- left - m
      kp <- kp + 1L
    }
    z <- do.call(rbind, parts)
    for (tr in unique(z$transport)) {
      zz <- z[z$transport == tr, , drop = FALSE]
      sel <- zz[zz$rule == "SEL", , drop = FALSE]
      sel2 <- e2pp_moment_from_sums(sel$sum_l2, sel$sum2_l2, sel$nmc)
      sel1 <- e2pp_moment_from_sums(sel$sum_l1, sel$sum2_l1, sel$nmc)
      for (rule in unique(zz$rule)) {
        q <- zz[zz$rule == rule, , drop = FALSE]
        q2 <- e2pp_moment_from_sums(q$sum_l2, q$sum2_l2, q$nmc)
        q1 <- e2pp_moment_from_sums(q$sum_l1, q$sum2_l1, q$nmc)

        # For paired gap MCSE we need chunk-level joint simulation. Re-run only the
        # algebra from the stored aggregate sums is insufficient, so calculate a
        # separate paired benchmark below in vectorised chunks.
        rows[[krow]] <- data.frame(
          s = s, transport = tr, rule = rule,
          limit_risk_stage2 = q2["mean"], mcse_limit_risk_stage2 = q2["se"],
          limit_risk_stage1 = q1["mean"], mcse_limit_risk_stage1 = q1["se"],
          stringsAsFactors = FALSE
        )
        krow <- krow + 1L
      }
    }
  }
  out <- do.call(rbind, rows)

  # Paired rule-minus-SEL limit gaps, with correct paired MCSEs.
  gap_rows <- list(); kg <- 1L
  for (s in control$s_grid) {
    d <- 3L; delta <- rep(s / sqrt(3), d)
    specs <- list(
      e2pp_transport_spec("multidirectional", control$rho),
      e2pp_transport_spec("collinear", control$rho)
    )
    left <- control$R_theory; chunk_id <- 1L
    accum <- list()
    for (sp in specs) for (rule in c("AIC_MA", "ST_a1", "ST_a2")) {
      accum[[paste(sp$type, rule, sep = "::")]] <- c(sum2 = 0, sum1 = 0, ss2 = 0, ss1 = 0, n = 0)
    }
    while (left > 0L) {
      m <- min(control$theory_chunk, left)
      set.seed(control$seed + 900000L + 10000L * round(s) + chunk_id)
      xi <- matrix(stats::rnorm(m * d), nrow = m, ncol = d)
      xi <- sweep(xi, 2L, delta, "+")
      F <- rowSums(xi^2)
      dm <- matrix(delta, nrow = m, ncol = d, byrow = TRUE)
      ww <- list(
        SEL = e2pp_selection_weight(F, d, 6),
        AIC_MA = e2pp_akaike_weight(F, d),
        ST_a1 = e2pp_stein_weight(F, 1, d, 6),
        ST_a2 = e2pp_stein_weight(F, 2, d, 6)
      )
      err <- lapply(ww, function(w) xi * w - dm)
      loss2 <- lapply(err, function(er) rowSums(er^2))
      for (sp in specs) {
        loss1 <- lapply(err, function(er) rowSums((er %*% sp$M) * er))
        for (rule in c("AIC_MA", "ST_a1", "ST_a2")) {
          g2 <- loss2[[rule]] - loss2$SEL
          g1 <- loss1[[rule]] - loss1$SEL
          nm <- paste(sp$type, rule, sep = "::")
          ac <- accum[[nm]]
          ac["sum2"] <- ac["sum2"] + sum(g2)
          ac["sum1"] <- ac["sum1"] + sum(g1)
          ac["ss2"] <- ac["ss2"] + sum(g2^2)
          ac["ss1"] <- ac["ss1"] + sum(g1^2)
          ac["n"] <- ac["n"] + m
          accum[[nm]] <- ac
        }
      }
      left <- left - m; chunk_id <- chunk_id + 1L
    }
    for (sp in specs) for (rule in c("AIC_MA", "ST_a1", "ST_a2")) {
      nm <- paste(sp$type, rule, sep = "::")
      ac <- accum[[nm]]; N <- ac["n"]
      mu2 <- ac["sum2"] / N; mu1 <- ac["sum1"] / N
      v2 <- max(0, (ac["ss2"] - N * mu2^2) / max(1, N - 1))
      v1 <- max(0, (ac["ss1"] - N * mu1^2) / max(1, N - 1))
      aval <- if (grepl("^ST_a", rule)) as.numeric(sub("^ST_a", "", rule)) else NA_real_
      exact_term <- if (is.finite(aval)) e2pp_exact_terminal_ST_minus_SEL(s, aval) else NA_real_
      exact_up <- if (sp$type == "multidirectional" && is.finite(exact_term)) control$rho^2 * exact_term else NA_real_
      gap_rows[[kg]] <- data.frame(
        s = s, transport = sp$type, rule = rule,
        limit_gap_rule_minus_SEL_stage2 = mu2,
        mcse_limit_gap_stage2 = sqrt(v2 / N),
        limit_gap_rule_minus_SEL_stage1 = mu1,
        mcse_limit_gap_stage1 = sqrt(v1 / N),
        exact_terminal_ST_minus_SEL = exact_term,
        exact_multidirectional_stage1_ST_minus_SEL = exact_up,
        stringsAsFactors = FALSE
      )
      kg <- kg + 1L
    }
  }
  gaps <- do.call(rbind, gap_rows)
  utils::write.csv(out, file.path(outdir, "E2pp_limit_risks.csv"), row.names = FALSE)
  utils::write.csv(gaps, file.path(outdir, "E2pp_limit_gaps.csv"), row.names = FALSE)
  invisible(list(risks = out, gaps = gaps))
}

# Restartable chunk engine

e2pp_rep_seed <- function(base, key, i) {
  raw <- utf8ToInt(key)
  h <- if (length(raw)) sum((seq_along(raw) * raw) %% 1000003L) else 1L
  as.integer((base + h * 1009 + i * 104729) %% 2147483000L + 1L)
}

e2pp_engine_md5 <- function(root) {
  fs <- c(file.path(root, "R", "e2pp_stein_functions.R"),
          file.path(root, "run_E2pp_stein_simulation.R"),
          file.path(root, "E2PP_VERSION"))
  fs <- fs[file.exists(fs)]
  unname(as.character(tools::md5sum(fs)))
}

e2pp_run_chunked_signal <- function(s, control, outdir, root, overwrite = FALSE) {
  key <- paste0("s", format(s, trim = TRUE, scientific = FALSE))
  rawdir <- e2pp_safe_dir(file.path(outdir, "raw", key))
  metafile <- file.path(rawdir, "metadata.rds")
  md5_now <- e2pp_engine_md5(root)
  chunks <- sort(list.files(rawdir, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = TRUE))

  if (overwrite) {
    if (length(chunks)) unlink(chunks, force = TRUE)
    if (file.exists(metafile)) unlink(metafile, force = TRUE)
    chunks <- character(0)
  } else if (length(chunks)) {
    if (!file.exists(metafile)) stop("Existing E2'' chunks have no metadata for ", key)
    old <- readRDS(metafile)
    if (!identical(as.integer(old$chunk_size), as.integer(control$chunk_size))) {
      stop("Chunk size changed for ", key, "; rerun with overwrite=TRUE")
    }
    if (!identical(old$engine_md5, md5_now)) {
      stop("E2'' source code changed since chunks were created for ", key,
           "; rerun with overwrite=TRUE")
    }
  }

  saveRDS(list(key = key, s = s, n = control$n, requested_R = control$R,
               chunk_size = control$chunk_size, engine_md5 = md5_now,
               timestamp = as.character(Sys.time()), R_version = R.version.string), metafile)

  starts <- seq.int(1L, control$R, by = control$chunk_size)
  for (st in starts) {
    en <- min(control$R, st + control$chunk_size - 1L)
    f <- file.path(rawdir, sprintf("chunk_%06d_%06d.rds", st, en))
    if (file.exists(f)) next
    ids <- st:en
    one <- function(i) {
      e2pp_simulate_replicate(i, s, control, e2pp_rep_seed(control$seed, key, i))
    }
    z <- if (control$cores > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(ids, one, mc.cores = control$cores,
                         mc.preschedule = TRUE, mc.set.seed = FALSE)
    } else {
      lapply(ids, one)
    }
    ans <- do.call(rbind, z)
    tmp <- paste0(f, ".tmp")
    saveRDS(ans, tmp, compress = "gzip")
    if (!file.rename(tmp, f)) stop("Could not finalise E2'' chunk: ", f)
    e2pp_log("completed ", key, " reps ", st, "-", en)
  }
  invisible(TRUE)
}

e2pp_read_signal_chunks <- function(s, control, outdir) {
  key <- paste0("s", format(s, trim = TRUE, scientific = FALSE))
  d <- file.path(outdir, "raw", key)
  fs <- sort(list.files(d, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = TRUE))
  if (!length(fs)) stop("No E2'' chunks found for ", key)
  z <- do.call(rbind, lapply(fs, readRDS))
  z <- z[z$rep_id <= control$R, , drop = FALSE]
  id <- paste(z$rep_id, z$transport, z$rule, sep = "::")
  if (anyDuplicated(id)) stop("Duplicate E2'' replicate/transport/rule rows for ", key)
  expected <- control$R * 2L * 4L
  if (nrow(z) != expected) {
    stop("Incomplete E2'' result for ", key, ": expected ", expected, " rows, found ", nrow(z))
  }
  z
}

e2pp_run_simulation_grid <- function(control, outdir, root, overwrite = FALSE) {
  e2pp_safe_dir(outdir)
  manifest <- do.call(rbind, lapply(control$s_grid, function(s) {
    p_rev <- if (s > 0) stats::pnorm(-control$tau0 * sqrt(control$n) / s) else 0
    data.frame(
      n = control$n, R = control$R, s = s, rho = control$rho, tau0 = control$tau0,
      beta1 = control$beta[1], beta2 = control$beta[2], beta3 = control$beta[3],
      d = 3L, AIC_threshold = 6,
      theoretical_true_sign_reversal_probability = p_rev,
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(manifest, file.path(outdir, "E2pp_manifest.csv"), row.names = FALSE)
  for (s in control$s_grid) {
    e2pp_log("running E2'' signal s=", s)
    e2pp_run_chunked_signal(s, control, outdir, root, overwrite = overwrite)
  }
  e2pp_summarise(control, outdir)
}

e2pp_summarise <- function(control, outdir) {
  risk_rows <- list(); gap_rows <- list(); kr <- 1L; kg <- 1L
  for (s in control$s_grid) {
    z <- e2pp_read_signal_chunks(s, control, outdir)
    for (tr in c("multidirectional", "collinear")) {
      zz <- z[z$transport == tr, , drop = FALSE]
      for (rule in c("SEL", "AIC_MA", "ST_a1", "ST_a2")) {
        q <- zz[zz$rule == rule, , drop = FALSE]
        r2 <- e2pp_mc_mean_se(q$risk_stage2)
        r1 <- e2pp_mc_mean_se(q$risk_stage1)
        risk_rows[[kr]] <- data.frame(
          n = control$n, R = control$R, s = s, transport = tr, rule = rule,
          risk_stage2 = r2["mean"], mcse_risk_stage2 = r2["se"], n_risk_stage2 = control$n * r2["mean"],
          risk_stage1 = r1["mean"], mcse_risk_stage1 = r1["se"], n_risk_stage1 = control$n * r1["mean"],
          mean_lambda2 = mean(q$lambda2), wide_selection_probability = mean(q$lambda2 > 6),
          mean_wide_weight = mean(q$wide_weight),
          mean_fitted_switch_rate = mean(q$fitted_switch_rate),
          mean_true_switch_rate = mean(q$true_switch_rate),
          reff_exact = q$reff_exact[1], stringsAsFactors = FALSE
        )
        kr <- kr + 1L
      }

      sel <- zz[zz$rule == "SEL", c("rep_id", "risk_stage2", "risk_stage1")]
      names(sel)[2:3] <- c("risk2_sel", "risk1_sel")
      for (rule in c("AIC_MA", "ST_a1", "ST_a2")) {
        alt <- zz[zz$rule == rule, c("rep_id", "risk_stage2", "risk_stage1")]
        names(alt)[2:3] <- c("risk2_alt", "risk1_alt")
        pp <- merge(sel, alt, by = "rep_id", sort = TRUE)
        d2 <- pp$risk2_alt - pp$risk2_sel
        d1 <- pp$risk1_alt - pp$risk1_sel
        q2 <- e2pp_mc_mean_se(d2); q1 <- e2pp_mc_mean_se(d1)
        gap_rows[[kg]] <- data.frame(
          n = control$n, R = control$R, s = s, transport = tr, rule = rule,
          gap_rule_minus_SEL_stage2 = q2["mean"], mcse_gap_stage2 = q2["se"],
          n_gap_rule_minus_SEL_stage2 = control$n * q2["mean"], n_mcse_gap_stage2 = control$n * q2["se"],
          gap_rule_minus_SEL_stage1 = q1["mean"], mcse_gap_stage1 = q1["se"],
          n_gap_rule_minus_SEL_stage1 = control$n * q1["mean"], n_mcse_gap_stage1 = control$n * q1["se"],
          stringsAsFactors = FALSE
        )
        kg <- kg + 1L
      }
    }
  }
  risks <- do.call(rbind, risk_rows)
  gaps <- do.call(rbind, gap_rows)
  utils::write.csv(risks, file.path(outdir, "E2pp_risk_summary.csv"), row.names = FALSE)
  utils::write.csv(gaps, file.path(outdir, "E2pp_pairwise_gaps.csv"), row.names = FALSE)

  # Publication-facing table matching the current main-text layout, but with
  # one MCSE column next to every estimate so no approximate SE needs to be guessed.
  getgap <- function(s, tr, rule, stage) {
    q <- gaps[gaps$s == s & gaps$transport == tr & gaps$rule == rule, , drop = FALSE]
    if (nrow(q) != 1L) stop("Could not identify unique E2'' main-table gap")
    if (stage == 2L) c(est = q$n_gap_rule_minus_SEL_stage2, se = q$n_mcse_gap_stage2)
    else c(est = q$n_gap_rule_minus_SEL_stage1, se = q$n_mcse_gap_stage1)
  }
  tab <- do.call(rbind, lapply(control$s_grid, function(s) {
    t1 <- getgap(s, "multidirectional", "ST_a1", 2L)
    m1 <- getgap(s, "multidirectional", "ST_a1", 1L)
    ma <- getgap(s, "multidirectional", "AIC_MA", 1L)
    c1 <- getgap(s, "collinear", "ST_a1", 1L)
    c2 <- getgap(s, "collinear", "ST_a2", 1L)
    data.frame(
      s = s,
      terminal_ST_a1 = t1["est"], mcse_terminal_ST_a1 = t1["se"],
      multidirectional_ST_a1 = m1["est"], mcse_multidirectional_ST_a1 = m1["se"],
      multidirectional_Akaike = ma["est"], mcse_multidirectional_Akaike = ma["se"],
      collinear_ST_a1 = c1["est"], mcse_collinear_ST_a1 = c1["se"],
      collinear_ST_a2 = c2["est"], mcse_collinear_ST_a2 = c2["se"],
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(tab, file.path(outdir, "E2pp_main_table.csv"), row.names = FALSE)

  specs <- list(e2pp_transport_spec("multidirectional", control$rho),
                e2pp_transport_spec("collinear", control$rho))
  validation <- do.call(rbind, lapply(specs, function(sp) {
    data.frame(
      transport = sp$type, rho = control$rho,
      reff_exact = sp$reff,
      M_eigen_1 = sp$eigenvalues_M[1], M_eigen_2 = sp$eigenvalues_M[2], M_eigen_3 = sp$eigenvalues_M[3],
      max_whitening_error = max(abs(sp$Sinv2 %*% sp$Sigma %*% sp$Sinv2 - diag(3))),
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(validation, file.path(outdir, "E2pp_transport_validation.csv"), row.names = FALSE)
  invisible(list(risks = risks, gaps = gaps, main_table = tab, validation = validation))
}

# Smoke checks

e2pp_smoke_checks <- function(control) {
  tol <- 1e-9
  md <- e2pp_transport_spec("multidirectional", control$rho)
  mc <- e2pp_transport_spec("collinear", control$rho)
  target_md <- control$rho^2 * diag(3)
  target_mc <- control$rho^2 / (1 + 2 * control$rho^2) * matrix(1, 3, 3)
  if (max(abs(md$M - target_md)) > tol) stop("Multidirectional M check failed")
  if (max(abs(mc$M - target_mc)) > tol) stop("Collinear M check failed")
  if (abs(md$reff - 3) > tol) stop("Multidirectional effective-rank check failed")
  if (abs(mc$reff - 1) > tol) stop("Collinear effective-rank check failed")
  if (max(abs(md$Sinv2 %*% md$Sigma %*% md$Sinv2 - diag(3))) > tol) stop("Whitening check failed")
  if (max(abs(mc$Sinv2 %*% mc$Sigma %*% mc$Sinv2 - diag(3))) > tol) stop("Collinear whitening check failed")

  Ftest <- c(0, 5.9, 6, 6.1, 20)
  w1 <- e2pp_stein_weight(Ftest, 1, 3, 6)
  w2 <- e2pp_stein_weight(Ftest, 2, 3, 6)
  if (any(w1 < 0 | w1 > 1) || any(w2 < 0 | w2 > 1)) stop("Stein convex-weight check failed")
  if (any(w1[Ftest <= 6] != 0) || any(w2[Ftest <= 6] != 0)) stop("Stein threshold check failed")

  ex <- vapply(c(0, 1, 2, 3, 4, 6), function(s) e2pp_exact_terminal_ST_minus_SEL(s, 1), numeric(1))
  if (any(!is.finite(ex)) || any(ex >= 0)) stop("Exact terminal Stein-risk check failed")

  small <- control
  small$n <- max(80L, min(120L, control$n))
  small$a_values <- c(1, 2)
  zz <- e2pp_simulate_replicate(1L, 1, small, 12345L)
  if (nrow(zz) != 8L || any(!is.finite(zz$risk_stage1)) || any(!is.finite(zz$risk_stage2))) {
    stop("Finite-sample replicate smoke check failed")
  }
  invisible(TRUE)
}
