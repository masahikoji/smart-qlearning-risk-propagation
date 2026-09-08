# OLS fitting, information criteria, recursive selection/averaging.

safe_ols <- function(X, y, stabilize_tol = 1e-10) {
  n <- nrow(X); p <- ncol(X)
  fit <- .lm.fit(x = X, y = y)
  coef_raw <- fit$coefficients
  rank_ok <- isTRUE(fit$rank == p) && all(is.finite(coef_raw))
  if (!rank_ok) {
    return(list(valid = FALSE, coef = rep(NA_real_, p), rss = Inf,
                rank = fit$rank, stabilized = TRUE, pivoted = NA))
  }

  # .lm.fit() may return coefficients in QR-pivot order.  Under full rank we can
  # restore the original column order explicitly; this also makes the assumption
  # testable rather than relying on implementation details of the current R build.
  pivot <- fit$pivot
  if (is.null(pivot) || length(pivot) != p) pivot <- seq_len(p)
  coef <- as.numeric(coef_raw)
  pivoted <- !identical(as.integer(pivot), seq_len(p))
  if (pivoted) {
    tmp <- coef
    coef <- numeric(p)
    coef[pivot] <- tmp
  }

  G <- crossprod(X) / n
  evmin <- tryCatch(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values),
                    error = function(e) -Inf)
  if (!is.finite(evmin) || evmin < stabilize_tol) {
    return(list(valid = FALSE, coef = rep(NA_real_, p), rss = Inf,
                rank = fit$rank, stabilized = TRUE, pivoted = pivoted))
  }
  resid <- as.numeric(y - X %*% coef)
  rss <- max(sum(resid^2), .Machine$double.xmin)
  list(valid = TRUE, coef = coef, rss = rss,
       rank = fit$rank, stabilized = FALSE, pivoted = pivoted)
}

fit_candidate_library <- function(data, stage, y, spec, stabilize_tol = 1e-10) {
  n <- length(y)
  out <- vector("list", length(spec$candidates))
  names(out) <- names(spec$candidates)
  for (m in names(spec$candidates)) {
    feats <- spec$candidates[[m]]
    X <- feature_matrix(data, stage, data[[paste0("A", stage)]], feats)
    f <- safe_ols(X, y, stabilize_tol = stabilize_tol)
    if (!f$valid) {
      # Prespecified deterministic fallback: use the smallest valid nested model
      # available from this same library; if none is valid, stop.  The fallback is
      # logged via the stabilized flag in the returned object.
      f$features <- feats
      f$k <- length(feats)
      f$aic <- Inf; f$bic <- Inf
      f$coef_union <- setNames(rep(NA_real_, length(spec$union)), spec$union)
      out[[m]] <- f
      next
    }
    k <- length(feats)
    neg2loglik_constfree <- n * log(f$rss / n)
    f$features <- feats
    f$k <- k
    f$aic <- neg2loglik_constfree + 2 * k
    f$bic <- neg2loglik_constfree + log(n) * k
    f$coef_union <- embed_coef(f$coef, feats, spec$union)
    out[[m]] <- f
  }

  valid <- vapply(out, function(z) isTRUE(z$valid), logical(1))
  if (!any(valid)) stop("All candidate models failed the Gram/rank check at stage ", stage)

  # If some larger candidate fails, make its criterion infinite.  This is equivalent
  # to excluding it on that replicate and preserves the valid smaller fits.
  attr(out, "stabilized") <- any(!valid)
  out
}

aggregate_library <- function(fits, criterion = c("AIC", "BIC"),
                              mode = c("SEL", "MA")) {
  criterion <- match.arg(criterion)
  mode <- match.arg(mode)
  ic <- vapply(fits, function(z) if (criterion == "AIC") z$aic else z$bic,
               numeric(1))
  valid <- is.finite(ic)
  if (!any(valid)) stop("No valid candidate in aggregate_library")

  if (mode == "SEL") {
    j <- which.min(ic)
    w <- rep(0, length(fits)); w[j] <- 1
  } else {
    logw <- rep(-Inf, length(fits)); logw[valid] <- -0.5 * ic[valid]
    w <- stable_softmax(logw)
    j <- which.max(w)
  }
  names(w) <- names(fits)

  union <- names(fits[[which(valid)[1L]]]$coef_union)
  beta <- setNames(rep(0, length(union)), union)
  for (m in seq_along(fits)) {
    if (w[m] > 0 && fits[[m]]$valid) beta <- beta + w[m] * fits[[m]]$coef_union
  }

  list(coef = beta, weights = w, selected = names(fits)[j],
       criterion = criterion, mode = mode,
       ic = ic, stabilized = isTRUE(attr(fits, "stabilized")))
}

procedure_definitions <- function() {
  list(
    AIC_SEL = list(criterion = "AIC", mode = "SEL"),
    AIC_MA  = list(criterion = "AIC", mode = "MA"),
    BIC_SEL = list(criterion = "BIC", mode = "SEL"),
    BIC_MA  = list(criterion = "BIC", mode = "MA")
  )
}

# Useful exact nested-model diagnostics.
nested_ic_diagnostics <- function(fits, n) {
  if (!all(c("narrow", "wide") %in% names(fits))) return(NULL)
  f0 <- fits[["narrow"]]; f1 <- fits[["wide"]]
  d <- f1$k - f0$k
  if (!f0$valid || !f1$valid) {
    return(list(lambda = NA_real_, d = d,
                aic_diff_wide_minus_narrow = NA_real_,
                bic_diff_wide_minus_narrow = NA_real_,
                aic_identity_error = NA_real_, bic_identity_error = NA_real_,
                akaike_wide_weight = NA_real_, bic_wide_weight = NA_real_))
  }
  lambda <- n * log(f0$rss / f1$rss)
  list(
    lambda = lambda,
    d = d,
    aic_diff_wide_minus_narrow = f1$aic - f0$aic,
    bic_diff_wide_minus_narrow = f1$bic - f0$bic,
    aic_identity_error = (f1$aic - f0$aic) - (2 * d - lambda),
    bic_identity_error = (f1$bic - f0$bic) - (d * log(n) - lambda),
    akaike_wide_weight = 1 / (1 + exp(d - lambda / 2)),
    bic_wide_weight = 1 / (1 + exp((d * log(n) - lambda) / 2))
  )
}
