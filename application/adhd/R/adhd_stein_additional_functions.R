model_matrix_new <- function(fit, newdata) {
  tt <- stats::delete.response(stats::terms(fit))
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = fit$xlevels)
  mm <- stats::model.matrix(tt, mf, contrasts.arg = fit$contrasts)
  train_names <- colnames(stats::model.matrix(fit))
  miss <- setdiff(train_names, colnames(mm))
  if (length(miss)) {
    add <- matrix(0, nrow = nrow(mm), ncol = length(miss), dimnames = list(NULL, miss))
    mm <- cbind(mm, add)
  }
  mm[, train_names, drop = FALSE]
}

fit_nested_block <- function(dat, narrow_formula, wide_formula) {
  f0 <- stats::lm(narrow_formula, data = dat, na.action = stats::na.fail, x = TRUE, y = TRUE, model = TRUE)
  f1 <- stats::lm(wide_formula, data = dat, na.action = stats::na.fail, x = TRUE, y = TRUE, model = TRUE)
  if (stats::nobs(f0) != stats::nobs(f1)) stop("Narrow and wide models use different samples")
  if (f0$rank != length(stats::coef(f0)) || f1$rank != length(stats::coef(f1))) stop("Rank-deficient nested block")
  X <- stats::model.matrix(f0)
  W <- stats::model.matrix(f1)
  qx <- qr(X)
  C <- qr.coef(qx, W)
  C[is.na(C)] <- 0
  R <- W - X %*% C
  sv <- svd(R, nu = 0)
  tol <- max(dim(R)) * max(sv$d) * .Machine$double.eps * 100
  keep <- which(sv$d > tol)
  dadd <- length(keep)
  expected <- length(stats::coef(f1)) - length(stats::coef(f0))
  if (dadd != expected) {
    stop("Residualized added-block rank mismatch: SVD rank=", dadd, ", coefficient difference=", expected)
  }
  Vd <- sv$v[, keep, drop = FALSE]
  Dd <- sv$d[keep]
  Bwhite <- Vd %*% diag(sqrt(nrow(R)) / Dd, nrow = dadd)
  Zstd <- R %*% Bwhite
  gram <- crossprod(Zstd) / nrow(Zstd)
  if (max(abs(gram - diag(dadd))) > 1e-7) stop("Whitening check failed")
  rss0 <- sum(stats::residuals(f0)^2)
  rss1 <- sum(stats::residuals(f1)^2)
  lambda <- stats::nobs(f1) * log(rss0 / rss1)
  aic_diff <- stats::AIC(f1) - stats::AIC(f0)
  aic_identity_error <- aic_diff - (2 * dadd - lambda)
  list(
    narrow = f0, wide = f1, X = X, W = W, projection_coef = C,
    Bwhite = Bwhite, Zstd_train = Zstd, d = dadd,
    lambda = max(0, lambda), c = 2 * dadd,
    aic_diff_wide_minus_narrow = aic_diff,
    aic_identity_error = aic_identity_error
  )
}

residualized_standardized_block_new <- function(block_fit, newdata) {
  Xn <- model_matrix_new(block_fit$narrow, newdata)
  Wn <- model_matrix_new(block_fit$wide, newdata)
  Rn <- Wn - Xn %*% block_fit$projection_coef
  Rn %*% block_fit$Bwhite
}

block_rule_weight <- function(lambda, d, rule = c("SEL", "AIC_MA", "ST"), a = NA_real_) {
  rule <- match.arg(rule)
  if (rule == "SEL") return(as.numeric(lambda > 2 * d))
  if (rule == "AIC_MA") return(1 / (1 + exp(d - lambda / 2)))
  if (d < 3L) stop("Pretest--Stein requires d >= 3")
  if (!is.finite(a) || !(a > 0 && a <= min(2 * d, 2 * (d - 2)))) {
    stop("For ST, a must satisfy 0<a<=min(2d,2(d-2))")
  }
  if (lambda <= 2 * d) 0 else 1 - a / lambda
}

predict_nested_aggregate <- function(block_fit, newdata, weight) {
  p0 <- as.numeric(stats::predict(block_fit$narrow, newdata = newdata))
  p1 <- as.numeric(stats::predict(block_fit$wide, newdata = newdata))
  (1 - weight) * p0 + weight * p1
}

stage2_actions_nested <- function(block_fit, data_c, weight) {
  dm <- data_c; dp <- data_c
  dm$a2 <- -1; dp$a2 <- 1
  qm <- predict_nested_aggregate(block_fit, dm, weight)
  qp <- predict_nested_aggregate(block_fit, dp, weight)
  contrast <- qp - qm
  data.frame(q_minus = qm, q_plus = qp, contrast = contrast,
             action = ifelse(contrast >= 0, 1, -1))
}

pseudo_from_nested_rule <- function(block_fit, data_c, weight) {
  nr <- stage2_eligible(data_c)
  out <- data_c$y
  acts <- data.frame(q_minus = NA_real_, q_plus = NA_real_, contrast = NA_real_, action = NA_integer_)
  acts <- acts[rep(1L, nrow(data_c)), , drop = FALSE]
  if (any(nr)) {
    aa <- stage2_actions_nested(block_fit, data_c[nr, , drop = FALSE], weight)
    out[nr] <- pmax(aa$q_minus, aa$q_plus)
    acts[nr, ] <- aa
  }
  list(pseudo = out, actions = acts)
}

fit_fixed_M11 <- function(data_c, pseudo) {
  dd <- data_c; dd$pseudo <- pseudo
  stats::lm(pseudo ~ o11c + o12c + o13 + a1 + o13:a1, data = dd, na.action = stats::na.fail)
}

stage1_actions_fit <- function(fit, data_c) {
  dm <- data_c; dp <- data_c
  dm$a1 <- -1; dp$a1 <- 1
  qm <- as.numeric(stats::predict(fit, newdata = dm))
  qp <- as.numeric(stats::predict(fit, newdata = dp))
  cc <- qp - qm
  data.frame(q_minus = qm, q_plus = qp, contrast = cc, action = ifelse(cc >= 0, 1, -1))
}

estimate_transport_matrix <- function(block_fit, data_c) {
  nr <- stage2_eligible(data_c)
  # Use the fitted wide Q model only to define an estimated designated terminal action.
  a_wide <- stage2_actions_nested(block_fit, data_c[nr, , drop = FALSE], weight = 1)$action
  dnew <- data_c[nr, , drop = FALSE]
  dnew$a2 <- a_wide
  zopt <- residualized_standardized_block_new(block_fit, dnew)

  # For responders the stage-1 pseudo-outcome is observed Y, so the stage-2
  # generated-response direction contributes zero. This embeds the response-triggered
  # SMART mechanism in the empirical transport diagnostic.
  Tobs <- matrix(0, nrow = nrow(data_c), ncol = block_fit$d)
  Tobs[nr, ] <- zopt
  X1 <- stats::model.matrix(~ o11c + o12c + o13 + a1 + o13:a1, data = data_c)
  qx <- qr(X1)
  muhat <- apply(Tobs, 2L, function(y) {
    cc <- qr.coef(qx, y)
    cc[is.na(cc)] <- 0
    as.numeric(X1 %*% cc)
  })
  if (block_fit$d == 1L) muhat <- matrix(muhat, ncol = 1L)
  Mhat <- crossprod(muhat) / nrow(muhat)
  ev <- eigen(Mhat, symmetric = TRUE, only.values = TRUE)$values
  ev[ev < 0 & ev > -1e-10] <- 0
  lam <- max(ev)
  reff <- if (lam > 0) sum(ev) / lam else NA_real_
  list(M = Mhat, eigenvalues = ev, reff = reff, muhat = muhat,
       designated_action = a_wide)
}

make_block_pairs <- function() {
  libs <- adhd_sensitivity_formula_libraries(response_stage2 = "y", response_stage1 = "pseudo")
  list(
    expanded_E22_E29 = list(
      narrow = libs$expanded_stage2_effect_modifiers$stage2$E22,
      wide = libs$expanded_stage2_effect_modifiers$stage2$E29,
      role = "near-threshold broad-block comparison",
      description = "Existing expanded stage-2 sensitivity: E22 versus the prespecified broad E29 block."
    ),
    expanded_E20_E29 = list(
      narrow = libs$expanded_stage2_effect_modifiers$stage2$E20,
      wide = libs$expanded_stage2_effect_modifiers$stage2$E29,
      role = "strong-separation broad-block control",
      description = "Existing expanded stage-2 sensitivity: E20 versus the prespecified broad E29 block."
    ),
    spline_N22_N26 = list(
      narrow = libs$spline_functional_form$stage2$N22,
      wide = libs$spline_functional_form$stage2$N26,
      role = "near-threshold spline-block comparison",
      description = "Existing spline sensitivity: N22 versus the prespecified broad N26 treatment-modifier block."
    ),
    spline_N20_N26 = list(
      narrow = libs$spline_functional_form$stage2$N20,
      wide = libs$spline_functional_form$stage2$N26,
      role = "strong-separation spline-block control",
      description = "Existing spline sensitivity: N20 versus the prespecified broad N26 treatment-modifier block."
    )
  )
}

run_one_adhd_block <- function(data_c, pair_name, pair_def, a_values = c(1)) {
  nr <- stage2_eligible(data_c)
  bf <- fit_nested_block(data_c[nr, , drop = FALSE], pair_def$narrow, pair_def$wide)
  d <- bf$d
  rules <- list(SEL = list(type = "SEL", a = NA_real_), AIC_MA = list(type = "AIC_MA", a = NA_real_))
  if (d >= 3L) {
    for (a in unique(c(a_values, d - 2))) {
      if (a > 0 && a <= min(2 * d, 2 * (d - 2))) {
        rules[[paste0("ST_a", format(a, trim = TRUE, scientific = FALSE))]] <- list(type = "ST", a = a)
      }
    }
  }
  tr <- estimate_transport_matrix(bf, data_c)
  result_rows <- list(); subject_rows <- list(); k <- 1L
  for (nm in names(rules)) {
    ww <- block_rule_weight(bf$lambda, d, rules[[nm]]$type, rules[[nm]]$a)
    ps <- pseudo_from_nested_rule(bf, data_c, ww)
    f1 <- fit_fixed_M11(data_c, ps$pseudo)
    a1 <- stage1_actions_fit(f1, data_c)
    cc <- stats::coef(f1)
    result_rows[[k]] <- data.frame(
      pair = pair_name, role = pair_def$role, description = pair_def$description,
      stage2_n = sum(nr), d = d,
      lambda = bf$lambda, AIC_threshold = 2 * d,
      AIC_diff_wide_minus_narrow = bf$aic_diff_wide_minus_narrow,
      AIC_identity_error = bf$aic_identity_error,
      rule = nm, a = rules[[nm]]$a,
      wide_weight = ww, stage2_wide_selected_by_AIC = bf$lambda > 2 * d,
      stage2_plus_action_rate = mean(ps$actions$action[nr] == 1),
      pseudo_mean = mean(ps$pseudo), pseudo_sd = stats::sd(ps$pseudo),
      stage1_beta_A1 = unname(cc["a1"]), stage1_beta_O13A1 = unname(cc["o13:a1"]),
      stage1_min_abs_contrast = min(abs(a1$contrast)),
      stage1_plus_action_rate = mean(a1$action == 1),
      transported_reff_hat = tr$reff,
      empirical_reff_condition_indicator = if (rules[[nm]]$type == "ST") {
        is.finite(tr$reff) && tr$reff > 2 && rules[[nm]]$a <= 2 * tr$reff - 4
      } else NA,
      stringsAsFactors = FALSE
    )
    subject_rows[[k]] <- data.frame(
      id = data_c$id, pair = pair_name, rule = nm,
      pseudo = ps$pseudo, stage1_contrast = a1$contrast, stage1_action = a1$action,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  ev <- data.frame(pair = pair_name, eigen_index = seq_along(tr$eigenvalues),
                   eigenvalue = tr$eigenvalues, transported_reff_hat = tr$reff,
                   stringsAsFactors = FALSE)
  list(summary = do.call(rbind, result_rows), subject = do.call(rbind, subject_rows),
       eigen = ev, block_fit = bf, transport = tr)
}

pairwise_rule_differences <- function(subject) {
  key <- unique(subject$pair)
  rules <- unique(subject$rule)
  sel <- subject[subject$rule == "SEL", c("id", "pair", "pseudo", "stage1_contrast", "stage1_action")]
  names(sel)[3:5] <- c("pseudo_SEL", "contrast_SEL", "action_SEL")
  rows <- list(); k <- 1L
  for (rr in setdiff(rules, "SEL")) {
    z <- subject[subject$rule == rr, c("id", "pair", "pseudo", "stage1_contrast", "stage1_action")]
    names(z)[3:5] <- c("pseudo_ALT", "contrast_ALT", "action_ALT")
    m <- merge(sel, z, by = c("id", "pair"), sort = FALSE)
    rows[[k]] <- data.frame(
      pair = key, alt_rule = rr,
      pseudo_RMS_SEL_minus_ALT = sqrt(mean((m$pseudo_SEL - m$pseudo_ALT)^2)),
      pseudo_max_abs_SEL_minus_ALT = max(abs(m$pseudo_SEL - m$pseudo_ALT)),
      contrast_RMS_SEL_minus_ALT = sqrt(mean((m$contrast_SEL - m$contrast_ALT)^2)),
      contrast_max_abs_SEL_minus_ALT = max(abs(m$contrast_SEL - m$contrast_ALT)),
      stage1_decision_discordance_n = sum(m$action_SEL != m$action_ALT),
      stage1_decision_discordance_rate = mean(m$action_SEL != m$action_ALT),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  do.call(rbind, rows)
}
