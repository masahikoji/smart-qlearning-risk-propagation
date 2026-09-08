# ADHD SMART: recursive model selection vs model averaging
# Base-R implementation designed to mirror the source-aligned Q-learning models.

stop_if_missing <- function(x, nm) {
  miss <- setdiff(nm, names(x))
  if (length(miss) > 0L) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }
}

stage2_eligible <- function(dat) {
  dat$r == 0 & !is.na(dat$a2)
}

load_adhd_data <- function() {
  if (!requireNamespace("DTRlearn2", quietly = TRUE)) {
    stop(
      "Package 'DTRlearn2' is required only to load the public ADHD data. ",
      "Install it with install.packages('DTRlearn2')."
    )
  }
  e <- new.env(parent = emptyenv())
  data("adhd", package = "DTRlearn2", envir = e)
  if (!exists("adhd", envir = e, inherits = FALSE)) {
    stop("Could not load DTRlearn2::adhd.")
  }
  dat <- get("adhd", envir = e, inherits = FALSE)
  dat <- as.data.frame(dat)
  req <- c("id", "o11", "o12", "o13", "o14", "a1", "r", "o21", "o22", "a2", "y")
  stop_if_missing(dat, req)

  if (nrow(dat) != 150L) {
    warning("Expected 150 rows from DTRlearn2::adhd; found ", nrow(dat), ".")
  }
  if (!all(na.omit(unique(dat$a1)) %in% c(-1, 1))) stop("a1 must be coded -1/+1.")
  if (!all(na.omit(unique(dat$a2)) %in% c(-1, 1))) stop("a2 must be coded -1/+1 when observed.")
  if (!all(na.omit(unique(dat$r)) %in% c(0, 1))) stop("r must be coded 0/1.")
  if (any(dat$r == 0 & is.na(dat$a2))) {
    stop("Source-aligned analysis requires observed a2 for every nonresponder.")
  }

  dat$s <- 1L - dat$r
  dat
}

compute_centers <- function(train) {
  nr <- stage2_eligible(train)
  if (sum(nr, na.rm = TRUE) < 20L) stop("Too few nonresponders in training data.")
  c(
    o11 = mean(train$o11, na.rm = TRUE),
    o12 = mean(train$o12, na.rm = TRUE),
    o13 = mean(train$o13, na.rm = TRUE),
    o14 = mean(train$o14, na.rm = TRUE),
    o21 = mean(train$o21[nr], na.rm = TRUE)
  )
}

apply_centers <- function(dat, centers) {
  out <- dat
  out$o11c <- out$o11 - centers[["o11"]]
  out$o12c <- out$o12 - centers[["o12"]]
  out$o13c <- out$o13 - centers[["o13"]]
  out$o14c <- out$o14 - centers[["o14"]]
  out$o21c <- out$o21 - centers[["o21"]]
  out
}

candidate_specs <- function(response = "y",
                            stage2_library = c("four_model", "nested_M22_M23", "strong_M20_M23")) {
  stage2_library <- match.arg(stage2_library)
  s2_all <- list(
    M20 = c("o11c", "o12c", "o13c", "o14c", "o21c", "a1", "o22", "a2"),
    M21 = c("o11c", "o12c", "o13c", "o14c", "o21c", "a1", "o22", "a2", "a1:a2"),
    M22 = c("o11c", "o12c", "o13c", "o14c", "o21c", "a1", "o22", "a2", "o22:a2"),
    M23 = c("o11c", "o12c", "o13c", "o14c", "o21c", "a1", "o22", "a2", "a1:a2", "o22:a2")
  )
  s1 <- list(
    M10 = c("o11c", "o12c", "o13", "a1"),
    M11 = c("o11c", "o12c", "o13", "a1", "o13:a1")
  )
  s2 <- switch(
    stage2_library,
    four_model = s2_all,
    nested_M22_M23 = s2_all[c("M22", "M23")],
    strong_M20_M23 = s2_all[c("M20", "M23")]
  )
  list(
    stage2 = lapply(s2, function(rhs) reformulate(rhs, response = response)),
    stage1 = lapply(s1, function(rhs) reformulate(rhs, response = response))
  )
}

aicc_lm <- function(fit) {
  ll <- logLik(fit)
  k <- attr(ll, "df")
  n <- nobs(fit)
  den <- n - k - 1
  if (den <= 0) return(Inf)
  AIC(fit) + 2 * k * (k + 1) / den
}

criterion_value <- function(fit, criterion) {
  criterion <- match.arg(criterion, c("AIC", "BIC", "AICc"))
  if (criterion == "AIC") return(AIC(fit))
  if (criterion == "BIC") return(BIC(fit))
  aicc_lm(fit)
}

ic_weights <- function(ic) {
  if (any(!is.finite(ic))) {
    finite <- is.finite(ic)
    if (!any(finite)) stop("No finite information criterion values.")
    w <- rep(0, length(ic))
    d <- ic[finite] - min(ic[finite])
    w[finite] <- exp(-0.5 * d)
  } else {
    d <- ic - min(ic)
    w <- exp(-0.5 * d)
  }
  w / sum(w)
}

fit_library <- function(dat, formulas, criterion) {
  fits <- lapply(formulas, function(f) lm(f, data = dat, na.action = na.omit))
  nobs_vec <- vapply(fits, nobs, numeric(1))
  if (length(unique(nobs_vec)) != 1L) {
    stop(
      "Information-criterion comparison requires identical analysis samples; nobs = ",
      paste(names(formulas), nobs_vec, sep = ":", collapse = ", ")
    )
  }
  rank_ok <- vapply(fits, function(f) f$rank == length(coef(f)), logical(1))
  if (any(!rank_ok)) {
    stop("Rank-deficient candidate model(s): ", paste(names(formulas)[!rank_ok], collapse = ", "))
  }

  ic <- vapply(fits, criterion_value, numeric(1), criterion = criterion)
  w <- ic_weights(ic)
  sel_idx <- which.min(ic)[1]
  names(ic) <- names(formulas)
  names(w) <- names(formulas)

  tab <- do.call(rbind, lapply(seq_along(fits), function(j) {
    fit <- fits[[j]]
    ll <- logLik(fit)
    k_all <- attr(ll, "df")
    data.frame(
      model = names(formulas)[j],
      n = nobs(fit),
      p_beta = length(coef(fit)),
      k_ic = k_all,
      residual_df = df.residual(fit),
      criterion = criterion,
      ic = ic[j],
      delta_ic = ic[j] - min(ic),
      weight = w[j],
      selected = (j == sel_idx),
      aicc_correction = aicc_lm(fit) - AIC(fit),
      stringsAsFactors = FALSE
    )
  }))

  list(
    fits = fits,
    ic = ic,
    weights = w,
    selected_index = sel_idx,
    selected_model = names(formulas)[sel_idx],
    table = tab
  )
}

predict_library <- function(lib, newdata, mode = c("SEL", "MA")) {
  mode <- match.arg(mode)
  if (mode == "SEL") {
    return(as.numeric(predict(lib$fits[[lib$selected_index]], newdata = newdata)))
  }
  pp <- sapply(lib$fits, function(fit) as.numeric(predict(fit, newdata = newdata)))
  if (is.null(dim(pp))) pp <- matrix(pp, ncol = length(lib$fits))
  as.numeric(pp %*% lib$weights)
}

coef_union <- function(lib, mode = c("SEL", "MA")) {
  mode <- match.arg(mode)
  all_names <- unique(unlist(lapply(lib$fits, function(f) names(coef(f)))))
  mat <- matrix(0, nrow = length(lib$fits), ncol = length(all_names),
                dimnames = list(names(lib$fits), all_names))
  for (j in seq_along(lib$fits)) {
    cc <- coef(lib$fits[[j]])
    cc[is.na(cc)] <- 0
    mat[j, names(cc)] <- cc
  }
  if (mode == "SEL") return(mat[lib$selected_index, ])
  out <- as.numeric(lib$weights %*% mat)
  setNames(out, colnames(mat))
}

make_q_actions <- function(lib, dat, treatment, mode = c("SEL", "MA")) {
  mode <- match.arg(mode)
  dm <- dat
  dp <- dat
  dm[[treatment]] <- -1
  dp[[treatment]] <- 1
  qm <- predict_library(lib, dm, mode = mode)
  qp <- predict_library(lib, dp, mode = mode)
  contrast <- qp - qm
  action <- ifelse(contrast >= 0, 1, -1)
  data.frame(q_minus = qm, q_plus = qp, contrast = contrast, action = action)
}

make_stage2_pseudo <- function(train_c, eval_c, lib2, mode = c("SEL", "MA")) {
  mode <- match.arg(mode)
  nr_tr <- stage2_eligible(train_c)
  nr_ev <- stage2_eligible(eval_c)

  pseudo_tr <- train_c$y
  pseudo_ev <- eval_c$y

  act_tr <- data.frame(q_minus = rep(NA_real_, nrow(train_c)),
                       q_plus = NA_real_, contrast = NA_real_, action = NA_integer_)
  act_ev <- data.frame(q_minus = rep(NA_real_, nrow(eval_c)),
                       q_plus = NA_real_, contrast = NA_real_, action = NA_integer_)

  if (any(nr_tr)) {
    tmp <- make_q_actions(lib2, train_c[nr_tr, , drop = FALSE], "a2", mode = mode)
    pseudo_tr[nr_tr] <- pmax(tmp$q_minus, tmp$q_plus)
    act_tr[nr_tr, ] <- tmp
  }
  if (any(nr_ev)) {
    tmp <- make_q_actions(lib2, eval_c[nr_ev, , drop = FALSE], "a2", mode = mode)
    pseudo_ev[nr_ev] <- pmax(tmp$q_minus, tmp$q_plus)
    act_ev[nr_ev, ] <- tmp
  }

  list(pseudo_train = pseudo_tr, pseudo_eval = pseudo_ev,
       stage2_train = act_tr, stage2_eval = act_ev)
}

fit_stage1_for_pseudo <- function(train_c, eval_c, pseudo_train, criterion, mode) {
  specs <- candidate_specs(response = "pseudo")$stage1
  d <- train_c
  d$pseudo <- pseudo_train
  lib1 <- fit_library(d, specs, criterion = criterion)
  act1 <- make_q_actions(lib1, eval_c, "a1", mode = mode)
  list(lib = lib1, actions = act1)
}

run_recursive_pair <- function(train, eval = train, criterion = "AIC",
                               stage2_library = c("four_model", "nested_M22_M23", "strong_M20_M23")) {
  stage2_library <- match.arg(stage2_library)
  criterion <- match.arg(criterion, c("AIC", "BIC", "AICc"))
  centers <- compute_centers(train)
  train_c <- apply_centers(train, centers)
  eval_c <- apply_centers(eval, centers)

  nr <- stage2_eligible(train_c)
  if (sum(nr) < 20L) stop("Too few stage-2 observations after restricting to nonresponders.")

  s2_specs <- candidate_specs(response = "y", stage2_library = stage2_library)$stage2
  lib2 <- fit_library(train_c[nr, , drop = FALSE], s2_specs, criterion = criterion)

  psel <- make_stage2_pseudo(train_c, eval_c, lib2, mode = "SEL")
  pma  <- make_stage2_pseudo(train_c, eval_c, lib2, mode = "MA")

  s1_sel <- fit_stage1_for_pseudo(train_c, eval_c, psel$pseudo_train, criterion, "SEL")
  s1_ma  <- fit_stage1_for_pseudo(train_c, eval_c, pma$pseudo_train, criterion, "MA")

  list(
    criterion = criterion,
    stage2_library_type = stage2_library,
    centers = centers,
    train_c = train_c,
    eval_c = eval_c,
    stage2_lib = lib2,
    stage2_sel = psel,
    stage2_ma = pma,
    stage1_sel = s1_sel,
    stage1_ma = s1_ma
  )
}

fixed_stage1_diagnostic <- function(res) {
  d <- res$train_c
  e <- res$eval_c
  f_wide <- reformulate(c("o11c", "o12c", "o13", "a1", "o13:a1"), response = "pseudo")
  f_narrow <- reformulate(c("o11c", "o12c", "o13", "a1"), response = "pseudo")

  ds <- d; ds$pseudo <- res$stage2_sel$pseudo_train
  dm <- d; dm$pseudo <- res$stage2_ma$pseudo_train

  wide_s <- lm(f_wide, data = ds)
  wide_m <- lm(f_wide, data = dm)
  nar_s <- lm(f_narrow, data = ds)
  nar_m <- lm(f_narrow, data = dm)

  pred_actions <- function(fit) {
    em <- e; ep <- e
    em$a1 <- -1; ep$a1 <- 1
    qm <- as.numeric(predict(fit, newdata = em))
    qp <- as.numeric(predict(fit, newdata = ep))
    cc <- qp - qm
    data.frame(q_minus = qm, q_plus = qp, contrast = cc,
               action = ifelse(cc >= 0, 1, -1))
  }

  asel <- pred_actions(wide_s)
  ama <- pred_actions(wide_m)

  u <- res$stage2_sel$pseudo_train - res$stage2_ma$pseudo_train
  X <- model.matrix(~ o11c + o12c + o13 + a1, data = d)
  z <- d$o13 * d$a1
  qx <- qr(X)
  z_res <- qr.resid(qx, z)
  gamma_projection <- sum(z_res * u) / sum(z_res^2)
  gamma_direct <- unname(coef(wide_s)["o13:a1"] - coef(wide_m)["o13:a1"])

  lr <- function(nar, wide) {
    n <- nobs(wide)
    rss0 <- sum(residuals(nar)^2)
    rss1 <- sum(residuals(wide)^2)
    n * log(rss0 / rss1)
  }
  lr_s <- lr(nar_s, wide_s)
  lr_m <- lr(nar_m, wide_m)

  coef_names <- union(names(coef(wide_s)), names(coef(wide_m)))
  cs <- setNames(rep(0, length(coef_names)), coef_names)
  cm <- cs
  cs[names(coef(wide_s))] <- coef(wide_s)
  cm[names(coef(wide_m))] <- coef(wide_m)

  list(
    action_sel = asel,
    action_ma = ama,
    coef_table = data.frame(
      term = coef_names,
      fixed_stage1_from_stage2_SEL = cs,
      fixed_stage1_from_stage2_MA = cm,
      difference = cs - cm,
      row.names = NULL
    ),
    diagnostics = data.frame(
      criterion = res$criterion,
      n = nrow(d),
      pseudo_diff_mean = mean(u),
      pseudo_diff_sd = sd(u),
      pseudo_diff_rms = sqrt(mean(u^2)),
      pseudo_diff_max_abs = max(abs(u)),
      residualized_score_over_sqrt_n = sum(z_res * u) / sqrt(nrow(d)),
      gamma_shift_direct = gamma_direct,
      gamma_shift_projection = gamma_projection,
      sqrt_n_gamma_shift_direct = sqrt(nrow(d)) * gamma_direct,
      sqrt_n_gamma_shift_projection = sqrt(nrow(d)) * gamma_projection,
      gamma_identity_error = gamma_direct - gamma_projection,
      lr_stage1_from_stage2_SEL = lr_s,
      lr_stage1_from_stage2_MA = lr_m,
      lr_shift = lr_s - lr_m,
      fixed_stage1_decision_discordance_n = sum(asel$action != ama$action),
      fixed_stage1_decision_discordance_rate = mean(asel$action != ama$action)
    )
  )
}

decision_summary <- function(a_sel, a_ma, ids = NULL, strata = NULL, criterion = NA_character_, stage = 1L) {
  if (is.null(ids)) ids <- seq_len(nrow(a_sel))
  ok <- !is.na(a_sel$action) & !is.na(a_ma$action)
  base <- data.frame(
    criterion = criterion,
    stage = stage,
    stratum = "all",
    n = sum(ok),
    sel_plus_rate = mean(a_sel$action[ok] == 1),
    ma_plus_rate = mean(a_ma$action[ok] == 1),
    discordant_n = sum(a_sel$action[ok] != a_ma$action[ok]),
    discordant_rate = mean(a_sel$action[ok] != a_ma$action[ok]),
    mean_abs_contrast_difference = mean(abs(a_sel$contrast[ok] - a_ma$contrast[ok])),
    max_abs_contrast_difference = max(abs(a_sel$contrast[ok] - a_ma$contrast[ok])),
    contrast_correlation = suppressWarnings(cor(a_sel$contrast[ok], a_ma$contrast[ok]))
  )
  if (is.null(strata)) return(base)
  lev <- unique(as.character(strata[ok]))
  extra <- lapply(lev, function(z) {
    ii <- ok & as.character(strata) == z
    data.frame(
      criterion = criterion,
      stage = stage,
      stratum = z,
      n = sum(ii),
      sel_plus_rate = mean(a_sel$action[ii] == 1),
      ma_plus_rate = mean(a_ma$action[ii] == 1),
      discordant_n = sum(a_sel$action[ii] != a_ma$action[ii]),
      discordant_rate = mean(a_sel$action[ii] != a_ma$action[ii]),
      mean_abs_contrast_difference = mean(abs(a_sel$contrast[ii] - a_ma$contrast[ii])),
      max_abs_contrast_difference = max(abs(a_sel$contrast[ii] - a_ma$contrast[ii])),
      contrast_correlation = suppressWarnings(cor(a_sel$contrast[ii], a_ma$contrast[ii]))
    )
  })
  do.call(rbind, c(list(base), extra))
}

decision_margin_diagnostics <- function(a_sel, a_ma, criterion = NA_character_,
                                        stage = 1L, analysis_variant = NA_character_) {
  ok <- is.finite(a_sel$contrast) & is.finite(a_ma$contrast)
  if (!any(ok)) stop("No finite decision contrasts available for margin diagnostics.")
  abs_sel <- abs(a_sel$contrast[ok])
  abs_ma <- abs(a_ma$contrast[ok])
  diff <- abs(a_sel$contrast[ok] - a_ma$contrast[ok])
  qfun <- function(x, p) unname(quantile(x, p, names = FALSE, type = 7))
  data.frame(
    analysis_variant = analysis_variant,
    criterion = criterion,
    stage = stage,
    n = sum(ok),
    discordant_n = sum(a_sel$action[ok] != a_ma$action[ok]),
    min_abs_contrast_SEL = min(abs_sel),
    q10_abs_contrast_SEL = qfun(abs_sel, 0.10),
    median_abs_contrast_SEL = median(abs_sel),
    min_abs_contrast_MA = min(abs_ma),
    q10_abs_contrast_MA = qfun(abs_ma, 0.10),
    median_abs_contrast_MA = median(abs_ma),
    mean_abs_SEL_MA_contrast_difference = mean(diff),
    q90_abs_SEL_MA_contrast_difference = qfun(diff, 0.90),
    max_abs_SEL_MA_contrast_difference = max(diff),
    max_difference_relative_to_min_margin = max(diff) / max(min(c(abs_sel, abs_ma)), .Machine$double.eps),
    stringsAsFactors = FALSE
  )
}

nested_stage2_diagnostic <- function(res, narrow_model, wide_model, analysis_variant = NA_character_) {
  nm <- names(res$stage2_lib$fits)
  if (!all(c(narrow_model, wide_model) %in% nm)) {
    stop("Requested nested models are not both present in the stage-2 library.")
  }
  f0 <- res$stage2_lib$fits[[narrow_model]]
  f1 <- res$stage2_lib$fits[[wide_model]]
  n <- nobs(f1)
  if (nobs(f0) != n) stop("Nested LR comparison requires identical nobs.")
  d <- length(coef(f1)) - length(coef(f0))
  if (d <= 0L) stop("wide_model must have more regression coefficients than narrow_model.")
  rss0 <- sum(residuals(f0)^2)
  rss1 <- sum(residuals(f1)^2)
  lambda <- n * log(rss0 / rss1)
  data.frame(
    analysis_variant = analysis_variant,
    criterion = res$criterion,
    narrow_model = narrow_model,
    wide_model = wide_model,
    n = n,
    added_parameters = d,
    lambda = lambda,
    AIC_threshold_2d = 2 * d,
    BIC_threshold_d_log_n = d * log(n),
    lambda_minus_AIC_threshold = lambda - 2 * d,
    lambda_minus_BIC_threshold = lambda - d * log(n),
    narrow_weight = unname(res$stage2_lib$weights[narrow_model]),
    wide_weight = unname(res$stage2_lib$weights[wide_model]),
    selected_model = res$stage2_lib$selected_model,
    stringsAsFactors = FALSE
  )
}

model_coefficients_table <- function(res) {
  c2_sel <- coef_union(res$stage2_lib, "SEL")
  c2_ma <- coef_union(res$stage2_lib, "MA")
  c1_sel <- coef_union(res$stage1_sel$lib, "SEL")
  c1_ma <- coef_union(res$stage1_ma$lib, "MA")

  make <- function(stage, cs, cm) {
    nm <- union(names(cs), names(cm))
    xs <- setNames(rep(0, length(nm)), nm)
    xm <- xs
    xs[names(cs)] <- cs
    xm[names(cm)] <- cm
    data.frame(
      criterion = res$criterion,
      stage = stage,
      term = nm,
      SEL = xs,
      MA = xm,
      SEL_minus_MA = xs - xm,
      row.names = NULL
    )
  }
  rbind(make(2, c2_sel, c2_ma), make(1, c1_sel, c1_ma))
}

reference_reproduction <- function(dat) {
  centers <- compute_centers(dat)
  d <- apply_centers(dat, centers)
  nr <- stage2_eligible(d)

  f2 <- reformulate(c("o11c", "o12c", "o13c", "o14c", "o21c", "a1", "o22", "a2", "a1:a2", "o22:a2"), response = "y")
  fit2 <- lm(f2, data = d[nr, , drop = FALSE])

  dm <- d[nr, , drop = FALSE]; dp <- dm
  dm$a2 <- -1; dp$a2 <- 1
  pseudo <- d$y
  pseudo[nr] <- pmax(predict(fit2, newdata = dm), predict(fit2, newdata = dp))

  d$pseudo <- pseudo
  f1 <- reformulate(c("o11c", "o12c", "o13", "a1", "o13:a1"), response = "pseudo")
  fit1 <- lm(f1, data = d)

  ref2 <- c(
    "(Intercept)" = 3.0039,
    "o11c" = -0.2462,
    "o12c" = -0.2961,
    "o13c" = 0.0391,
    "o14c" = 0.4868,
    "o21c" = -0.0097,
    "a1" = 0.0758,
    "o22" = -0.0980,
    "a2" = -0.8640,
    "a1:a2" = -0.1934,
    "o22:a2" = 1.1826
  )
  ref1 <- c(
    "(Intercept)" = 3.4497,
    "o11c" = -0.4556,
    "o12c" = -0.3458,
    "o13" = -0.0236,
    "a1" = 0.2934,
    "o13:a1" = -0.5254
  )

  compare <- function(stage, est, ref) {
    nm <- names(ref)
    ee <- est[nm]
    data.frame(
      stage = stage,
      term = nm,
      estimate = unname(ee),
      published_reference = unname(ref),
      absolute_difference = abs(unname(ee - ref)),
      reference_source = "d3center PROC QLEARN published output",
      row.names = NULL
    )
  }
  tab <- rbind(compare(2, coef(fit2), ref2), compare(1, coef(fit1), ref1))
  attr(tab, "max_abs_difference") <- max(tab$absolute_difference, na.rm = TRUE)
  tab
}

write_base_plots <- function(res_list, dat, outdir) {
  for (criterion in names(res_list)) {
    res <- res_list[[criterion]]
    a1s <- res$stage1_sel$actions
    a1m <- res$stage1_ma$actions

    png(file.path(outdir, paste0("stage1_contrast_SEL_vs_MA_", criterion, ".png")), width = 900, height = 800)
    plot(a1s$contrast, a1m$contrast,
         xlab = paste0(criterion, " selection: stage-1 Q contrast (+1 minus -1)"),
         ylab = paste0(criterion, " averaging: stage-1 Q contrast (+1 minus -1)"),
         main = paste0("Stage-1 decision contrast: ", criterion))
    abline(h = 0, v = 0, lty = 3)
    abline(0, 1, lty = 2)
    dev.off()

    nr <- stage2_eligible(dat)
    u <- res$stage2_sel$pseudo_eval - res$stage2_ma$pseudo_eval
    png(file.path(outdir, paste0("stage1_pseudooutcome_shift_", criterion, ".png")), width = 900, height = 700)
    boxplot(u ~ factor(dat$r, levels = c(0, 1), labels = c("Nonresponder", "Responder")),
            ylab = "Stage-1 pseudo-outcome: stage-2 SEL minus MA",
            xlab = "Stage-1 response status",
            main = paste0("Generated-response perturbation: ", criterion))
    abline(h = 0, lty = 3)
    dev.off()

    a2s <- res$stage2_sel$stage2_eval[nr, ]
    a2m <- res$stage2_ma$stage2_eval[nr, ]
    png(file.path(outdir, paste0("stage2_contrast_SEL_vs_MA_", criterion, ".png")), width = 900, height = 800)
    plot(a2s$contrast, a2m$contrast,
         xlab = paste0(criterion, " selection: stage-2 Q contrast (+1 minus -1)"),
         ylab = paste0(criterion, " averaging: stage-2 Q contrast (+1 minus -1)"),
         main = paste0("Stage-2 decision contrast: ", criterion))
    abline(h = 0, v = 0, lty = 3)
    abline(0, 1, lty = 2)
    dev.off()
  }
}

bootstrap_stability <- function(dat, B = 2000L, criterion = "AIC", cores = 1L, seed = 20260907L,
                                stage2_library = "four_model") {
  if (B <= 0L) return(NULL)
  criterion <- match.arg(criterion, c("AIC", "BIC", "AICc"))
  n <- nrow(dat)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)

  one <- function(b) {
    set.seed(seeds[b])
    idx <- sample.int(n, n, replace = TRUE)
    tr <- dat[idx, , drop = FALSE]
    z <- try(run_recursive_pair(tr, dat, criterion = criterion, stage2_library = stage2_library), silent = TRUE)
    if (inherits(z, "try-error")) {
      return(list(ok = FALSE, error = as.character(z)))
    }
    list(
      ok = TRUE,
      a_sel = z$stage1_sel$actions$action,
      a_ma = z$stage1_ma$actions$action,
      c_sel = z$stage1_sel$actions$contrast,
      c_ma = z$stage1_ma$actions$contrast,
      s2_model = z$stage2_lib$selected_model,
      s1_sel_model = z$stage1_sel$lib$selected_model,
      s1_ma_top_model = names(which.max(z$stage1_ma$lib$weights))[1]
    )
  }

  if (.Platform$OS.type == "unix" && cores > 1L) {
    rr <- parallel::mclapply(seq_len(B), one, mc.cores = cores, mc.set.seed = FALSE)
  } else {
    rr <- lapply(seq_len(B), one)
  }

  ok <- vapply(rr, function(x) isTRUE(x$ok), logical(1))
  if (!any(ok)) stop("All bootstrap replicates failed.")
  good <- rr[ok]

  amat_s <- do.call(rbind, lapply(good, `[[`, "a_sel"))
  amat_m <- do.call(rbind, lapply(good, `[[`, "a_ma"))
  cmat_s <- do.call(rbind, lapply(good, `[[`, "c_sel"))
  cmat_m <- do.call(rbind, lapply(good, `[[`, "c_ma"))

  subject <- data.frame(
    id = dat$id,
    criterion = criterion,
    B_requested = B,
    B_success = sum(ok),
    SEL_prob_plus1 = colMeans(amat_s == 1),
    MA_prob_plus1 = colMeans(amat_m == 1),
    SEL_MA_discordance_prob = colMeans(amat_s != amat_m),
    SEL_contrast_mean = colMeans(cmat_s),
    MA_contrast_mean = colMeans(cmat_m),
    stringsAsFactors = FALSE
  )

  rep_summary <- data.frame(
    replicate = which(ok),
    SEL_plus1_rate = rowMeans(amat_s == 1),
    MA_plus1_rate = rowMeans(amat_m == 1),
    SEL_MA_discordance_rate = rowMeans(amat_s != amat_m),
    stage2_selected_model = vapply(good, `[[`, character(1), "s2_model"),
    stage1_SEL_selected_model = vapply(good, `[[`, character(1), "s1_sel_model"),
    stage1_MA_top_weight_model = vapply(good, `[[`, character(1), "s1_ma_top_model"),
    stringsAsFactors = FALSE
  )

  list(
    subject = subject,
    replicate = rep_summary,
    failures = data.frame(
      replicate = which(!ok),
      error = vapply(rr[!ok], function(x) x$error, character(1)),
      stringsAsFactors = FALSE
    )
  )
}
