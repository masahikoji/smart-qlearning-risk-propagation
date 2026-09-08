# Recursive Q-learning procedures and generated-response diagnostics.

predicted_value <- function(beta_union, data, stage) {
  pmax(predict_embedded_q(beta_union, data, stage, 1),
       predict_embedded_q(beta_union, data, stage, -1))
}

oracle_pseudo_outcome <- function(scenario, n, data, stage, quad) {
  stopifnot(stage < scenario$K)
  data[[paste0("Y", stage)]] + true_v(scenario, n, stage + 1L, data, quad)
}

wide_added_coef <- function(fits, added_feature) {
  if (is.na(added_feature) || !("wide" %in% names(fits))) return(NA_real_)
  f <- fits[["wide"]]
  if (!f$valid || !(added_feature %in% names(f$coef_union))) return(NA_real_)
  unname(f$coef_union[added_feature])
}

residualized_score_shift <- function(data, stage, spec, u) {
  if (is.na(spec$added_feature) || !all(c("narrow", "wide") %in% names(spec$candidates))) return(NA_real_)
  narrow_feats <- spec$candidates[["narrow"]]
  zfeat <- spec$added_feature
  X <- feature_matrix(data, stage, data[[paste0("A", stage)]], narrow_feats)
  Z <- feature_matrix(data, stage, data[[paste0("A", stage)]], zfeat)
  zres <- .lm.fit(X, Z[, 1L])$residuals
  sum(zres * u)
}

fit_recursive_procedures <- function(scenario, n, data, eval_cache,
                                     stabilize_tol = 1e-10,
                                     check_ic_identity = TRUE) {
  K <- scenario$K
  pdefs <- procedure_definitions()
  proc_names <- names(pdefs)
  quad <- eval_cache$quad

  # Storage: aggregate[[proc]][[stage]], fits[[proc]][[stage]], pseudo[[proc]][[stage]].
  aggregate <- setNames(lapply(proc_names, function(.) vector("list", K)), proc_names)
  fits_store <- setNames(lapply(proc_names, function(.) vector("list", K)), proc_names)
  pseudo_store <- setNames(lapply(proc_names, function(.) vector("list", K)), proc_names)
  oracle_fits <- vector("list", K)
  oracle_aggregate <- setNames(lapply(proc_names, function(.) vector("list", K)), proc_names)
  oracle_pseudo <- vector("list", K)

  # Final stage: same response and candidate fits for every recursive procedure.
  s <- K
  specK <- stage_model_spec(scenario, s)
  fK <- fit_candidate_library(data, s, data[[paste0("Y", s)]], specK, stabilize_tol)
  if (check_ic_identity && scenario$candidate_count == 2L) validate_nested_ic_identity(fK, n)
  for (r in proc_names) {
    aggregate[[r]][[s]] <- aggregate_library(fK, pdefs[[r]]$criterion, pdefs[[r]]$mode)
    fits_store[[r]][[s]] <- fK
    pseudo_store[[r]][[s]] <- data[[paste0("Y", s)]]
  }

  if (K > 1L) {
    for (s in (K - 1L):1L) {
      spec <- stage_model_spec(scenario, s)
      # Oracle fit is common because the candidate specifications are common.
      y0 <- oracle_pseudo_outcome(scenario, n, data, s, quad)
      oracle_pseudo[[s]] <- y0
      f0 <- fit_candidate_library(data, s, y0, spec, stabilize_tol)
      oracle_fits[[s]] <- f0
      if (check_ic_identity && length(f0) == 2L) validate_nested_ic_identity(f0, n)
      for (r in proc_names) {
        oracle_aggregate[[r]][[s]] <- aggregate_library(f0, pdefs[[r]]$criterion, pdefs[[r]]$mode)
      }

      for (r in proc_names) {
        vnext <- predicted_value(aggregate[[r]][[s + 1L]]$coef, data, s + 1L)
        yr <- data[[paste0("Y", s)]] + vnext
        fr <- fit_candidate_library(data, s, yr, spec, stabilize_tol)
        if (check_ic_identity && length(fr) == 2L) validate_nested_ic_identity(fr, n)
        ar <- aggregate_library(fr, pdefs[[r]]$criterion, pdefs[[r]]$mode)
        aggregate[[r]][[s]] <- ar
        fits_store[[r]][[s]] <- fr
        pseudo_store[[r]][[s]] <- yr
      }
    }
  }

  list(aggregate = aggregate, fits = fits_store, pseudo = pseudo_store,
       oracle_fits = oracle_fits, oracle_aggregate = oracle_aggregate,
       oracle_pseudo = oracle_pseudo)
}

flatten_replicate <- function(scenario, n, fit_obj, eval_cache) {
  K <- scenario$K
  pdefs <- procedure_definitions()
  proc_names <- names(pdefs)
  out <- c()

  for (r in proc_names) {
    for (s in seq_len(K)) {
      ar <- fit_obj$aggregate[[r]][[s]]
      fr <- fit_obj$fits[[r]][[s]]
      cache <- eval_cache$risk[[paste0("stage", s)]]
      prefix <- paste0(r, ".s", s, ".")
      out[paste0(prefix, "risk")] <- risk_from_cache(ar$coef, cache)
      out[paste0(prefix, "stabilized")] <- as.numeric(ar$stabilized)
      # Candidate selection/weights.
      cand_names <- names(fr)
      out[paste0(prefix, "selected_index")] <- match(ar$selected, cand_names)
      for (m in seq_along(ar$weights)) {
        out[paste0(prefix, "w_", names(ar$weights)[m])] <- ar$weights[m]
      }
      # Save aggregate coefficients so secondary metrics can be computed without refitting.
      for (nm in names(ar$coef)) out[paste0(prefix, "coef_", nm)] <- ar$coef[nm]

      nd <- nested_ic_diagnostics(fr, n)
      if (!is.null(nd)) {
        out[paste0(prefix, "lambda")] <- nd$lambda
        out[paste0(prefix, "aic_identity_error")] <- nd$aic_identity_error
        out[paste0(prefix, "bic_identity_error")] <- nd$bic_identity_error
      }

      if (s < K && !is.null(fit_obj$oracle_fits[[s]])) {
        spec <- stage_model_spec(scenario, s)
        f0 <- fit_obj$oracle_fits[[s]]
        a0 <- fit_obj$oracle_aggregate[[r]][[s]]
        y0 <- fit_obj$oracle_pseudo[[s]]
        yr <- fit_obj$pseudo[[r]][[s]]
        u <- yr - y0
        nd0 <- nested_ic_diagnostics(f0, n)
        out[paste0(prefix, "pseudo_rms")] <- sqrt(mean(u^2))
        score_shift <- residualized_score_shift(
          data = attr(fit_obj, "data"), stage = s, spec = spec, u = u)
        out[paste0(prefix, "score_shift")] <- score_shift
        out[paste0(prefix, "score_shift_over_sqrt_n")] <- score_shift / sqrt(n)
        wide_shift <- wide_added_coef(fr, spec$added_feature) -
          wide_added_coef(f0, spec$added_feature)
        out[paste0(prefix, "wide_coef_shift")] <- wide_shift
        out[paste0(prefix, "sqrt_n_wide_coef_shift")] <- sqrt(n) * wide_shift
        if (!is.null(nd) && !is.null(nd0)) {
          out[paste0(prefix, "lambda_shift")] <- nd$lambda - nd0$lambda
        }
        out[paste0(prefix, "selected_disagree_oracle")] <-
          if (pdefs[[r]]$mode == "SEL") as.numeric(ar$selected != a0$selected) else NA_real_
        if (pdefs[[r]]$mode == "MA") {
          common_w <- intersect(names(ar$weights), names(a0$weights))
          for (wm in common_w) {
            out[paste0(prefix, "weight_shift_oracle_", wm)] <- ar$weights[wm] - a0$weights[wm]
          }
        }
        if (pdefs[[r]]$mode == "MA" && "wide" %in% names(ar$weights) && "wide" %in% names(a0$weights)) {
          out[paste0(prefix, "wide_weight_shift_oracle")] <- ar$weights["wide"] - a0$weights["wide"]
        } else {
          out[paste0(prefix, "wide_weight_shift_oracle")] <- NA_real_
        }
        if ("narrow" %in% names(fr) && "narrow" %in% names(f0)) {
          bn1 <- fr[["narrow"]]$coef_union
          bn0 <- f0[["narrow"]]$coef_union
          out[paste0(prefix, "narrow_q_rms_shift")] <- qdiff_rms_from_cache(bn1, bn0, cache)
          out[paste0(prefix, "narrow_rss_shift")] <- fr[["narrow"]]$rss - f0[["narrow"]]$rss
        }
        if ("wide" %in% names(fr) && "wide" %in% names(f0)) {
          # Integrated RMS change of the fitted wide candidate on the evaluation law.
          b1 <- fr[["wide"]]$coef_union
          b0 <- f0[["wide"]]$coef_union
          out[paste0(prefix, "wide_q_rms_shift")] <- qdiff_rms_from_cache(b1, b0, cache)
          out[paste0(prefix, "wide_rss_shift")] <- fr[["wide"]]$rss - f0[["wide"]]$rss
        }
      }
    }
  }
  out
}

simulate_one_replicate <- function(scenario, n, rep_seed, eval_cache) {
  data <- simulate_smart_data(scenario, n, seed = rep_seed, include_rewards = TRUE)
  fit_obj <- fit_recursive_procedures(scenario, n, data, eval_cache)
  # Avoid carrying the full training dataset in the returned object, but make it
  # temporarily available to the diagnostic flattener.
  attr(fit_obj, "data") <- data
  ans <- flatten_replicate(scenario, n, fit_obj, eval_cache)
  attr(fit_obj, "data") <- NULL
  ans
}
