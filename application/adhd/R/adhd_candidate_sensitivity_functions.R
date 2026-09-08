# Candidate-library sensitivity extensions for the ADHD SMART SEL/MA analysis.
# This file is sourced AFTER R/adhd_sel_ma_functions.R.
# It reuses the original data loading, centering, IC, prediction, and diagnostic
# functions so that the sensitivity analysis differs only in prespecified
# candidate-library definitions.

make_formula_from_terms <- function(response, rhs_terms) {
  if (length(rhs_terms) == 0L) {
    return(stats::as.formula(paste(response, "~ 1")))
  }
  stats::as.formula(
    paste(response, "~", paste(rhs_terms, collapse = " + ")),
    env = parent.frame()
  )
}

make_formula_library <- function(response, rhs_list) {
  out <- lapply(rhs_list, function(rhs) make_formula_from_terms(response, rhs))
  names(out) <- names(rhs_list)
  out
}

adhd_sensitivity_library_terms <- function() {
  # Primary/source-aligned prognostic blocks from the original analysis.
  b2 <- c("o11c", "o12c", "o13c", "o14c", "o21c", "a1", "o22", "a2")
  b1 <- c("o11c", "o12c", "o13", "a1")

  primary_s2 <- list(
    M20 = b2,
    M21 = c(b2, "a1:a2"),
    M22 = c(b2, "o22:a2"),
    M23 = c(b2, "a1:a2", "o22:a2")
  )
  primary_s1 <- list(
    M10 = b1,
    M11 = c(b1, "o13:a1")
  )

  # Sensitivity A: broaden the set of plausible stage-2 treatment modifiers.
  # The common prognostic block is unchanged. We include the four original
  # candidates, one-covariate alternatives for every other observed stage-2
  # history variable already present in B2, and one deliberately broad model.
  # This avoids selecting a large collection of nearly duplicate subset models
  # after looking at the data, while still challenging the primary library.
  expanded_s2 <- list(
    E20 = b2,
    E21 = c(b2, "a1:a2"),
    E22 = c(b2, "o22:a2"),
    E23 = c(b2, "a1:a2", "o22:a2"),
    E24 = c(b2, "o11c:a2"),
    E25 = c(b2, "o12c:a2"),
    E26 = c(b2, "o13c:a2"),
    E27 = c(b2, "o14c:a2"),
    E28 = c(b2, "o21c:a2"),
    E29 = c(
      b2,
      "a1:a2", "o22:a2", "o11c:a2", "o12c:a2",
      "o13c:a2", "o14c:a2", "o21c:a2"
    )
  )

  # Sensitivity B: broaden stage-1 treatment-effect modification without
  # altering the source-aligned prognostic block. Because o11c, o12c, and o13
  # are already in B1, all 2^3 interaction subsets satisfy strong hierarchy.
  expanded_s1 <- list(
    E10 = b1,
    E11 = c(b1, "o13:a1"),
    E12 = c(b1, "o11c:a1"),
    E13 = c(b1, "o12c:a1"),
    E14 = c(b1, "o11c:a1", "o12c:a1"),
    E15 = c(b1, "o11c:a1", "o13:a1"),
    E16 = c(b1, "o12c:a1", "o13:a1"),
    E17 = c(b1, "o11c:a1", "o12c:a1", "o13:a1")
  )

  # Sensitivity C: change functional form using natural cubic splines for the
  # continuous variables o12 and o21 while retaining Gaussian linear-model
  # likelihoods, so AIC/AICc/BIC and normalized exp(-IC/2) weighting remain
  # directly comparable to the primary analysis.
  b2_spline <- c(
    "o11c", "splines::ns(o12c, df = 3)", "o13c", "o14c",
    "splines::ns(o21c, df = 3)", "a1", "o22", "a2"
  )
  spline_s2 <- list(
    N20 = b2_spline,
    N21 = c(b2_spline, "a1:a2"),
    N22 = c(b2_spline, "o22:a2"),
    N23 = c(b2_spline, "a1:a2", "o22:a2"),
    N24 = c(b2_spline, "a2:splines::ns(o12c, df = 3)"),
    N25 = c(b2_spline, "a2:splines::ns(o21c, df = 3)"),
    N26 = c(
      b2_spline,
      "a1:a2", "o22:a2",
      "a2:splines::ns(o12c, df = 3)",
      "a2:splines::ns(o21c, df = 3)"
    )
  )

  b1_spline <- c("o11c", "splines::ns(o12c, df = 3)", "o13", "a1")
  spline_s1 <- list(
    N10 = b1_spline,
    N11 = c(b1_spline, "o13:a1"),
    N12 = c(b1_spline, "a1:splines::ns(o12c, df = 3)"),
    N13 = c(
      b1_spline,
      "o13:a1",
      "a1:splines::ns(o12c, df = 3)"
    )
  )

  list(
    primary_reference = list(
      stage2 = primary_s2,
      stage1 = primary_s1,
      description = paste(
        "Original source-aligned application library:",
        "four stage-2 candidates and two stage-1 candidates."
      )
    ),
    expanded_stage2_effect_modifiers = list(
      stage2 = expanded_s2,
      stage1 = primary_s1,
      description = paste(
        "Stage-2 treatment-modifier library expanded to alternative baseline/intermediate",
        "modifiers; original stage-1 library retained."
      )
    ),
    expanded_stage1_effect_modifiers = list(
      stage2 = primary_s2,
      stage1 = expanded_s1,
      description = paste(
        "Original stage-2 library retained; stage-1 library expanded to all interaction",
        "subsets of o11c, o12c, and o13 with A1."
      )
    ),
    expanded_both_effect_modifiers = list(
      stage2 = expanded_s2,
      stage1 = expanded_s1,
      description = paste(
        "Both stages use the expanded effect-modifier libraries."
      )
    ),
    spline_functional_form = list(
      stage2 = spline_s2,
      stage1 = spline_s1,
      description = paste(
        "Natural-cubic-spline sensitivity for continuous history variables,",
        "with spline treatment interactions included as prespecified alternatives."
      )
    )
  )
}

adhd_sensitivity_formula_libraries <- function(response_stage2 = "y", response_stage1 = "pseudo") {
  defs <- adhd_sensitivity_library_terms()
  lapply(defs, function(v) {
    list(
      stage2 = make_formula_library(response_stage2, v$stage2),
      stage1 = make_formula_library(response_stage1, v$stage1),
      description = v$description
    )
  })
}

library_definition_table <- function() {
  defs <- adhd_sensitivity_library_terms()
  rows <- list()
  k <- 1L
  for (variant in names(defs)) {
    for (stage in c("stage2", "stage1")) {
      ll <- defs[[variant]][[stage]]
      for (model in names(ll)) {
        rows[[k]] <- data.frame(
          analysis_variant = variant,
          description = defs[[variant]]$description,
          stage = ifelse(stage == "stage2", 2L, 1L),
          model = model,
          rhs = paste(ll[[model]], collapse = " + "),
          number_of_candidates_at_stage = length(ll),
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  do.call(rbind, rows)
}

fit_stage1_for_custom_pseudo <- function(train_c, eval_c, pseudo_train,
                                         stage1_formulas, criterion, mode) {
  d <- train_c
  d$pseudo <- pseudo_train
  lib1 <- fit_library(d, stage1_formulas, criterion = criterion)
  act1 <- make_q_actions(lib1, eval_c, "a1", mode = mode)
  list(lib = lib1, actions = act1)
}

run_recursive_pair_custom <- function(train, eval = train, criterion = "AIC",
                                      stage2_formulas, stage1_formulas,
                                      analysis_variant = NA_character_) {
  criterion <- match.arg(criterion, c("AIC", "BIC", "AICc"))
  centers <- compute_centers(train)
  train_c <- apply_centers(train, centers)
  eval_c <- apply_centers(eval, centers)

  nr <- stage2_eligible(train_c)
  if (sum(nr) < 20L) stop("Too few stage-2 observations after restricting to nonresponders.")

  lib2 <- fit_library(
    train_c[nr, , drop = FALSE],
    stage2_formulas,
    criterion = criterion
  )

  psel <- make_stage2_pseudo(train_c, eval_c, lib2, mode = "SEL")
  pma  <- make_stage2_pseudo(train_c, eval_c, lib2, mode = "MA")

  s1_sel <- fit_stage1_for_custom_pseudo(
    train_c, eval_c, psel$pseudo_train,
    stage1_formulas = stage1_formulas,
    criterion = criterion, mode = "SEL"
  )
  s1_ma <- fit_stage1_for_custom_pseudo(
    train_c, eval_c, pma$pseudo_train,
    stage1_formulas = stage1_formulas,
    criterion = criterion, mode = "MA"
  )

  list(
    criterion = criterion,
    analysis_variant = analysis_variant,
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

weight_uncertainty_metrics <- function(w) {
  w <- as.numeric(w)
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) {
    return(c(top_weight = NA_real_, effective_models = NA_real_, entropy = NA_real_))
  }
  c(
    top_weight = max(w),
    effective_models = 1 / sum(w^2),
    entropy = -sum(w * log(w))
  )
}

sensitivity_propagation_summary <- function(res, dat) {
  nr <- stage2_eligible(dat)
  fd <- fixed_stage1_diagnostic(res)$diagnostics
  md1 <- decision_margin_diagnostics(
    res$stage1_sel$actions, res$stage1_ma$actions,
    criterion = res$criterion, stage = 1L,
    analysis_variant = res$analysis_variant
  )
  ds1 <- decision_summary(
    res$stage1_sel$actions, res$stage1_ma$actions,
    criterion = res$criterion, stage = 1L
  )
  ds2 <- decision_summary(
    res$stage2_sel$stage2_eval[nr, , drop = FALSE],
    res$stage2_ma$stage2_eval[nr, , drop = FALSE],
    criterion = res$criterion, stage = 2L
  )
  u <- res$stage2_sel$pseudo_eval - res$stage2_ma$pseudo_eval

  w2 <- weight_uncertainty_metrics(res$stage2_lib$weights)
  w1 <- weight_uncertainty_metrics(res$stage1_ma$lib$weights)

  data.frame(
    analysis_variant = res$analysis_variant,
    criterion = res$criterion,
    stage2_n_candidates = length(res$stage2_lib$fits),
    stage2_selected_model = res$stage2_lib$selected_model,
    stage2_top_MA_model = names(which.max(res$stage2_lib$weights))[1],
    stage2_top_weight = unname(w2["top_weight"]),
    stage2_effective_models = unname(w2["effective_models"]),
    stage2_weight_entropy = unname(w2["entropy"]),
    stage2_SEL_MA_decision_discordance_n = ds2$discordant_n[1],
    stage2_SEL_MA_decision_discordance_rate = ds2$discordant_rate[1],
    pseudo_diff_mean = mean(u),
    pseudo_diff_sd = stats::sd(u),
    pseudo_diff_rms = sqrt(mean(u^2)),
    pseudo_diff_max_abs = max(abs(u)),
    fixed_M11_gamma_shift = fd$gamma_shift_direct,
    fixed_M11_sqrt_n_gamma_shift = fd$sqrt_n_gamma_shift_direct,
    fixed_M11_LR_SEL = fd$lr_stage1_from_stage2_SEL,
    fixed_M11_LR_MA = fd$lr_stage1_from_stage2_MA,
    fixed_M11_LR_shift = fd$lr_shift,
    fixed_M11_decision_discordance_n = fd$fixed_stage1_decision_discordance_n,
    stage1_n_candidates = length(res$stage1_sel$lib$fits),
    stage1_SEL_selected_model = res$stage1_sel$lib$selected_model,
    stage1_MA_top_model = names(which.max(res$stage1_ma$lib$weights))[1],
    stage1_MA_top_weight = unname(w1["top_weight"]),
    stage1_MA_effective_models = unname(w1["effective_models"]),
    stage1_MA_weight_entropy = unname(w1["entropy"]),
    stage1_recursive_discordance_n = ds1$discordant_n[1],
    stage1_recursive_discordance_rate = ds1$discordant_rate[1],
    stage1_min_abs_contrast_SEL = md1$min_abs_contrast_SEL,
    stage1_min_abs_contrast_MA = md1$min_abs_contrast_MA,
    stage1_max_abs_SEL_MA_contrast_difference = md1$max_abs_SEL_MA_contrast_difference,
    stringsAsFactors = FALSE
  )
}

sensitivity_ic_table <- function(res) {
  bind_one <- function(tab, stage, target) {
    out <- tab
    out$analysis_variant <- res$analysis_variant
    out$stage <- stage
    out$procedure_target <- target
    out[, c("analysis_variant", "stage", "procedure_target", setdiff(names(out), c("analysis_variant", "stage", "procedure_target")))]
  }
  rbind(
    bind_one(res$stage2_lib$table, 2L, "common_terminal_fit"),
    bind_one(res$stage1_sel$lib$table, 1L, paste0(res$criterion, "_SEL_pseudooutcome")),
    bind_one(res$stage1_ma$lib$table, 1L, paste0(res$criterion, "_MA_pseudooutcome"))
  )
}

sensitivity_subject_table <- function(res, dat) {
  fd <- fixed_stage1_diagnostic(res)
  data.frame(
    id = dat$id,
    analysis_variant = res$analysis_variant,
    criterion = res$criterion,
    o11 = dat$o11,
    o12 = dat$o12,
    o13 = dat$o13,
    o14 = dat$o14,
    a1_observed = dat$a1,
    response = dat$r,
    stage1_contrast_SEL = res$stage1_sel$actions$contrast,
    stage1_contrast_MA = res$stage1_ma$actions$contrast,
    stage1_recommendation_SEL = ifelse(res$stage1_sel$actions$action == 1, "BMOD", "MED"),
    stage1_recommendation_MA = ifelse(res$stage1_ma$actions$action == 1, "BMOD", "MED"),
    stage1_recursive_discordant = res$stage1_sel$actions$action != res$stage1_ma$actions$action,
    stage1_pseudo_SEL = res$stage2_sel$pseudo_eval,
    stage1_pseudo_MA = res$stage2_ma$pseudo_eval,
    pseudo_SEL_minus_MA = res$stage2_sel$pseudo_eval - res$stage2_ma$pseudo_eval,
    fixed_M11_contrast_from_stage2_SEL = fd$action_sel$contrast,
    fixed_M11_contrast_from_stage2_MA = fd$action_ma$contrast,
    fixed_M11_downstream_only_discordant = fd$action_sel$action != fd$action_ma$action,
    stringsAsFactors = FALSE
  )
}

primary_reference_validation <- function(res) {
  if (res$analysis_variant != "primary_reference" || res$criterion != "AIC") {
    stop("primary_reference_validation requires primary_reference/AIC result.")
  }
  fd <- fixed_stage1_diagnostic(res)$diagnostics
  w <- res$stage2_lib$weights
  getw <- function(nm) if (nm %in% names(w)) unname(w[nm]) else NA_real_

  observed <- c(
    w_M22 = getw("M22"),
    w_M23 = getw("M23"),
    pseudo_rms = fd$pseudo_diff_rms,
    gamma_shift = fd$gamma_shift_direct,
    sqrt_n_gamma_shift = fd$sqrt_n_gamma_shift_direct,
    lr_shift = fd$lr_shift,
    recursive_stage1_discordance_n = sum(res$stage1_sel$actions$action != res$stage1_ma$actions$action)
  )
  expected <- c(
    w_M22 = 0.245,
    w_M23 = 0.755,
    pseudo_rms = 0.04186,
    gamma_shift = -0.00808,
    sqrt_n_gamma_shift = -0.09895,
    lr_shift = 0.24019,
    recursive_stage1_discordance_n = 0
  )
  tolerance <- c(
    w_M22 = 0.015,
    w_M23 = 0.015,
    pseudo_rms = 0.003,
    gamma_shift = 0.0015,
    sqrt_n_gamma_shift = 0.02,
    lr_shift = 0.04,
    recursive_stage1_discordance_n = 0
  )
  pass <- abs(observed - expected) <= tolerance
  pass["recursive_stage1_discordance_n"] <- observed["recursive_stage1_discordance_n"] == expected["recursive_stage1_discordance_n"]

  data.frame(
    quantity = names(observed),
    expected_from_current_manuscript = as.numeric(expected),
    observed = as.numeric(observed),
    absolute_difference = abs(as.numeric(observed - expected)),
    tolerance = as.numeric(tolerance),
    pass = as.logical(pass),
    stringsAsFactors = FALSE
  )
}

write_sensitivity_plots <- function(all_results, outdir) {
  variants <- names(all_results)
  # One plot per variant for AIC, the theory-aligned primary criterion.
  for (variant in variants) {
    if (!"AIC" %in% names(all_results[[variant]])) next
    res <- all_results[[variant]][["AIC"]]
    f <- file.path(outdir, paste0("AIC_stage1_contrast_SEL_vs_MA__", variant, ".png"))
    grDevices::png(f, width = 1000, height = 850)
    graphics::plot(
      res$stage1_sel$actions$contrast,
      res$stage1_ma$actions$contrast,
      xlab = "AIC selection: stage-1 Q contrast (+1 minus -1)",
      ylab = "Akaike weighting: stage-1 Q contrast (+1 minus -1)",
      main = paste("Candidate-library sensitivity:", variant)
    )
    graphics::abline(h = 0, v = 0, lty = 3)
    graphics::abline(0, 1, lty = 2)
    grDevices::dev.off()
  }
}
