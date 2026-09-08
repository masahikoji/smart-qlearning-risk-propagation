# Summaries for primary risk and generated-response diagnostics.

mc_summary <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(c(mean = NA_real_, se = NA_real_, sd = NA_real_, n = 0))
  s <- if (n > 1L) sd(x) else NA_real_
  c(mean = mean(x), se = s / sqrt(n), sd = s, n = n)
}

summarize_scenario_n <- function(scenario, n) {
  mat <- read_scenario_matrix(scenario$id, n)
  procs <- names(procedure_definitions())
  risk_rows <- list(); gap_rows <- list(); diag_rows <- list(); model_rows <- list()

  rr <- 0L; mm <- 0L
  for (r in procs) {
    for (s in seq_len(scenario$K)) {
      col <- paste0(r, ".s", s, ".risk")
      z <- mc_summary(mat[, col])
      rr <- rr + 1L
      risk_rows[[rr]] <- data.frame(
        scenario = scenario$id, K = scenario$K, n = n,
        procedure = r, stage = s, mean_risk = z["mean"],
        mc_se = z["se"], mc_sd = z["sd"], R_used = z["n"],
        stringsAsFactors = FALSE
      )
      wprefix <- paste0(r, ".s", s, ".w_")
      wcols <- colnames(mat)[startsWith(colnames(mat), wprefix)]
      for (wc in wcols) {
        wz <- mc_summary(mat[, wc])
        mm <- mm + 1L
        model_rows[[mm]] <- data.frame(
          scenario = scenario$id, K = scenario$K, n = n, procedure = r, stage = s,
          candidate = sub(wprefix, "", wc, fixed = TRUE), mean_weight = wz["mean"],
          mc_se_weight = wz["se"], R_used = wz["n"], stringsAsFactors = FALSE
        )
      }
    }
  }

  gg <- 0L
  for (s in seq_len(scenario$K)) {
    for (crit in c("AIC", "BIC")) {
      sel <- paste0(crit, "_SEL.s", s, ".risk")
      ma  <- paste0(crit, "_MA.s", s, ".risk")
      d <- mat[, sel] - mat[, ma]
      z <- mc_summary(d)
      gg <- gg + 1L
      th <- if (s == scenario$K) terminal_theory_n_gap(scenario, crit) else NA_real_
      delt <- if (s == scenario$K) terminal_delta_exact(scenario) else NA_real_
      gap_rows[[gg]] <- data.frame(
        scenario = scenario$id, K = scenario$K, n = n,
        criterion = crit, stage = s,
        mean_delta = z["mean"], mc_se = z["se"], mc_sd = z["sd"],
        n_times_delta = n * z["mean"],
        n_times_mc_se = n * z["se"],
        n32_times_delta = n^(3/2) * z["mean"],
        n32_times_mc_se = n^(3/2) * z["se"],
        delta_terminal = delt,
        theory_n_times_delta = th,
        theory_minus_mc = if (is.finite(th)) th - n * z["mean"] else NA_real_,
        R_used = z["n"], stringsAsFactors = FALSE
      )
    }
  }

  # Summarize diagnostic columns generically.
  dcols <- grep("(pseudo_rms|score_shift|score_shift_over_sqrt_n|wide_coef_shift|sqrt_n_wide_coef_shift|lambda_shift|selected_disagree_oracle|wide_weight_shift_oracle|weight_shift_oracle_[A-Za-z0-9_]+|narrow_q_rms_shift|wide_q_rms_shift|narrow_rss_shift|wide_rss_shift|stabilized|identity_error)$",
                colnames(mat), value = TRUE)
  if (length(dcols)) {
    diag_rows <- lapply(dcols, function(nm) {
      z <- mc_summary(mat[, nm])
      data.frame(scenario = scenario$id, K = scenario$K, n = n,
                 metric = nm, mean = z["mean"], mc_se = z["se"],
                 mc_sd = z["sd"], R_used = z["n"], stringsAsFactors = FALSE)
    })
  }

  list(risk = do.call(rbind, risk_rows), gap = do.call(rbind, gap_rows),
       model = if (length(model_rows)) do.call(rbind, model_rows) else data.frame(),
       diagnostics = if (length(diag_rows)) do.call(rbind, diag_rows) else data.frame())
}

summarize_all <- function(K = NULL) {
  scens <- if (is.null(K)) c(smart_scenarios(2L), smart_scenarios(3L)) else smart_scenarios(K)
  risk <- list(); gap <- list(); model <- list(); diag <- list(); ir <- ig <- im <- id <- 0L
  for (s in scens) {
    base <- file.path(smart_results_root(), "raw", s$id)
    if (!dir.exists(base)) next
    ndirs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
    ns <- suppressWarnings(as.integer(sub("^n", "", ndirs)))
    ns <- ns[is.finite(ns)]
    for (n in sort(ns)) {
      z <- summarize_scenario_n(s, n)
      ir <- ir + 1L; risk[[ir]] <- z$risk
      ig <- ig + 1L; gap[[ig]] <- z$gap
      if (nrow(z$model)) { im <- im + 1L; model[[im]] <- z$model }
      if (nrow(z$diagnostics)) { id <- id + 1L; diag[[id]] <- z$diagnostics }
    }
  }
  outdir <- smart_dir("results", "summary")
  riskdf <- if (length(risk)) do.call(rbind, risk) else data.frame()
  gapdf <- if (length(gap)) do.call(rbind, gap) else data.frame()
  modeldf <- if (length(model)) do.call(rbind, model) else data.frame()
  diagdf <- if (length(diag)) do.call(rbind, diag) else data.frame()
  write.csv(riskdf, file.path(outdir, "risk_summary.csv"), row.names = FALSE)
  write.csv(gapdf, file.path(outdir, "risk_gap_summary.csv"), row.names = FALSE)
  write.csv(modeldf, file.path(outdir, "model_weight_summary.csv"), row.names = FALSE)
  write.csv(diagdf, file.path(outdir, "diagnostic_summary.csv"), row.names = FALSE)
  invisible(list(risk = riskdf, gap = gapdf, model = modeldf, diagnostics = diagdf))
}

plot_scaled_gaps <- function(summary_obj = summarize_all()) {
  d <- summary_obj$gap
  if (!nrow(d)) return(invisible(NULL))
  f <- file.path(smart_dir("results", "figures"), "scaled_risk_gaps.pdf")
  pdf(f, width = 8.5, height = 6.5)
  on.exit(dev.off(), add = TRUE)
  scenarios <- unique(d$scenario)
  for (sid in scenarios) {
    z <- d[d$scenario == sid, , drop = FALSE]
    stages <- sort(unique(z$stage))
    for (s in stages) {
      zz <- z[z$stage == s, , drop = FALSE]
      lo <- zz$n_times_delta - 2 * zz$n_times_mc_se
      hi <- zz$n_times_delta + 2 * zz$n_times_mc_se
      ylim <- range(c(lo, hi), finite = TRUE)
      plot(range(zz$n), ylim, type = "n", xlab = "n", ylab = expression(n * Delta[t]),
           main = paste(sid, "stage", s))
      abline(h = 0, lty = 2)
      thv <- unique(zz$theory_n_times_delta[is.finite(zz$theory_n_times_delta)])
      if (length(thv) == 1L) abline(h = thv, lty = 3)
      for (crit in unique(zz$criterion)) {
        q <- zz[zz$criterion == crit, , drop = FALSE]
        o <- order(q$n)
        xx <- q$n[o]; yy <- q$n_times_delta[o]; ee <- 2 * q$n_times_mc_se[o]
        lines(xx, yy, type = "b", pch = if (crit == "AIC") 1 else 2)
        segments(xx, yy - ee, xx, yy + ee)
      }
      legend("topright", legend = unique(zz$criterion), pch = c(1, 2)[seq_along(unique(zz$criterion))],
             lty = 1, bty = "n")
    }
  }
  invisible(f)
}

summarize_secondary_all <- function(K = NULL) {
  scens <- if (is.null(K)) c(smart_scenarios(2L), smart_scenarios(3L)) else smart_scenarios(K)
  rows <- list(); k <- 0L
  for (s in scens) {
    base <- file.path(smart_results_root(), "raw", s$id)
    if (!dir.exists(base)) next
    ndirs <- list.dirs(base, full.names = TRUE, recursive = FALSE)
    for (d in ndirs) {
      n <- suppressWarnings(as.integer(sub("^n", "", basename(d))))
      f <- file.path(d, "secondary_metrics.rds")
      if (!is.finite(n) || !file.exists(f)) next
      z <- readRDS(f)
      mats <- list(z$switching)
      if (!is.null(z$regret)) mats <- c(mats, list(z$regret))
      M <- do.call(cbind, mats)
      rr_meta <- z$regret_meta
      if (is.null(rr_meta)) {
        rr_meta <- list(N_value = if (!is.null(z$regret)) attr(z$regret, "N_value") else NA_integer_,
                        regret_method = if (!is.null(z$regret)) attr(z$regret, "regret_method") else NA_character_,
                        evaluation_design = if (!is.null(z$regret)) attr(z$regret, "evaluation_design") else NA_character_)
      }
      nv_meta <- rr_meta$N_value; if (is.null(nv_meta) || !length(nv_meta)) nv_meta <- NA_integer_
      method_meta <- rr_meta$regret_method; if (is.null(method_meta) || !length(method_meta)) method_meta <- NA_character_
      design_meta <- rr_meta$evaluation_design; if (is.null(design_meta) || !length(design_meta)) design_meta <- NA_character_
      for (nm in colnames(M)) {
        q <- mc_summary(M[, nm])
        is_regret <- grepl("\\.regime_regret$", nm)
        k <- k + 1L
        rows[[k]] <- data.frame(scenario=s$id, K=s$K, n=n, metric=nm,
                                mean=q["mean"], mc_se=q["se"], mc_sd=q["sd"],
                                R_used=q["n"],
                                N_value=if (is_regret) nv_meta else NA_integer_,
                                regret_method=if (is_regret) method_meta else NA_character_,
                                evaluation_design=if (is_regret) design_meta else NA_character_,
                                stringsAsFactors=FALSE)
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write.csv(out, file.path(smart_dir("results", "summary"), "secondary_summary.csv"), row.names=FALSE)
  invisible(out)
}

summarize_decomposition_all <- function(K = NULL) {
  scens <- if (is.null(K)) c(smart_scenarios(2L), smart_scenarios(3L)) else smart_scenarios(K)
  rows <- list(); k <- 0L
  for (s in scens) {
    base <- file.path(smart_results_root(), "raw", s$id)
    if (!dir.exists(base)) next
    ndirs <- list.dirs(base, full.names = TRUE, recursive = FALSE)
    for (d in ndirs) {
      n <- suppressWarnings(as.integer(sub("^n", "", basename(d))))
      f <- file.path(d, "decomposition_metrics.rds")
      if (!is.finite(n) || !file.exists(f)) next
      z <- readRDS(f)
      for (entry in z$entries) {
        M <- entry$metrics
        for (nm in colnames(M)) {
          q <- mc_summary(M[, nm])
          k <- k + 1L
          rows[[k]] <- data.frame(scenario=s$id, K=s$K, n=n,
                                  target_stage=entry$target_stage,
                                  procedure=entry$procedure, metric=nm,
                                  mean=q["mean"], mc_se=q["se"], mc_sd=q["sd"],
                                  R_used=q["n"], N_decomp=z$N_decomp,
                                  N_pairs=if (!is.null(z$N_pairs)) z$N_pairs else NA_integer_,
                                  cores_decomp=if (!is.null(z$cores_decomp)) z$cores_decomp else NA_integer_,
                                  decomposition_method=if (!is.null(z$method)) z$method else "legacy",
                                  quadrature_order=if (!is.null(z$quadrature_order)) z$quadrature_order else NA_integer_,
                                  stringsAsFactors=FALSE)
        }
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write.csv(out, file.path(smart_dir("results", "summary"), "decomposition_summary.csv"), row.names=FALSE)
  invisible(out)
}
