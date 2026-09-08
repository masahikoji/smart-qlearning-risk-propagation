# Chunked, restartable Monte Carlo runner.

scenario_index_global <- function(id) {
  all <- c(smart_scenarios(2L), smart_scenarios(3L))
  ids <- vapply(all, `[[`, character(1), "id")
  match(id, ids)
}

scenario_result_dir <- function(scenario, n) {
  smart_dir("results", "raw", scenario$id, paste0("n", n))
}

run_one_chunk <- function(scenario, n, rep_ids, control, eval_cache) {
  f <- function(i) {
    seed <- replicate_seed(control$seed, scenario$id, n, i)
    z <- simulate_one_replicate(scenario, n, seed, eval_cache)
    c(rep_id = i, rep_seed = seed, z)
  }
  if (control$cores > 1L && .Platform$OS.type != "windows") {
    ans <- parallel::mclapply(rep_ids, f, mc.cores = control$cores,
                             mc.preschedule = TRUE, mc.set.seed = FALSE)
  } else {
    ans <- lapply(rep_ids, f)
  }
  mat <- do.call(rbind, ans)
  storage.mode(mat) <- "double"
  mat
}

run_scenario_n <- function(scenario, n, control, overwrite = FALSE) {
  if (is.na(scenario_index_global(scenario$id))) stop("Scenario not registered: ", scenario$id)
  outdir <- scenario_result_dir(scenario, n)
  metafile <- file.path(outdir, "metadata.rds")
  chunk_files <- list.files(outdir, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = TRUE)

  if (overwrite) {
    if (length(chunk_files)) unlink(chunk_files, force = TRUE)
    unlink(file.path(outdir, c("secondary_metrics.rds", "decomposition_metrics.rds")), force = TRUE)
    chunk_files <- character(0)
  } else if (length(chunk_files)) {
    if (!file.exists(metafile)) {
      stop("Existing result chunks have no metadata.rds in ", outdir,
           ". Move/delete them or rerun with overwrite=TRUE.")
    }
    old <- readRDS(metafile)
    if (!identical(scenario_result_signature(old$scenario), scenario_result_signature(scenario))) {
      stop("Scenario definition changed since existing chunks were created for ",
           scenario$id, ", n=", n, ". Use a new scenario id or overwrite=TRUE.")
    }
    if (!identical(control_result_signature(old$control), control_result_signature(control))) {
      stop("Result-defining run control changed since existing chunks were created for ",
           scenario$id, ", n=", n,
           ". R and worker cores may change on restart; N_eval, chunk size, quadrature order and base seed may not. Use overwrite=TRUE otherwise.")
    }
    now_md5 <- smart_engine_checksums()
    if (!is.null(old$engine_checksums) && !identical(old$engine_checksums, now_md5)) {
      stop("Simulation-engine source files changed since existing chunks were created for ",
           scenario$id, ", n=", n, ". Use overwrite=TRUE to avoid mixing code versions.")
    }
  }

  quad <- gauss_legendre(control$quadrature_order)
  validate_scenario(scenario, n, control, quad = quad, verbose = TRUE)
  log_msg("Building evaluation cache for ", scenario$id, ", n=", n,
          ", N_eval=", control$N_eval)
  eval_cache <- build_evaluation_cache(scenario, n, control)

  # Metadata is rewritten on every compatible restart so requested_R reflects
  # the current target. This permits a pilot R=1000 to be extended to R=20000.
  meta <- list(scenario = scenario, n = n, control = control,
               requested_R = control$R,
               evaluation_seed = eval_cache$seed,
               engine_checksums = smart_engine_checksums(),
               R_version = R.version.string,
               platform = R.version$platform,
               timestamp = as.character(Sys.time()))
  saveRDS(meta, metafile)

  starts <- seq.int(1L, control$R, by = control$chunk_size)
  for (st in starts) {
    en <- min(control$R, st + control$chunk_size - 1L)
    outfile <- file.path(outdir, sprintf("chunk_%06d_%06d.rds", st, en))
    if (file.exists(outfile) && !overwrite) {
      log_msg("Skipping existing ", basename(outfile))
      next
    }
    # If a shorter last pilot chunk overlaps the requested block, do not silently
    # create duplicate replicate ids. Such a case requires a clean overwrite.
    existing <- list.files(outdir, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = FALSE)
    if (length(existing) && !file.exists(outfile)) {
      rg <- strcapture("^chunk_([0-9]+)_([0-9]+)\\.rds$", existing,
                       proto = list(lo=integer(), hi=integer()))
      if (nrow(rg) && any(rg[,"lo"] <= en & rg[,"hi"] >= st)) {
        stop("Existing chunk boundaries overlap requested block ", st, "-", en,
             " for ", scenario$id, ". Keep the same SMART_CHUNK or rerun with overwrite=TRUE.")
      }
    }
    log_msg("Running ", scenario$id, " n=", n, " reps ", st, "-", en)
    mat <- run_one_chunk(scenario, n, st:en, control, eval_cache)
    tmp <- paste0(outfile, ".tmp")
    saveRDS(mat, tmp, compress = "gzip")
    if (file.exists(outfile)) file.remove(outfile)
    file.rename(tmp, outfile)
  }
  invisible(TRUE)
}

run_scenario <- function(scenario, control, overwrite = FALSE) {
  for (n in control$n_values) run_scenario_n(scenario, as.integer(n), control, overwrite)
  invisible(TRUE)
}

run_scenario_set <- function(K, profile = "paper", ids = NULL, overwrite = FALSE) {
  control <- smart_run_control(profile)
  scens <- smart_scenarios(K)
  if (!is.null(ids)) {
    keep <- vapply(scens, function(s) s$id %in% ids, logical(1))
    scens <- scens[keep]
    missing <- setdiff(ids, vapply(scens, `[[`, character(1), "id"))
    if (length(missing)) stop("Unknown scenario ids: ", paste(missing, collapse = ", "))
  }
  log_msg("Run profile=", profile, ", K=", K, ", scenarios=", length(scens),
          ", R=", control$R, ", N_eval=", control$N_eval,
          ", cores=", control$cores, ", results=", smart_results_root())
  for (s in scens) run_scenario(s, control, overwrite)
  invisible(TRUE)
}

read_scenario_matrix <- function(scenario_id, n) {
  d <- smart_dir("results", "raw", scenario_id, paste0("n", n))
  fs <- sort(list.files(d, pattern = "^chunk_[0-9]+_[0-9]+\\.rds$", full.names = TRUE))
  if (!length(fs)) stop("No chunks found in ", d)
  mats <- lapply(fs, readRDS)
  out <- do.call(rbind, mats)
  if (anyDuplicated(out[, "rep_id"])) {
    stop("Duplicate replicate ids detected in ", d,
         ". This usually means incompatible chunk files were mixed.")
  }
  out <- out[order(out[, "rep_id"]), , drop = FALSE]
  mf <- file.path(d, "metadata.rds")
  if (file.exists(mf)) {
    meta <- readRDS(mf)
    target_R <- if (!is.null(meta$requested_R)) meta$requested_R else meta$control$R
    if (is.finite(target_R)) out <- out[out[, "rep_id"] <= target_R, , drop = FALSE]
  }
  out
}
