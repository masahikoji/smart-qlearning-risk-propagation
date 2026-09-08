# General utilities: paths, quadrature, stable softmax, deterministic seeds.

smart_project_root <- function() {
  root <- Sys.getenv("SMART_PROGRAM_DIR", unset = "")
  if (nzchar(root)) return(normalizePath(root, mustWork = FALSE))
  normalizePath(getwd(), mustWork = FALSE)
}

# Revised results are written separately from the first development version so
# that old pilot chunks can never be mixed with the corrected simulation engine.
smart_results_root <- function() {
  p <- Sys.getenv("SMART_RESULTS_DIR", unset = "")
  if (!nzchar(p)) p <- file.path(smart_project_root(), "results_v2")
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  normalizePath(p, mustWork = FALSE)
}

smart_dir <- function(...) {
  parts <- list(...)
  if (length(parts) && identical(as.character(parts[[1L]]), "results")) {
    tailparts <- if (length(parts) > 1L) parts[-1L] else list()
    p <- do.call(file.path, c(list(smart_results_root()), tailparts))
  } else {
    p <- do.call(file.path, c(list(smart_project_root()), parts))
  }
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

stable_softmax <- function(logw) {
  z <- logw - max(logw)
  ez <- exp(z)
  ez / sum(ez)
}

# n-point Gauss-Legendre nodes and weights on [-1,1]. No external package.
gauss_legendre <- function(n = 32L) {
  n <- as.integer(n)
  if (n < 2L) stop("quadrature order must be >= 2")
  i <- seq_len(n - 1L)
  beta <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- beta
  J[cbind(i + 1L, i)] <- beta
  eg <- eigen(J, symmetric = TRUE)
  ord <- order(eg$values)
  nodes <- eg$values[ord]
  V <- eg$vectors[, ord, drop = FALSE]
  weights <- 2 * V[1L, ]^2
  list(nodes = nodes, weights = weights, prob_weights = weights / 2)
}

# Stable scenario-id hash. Seeds therefore do not change when another scenario
# is inserted earlier in the configuration file.
scenario_hash <- function(id) {
  z <- utf8ToInt(enc2utf8(as.character(id)))
  h <- 104729
  if (length(z)) {
    for (v in z) h <- (h * 1009 + v + 97) %% 2147483000
  }
  as.double(h)
}

replicate_seed <- function(base_seed, scenario_id, n, rep_id) {
  z <- as.double(base_seed) + 1000003 * scenario_hash(scenario_id) +
    1009 * as.integer(n) + as.integer(rep_id)
  as.integer(z %% 2147483000)
}

evaluation_seed <- function(base_seed, scenario_id, n) {
  as.integer((as.double(base_seed) + 700000003 +
                500009 * scenario_hash(scenario_id) + 1009 * as.integer(n)) %% 2147483000)
}

safe_mean <- function(x) if (length(x)) mean(x) else NA_real_

prefixed <- function(x, prefix) {
  x <- as.numeric(x)
  names(x) <- paste0(prefix, names(x))
  x
}

past_action_grid <- function(s) {
  if (s <= 1L) return(matrix(numeric(0), nrow = 1L, ncol = 0L))
  g <- expand.grid(rep(list(c(-1, 1)), s - 1L), KEEP.OUT.ATTRS = FALSE,
                   stringsAsFactors = FALSE)
  as.matrix(g)
}

scenario_result_signature <- function(scenario) {
  z <- unclass(scenario)
  z$notes <- NULL
  z
}

# R is deliberately excluded: a pilot with R=1000 may be extended to R=20000
# or R=100000 without invalidating the first replicates. Cores are also excluded.
control_result_signature <- function(control) {
  keep <- c("N_eval", "chunk_size", "quadrature_order", "seed")
  control[intersect(keep, names(control))]
}

# Configuration is excluded from the engine checksum. Scenario numerical values
# are protected by scenario_result_signature(), and run-control values by
# control_result_signature(). This prevents an unrelated scenario edit from
# invalidating every existing result set.
smart_engine_checksums <- function() {
  fs <- file.path(smart_project_root(), "R", c(
    "01_utils.R", "02_features.R", "03_fit.R", "04_truth_dgp.R",
    "05_validation.R", "06_eval_cache.R", "07_recursive_fit.R", "08_runner.R"
  ))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) return(character(0))
  setNames(unname(tools::md5sum(fs)), basename(fs))
}
