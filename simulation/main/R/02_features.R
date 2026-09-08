# Candidate-model definitions and feature matrices.

past_feature_names <- function(stage) {
  if (stage <= 1L) character(0) else paste0("pastA", seq_len(stage - 1L))
}

# X^2 is common in the accumulation study so transported quadratic prognostic
# components are absorbed upstream. X^3 is included only for the smooth-boundary
# construction (xi_X=0), where maximization generates a cubic prognostic term.
standard_common_features <- function(scenario, stage) {
  out <- c("int", "X", "X_sq")
  if (identical(scenario$boundary, "smooth")) out <- c(out, "X_cu")
  c(out, past_feature_names(stage))
}

standard_union_features <- function(scenario, stage) {
  c(standard_common_features(scenario, stage), "A", "AX", "AX2")
}

standard_candidate_features <- function(scenario, stage, candidate_count = 2L) {
  common <- standard_common_features(scenario, stage)
  if (candidate_count == 2L) {
    list(
      narrow = c(common, "A", "AX"),
      wide   = c(common, "A", "AX", "AX2")
    )
  } else if (candidate_count == 4L) {
    list(
      M0 = c(common, "A"),
      M1 = c(common, "A", "AX"),
      M2 = c(common, "A", "AX2"),
      M3 = c(common, "A", "AX", "AX2")
    )
  } else stop("candidate_count must be 2 or 4")
}

criterion_union_features <- function(stage, variant = c("aligned", "control")) {
  variant <- match.arg(variant)
  wide_only <- if (variant == "aligned") "X_sq" else "AX2"
  unique(c("int", "X", past_feature_names(stage), "A", "AX", wide_only))
}

criterion_candidate_features <- function(stage, variant = c("aligned", "control")) {
  variant <- match.arg(variant)
  common <- c("int", "X", past_feature_names(stage), "A", "AX")
  wide_only <- if (variant == "aligned") "X_sq" else "AX2"
  list(narrow = common, wide = c(common, wide_only))
}

stage_model_spec <- function(scenario, stage) {
  if (scenario$family == "criterion" && stage < scenario$K) {
    cand <- criterion_candidate_features(stage, scenario$criterion_variant)
    union <- criterion_union_features(stage, scenario$criterion_variant)
    added <- if (scenario$criterion_variant == "aligned") "X_sq" else "AX2"
  } else {
    cand <- standard_candidate_features(scenario, stage, scenario$candidate_count)
    union <- standard_union_features(scenario, stage)
    added <- if (scenario$candidate_count == 2L) "AX2" else NA_character_
  }
  list(union = union, candidates = cand, added_feature = added)
}

# data must contain X1,...,XK and past A columns; action may be scalar or vector.
feature_matrix <- function(data, stage, action, features) {
  x <- data[[paste0("X", stage)]]
  n <- length(x)
  if (length(action) == 1L) action <- rep(action, n)
  if (length(action) != n) stop("action length mismatch")
  ans <- matrix(NA_real_, nrow = n, ncol = length(features),
                dimnames = list(NULL, features))
  for (j in seq_along(features)) {
    f <- features[j]
    ans[, j] <- switch(
      f,
      int = 1,
      X = x,
      X_sq = x^2,
      X_cu = x^3,
      A = action,
      AX = action * x,
      AX2 = action * x^2,
      {
        if (grepl("^pastA[0-9]+$", f)) {
          k <- as.integer(sub("pastA", "", f, fixed = TRUE))
          data[[paste0("A", k)]]
        } else stop("Unknown feature: ", f)
      }
    )
  }
  ans
}

embed_coef <- function(beta, candidate_features, union_features) {
  out <- setNames(rep(0, length(union_features)), union_features)
  out[candidate_features] <- beta
  out
}

predict_embedded_q <- function(beta_union, data, stage, action) {
  X <- feature_matrix(data, stage, action, names(beta_union))
  as.numeric(X %*% beta_union)
}

contrast_matrix <- function(data, stage, union_features) {
  xp <- feature_matrix(data, stage, 1, union_features)
  xm <- feature_matrix(data, stage, -1, union_features)
  xp - xm
}
