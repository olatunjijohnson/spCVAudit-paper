# Cross-validation scheme construction and a model-agnostic CV runner, used by
# estimate_optimism(). See docs/METHOD-EMULATOR.md section 1.

#' Spatial block cross-validation folds
#'
#' Assigns observations to `k` spatially contiguous folds by k-means clustering of
#' the coordinates. Spatial separation alone does not establish target-risk
#' validity. This is coordinate k-means: project longitude/latitude data first
#' when approximately distance-based compact blocks are required.
#'
#' @param data An `sf` object, numeric coordinate matrix, or `data.frame`.
#' @param k Number of spatial folds.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @return An integer fold-id vector of length `nrow(data)`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(120), y = runif(120), z = rnorm(120))
#' folds <- spatial_block_cv(d, k = 5, coords = c("x", "y"))
#' table(folds)
#' @export
spatial_block_cv <- function(data, k = 10L, coords = NULL) {
  xy <- .extract_coords(data, coords)
  k <- min(k, nrow(unique(xy$coords)))
  km <- suppressWarnings(stats::kmeans(xy$coords, centers = k, nstart = 5L, iter.max = 50L))
  as.integer(km$cluster)
}

# Buffered folds: from a base fold structure, drop training points within `buffer`
# of any test point (buffered/dead-zone CV; Ploton et al. 2020).
.buffer_folds <- function(coords, folds, buffer, info) {
  lapply(folds, function(f) {
    if (!length(f$test) || !length(f$train)) return(f)
    d <- .nn_dist(coords[f$train, , drop = FALSE], coords[f$test, , drop = FALSE], info)
    list(test = f$test, train = f$train[d > buffer])
  })
}

# Model-agnostic out-of-fold prediction. `predict_fun(train_df, test_df)` returns
# numeric predictions for the test rows.
.cv_predict <- function(df, folds, predict_fun) {
  pred <- rep(NA_real_, nrow(df))
  for (f in folds) {
    if (!length(f$test) || !length(f$train)) next
    pred[f$test] <- predict_fun(df[f$train, , drop = FALSE], df[f$test, , drop = FALSE])
  }
  pred
}

# Inverse-distance-weighting predictor factory. A score contrast for this learner
# does not establish a contrast for a different prediction algorithm.
.idw_predictor <- function(coord_cols, response, power = 2, info = NULL) {
  function(train, test) {
    tr <- as.matrix(train[, coord_cols, drop = FALSE])
    te <- as.matrix(test[, coord_cols, drop = FALSE])
    y <- train[[response]]
    geographic <- isTRUE(info$geographic)
    if (geographic) {
      q <- sf::st_as_sf(as.data.frame(te), coords = c(1L, 2L), crs = info$crs)
      r <- sf::st_as_sf(as.data.frame(tr), coords = c(1L, 2L), crs = info$crs)
    }
    vapply(seq_len(nrow(te)), function(i) {
      d <- if (geographic) as.numeric(sf::st_distance(q[i, ], r)) else
        sqrt(colSums((t(tr) - te[i, ])^2))
      hit <- which(d == 0)
      if (length(hit)) return(y[hit[1L]])
      w <- 1 / d^power
      sum(w * y) / sum(w)
    }, numeric(1))
  }
}

# Error metrics, including proper scoring rules for non-Gaussian responses. Each
# spec gives a per-point `loss` (so the fold/cluster bootstrap can resample it) and
# an `agg` aggregator. All are oriented so lower = better.
#   rmse/mae  : continuous;  brier/logloss : binary (p in [0,1]);  poisson : counts.
.metric_spec <- function(metric) {
  switch(
    metric,
    rmse    = list(loss = function(o, p) (o - p)^2, agg = function(x) sqrt(mean(x))),
    mae     = list(loss = function(o, p) abs(o - p), agg = mean),
    brier   = list(loss = function(o, p) (o - p)^2, agg = mean),
    logloss = list(loss = function(o, p) {
      p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -(o * log(p) + (1 - o) * log(1 - p))
    }, agg = mean),
    poisson = list(loss = function(o, p) {
      p <- pmax(p, 1e-9); 2 * (p - o + ifelse(o > 0, o * log(o / p), 0))
    }, agg = mean),
    stop(sprintf("Unknown metric '%s'.", metric)))
}

# Full metric as a single function(observed, predicted).
.metric_fun <- function(metric) {
  s <- .metric_spec(metric); function(o, p) s$agg(s$loss(o, p))
}
