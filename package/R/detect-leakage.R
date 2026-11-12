# Core diagnostic: detect_leakage(). Computes the observed (test->train) and target
# (prediction->sample) nearest-neighbour distances, then the signed SLI in both the
# distance form (SLI_d) and the dependence form (SLI_rho). See docs/METHOD-SLI.md.

# Parse a split specification into a list of folds, each list(test=, train=).
# Accepts: a list with $test (+ optional $train); a list of test-index vectors
# (k folds, train = complement); or an integer/factor fold-id vector of length n.
.parse_split <- function(split, n) {
  # Adaptor: tidymodels rsample rset (e.g. vfold_cv, spatial_block_cv) -> folds.
  if (inherits(split, "rset") && requireNamespace("rsample", quietly = TRUE)) {
    return(lapply(split$splits, function(s) list(
      test = as.integer(s, data = "assessment"),
      train = as.integer(s, data = "analysis"))))
  }
  if (is.list(split) && !is.null(split$test)) {
    test <- split$test
    train <- if (!is.null(split$train)) split$train else setdiff(seq_len(n), test)
    return(list(list(test = test, train = train)))
  }
  # A list of pre-built fold objects, each list(test=, train=) (e.g. from
  # temporal_kfold() or a workflow adaptor).
  if (is.list(split) && length(split) && is.list(split[[1L]]) && !is.null(split[[1L]]$test)) {
    return(lapply(split, function(f) list(
      test = f$test,
      train = if (!is.null(f$train)) f$train else setdiff(seq_len(n), f$test))))
  }
  if (is.list(split)) {
    folds <- lapply(split, function(test) {
      list(test = test, train = setdiff(seq_len(n), test))
    })
    return(folds)
  }
  if (length(split) == n) {
    ids <- unique(split)
    folds <- lapply(ids, function(k) {
      test <- which(split == k)
      list(test = test, train = setdiff(seq_len(n), test))
    })
    return(folds)
  }
  stop("Unrecognised `split`: supply a list(test=, train=), a list of test-index ",
       "vectors, or a fold-id vector of length nrow(data).")
}

# Signed areas between the observed and target NN-distance ECDFs. Uses exact
# step-function (left-Riemann) integration over the union of jump points, so the
# signed area A equals mean(fpred) - mean(gobs) exactly (METHOD-SLI section 1).
.sli_areas <- function(gobs, fpred) {
  rgrid <- sort(unique(c(0, gobs, fpred)))
  Gobs <- stats::ecdf(gobs)(rgrid)
  Gpred <- stats::ecdf(fpred)(rgrid)
  m <- length(rgrid)
  dr <- diff(rgrid)
  dleft <- (Gobs - Gpred)[-m]    # constant value of (Gobs - Gpred) on each interval
  Aplus <- sum(dr * pmax(dleft, 0))
  Aminus <- sum(dr * pmax(-dleft, 0))
  A <- Aplus - Aminus            # signed area = mean(fpred) - mean(gobs)
  W <- Aplus + Aminus            # unsigned = NNDM/kNNDM Wasserstein-W
  list(A = A, Aplus = Aplus, Aminus = Aminus, W = W,
       delta = if (W > 0) A / W else 0,
       crosses = Aplus > 1e-12 * max(W, .Machine$double.eps) &&
                 Aminus > 1e-12 * max(W, .Machine$double.eps),
       rgrid = rgrid, Gobs = Gobs, Gpred = Gpred)
}

#' Audit spatial distance mismatch between validation and a prediction target
#'
#' Audits an existing split or cross-validation scheme: computes the Spatial
#' Leakage Index (signed) in a distance form (`SLI_d`) and a dependence-weighted
#' form (`SLI_rho`), with per-point contributions. Positive `SLI_rho` means that
#' validation observations have larger mean training-neighbour correlation under
#' the supplied dependence function than target observations. This is a diagnostic
#' of mismatch, not a general estimate of predictive bias. Inspect the full distance
#' distributions and `W`: a signed value near zero can hide cancellation.
#' For `k = 1`, `SLI_rho2` reports the difference in mean squared correlations.
#' Multiplied by the observation variance, this equals the expected MSE gap for
#' the known-mean, known-covariance Gaussian predictor using only its nearest
#' training observation. It is not an exact expression for full kriging or a
#' general machine-learning model. `W_rho` is the unsigned Wasserstein distance
#' between the retained-correlation distributions.
#'
#' @param data An `sf` object, numeric coordinate matrix, or `data.frame`.
#' @param split A split specification: a `list(test=, train=)`, a list of
#'   test-index vectors (k folds), or a fold-id vector of length `nrow(data)`.
#' @param target A [prediction_target()] (or coordinates coercible to one).
#' @param response Name of the numeric response column, used to estimate spatial
#'   dependence (required unless `dependence` is supplied).
#' @param dependence A precomputed [estimate_dependence()] object (optional;
#'   overrides `response`). Carries any anisotropy used.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param n_boot Number of Monte-Carlo draws from the approximate covariance of
#'   the variogram parameters. `0` (default) skips it. The resulting 90% parameter
#'   propagation interval conditions on the geometry and fitted covariance family;
#'   it is not a confidence interval for predictive optimism. If parameter
#'   covariance is unavailable, no interval is returned. Supported for `k = 1`.
#' @param k Number of nearest training neighbours for the dependence-form `SLI_rho`
#'   (noisy-OR combination). `1` (default) is the single-nearest-neighbour index.
#'   For `k > 1` the combination is heuristic and has no kriging-variance identity.
#' @return An object of class `leakage_diagnosis`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(80), y = runif(80), z = rnorm(80))
#' folds <- sample(rep_len(1:5, 80))
#' grid <- as.matrix(expand.grid(x = seq(0, 1, 0.2), y = seq(0, 1, 0.2)))
#' tgt <- prediction_target(grid = grid, type = "grid")
#' detect_leakage(d, folds, tgt, response = "z", coords = c("x", "y"))
#' @export
detect_leakage <- function(data, split, target, response = NULL,
                           dependence = NULL, coords = NULL, n_boot = 0L, k = 1L) {
  xy <- .extract_coords(data, coords)
  n <- nrow(xy$coords)
  if (length(k) != 1L || !is.finite(k) || k < 1 || k != as.integer(k))
    stop("`k` must be a positive integer.")
  if (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 0 || n_boot != as.integer(n_boot))
    stop("`n_boot` must be a nonnegative integer.")
  if (k > 1L && n_boot > 0L)
    stop("Parameter propagation is supported only for `k = 1`.")
  info <- list(crs = xy$crs, geographic = xy$geographic)

  tgt <- .as_prediction_target(target, data = data, coords = coords)
  has_crs <- function(z) inherits(z, "crs") && !is.na(z)
  if (has_crs(xy$crs) && has_crs(tgt$crs) && xy$crs != tgt$crs)
    stop("Training and target CRS differ. Transform target locations to the training CRS before auditing.")
  folds <- .parse_split(split, n)
  valid_index <- function(z) is.numeric(z) && length(z) > 0L &&
    all(is.finite(z)) && all(z == as.integer(z)) && all(z >= 1L & z <= n) &&
    !anyDuplicated(z)
  for (f in folds) {
    if (!valid_index(f$test) || !valid_index(f$train))
      stop("Each fold must contain nonempty, unique, valid training and test row indices.")
    if (length(intersect(f$train, f$test)))
      stop("Training and test rows must not overlap within a fold.")
  }
  if (anyDuplicated(unlist(lapply(folds, `[[`, "test"))))
    stop("A row is tested more than once. Audit each cross-validation repeat separately.")

  # Dependence model (for SLI_rho and the phi normalisation of SLI_d).
  if (is.null(dependence)) {
    if (is.null(response)) {
      stop("Supply either `response` or a precomputed `dependence` object.")
    }
    dependence <- estimate_dependence(data, response = response, coords = coords)
  }
  rho <- dependence$rho
  phi <- dependence$practical_range

  # Measure distances in the dependence metric: under (geometric) anisotropy this
  # is Euclidean distance in the transformed space; under isotropy it is unchanged.
  aniso <- dependence$anisotropy
  sxy <- .apply_aniso(xy$coords, aniso)
  txy <- .apply_aniso(tgt$coords, aniso)

  # Observed test->train nearest-neighbour distance (scalar, for SLI_d/ECDF) and the
  # k-neighbour retained correlation (for SLI_rho), tracked back to original rows.
  gobs <- rep(NA_real_, n); robs <- rep(NA_real_, n)
  fold_id <- rep(NA_integer_, n)
  for (fk in seq_along(folds)) {
    f <- folds[[fk]]
    if (length(f$test) == 0L || length(f$train) == 0L) next
    qy <- sxy[f$test, , drop = FALSE]; rf <- sxy[f$train, , drop = FALSE]
    gobs[f$test] <- .nn_dist(qy, rf, info)
    robs[f$test] <- .retained_corr(qy, rf, info, rho, k)
    fold_id[f$test] <- fk
  }
  tested <- which(!is.na(gobs))
  if (length(tested) == 0L) stop("No test points found in `split`.")

  # Target prediction->sample NN distances.
  fpred <- .nn_dist(txy, sxy, info)

  # Distance-form SLI (single nearest neighbour).
  ar <- .sli_areas(gobs[tested], fpred)
  SLI_d <- ar$A / phi

  # Dependence-form SLI and per-point leakage map (k-neighbour retained correlation).
  rpred <- .retained_corr(txy, sxy, info, rho, k)
  c_pred <- mean(rpred)
  c_obs <- mean(robs[tested])
  SLI_rho <- c_obs - c_pred
  SLI_rho2 <- if (k == 1L) mean(robs[tested]^2) - mean(rpred^2) else NA_real_
  W_rho <- .sli_areas(robs[tested], rpred)$W
  leak_point <- robs - c_pred                # length n; NA for untested rows
  fold_sli <- tapply(leak_point[tested], fold_id[tested], mean)

  # Monte-Carlo uncertainty: geometry fixed, resample variogram parameters.
  SLI_rho_ci <- NULL; SLI_d_ci <- NULL; boot <- NULL
  interval_method <- "none"
  if (n_boot > 0L) {
    pars <- .sample_params(dependence, n_boot)
    if (is.null(pars)) {
      warning("Variogram parameter covariance is unavailable or invalid; no interval returned.",
              call. = FALSE)
    } else {
    g <- gobs[tested]; f <- fpred
    sl_rho <- numeric(n_boot); sl_d <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      sp <- pars$signal_prop[b]; rg <- pars$range[b]
      rb <- function(h) sp * exp(-h / rg)
      sl_rho[b] <- mean(rb(g)) - mean(rb(f))
      sl_d[b] <- ar$A / (3 * rg)
    }
    SLI_rho_ci <- unname(stats::quantile(sl_rho, c(0.05, 0.95)))
    SLI_d_ci <- unname(stats::quantile(sl_d, c(0.05, 0.95)))
    boot <- list(SLI_rho = sl_rho, SLI_d = sl_d)
    interval_method <- "approximate variogram parameter propagation (geometry and family fixed)"
    }
  }

  structure(
    list(
      SLI_rho = SLI_rho, SLI_rho2 = SLI_rho2, SLI_d = SLI_d, W_rho = W_rho,
      SLI_rho_ci = SLI_rho_ci, SLI_d_ci = SLI_d_ci, boot = boot,
      interval_method = interval_method,
      c_obs = c_obs, c_pred = c_pred,
      A = ar$A, Aplus = ar$Aplus, Aminus = ar$Aminus, W = ar$W, delta = ar$delta,
      crosses = ar$crosses, k = k,
      phi = phi,
      gobs = gobs, fpred = fpred, leak_point = leak_point,
      fold_id = fold_id, fold_sli = fold_sli, tested = tested,
      ecdf = list(r = ar$rgrid, Gobs = ar$Gobs, Gpred = ar$Gpred),
      coords = xy$coords, dependence = dependence, target_type = tgt$type,
      anisotropic = !is.null(aniso),
      n = n, n_test = length(tested), n_folds = length(folds)
    ),
    class = "leakage_diagnosis")
}
