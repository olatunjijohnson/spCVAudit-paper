# Legacy API: compare CV schemes. Historical field names are retained, but neither
# the comparison score nor a rescaled published score is an identified target risk.

# Per-point loss contributions per fold (so the cluster bootstrap can resample them).
.fold_errors <- function(df, folds, predict_fun, response, metric) {
  y <- df[[response]]; loss <- .metric_spec(metric)$loss
  out <- lapply(folds, function(f) {
    if (!length(f$test) || !length(f$train)) return(numeric(0))
    loss(y[f$test], predict_fun(df[f$train, , drop = FALSE], df[f$test, , drop = FALSE]))
  })
  out[lengths(out) > 0L]
}

.agg_err <- function(err_list, metric) {
  .metric_spec(metric)$agg(unlist(err_list, use.names = FALSE))
}

# Cluster bootstrap over folds: resample whole folds with replacement.
.boot_metric <- function(err_list, metric, B) {
  K <- length(err_list)
  vapply(seq_len(B), function(b) .agg_err(err_list[sample.int(K, K, replace = TRUE)], metric),
         numeric(1))
}

# Target-matched control: among candidate schemes (random, block, and block buffered
# by increasing distances), pick the smallest unsigned distance mismatch. This is
# a geometric selection criterion; it does not establish unbiased target risk.
.matched_control <- function(data, target, response, coords, k, xy, info, dependence) {
  n <- nrow(xy$coords)
  if (is.null(dependence)) dependence <- estimate_dependence(data, response = response, coords = coords)
  tgt <- .as_prediction_target(target, data = data, coords = coords)
  diam <- .domain_diameter(xy$coords, info)
  base_block <- .parse_split(spatial_block_cv(xy$coords, k), n)
  cand <- list(random = .parse_split(sample(rep_len(seq_len(k), n)), n), block = base_block)
  for (frac in c(0.05, 0.1, 0.15, 0.2, 0.25)) {
    cand[[sprintf("buffer%.2f", frac)]] <- .buffer_folds(xy$coords, base_block, frac * diam, info)
  }
  diagnostics <- lapply(cand, function(f) tryCatch(
    detect_leakage(data, f, tgt, dependence = dependence, coords = coords),
    error = function(e) NULL))
  W <- vapply(diagnostics, function(z) if (is.null(z)) Inf else z$W, numeric(1))
  ok <- which(is.finite(W))
  if (!length(ok)) stop("No candidate scheme has complete, valid folds.")
  best <- ok[which.min(W[ok])]
  list(folds = cand[[best]], sli = diagnostics[[best]]$SLI_rho,
       W = W[best], label = names(cand)[best])
}

#' Compare a validation score with an alternative scheme (legacy API)
#'
#' Re-evaluates the supplied learner under two validation schemes. The historical
#' name and result fields are retained for compatibility. `deleaked` is the score
#' from the comparison scheme, not an identified bias-corrected target score.
#' Fold resampling gives a descriptive interval conditional on the fitted fold
#' predictions; overlapping training sets and scheme selection are not accounted for.
#'
#' @inheritParams estimate_optimism
#' @param reported Optional metric value to multiply by the comparison/user score
#'   ratio. This legacy calculation is a sensitivity scenario, not a statistically
#'   identified correction of another model's reported performance.
#' @param target Optional [prediction_target()]. When supplied, the controlled
#'   scheme is selected by the smallest unsigned nearest-neighbour distance
#'   mismatch among the candidates, instead of a fixed spatial block. This does
#'   not guarantee a match in predictive risk.
#' @param n_boot Number of fold-bootstrap resamples for the interval.
#' @param level Central probability for the descriptive fold-resampling interval.
#' @return An object of class `deleak_estimate`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(150), y = runif(150), z = rnorm(150))
#' deleak_estimate(d, split = sample(rep_len(1:10, 150)), response = "z",
#'                 coords = c("x", "y"), n_boot = 100)
#' @export
deleak_estimate <- function(data, split, response, reported = NULL, target = NULL,
                            predict_fun = NULL, metric = c("rmse", "mae", "brier", "logloss", "poisson"),
                            control = c("block", "buffer"), coords = NULL,
                            dependence = NULL, k_control = NULL, buffer = NULL,
                            n_boot = 500L, level = 0.90) {
  metric <- match.arg(metric); control <- match.arg(control)
  xy <- .extract_coords(data, coords); n <- nrow(xy$coords)
  info <- list(crs = xy$crs, geographic = xy$geographic)
  base <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
  df <- base; df$.x <- xy$coords[, 1]; df$.y <- xy$coords[, 2]
  if (!response %in% names(df)) stop(sprintf("Response column '%s' not found.", response))
  if (is.null(predict_fun)) predict_fun <- .idw_predictor(c(".x", ".y"), response, info = info)

  folds_user <- .parse_split(split, n)
  k <- if (!is.null(k_control)) k_control else max(length(folds_user), 5L)
  control_label <- control; matched_sli <- NA_real_
  if (!is.null(target)) {
    mc <- .matched_control(data, target, response, coords, k, xy, info, dependence)
    ctrl_folds <- mc$folds
    control_label <- paste0("distance-matched (", mc$label, ")"); matched_sli <- mc$sli
  } else {
    ctrl_folds <- switch(
      control,
      block = .parse_split(spatial_block_cv(xy$coords, k), n),
      buffer = {
        buf <- buffer
        if (is.null(buf)) {
          if (is.null(dependence)) dependence <- estimate_dependence(data, response = response, coords = coords)
          buf <- dependence$practical_range
        }
        diameter <- .domain_diameter(xy$coords, info)
        if (buf > 0.25 * diameter) buf <- 0.25 * diameter
        .buffer_folds(xy$coords, .parse_split(spatial_block_cv(xy$coords, k), n), buf, info)
      })
  }

  eu <- .fold_errors(df, folds_user, predict_fun, response, metric)
  ec <- .fold_errors(df, ctrl_folds, predict_fun, response, metric)
  E_cv <- .agg_err(eu, metric); E_dl <- .agg_err(ec, metric)
  ratio <- if (E_cv > 0) E_dl / E_cv else NA_real_

  a <- (1 - level) / 2
  bc <- .boot_metric(ec, metric, n_boot)
  ratio_b <- bc / .boot_metric(eu, metric, n_boot)
  dl_ci <- unname(stats::quantile(bc, c(a, 1 - a)))
  ratio_ci <- unname(stats::quantile(ratio_b, c(a, 1 - a)))

  anchored <- anchored_ci <- NULL
  if (!is.null(reported)) { anchored <- reported * ratio; anchored_ci <- reported * ratio_ci }

  structure(
    list(metric = metric, control = control_label, matched_sli = matched_sli,
         reported_cv = E_cv, deleaked = E_dl, deleaked_ci = dl_ci,
         optimism = E_dl - E_cv, optimism_rel = if (E_dl > 0) (E_dl - E_cv) / E_dl else NA_real_,
         ratio = ratio, ratio_ci = ratio_ci,
         reported = reported, anchored = anchored, anchored_ci = anchored_ci,
         level = level, n = n),
    class = "deleak_estimate")
}

#' @export
print.deleak_estimate <- function(x, ...) {
  pct <- round(100 * (x$level), 0)
  cat("<deleak_estimate>\n")
  cat(sprintf("  metric / control     : %s / %s\n", toupper(x$metric), x$control))
  if (is.finite(x$matched_sli)) {
    cat(sprintf("  comparison SLI      : %+.3f (inspect unsigned mismatch as well)\n", x$matched_sli))
  }
  cat(sprintf("  reported (your CV)   : %.4g\n", x$reported_cv))
  cat(sprintf("  comparison CV score : %.4g   (%d%% fold-resampling interval %.4g, %.4g)\n",
              x$deleaked, pct, x$deleaked_ci[1], x$deleaked_ci[2]))
  cat(sprintf("  relative score gap  : %+.0f%% of comparison score (ratio %.2f)\n",
              100 * x$optimism_rel, x$ratio))
  if (!is.null(x$reported)) {
    cat(sprintf("  rescaling scenario  : %.4g -> %.4g (%d%% interval %.4g, %.4g)\n",
                x$reported, x$anchored, pct, x$anchored_ci[1], x$anchored_ci[2]))
  }
  cat("  This comparison does not establish or correct deployment accuracy.\n")
  invisible(x)
}
