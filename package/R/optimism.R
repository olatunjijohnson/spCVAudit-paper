# Legacy name for a refitted comparison between validation schemes. The observed
# score gap equals deployment optimism only if the comparison score targets the
# deployment risk, which this function cannot establish.

#' Compare prediction losses under two cross-validation schemes
#'
#' Re-evaluates a model under the user's split and a spatial comparison scheme.
#' The result field `optimism` is retained for compatibility and contains the
#' comparison loss minus the user's CV loss. Neither blocking nor buffering
#' establishes that the comparison loss equals deployment error. Interpret this
#' quantity as a difference between validation schemes unless separate evidence
#' justifies the comparison scheme for the declared target.
#'
#' @param data An `sf` object, numeric coordinate matrix, or `data.frame`.
#' @param split The user's split: a `list(test=, train=)`, a list of test-index
#'   vectors, or a fold-id vector of length `nrow(data)`.
#' @param response Name of the numeric response column.
#' @param predict_fun A function `(train_df, test_df)` returning numeric predictions
#'   for the test rows. Defaults to inverse-distance weighting on the coordinates.
#'   Geographic inputs use geodesic distances; other inputs use coordinate units.
#'   Supply the learner of interest; a score contrast from another learner need
#'   not describe its performance.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param metric Scoring rule (lower = better): `"rmse"` (default) or `"mae"` for
#'   continuous responses; `"brier"` or `"logloss"` for binary responses (predictions
#'   in `[0, 1]`); `"poisson"` deviance for counts. Use a proper rule for the response
#'   type so optimism is meaningful (e.g. Brier/log-loss for disease prevalence).
#' @param control The comparison scheme: `"block"` (spatial
#'   k-means blocks, default) or `"buffer"` (buffer spatial block folds by the
#'   practical range, subject to a heuristic extent cap).
#' @param dependence Optional [estimate_dependence()] object (used to set the buffer
#'   for `control = "buffer"`).
#' @param k_control Number of control folds (defaults to the number of user folds).
#' @param buffer Buffer distance for `control = "buffer"` (defaults to the practical
#'   range from `dependence`). Uses metres for geographic data and coordinate
#'   units otherwise. The automatic cap does not guarantee feasible folds.
#' @return An object of class `optimism_estimate`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(120), y = runif(120), z = rnorm(120))
#' estimate_optimism(d, split = sample(rep_len(1:10, 120)), response = "z",
#'                   coords = c("x", "y"))
#' @export
estimate_optimism <- function(data, split, response, predict_fun = NULL,
                              coords = NULL, metric = c("rmse", "mae", "brier", "logloss", "poisson"),
                              control = c("block", "buffer"),
                              dependence = NULL, k_control = NULL, buffer = NULL) {
  metric <- match.arg(metric); control <- match.arg(control)
  xy <- .extract_coords(data, coords)
  n <- nrow(xy$coords)
  info <- list(crs = xy$crs, geographic = xy$geographic)

  base <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
  df <- base
  df$.x <- xy$coords[, 1]; df$.y <- xy$coords[, 2]
  if (!response %in% names(df)) stop(sprintf("Response column '%s' not found.", response))
  if (is.null(predict_fun)) predict_fun <- .idw_predictor(c(".x", ".y"), response, info = info)
  y <- df[[response]]

  folds_user <- .parse_split(split, n)
  k <- if (!is.null(k_control)) k_control else max(length(folds_user), 5L)

  ctrl_folds <- switch(
    control,
    block = .parse_split(spatial_block_cv(xy$coords, k), n),
    buffer = {
      buf <- buffer
      if (is.null(buf)) {
        if (is.null(dependence)) {
          dependence <- estimate_dependence(data, response = response, coords = coords)
        }
        buf <- dependence$practical_range
      }
      # Heuristic cap for very large buffers; feasibility still depends on folds.
      diameter <- .domain_diameter(xy$coords, info)
      if (buf > 0.25 * diameter) {
        warning(sprintf("Buffer (%.3g) exceeds 25%% of the domain; capping at %.3g.",
                        buf, 0.25 * diameter), call. = FALSE)
        buf <- 0.25 * diameter
      }
      # Buffered CV requires spatially contiguous test folds (buffering a dispersed
      # random fold would remove a ring around the whole domain). Buffer blocks.
      base <- .parse_split(spatial_block_cv(xy$coords, k), n)
      .buffer_folds(xy$coords, base, buf, info)
    })

  p_user <- .cv_predict(df, folds_user, predict_fun)
  p_ctrl <- .cv_predict(df, ctrl_folds, predict_fun)
  mf <- .metric_fun(metric)
  tu <- !is.na(p_user); tc <- !is.na(p_ctrl)
  E_cv <- if (any(tu)) mf(y[tu], p_user[tu]) else NA_real_
  E_ctrl <- if (any(tc)) mf(y[tc], p_ctrl[tc]) else NA_real_
  if (is.na(E_ctrl)) {
    warning("Controlled scheme left no usable training data; optimism is undefined. ",
            "Reduce the buffer or use control = 'block'.", call. = FALSE)
  }
  optimism <- E_ctrl - E_cv
  optimism_rel <- if (!is.na(E_ctrl) && E_ctrl > 0) optimism / E_ctrl else NA_real_

  structure(
    list(metric = metric, control = control,
         E_cv = E_cv, E_control = E_ctrl,
         optimism = optimism, optimism_rel = optimism_rel,
         n = n, n_user_folds = length(folds_user), k_control = k),
    class = "optimism_estimate")
}

#' @export
print.optimism_estimate <- function(x, ...) {
  dir <- if (is.na(x$optimism_rel)) "undefined"
         else if (x$optimism_rel > 0) "comparison loss is higher"
         else if (x$optimism_rel < 0) "comparison loss is lower"
         else "equal observed losses"
  cat("<optimism_estimate>\n")
  cat(sprintf("  metric / control  : %s / %s\n", toupper(x$metric), x$control))
  cat(sprintf("  user CV error     : %.4g\n", x$E_cv))
  cat(sprintf("  comparison error  : %.4g\n", x$E_control))
  cat(sprintf("  CV score gap      : %+.4g  (%+.1f%% of comparison)  [%s]\n",
              x$optimism, 100 * x$optimism_rel, dir))
  cat("  A CV-scheme difference does not by itself establish deployment bias.\n")
  invisible(x)
}
