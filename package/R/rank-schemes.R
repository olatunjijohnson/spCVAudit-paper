# Rank spatial mismatch diagnostics. Ranking does not establish target-risk validity.

#' Rank candidate cross-validation schemes by deployment match
#'
#' Builds each candidate CV scheme, measures its signed contrast against the
#' declared prediction target, and ranks them by a chosen mismatch summary.
#' The default unsigned Wasserstein distance avoids cancellation between opposing
#' CDF differences. A small value describes distance matching; it does not establish
#' unbiased predictive performance. The legacy signed criterion remains available
#' for comparison and should be interpreted alongside the full distributions.
#' Buffers use metres for geographic inputs and coordinate units otherwise.
#' Their automatic cap is a heuristic fraction of the observed extent, not a
#' guarantee of independent observations or feasible folds. Block construction
#' remains coordinate k-means; project geographic inputs for metric blocks.
#'
#' @param data An `sf` object, numeric matrix, or `data.frame`.
#' @param target A [prediction_target()] (or coordinates coercible to one).
#' @param response Name of the numeric response column (for the dependence model).
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param k Number of folds for each candidate scheme.
#' @param schemes Candidate schemes: any of `"random"`, `"block"`, `"buffered"`.
#' @param dependence Optional precomputed [estimate_dependence()] object.
#' @param seed RNG seed for the random scheme.
#' @param criterion Ranking criterion: `"W"` (default, distance-distribution
#'   mismatch), `"W_rho"` (retained-correlation distribution mismatch), or
#'   `"abs_sli"` (absolute signed contrast, which can hide cancellation).
#' @return An object of class `scheme_ranking`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(150), y = runif(150), z = rnorm(150))
#' grid <- as.matrix(expand.grid(x = seq(0, 1, .2), y = seq(0, 1, .2)))
#' rank_cv_schemes(d, prediction_target(grid = grid, type = "grid"),
#'                 response = "z", coords = c("x", "y"))
#' @export
rank_cv_schemes <- function(data, target, response, coords = NULL, k = 10L,
                            schemes = c("random", "block", "buffered"),
                            dependence = NULL, seed = 1L,
                            criterion = c("W", "W_rho", "abs_sli")) {
  schemes <- match.arg(schemes, several.ok = TRUE)
  criterion <- match.arg(criterion)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  xy <- .extract_coords(data, coords); n <- nrow(xy$coords)
  info <- list(crs = xy$crs, geographic = xy$geographic)
  if (is.null(dependence)) dependence <- estimate_dependence(data, response = response, coords = coords)
  tgt <- .as_prediction_target(target, data = data, coords = coords)

  build <- function(s) {
    set.seed(seed)
    switch(
      s,
      random   = .parse_split(sample(rep_len(seq_len(k), n)), n),
      block    = .parse_split(spatial_block_cv(xy$coords, k), n),
      buffered = {
        diam <- .domain_diameter(xy$coords, info)
        buf <- min(dependence$practical_range, 0.25 * diam)
        .buffer_folds(xy$coords, .parse_split(spatial_block_cv(xy$coords, k), n), buf, info)
      })
  }
  rows <- lapply(schemes, function(s) {
    lk <- detect_leakage(data, build(s), tgt, dependence = dependence, coords = coords)
    data.frame(scheme = s, sli_rho = lk$SLI_rho, abs_sli = abs(lk$SLI_rho),
               W = lk$W, W_rho = lk$W_rho, crosses = lk$crosses,
               verdict = .leak_verdict(lk$SLI_rho, lk$delta),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab[[criterion]]), , drop = FALSE]
  rownames(tab) <- NULL
  structure(list(ranking = tab, best = tab$scheme[1], target = tgt$type, criterion = criterion),
            class = "scheme_ranking")
}

#' @export
print.scheme_ranking <- function(x, ...) {
  cat(sprintf("<scheme_ranking>  target = '%s'\n", x$target))
  criterion <- if (is.null(x$criterion)) "abs_sli (legacy)" else x$criterion
  cat(sprintf("  ranked by %s (smaller mismatch first):\n", criterion))
  r <- x$ranking
  for (i in seq_len(nrow(r))) {
    cat(sprintf("    %s%-9s SLI_rho %+.3f  W %.3g  [%s]\n",
                if (i == 1L) "* " else "  ", r$scheme[i], r$sli_rho[i], r$W[i], r$verdict[i]))
  }
  cat(sprintf("  smallest diagnostic: %s; target-risk validity requires separate assessment.\n", x$best))
  invisible(x)
}
