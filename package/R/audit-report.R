# audit_workflow() and report_leakage(): assemble the diagnostics into a
# descriptive report. No score certifies the validity of a predictive error estimate.

# Legacy grade retained in returned objects for compatibility, not used for interpretation.
.grade <- function(sli) {
  a <- abs(sli)
  if (a < 0.05) "A" else if (a < 0.10) "B" else if (a < 0.20) "C" else
    if (a < 0.35) "D" else "F"
}

#' Audit a cross-validation workflow for spatial leakage
#'
#' Runs [detect_leakage()] and adds data-hygiene checks (duplicated coordinates,
#' documented CRS) to produce a scorecard.
#'
#' @inheritParams detect_leakage
#' @param group Optional grouping column for [detect_group_leakage()] (defaults to
#'   grouping by coordinates).
#' @param time Optional name of a time column. When supplied, the audit also runs the
#'   temporal ([detect_temporal_leakage()]) and joint space-time
#'   ([detect_st_leakage()]) channels, so that space-time data are not audited as if
#'   they were purely spatial.
#' @return An object of class `workflow_audit`. The legacy `leakage_grade` field
#'   is an arbitrary binning of the signed contrast, not an accuracy grade.
#' @examples
#' set.seed(1)
#' n <- 150
#' d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
#' g <- as.matrix(expand.grid(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)))
#' tgt <- prediction_target(grid = g, type = "grid")
#' folds <- sample(rep_len(1:10, n))
#' # every channel the supplied data support, in one call
#' audit_workflow(d, split = folds, target = tgt, response = "z",
#'                coords = c("x", "y"))
#'
#' # pass `time` and the temporal and joint space-time channels run too, so
#' # space-time data are never audited as if they were purely spatial
#' d$year <- sample(2000:2015, n, replace = TRUE)
#' audit_workflow(d, split = folds, target = tgt, response = "z",
#'                coords = c("x", "y"), time = "year")
#' @export
audit_workflow <- function(data, split, target, response = NULL,
                           dependence = NULL, coords = NULL, group = NULL,
                           time = NULL) {
  lk <- detect_leakage(data, split, target, response = response,
                       dependence = dependence, coords = coords)
  gl <- detect_group_leakage(data, split, group = group, coords = coords)
  xy <- .extract_coords(data, coords)
  trend <- if (!is.null(response)) tryCatch(trend_strength(data, response, coords),
                                            error = function(e) NA_real_) else NA_real_
  tl <- st <- NULL
  if (!is.null(time)) {
    tl <- tryCatch(detect_temporal_leakage(data, split, time = time),
                   error = function(e) NULL)
    st <- tryCatch(detect_st_leakage(data, split, time = time, coords = coords),
                   error = function(e) NULL)
  }
  flags <- list(
    leakage_grade = .grade(lk$SLI_rho),
    sli_rho = lk$SLI_rho,
    unsigned_mismatch = lk$W,
    distance_cdfs_cross = lk$crosses,
    trend_strength = trend,
    duplicated_coords = sum(duplicated(xy$coords)),
    group_leak_frac = gl$frac_leaked,
    lookahead_frac = if (!is.null(tl)) tl$lookahead_frac else NA_real_,
    sli_st = if (!is.null(st)) st$sli_st else NA_real_,
    crs_documented = inherits(xy$crs, "crs") && !is.na(xy$crs))
  structure(list(leakage = lk, group_leakage = gl, temporal = tl, st = st,
                 flags = flags),
            class = "workflow_audit")
}

#' @export
print.workflow_audit <- function(x, ...) {
  f <- x$flags
  cat("<workflow_audit>\n")
  cat(sprintf("  signed correlation contrast: %+.3f\n", f$sli_rho))
  cat(sprintf("  unsigned distance mismatch : %.4g\n", x$leakage$W))
  cat(sprintf("  distance ECDFs cross       : %s\n", isTRUE(x$leakage$crosses)))
  cat("  These diagnostics do not certify unbiased validation or corrected accuracy.\n")
  if (is.finite(f$trend_strength)) {
    cat(sprintf("  %s trend strength: %.2f (inspect the fitted mean and prediction target)\n",
                if (f$trend_strength > 0.5) "[!]" else "[ok]", f$trend_strength))
  }
  cat(sprintf("  %s duplicated coordinates: %d\n",
              if (f$duplicated_coords > 0) "[!]" else "[ok]", f$duplicated_coords))
  cat(sprintf("  %s groups/locations shared across training and test: %.1f%% of test points\n",
              if (f$group_leak_frac > 0) "[!]" else "[ok]", 100 * f$group_leak_frac))
  if (is.finite(f$lookahead_frac)) {
    cat(sprintf("  %s temporal lookahead: %.1f%% of test points trained on the future\n",
                if (f$lookahead_frac > 0) "[!]" else "[ok]", 100 * f$lookahead_frac))
    if (is.finite(f$sli_st))
      cat(sprintf("  joint space-time contrast: SLI = %+.3f (experimental diagnostic)\n", f$sli_st))
  }
  cat(sprintf("  %s CRS documented: %s\n",
              if (f$crs_documented) "[ok]" else "[!]", f$crs_documented))
  invisible(x)
}

#' Assemble a leakage report (scorecard)
#'
#' Combines a leakage diagnosis with optional optimism and validation
#' recommendation into a single descriptive report. The legacy `grade` field is
#' retained for compatibility and is not a measure of predictive accuracy.
#'
#' @param diagnosis A [detect_leakage()] result.
#' @param optimism An optional [estimate_optimism()] result.
#' @param recommendation An optional [recommend_validation()] result.
#' @param deleak An optional [deleak_estimate()] legacy score-comparison result.
#' @return An object of class `leakage_report`.
#' @examples
#' set.seed(1)
#' n <- 150
#' d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
#' g <- as.matrix(expand.grid(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)))
#' tgt <- prediction_target(grid = g, type = "grid")
#' folds <- sample(rep_len(1:10, n))
#' lk <- detect_leakage(d, split = folds, target = tgt, response = "z",
#'                      coords = c("x", "y"))
#' op <- estimate_optimism(d, folds, response = "z", coords = c("x", "y"))
#' report_leakage(lk, optimism = op)
#' @export
report_leakage <- function(diagnosis, optimism = NULL, recommendation = NULL,
                           deleak = NULL) {
  stopifnot(inherits(diagnosis, "leakage_diagnosis"))
  structure(
    list(diagnosis = diagnosis, optimism = optimism, deleak = deleak,
         recommendation = recommendation, grade = .grade(diagnosis$SLI_rho)),
    class = "leakage_report")
}

#' @export
print.leakage_report <- function(x, ...) {
  cat("================ spLeakage report ================\n")
  cat(" Validation-target mismatch and score comparisons\n")
  cat("--------------------------------------------------\n")
  print(x$diagnosis)
  if (!is.null(x$optimism)) { cat("\n"); print(x$optimism) }
  if (!is.null(x$deleak)) { cat("\n"); print(x$deleak) }
  if (!is.null(x$recommendation)) { cat("\n"); print(x$recommendation) }
  cat("==================================================\n")
  invisible(x)
}
