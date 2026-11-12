# Recommendations identify candidate assessments and their assumptions.
# Sampling design and clustering are distinct attributes; neither determines CV bias.

# Clark-Evans nearest-neighbour index (R < 1 clustered, ~1 random, > 1 regular).
# A risk flag only -- it detects spatial pattern, not whether the sample is a
# probability sample (which geometry cannot reveal).
.nn_index <- function(coords, info) {
  n <- nrow(coords)
  if (n < 2L || isTRUE(info$geographic)) return(NA_real_)
  if (requireNamespace("FNN", quietly = TRUE) && !isTRUE(info$geographic)) {
    nn <- FNN::get.knn(coords, k = 1L)$nn.dist[, 1]
  } else {
    nn <- vapply(seq_len(n), function(i) {
      d <- sqrt(colSums((t(coords[-i, , drop = FALSE]) - coords[i, ])^2))
      min(d)
    }, numeric(1))
  }
  area <- prod(apply(coords, 2L, function(z) diff(range(z))))
  if (area <= 0) return(NA_real_)
  expected <- 0.5 / sqrt(n / area)
  mean(nn) / expected
}

# Assessment of a fixed map differs from CV of a fitted learning procedure.
.decision <- function(estimand, design, target) {
  if (estimand == "population") {
    # Design-based map-accuracy estimand.
    if (design == "probability") {
      return(list(
        recommended = "Design-based assessment of the fixed map using a probability validation sample",
        avoid = "Treating CV of refitted models as a design-based assessment of a fixed map",
        spatial = FALSE,
        rationale = paste("Use validation observations excluded from model fitting, known inclusion",
                          "probabilities, and an estimator appropriate to the survey design.",
                          "A probability sample can be clustered. Probability sampling alone",
                          "does not make ordinary K-fold CV unbiased.")))
    }
    return(list(
      recommended = c("Collect a probability validation sample", "Model-based assessment with explicit sampling assumptions"),
      avoid = "Claiming design-unbiased map accuracy from a non-probability validation sample",
      spatial = NA,
      rationale = paste("Population-mean accuracy from a non-probability sample lacks a",
                        "design basis; results are model-dependent. Treat with caution.")))
  }

  # A declaration of probability sampling does not specify the prediction task,
  # the learning curve, or the learner's use of spatial dependence.
  if (design == "probability") {
    return(list(
      recommended = c("Compare random CV geometry with the declared prediction target",
                      "Compare target-aware CV when prediction conditions differ"),
      avoid = "Assuming random or spatial CV is unbiased solely because sampling was probabilistic",
      spatial = NA,
      rationale = paste("CV assesses refitted models with smaller training samples.",
                        "Its suitability depends on the target, learner, loss and data-generating process.",
                        "Probability sampling and clustering are separate attributes.")))
  }
  # clustered / convenience / unknown -> spatial CV matched to the target geometry.
  rec <- switch(target,
    interpolation = c("Compare random CV and target-aware NNDM / kNNDM CV"),
    grid          = c("Compare random, spatial block and NNDM / kNNDM CV"),
    newdata       = c("Compare target-aware NNDM / kNNDM and other feasible CV schemes"))
  if (target == "newdata") rec <- c(rec, "Inspect covariate support and area of applicability")
  list(
    recommended = rec,
    avoid = "Interpreting a small signed mismatch or a spatial split as proof of unbiased prediction error",
    spatial = NA,
    rationale = paste0(
      "For a '", target, "' prediction target and declared '", design,
      "' sampling information, inspect complete validation and target distance distributions,",
      " dependence assumptions, training sizes and covariate support. Distance matching is",
      " a diagnostic criterion; its relationship to error depends on the prediction model."))
}

#' Recommend a validation strategy for a sampling design and prediction target
#'
#' Operationalises the design x estimand x target framework: it asks for the
#' estimand and the sampling design (which cannot be inferred from coordinates) and
#' returns candidate assessments and their assumptions. Any supplied data is used
#' only to compute a projected-coordinate clustering flag, never to infer the design.
#' Spatial CV suitability is indeterminate for a prediction task without further
#' information about the learner and target. Probability sampling can be clustered.
#'
#' @param data Optional `sf`/matrix/`data.frame`, used only for the clustering flag.
#' @param estimand `"prediction"` (conditional predictive skill at locations,
#'   default) or `"population"` (population-mean map accuracy over a region).
#' @param design `"unknown"` (default), `"probability"`, `"clustered"`, or
#'   `"convenience"`. These legacy labels record supplied sampling information;
#'   clustering alone does not establish whether sampling was probabilistic.
#' @param target `"grid"` (wall-to-wall, default), `"interpolation"`, or `"newdata"`.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @return An object of class `validation_recommendation`.
#' @examples
#' # Candidate assessments for a clustered sample to be mapped wall-to-wall.
#' recommend_validation(estimand = "prediction", design = "clustered", target = "grid")
#' # Probability sampling does not by itself determine the best CV scheme.
#' recommend_validation(estimand = "prediction", design = "probability", target = "grid")
#' @export
recommend_validation <- function(data = NULL,
                                 estimand = c("prediction", "population"),
                                 design = c("unknown", "probability", "clustered", "convenience"),
                                 target = c("grid", "interpolation", "newdata"),
                                 coords = NULL) {
  estimand <- match.arg(estimand); design <- match.arg(design); target <- match.arg(target)
  flag <- NULL
  if (!is.null(data)) {
    xy <- .extract_coords(data, coords)
    R <- .nn_index(xy$coords, list(geographic = xy$geographic))
    flag <- list(nn_index = R, clustered = isTRUE(R < 0.85))
  }
  eff_design <- design
  rec <- .decision(estimand, eff_design, target)

  structure(
    list(estimand = estimand, design = design, target = target,
         recommended = rec$recommended, avoid = rec$avoid,
         spatial_cv_appropriate = rec$spatial, rationale = rec$rationale,
         clustering = flag, assumed_design = eff_design),
    class = "validation_recommendation")
}

#' @export
print.validation_recommendation <- function(x, ...) {
  cat("<validation_recommendation>\n")
  cat(sprintf("  estimand / design / target : %s / %s / %s\n",
              x$estimand, x$design, x$target))
  if (x$design == "unknown") {
    cat("  NOTE: sampling design is unknown; it has not been inferred from geometry.\n")
  }
  sp <- if (is.na(x$spatial_cv_appropriate)) "depends (see rationale)"
        else if (x$spatial_cv_appropriate) "YES" else "NO"
  cat(sprintf("  spatial CV appropriate     : %s\n", sp))
  cat("  recommended:\n"); for (r in x$recommended) cat(sprintf("    - %s\n", r))
  cat(sprintf("  avoid: %s\n", x$avoid))
  if (!is.null(x$clustering) && is.finite(x$clustering$nn_index)) {
    cat(sprintf("  clustering flag (risk only): NN index = %.2f%s\n",
                x$clustering$nn_index,
                if (x$clustering$clustered) " (clustered)" else ""))
  }
  cat(sprintf("  rationale: %s\n", x$rationale))
  invisible(x)
}
