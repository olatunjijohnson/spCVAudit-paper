# Spatial dependence estimation: empirical variogram + exponential fit, giving the
# correlation function rho(h), the practical range phi, optional geometric
# anisotropy, and the parameter covariance used to propagate uncertainty into SLI.
# See docs/METHOD-SLI.md sections 1-2, 6 (anisotropy) and 7 (uncertainty).

# Pairwise distances + directional angles for a (sub)set of pairs.
.pair_geometry <- function(coords, info, max_pairs = 2e6) {
  n <- nrow(coords)
  if (n < 2L) stop("At least two locations are required for pair geometry.")
  if (length(max_pairs) != 1L || !is.finite(max_pairs) || max_pairs < 1 ||
      max_pairs != floor(max_pairs) || max_pairs > .Machine$integer.max)
    stop("`max_pairs` must be a positive integer within R's indexing limit.")
  total <- as.double(n) * (n - 1) / 2
  if (total <= max_pairs) {
    pairs <- utils::combn(n, 2L)
  } else {
    if (total > 2^53) stop("Too many location pairs for exact pair indexing.")
    # Sample ranks in combn's lexicographic ordering before allocating pairs.
    # Constructing all n*(n-1)/2 pairs first defeats the max_pairs memory limit.
    ranks <- sample.int(total, max_pairs)
    i <- floor(n - .5 - sqrt(pmax((n - .5)^2 - 2 * (ranks - 1), 0))) + 1
    i <- pmin(pmax(i, 1), n - 1)
    start <- (i - 1) * (2 * n - i) / 2
    i[ranks <= start] <- i[ranks <= start] - 1
    start <- (i - 1) * (2 * n - i) / 2
    i[ranks > start + n - i] <- i[ranks > start + n - i] + 1
    start <- (i - 1) * (2 * n - i) / 2
    j <- ranks - start + i
    stopifnot(all(i >= 1), all(i < j), all(j <= n))
    pairs <- rbind(as.integer(i), as.integer(j))
  }
  i <- pairs[1L, ]; j <- pairs[2L, ]
  dxy <- coords[i, , drop = FALSE] - coords[j, , drop = FALSE]
  if (isTRUE(info$geographic)) {
    qi <- sf::st_as_sf(as.data.frame(coords[i, , drop = FALSE]), coords = c(1L, 2L), crs = info$crs)
    qj <- sf::st_as_sf(as.data.frame(coords[j, , drop = FALSE]), coords = c(1L, 2L), crs = info$crs)
    h <- as.numeric(sf::st_distance(qi, qj, by_element = TRUE))
    angle <- NA_real_
  } else {
    h <- sqrt(rowSums(dxy^2))
    angle <- atan2(dxy[, 2], dxy[, 1]) %% pi   # undirected, in [0, pi)
  }
  list(i = i, j = j, h = h, angle = angle)
}

# Bin pairs into an empirical semivariogram. `keep` optionally restricts pairs
# (e.g. to a direction).
.bin_variogram <- function(h, gamma2, n_bins, cutoff, keep = NULL) {
  if (!is.null(keep)) { h <- h[keep]; gamma2 <- gamma2[keep] }
  if (is.null(cutoff)) cutoff <- max(h) / 2
  sel <- h <= cutoff & h > 0
  h <- h[sel]; gamma2 <- gamma2[sel]
  if (length(h) < 3L) return(NULL)
  br <- seq(0, cutoff, length.out = n_bins + 1L)
  bin <- cut(h, breaks = br, include.lowest = TRUE)
  dist <- tapply(h, bin, mean); gamma <- tapply(gamma2, bin, mean)
  npb <- tapply(h, bin, length)
  ok <- !is.na(dist) & !is.na(gamma)
  data.frame(dist = as.numeric(dist[ok]), gamma = as.numeric(gamma[ok]),
             np = as.numeric(npb[ok]))
}

# Empirical (omnidirectional) semivariogram.
.empirical_variogram <- function(coords, z, info, n_bins = 15L, cutoff = NULL,
                                 max_pairs = 2e6) {
  if (nrow(coords) < 5L) stop("Need at least 5 locations to estimate a variogram.")
  pg <- .pair_geometry(coords, info, max_pairs)
  gamma2 <- 0.5 * (z[pg$i] - z[pg$j])^2
  vg <- .bin_variogram(pg$h, gamma2, n_bins, cutoff)
  if (is.null(vg)) stop("Could not form variogram bins (check cutoff / data).")
  vg
}

# Fit exponential model gamma(h) = nugget + psill * (1 - exp(-h / range)).
# Returns coefficients and (when available) the parameter covariance.
.fit_exponential <- function(vg) {
  if (nrow(vg) < 3L) stop("Too few variogram bins to fit a model.")
  g_max <- max(vg$gamma); g_min <- min(vg$gamma)
  if (!all(is.finite(c(vg$gamma, vg$dist, vg$np))) || g_min < 0 || max(vg$dist) <= 0)
    stop("Variogram bins must have finite, nonnegative semivariances and positive distances.")
  if (g_max == 0) return(list(coef = list(nugget = 0, psill = 0,
    range = max(vg$dist) / 3), vcov = NULL))
  start <- list(nugget = max(g_min, 1e-8 * g_max),
                psill = max(g_max - g_min, 1e-8 * g_max),
                range = max(vg$dist) / 3)
  fit <- tryCatch(
    suppressWarnings(stats::nls(
      gamma ~ nugget + psill * (1 - exp(-dist / range)),
      data = vg, start = start, weights = vg$np,
      control = stats::nls.control(maxiter = 200, warnOnly = TRUE))),
    error = function(e) NULL)
  V <- NULL
  valid_fit <- !is.null(fit) && isTRUE(fit$convInfo$isConv)
  if (valid_fit) {
    candidate <- stats::coef(fit)
    valid_fit <- all(is.finite(candidate)) && candidate[["nugget"]] >= 0 &&
      candidate[["psill"]] >= 0 && candidate[["range"]] > 0
  }
  # Clipping a nonphysical unconstrained fit is not constrained least squares.
  # Refit within the existing physical bounds and do not retain its old vcov.
  if (valid_fit) {
    co <- as.list(stats::coef(fit))
    V <- tryCatch(stats::vcov(fit)[c("nugget", "psill", "range"),
                                   c("nugget", "psill", "range")],
                  error = function(e) NULL)
  } else {
    # Normalize both axes so metres versus kilometres, or a response-unit
    # change, cannot alter optimizer step sizes and convergence behavior.
    h_scale <- max(vg$dist)
    h <- vg$dist / h_scale; gamma <- vg$gamma / g_max
    weight <- vg$np / sum(vg$np)
    obj <- function(p) {
      residual <- p[1] + p[2] * (1 - exp(-h / p[3])) - gamma
      sum(weight * residual^2)
    }
    gradient <- function(p) {
      e <- exp(-h / p[3]); residual <- p[1] + p[2] * (1 - e) - gamma
      2 * c(sum(weight * residual), sum(weight * residual * (1 - e)),
            sum(weight * residual * (-p[2] * e * h / p[3]^2)))
    }
    initial <- c(start$nugget / g_max, start$psill / g_max, start$range / h_scale)
    lower <- c(0, 1e-8, 1e-6); upper <- c(1, 2, 5)
    op <- stats::optim(initial, obj, gr = gradient, method = "L-BFGS-B",
                       lower = lower, upper = upper,
                       control = list(maxit = 1000L, factr = 1e7))
    if (op$convergence != 0L || !is.finite(op$value)) {
      alternate <- stats::nlminb(initial, obj, gradient = gradient,
        lower = lower, upper = upper, control = list(iter.max = 1000L, eval.max = 2000L))
      if (alternate$convergence != 0L || !is.finite(alternate$objective))
        stop("Constrained variogram fit did not converge.")
      op$par <- alternate$par
    }
    # A bounded optimizer can return a value a few ulps outside a bound (for
    # example a nugget of -1e-18). Reject material violations, then restore exact
    # feasibility only for this numerical roundoff from the constrained fit.
    tolerance <- 64 * .Machine$double.eps * pmax(1, abs(lower), abs(upper))
    if (any(!is.finite(op$par)) || any(op$par < lower - tolerance) || any(op$par > upper + tolerance))
      stop("Constrained variogram optimizer returned parameters outside its bounds.")
    op$par <- pmin(pmax(op$par, lower), upper)
    co <- list(nugget = op$par[1] * g_max, psill = op$par[2] * g_max,
               range = op$par[3] * h_scale)
  }
  list(coef = co, vcov = V)
}

# Estimate geometric anisotropy from directional ranges (projected data only).
.estimate_anisotropy <- function(coords, z, info, n_bins, cutoff) {
  if (isTRUE(info$geographic)) return(NULL)
  pg <- .pair_geometry(coords, info)
  gamma2 <- 0.5 * (z[pg$i] - z[pg$j])^2
  dirs <- c(0, pi / 4, pi / 2, 3 * pi / 4); tol <- pi / 8
  ranges <- vapply(dirs, function(a) {
    keep <- pmin(abs(pg$angle - a), pi - abs(pg$angle - a)) <= tol
    vg <- .bin_variogram(pg$h, gamma2, n_bins, cutoff, keep = keep)
    if (is.null(vg)) return(NA_real_)
    tryCatch(.fit_exponential(vg)$coef$range, error = function(e) NA_real_)
  }, numeric(1))
  if (sum(is.finite(ranges)) < 2L) return(NULL)
  imax <- which.max(ranges); imin <- which.min(ranges)
  ratio <- max(min(ranges[imin] / ranges[imax], 1), 0.05)
  list(angle = dirs[imax], ratio = ratio, range_major = ranges[imax],
       dirs = dirs, ranges = ranges)
}

# Linear transform implementing geometric anisotropy: rotate so the major axis is
# x, then stretch the minor axis by 1/ratio. Euclidean distance in the transformed
# space equals the anisotropic (dependence) distance in major-range units. Under
# isotropy (ratio = 1) this is a pure rotation and preserves all distances.
.apply_aniso <- function(coords, aniso) {
  if (is.null(aniso) || isTRUE(aniso$ratio == 1)) return(coords)  # isotropic no-op
  th <- aniso$angle
  Rt <- matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2L, 2L)  # rotate by -theta
  rot <- coords %*% Rt
  rot[, 2] <- rot[, 2] / aniso$ratio
  rot
}

#' Estimate spatial dependence (variogram, range, correlation, anisotropy)
#'
#' Fits an empirical semivariogram and an exponential model, returning the
#' correlation function `rho(h)` between distinct observations, the variogram range, the practical range `phi`,
#' optional geometric anisotropy, and the parameter covariance used to propagate
#' uncertainty into the Spatial Leakage Index.
#' The nugget is treated as independent observation noise. Thus `rho(0)` is the
#' signal proportion for two distinct observations at the same location; it is
#' not the self-correlation of one observation. Shared site effects or correlated
#' measurement errors require a different covariance model.
#'
#' @param data An `sf` object, numeric coordinate matrix, or `data.frame`.
#' @param response Name of the (numeric) response column used for the variogram.
#' @param coords For non-`sf` input, column names/indices of the coordinates.
#' @param n_bins Number of variogram bins.
#' @param cutoff Maximum lag distance; defaults to half the maximum pairwise distance.
#' @param anisotropy `NULL` (isotropic, default), `"auto"` (estimate geometric
#'   anisotropy from directional variograms; projected data only), or a
#'   `list(angle =, ratio =)` giving the major-axis angle (radians) and the
#'   minor/major range ratio in `(0, 1]`.
#' @return An object of class `sp_dependence` with elements `rho` (a function),
#'   `range`, `practical_range` (`phi`), `psill`, `nugget`, `signal_prop`,
#'   `anisotropy` (or `NULL`), `coef`, `vcov`, and `variogram`.
#' @examples
#' set.seed(1)
#' n <- 200
#' d <- data.frame(x = runif(n), y = runif(n))
#' d$z <- sin(6 * d$x) + cos(6 * d$y) + rnorm(n, 0, 0.3)
#' estimate_dependence(d, response = "z", coords = c("x", "y"))
#' @export
estimate_dependence <- function(data, response, coords = NULL,
                                n_bins = 15L, cutoff = NULL, anisotropy = NULL) {
  xy <- .extract_coords(data, coords)
  z <- .get_response(data, response)
  if (is.null(z)) stop("`response` is required to estimate spatial dependence.")
  info <- list(crs = xy$crs, geographic = xy$geographic)

  # Geometric anisotropy (optional). Work in transformed coordinates so the
  # exponential range and rho are defined on the dependence distance.
  aniso <- NULL
  if (!is.null(anisotropy)) {
    if (isTRUE(info$geographic)) {
      warning("Anisotropy is not supported for geographic CRS; ignoring.")
    } else if (identical(anisotropy, "auto")) {
      aniso <- .estimate_anisotropy(xy$coords, z, info, n_bins, cutoff)
    } else if (is.list(anisotropy)) {
      aniso <- list(angle = anisotropy$angle, ratio = anisotropy$ratio)
    } else {
      stop("`anisotropy` must be NULL, 'auto', or list(angle=, ratio=).")
    }
  }
  work <- .apply_aniso(xy$coords, aniso)

  vg <- .empirical_variogram(work, z, info, n_bins = n_bins, cutoff = cutoff)
  fit <- .fit_exponential(vg)
  co <- fit$coef
  total <- co$psill + co$nugget
  signal_prop <- if (total > 0) co$psill / total else 0
  range_par <- co$range
  rho <- function(h) signal_prop * exp(-h / range_par)

  structure(
    list(rho = rho, range = range_par, practical_range = 3 * range_par,
         psill = co$psill, nugget = co$nugget, signal_prop = signal_prop,
         anisotropy = aniso, coef = co, vcov = fit$vcov,
         variogram = vg, geographic = info$geographic,
         nugget_convention = "independent observation noise"),
    class = "sp_dependence")
}

# Sample (range, signal_prop) from the fit uncertainty for SLI error propagation.
.sample_params <- function(dep, n) {
  co <- dep$coef; V <- dep$vcov; mu <- c(co$nugget, co$psill, co$range)
  draws <- NULL
  if (!is.null(V) && all(is.finite(V))) {
    U <- tryCatch(chol(V), error = function(e) NULL)
    if (!is.null(U)) draws <- sweep(matrix(stats::rnorm(n * 3L), n, 3L) %*% U, 2L, mu, "+")
  }
  if (is.null(draws)) return(NULL)
  draws[, 1] <- pmax(draws[, 1], 0)
  draws[, 2] <- pmax(draws[, 2], 1e-10)
  draws[, 3] <- pmax(draws[, 3], 1e-8)
  data.frame(range = draws[, 3], signal_prop = draws[, 2] / (draws[, 2] + draws[, 1]))
}

#' @export
print.sp_dependence <- function(x, ...) {
  cat("<sp_dependence> exponential variogram\n")
  cat(sprintf("  range parameter : %.4g\n", x$range))
  cat(sprintf("  practical range : %.4g (phi)\n", x$practical_range))
  cat(sprintf("  partial sill    : %.4g\n", x$psill))
  cat(sprintf("  nugget          : %.4g\n", x$nugget))
  cat(sprintf("  signal prop.    : %.3f\n", x$signal_prop))
  if (!is.null(x$anisotropy)) {
    cat(sprintf("  anisotropy      : angle %.1f deg, ratio %.2f\n",
                x$anisotropy$angle * 180 / pi, x$anisotropy$ratio))
  }
  cat(sprintf("  param cov.      : %s\n", if (is.null(x$vcov)) "unavailable" else "available (approximate)"))
  invisible(x)
}
