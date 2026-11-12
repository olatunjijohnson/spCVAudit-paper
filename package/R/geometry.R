# Geometry engine: coordinate extraction, CRS handling, nearest-neighbour distances.
# Internal helpers (not exported). See docs/METHOD-SLI.md notation.

# Geometry/CRS info from an sf crs (or NULL).
.geom_info <- function(crs = NA) {
  geographic <- FALSE
  if (inherits(crs, "crs")) {
    lonlat <- sf::st_is_longlat(crs)
    geographic <- isTRUE(lonlat)
  }
  list(crs = crs, geographic = geographic)
}

# Extent used by the heuristic buffer cap, in the same units as .nn_dist.
# Preserve the historical bounding-box diagonal for planar coordinates. For a
# geographic CRS use maximum observed geodesic separation (metres), calculated
# in blocks to avoid storing the complete pairwise distance matrix.
.domain_diameter <- function(coords, info) {
  if (nrow(coords) < 2L) return(0)
  if (!isTRUE(info$geographic))
    return(sqrt(sum(apply(coords, 2L, function(z) diff(range(z)))^2)))
  points <- sf::st_as_sf(as.data.frame(coords), coords = c(1L, 2L), crs = info$crs)
  blocks <- split(seq_len(nrow(coords)), ceiling(seq_len(nrow(coords)) / 512L))
  largest <- 0
  for (idx in blocks) largest <- max(largest, as.numeric(sf::st_distance(points[idx, ], points)))
  largest
}

# Extract an n x 2 coordinate matrix + geometry info from sf, matrix, or data.frame.
# `coords` names/indices are required for non-sf input.
.extract_coords <- function(x, coords = NULL) {
  if (inherits(x, "sf") || inherits(x, "sfc")) {
    xy <- sf::st_coordinates(x)
    if (ncol(xy) < 2) stop("Could not extract 2D coordinates from sf object.")
    info <- .geom_info(sf::st_crs(x))
    return(list(coords = unname(xy[, 1:2, drop = FALSE]),
                crs = info$crs, geographic = info$geographic))
  }
  if (is.matrix(x) && is.numeric(x) && ncol(x) >= 2) {
    return(list(coords = unname(x[, 1:2, drop = FALSE]),
                crs = NA, geographic = FALSE))
  }
  if (is.data.frame(x)) {
    if (is.null(coords)) {
      stop("`coords` (column names or indices) is required for data.frame input.")
    }
    m <- as.matrix(x[, coords, drop = FALSE])
    storage.mode(m) <- "double"
    return(list(coords = unname(m), crs = NA, geographic = FALSE))
  }
  stop("Unsupported spatial input: provide an sf object, numeric matrix, or data.frame + `coords`.")
}

# Pull a response vector from sf/data.frame by column name.
.get_response <- function(data, response) {
  if (is.null(response)) return(NULL)
  if (inherits(data, "sf")) data <- sf::st_drop_geometry(data)
  if (!response %in% names(data)) {
    stop(sprintf("Response column '%s' not found in data.", response))
  }
  out <- data[[response]]
  if (!is.numeric(out)) {
    stop(sprintf("Response column '%s' must be numeric for variogram estimation.", response))
  }
  out
}

# 1-NN Euclidean distance from each query row to the reference set.
.knn1_euclid <- function(query, ref) {
  if (nrow(ref) == 0L) return(rep(NA_real_, nrow(query)))
  if (requireNamespace("FNN", quietly = TRUE)) {
    return(FNN::get.knnx(ref, query, k = 1)$nn.dist[, 1])
  }
  # Base fallback (O(n*m)); fine for MVP-scale data.
  apply(query, 1L, function(p) sqrt(min(colSums((t(ref) - p)^2))))
}

# Nearest-neighbour distance from each query point to a reference set, honouring
# the geometry: geodesic for geographic CRS, Euclidean otherwise.
.nn_dist <- function(query, ref, info) {
  if (nrow(ref) == 0L) return(rep(NA_real_, nrow(query)))
  if (!isTRUE(info$geographic)) {
    return(.knn1_euclid(query, ref))
  }
  q <- sf::st_as_sf(as.data.frame(query), coords = c(1L, 2L), crs = info$crs)
  r <- sf::st_as_sf(as.data.frame(ref), coords = c(1L, 2L), crs = info$crs)
  # Selecting a planar nearest feature and measuring its geodesic distance can
  # choose the wrong neighbour near a pole or the date line when s2 is disabled.
  # Minimize actual geographic distances, with bounded intermediate memory.
  blocks <- split(seq_len(nrow(query)), ceiling(seq_len(nrow(query)) / 512L))
  out <- numeric(nrow(query))
  for (idx in blocks) {
    dd <- matrix(as.numeric(sf::st_distance(q[idx, ], r)), length(idx), nrow(ref))
    out[idx] <- apply(dd, 1L, min)
  }
  out
}

# Retained spatial correlation at each query point from its k nearest training
# points (noisy-OR combination 1 - prod(1 - rho(d_j))). Reduces exactly to
# rho(nearest) when k = 1, so the single-neighbour SLI is the k = 1 special case.
# This combination is a heuristic, not the conditional variance explained by
# several correlated training observations. The one-neighbour theory applies
# only to k = 1. Geographic distances must be handled geodesically for every k.
.retained_corr <- function(query, ref, info, rho, k = 1L) {
  if (nrow(ref) == 0L) return(rep(NA_real_, nrow(query)))
  if (k <= 1L || nrow(ref) == 1L) return(rho(.nn_dist(query, ref, info)))
  k <- min(k, nrow(ref))
  if (isTRUE(info$geographic)) {
    q <- sf::st_as_sf(as.data.frame(query), coords = c(1L, 2L), crs = info$crs)
    r <- sf::st_as_sf(as.data.frame(ref), coords = c(1L, 2L), crs = info$crs)
    dd <- matrix(as.numeric(sf::st_distance(q, r)), nrow(query), nrow(ref))
    D <- t(apply(dd, 1L, function(z) sort(z)[seq_len(k)]))
  } else if (requireNamespace("FNN", quietly = TRUE)) {
    D <- FNN::get.knnx(ref, query, k = k)$nn.dist
  } else {
    D <- t(apply(query, 1L, function(p) sort(sqrt(colSums((t(ref) - p)^2)))[seq_len(k)]))
  }
  1 - apply(1 - rho(D), 1L, prod)
}

# Trapezoidal integration.
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(0)
  sum(diff(x) * (y[-n] + y[-1L]) / 2)
}
