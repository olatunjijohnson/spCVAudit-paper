# Shared reference calculations. Run scripts from the repository root.
# These routines do not depend on the historical package's fitted covariance code.

revision_root <- normalizePath("revision", mustWork = TRUE)

pair_distance <- function(x, y = x) {
  x <- as.matrix(x); y <- as.matrix(y)
  sqrt(pmax(outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y), 0))
}

correlation_matrix <- function(x, y = x, range = 0.15, family = "exponential",
                               axis_ratio = 1, angle = 0) {
  transform <- function(z) {
    rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2)
    z <- as.matrix(z) %*% rotation
    z[, 2] <- z[, 2] / axis_ratio
    z
  }
  h <- pair_distance(transform(x), transform(y)) / range
  switch(family,
    exponential = exp(-h),
    matern32 = (1 + sqrt(3) * h) * exp(-sqrt(3) * h),
    stop("Unknown covariance family"))
}

correlation_at <- function(h, range, signal = 1, family = "exponential") {
  z <- h / range
  signal * switch(family, exponential = exp(-z),
                  matern32 = (1 + sqrt(3) * z) * exp(-sqrt(3) * z),
                  stop("Unknown covariance family"))
}

variance_reduction <- function(training_covariance, cross_covariance) {
  U <- chol(training_covariance)
  z <- forwardsolve(t(U), t(cross_covariance))
  as.numeric(colSums(z^2))
}

kriging_prediction <- function(y, training_covariance, cross_covariance) {
  U <- chol(training_covariance)
  as.numeric(cross_covariance %*% backsolve(U, forwardsolve(t(U), y)))
}

distance_summary <- function(cv, target, cfun) {
  stopifnot(length(cv) > 0, length(target) > 0,
            all(is.finite(cv)), all(is.finite(target)), min(cv) >= 0, min(target) >= 0)
  areas <- function(x, y) {
    cuts <- sort(unique(c(x, y)))
    if (length(cuts) < 2) return(c(positive = 0, negative = 0, W = 0))
    diffcdf <- (ecdf(x)(cuts) - ecdf(y)(cuts))[-length(cuts)]
    pos <- sum(diff(cuts) * pmax(diffcdf, 0))
    neg <- sum(diff(cuts) * pmax(-diffcdf, 0))
    c(positive = pos, negative = neg, W = pos + neg)
  }
  ar <- areas(cv, target)
  rc <- cfun(cv); rt <- cfun(target)
  c(A = mean(target) - mean(cv), W = unname(ar["W"]),
    J1 = mean(rc) - mean(rt), J2 = mean(rc^2) - mean(rt^2),
    U1 = unname(areas(rc, rt)["W"]), U2 = unname(areas(rc^2, rt^2)["W"]),
    Aplus = unname(ar["positive"]), Aminus = unname(ar["negative"]),
    crosses = as.numeric(all(ar[c("positive", "negative")] > 1e-12 * max(ar["W"], .Machine$double.eps))))
}

sampling_locations <- function(n, design) {
  if (design == "random") return(matrix(runif(2 * n), n, 2))
  if (design == "clustered") {
    centers <- matrix(runif(12, 0.12, 0.88), 6, 2)
    id <- sample.int(6, n, replace = TRUE)
    z <- centers[id, ] + matrix(rnorm(2 * n, sd = 0.04), n, 2)
    # Reflect at the boundary; unlike clipping, this does not create coincident sites.
    z <- abs(z %% 2)
    z[z > 1] <- 2 - z[z > 1]
    return(z)
  }
  stop("Unknown sampling design")
}

target_locations <- function(type, side = 15L) {
  mid <- (seq_len(side) - 0.5) / side
  out <- as.matrix(expand.grid(x = mid, y = mid))
  if (type == "shifted") out[, 1] <- out[, 1] + 0.5
  else if (type != "domain") stop("Unknown target")
  out
}

folds_from_ids <- function(ids) {
  lapply(sort(unique(ids)), function(f) list(test = which(ids == f), train = which(ids != f)))
}

make_candidates <- function(xy, k = 5L, buffer = 0.05) {
  n <- nrow(xy)
  random <- folds_from_ids(sample(rep_len(seq_len(k), n)))
  blocks <- folds_from_ids(kmeans(xy, k, nstart = 5, iter.max = 100)$cluster)
  dd <- pair_distance(xy)
  buffered <- lapply(blocks, function(f) {
    keep <- apply(dd[f$train, f$test, drop = FALSE], 1, min) > buffer
    list(test = f$test, train = f$train[keep])
  })
  list(random = random, block = blocks, buffered = buffered)
}

split_fingerprint <- function(folds) {
  # Exact textual indices, retained in provenance, avoid relying on a hash collision.
  paste(vapply(folds, function(f) paste(paste(f$test, collapse = ","),
                                       paste(f$train, collapse = ","), sep = "/"), ""), collapse = "|")
}

gaussian_scheme <- function(xy, target, folds, range, signal,
                             family = "exponential", axis_ratio = 1, angle = 0,
                             diagnostic_axis_ratio = axis_ratio) {
  if (any(vapply(folds, function(f) length(f$train) < 2 || length(f$test) < 1, logical(1))))
    stop("Infeasible fold: fewer than two training observations or no test observations")
  n <- nrow(xy)
  R <- correlation_matrix(xy, range = range, family = family, axis_ratio = axis_ratio, angle = angle)
  Rq <- correlation_matrix(target, xy, range, family, axis_ratio, angle)
  S <- signal * R + diag(1 - signal, n)
  diag(S) <- 1
  ep <- variance_reduction(S, signal * Rq)
  diagnostic_xy <- xy; diagnostic_target <- target
  if (diagnostic_axis_ratio != 1) {
    rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2)
    diagnostic_xy <- xy %*% rotation; diagnostic_target <- target %*% rotation
    diagnostic_xy[, 2] <- diagnostic_xy[, 2] / diagnostic_axis_ratio
    diagnostic_target[, 2] <- diagnostic_target[, 2] / diagnostic_axis_ratio
  }
  D <- pair_distance(diagnostic_xy)
  Dq <- pair_distance(diagnostic_target, diagnostic_xy)
  dp <- apply(Dq, 1, min)
  ev <- dist <- rep(NA_real_, n)
  for (f in folds) {
    ev[f$test] <- variance_reduction(S[f$train, f$train, drop = FALSE],
                                     S[f$test, f$train, drop = FALSE])
    dist[f$test] <- apply(D[f$test, f$train, drop = FALSE], 1, min)
  }
  stopifnot(all(is.finite(ev)), all(is.finite(dist)), all(ep >= -1e-9), all(ep <= 1 + 1e-9),
            all(ev >= -1e-9), all(ev <= 1 + 1e-9))
  cfun <- function(h) correlation_at(h, range, signal, family)
  di <- distance_summary(dist, dp, cfun)
  c(di, exact = mean(ev) - mean(ep), risk_target = 1 - mean(ep), risk_cv = 1 - mean(ev),
    remainder = mean(ev) - mean(ep) - di[["J2"]],
    min_train = min(vapply(folds, function(f) length(f$train), integer(1))),
    mean_train = mean(vapply(folds, function(f) length(f$train), integer(1))))
}

save_run <- function(object, stem, config = NULL) {
  src <- list.files(file.path(revision_root, "R"), pattern = "[.]R$", full.names = TRUE)
  metadata <- list(time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                   config = config, source_md5 = tools::md5sum(src), session = sessionInfo())
  saveRDS(list(results = object, metadata = metadata), file.path(revision_root, "results", paste0(stem, ".rds")))
  if (is.data.frame(object)) write.csv(object, file.path(revision_root, "results", paste0(stem, ".csv")), row.names = FALSE)
  invisible(metadata)
}
