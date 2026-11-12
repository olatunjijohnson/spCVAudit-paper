# This sensitivity study is separate from the pre-existing main Gaussian grid.
source("revision/R/local_libraries.R")
source("revision/R/common.R")
stopifnot(requireNamespace("CAST", quietly = TRUE), requireNamespace("spLeakage", quietly = TRUE))
pilot <- "--pilot" %in% commandArgs(trailingOnly = TRUE)
grid <- expand.grid(range = c(.05, .2), signal = c(.3, .7, .95),
  design = c("random", "clustered"), target = c("domain", "shifted"),
  axis_ratio = c(1, .25), replicate = seq_len(if (pilot) 1L else 30L), stringsAsFactors = FALSE)
grid$n <- if (pilot) 80L else 160L; grid$id <- seq_len(nrow(grid))
grid$seed <- (if (pilot) 710000L else 810000L) + grid$id
grid$response_seed <- grid$seed + 10000000L
stem <- if (pilot) "sensitivity-pilot" else "sensitivity"
stem <- paste0(stem, Sys.getenv("REVISION_RUN_SUFFIX", ""))
write.csv(grid, file.path("revision/results", paste0(stem, "-configuration.csv")), row.names = FALSE)
workers <- if (pilot) 1L else as.integer(Sys.getenv("SLURM_NTASKS", "2"))
source_files <- c("revision/R/run_sensitivity.R", "revision/R/common.R", "revision/R/local_libraries.R")
source_snapshot <- tools::md5sum(source_files)
package_snapshot <- tools::md5sum(list.files("package/R", full.names = TRUE, pattern = "[.]R$"))
dependence_implementation <- paste(deparse(body(get(".fit_exponential", asNamespace("spLeakage")))), collapse = "\n")

run_one <- function(i) {
  cfg <- grid[i, ]; set.seed(cfg$seed)
  xy <- sampling_locations(cfg$n, cfg$design)
  target <- target_locations(cfg$target, if (pilot) 8L else 15L)
  candidates <- make_candidates(xy)
  points <- function(x) sf::st_as_sf(as.data.frame(x), coords = 1:2, crs = 3857)
  km <- suppressMessages(CAST::knndm(points(xy), predpoints = points(target), k = 5))
  candidates$knndm <- Map(function(te, tr) list(test = te, train = tr), km$indx_test, km$indx_train)
  angle <- pi / 6
  R <- correlation_matrix(xy, range = cfg$range, axis_ratio = cfg$axis_ratio, angle = angle)
  Bq <- cfg$signal * correlation_matrix(target, xy, cfg$range, axis_ratio = cfg$axis_ratio, angle = angle)
  S <- cfg$signal * R + diag(1 - cfg$signal, cfg$n)
  # CAST can reset the global RNG internally. Response innovations use their
  # own declared seed, independently of candidate-construction random draws.
  set.seed(cfg$response_seed)
  y <- as.numeric(t(chol(S)) %*% rnorm(cfg$n))
  fit_warnings <- character()
  fit <- tryCatch(withCallingHandlers(spLeakage::estimate_dependence(
    data.frame(x = xy[, 1], y = xy[, 2], z = y), "z", coords = c("x", "y")),
    warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) list(error = conditionMessage(e)))
  fitted_ok <- is.null(fit$error) && all(is.finite(c(fit$psill, fit$nugget, fit$range))) &&
    fit$range > 0 && fit$psill + fit$nugget > 0
  distances <- pair_distance(xy); target_distances <- pair_distance(target, xy)
  eq <- variance_reduction(S, Bq)
  nearest_and_bound <- function(training_covariance, cross_covariance) {
    b1 <- apply(cross_covariance, 1, max)
    nn <- max.col(cross_covariance, ties.method = "first")
    bound <- vapply(seq_len(nrow(cross_covariance)), function(j) {
      remaining <- setdiff(seq_len(ncol(cross_covariance)), nn[j])
      residual <- cross_covariance[j, remaining] - training_covariance[remaining, nn[j]] * b1[j]
      min(sum(residual^2) / (1 - cfg$signal), cfg$signal - b1[j]^2)
    }, 0.0)
    list(one = b1^2, bound = pmax(bound, 0))
  }
  target_extra <- nearest_and_bound(S, Bq)
  dp <- apply(target_distances, 1, min)
  rows <- list()
  for (scheme in names(candidates)) {
    folds <- candidates[[scheme]]
    if (any(vapply(folds, function(f) length(f$train) < 2L, logical(1)))) stop("Infeasible fold")
    ev <- dc <- one_cv <- bound_cv <- numeric(cfg$n)
    matched_target <- 0
    for (f in folds) {
      ss <- S[f$train, f$train, drop = FALSE]
      bb <- S[f$test, f$train, drop = FALSE]
      ev[f$test] <- variance_reduction(ss, bb)
      matched_target <- matched_target + length(f$test) / cfg$n * mean(variance_reduction(ss, Bq[, f$train, drop = FALSE]))
      dc[f$test] <- apply(distances[f$test, f$train, drop = FALSE], 1, min)
      extra <- nearest_and_bound(ss, bb)
      one_cv[f$test] <- extra$one; bound_cv[f$test] <- extra$bound
    }
    exact <- mean(ev) - mean(eq)
    matched <- mean(ev) - matched_target
    training_effect <- matched_target - mean(eq)
    true_J2 <- mean(one_cv) - mean(target_extra$one)
    lower <- true_J2 - mean(target_extra$bound)
    upper <- true_J2 + mean(bound_cv)
    stopifnot(abs(exact - matched - training_effect) < 1e-8, training_effect <= 1e-8,
              exact >= lower - 1e-8, exact <= upper + 1e-8)
    for (model in c("known_isotropic", "estimated_isotropic")) {
      available <- model == "known_isotropic" || fitted_ok
      cfun <- if (model == "known_isotropic") function(h) cfg$signal * exp(-h / cfg$range) else if (fitted_ok) fit$rho else NULL
      di <- if (available) distance_summary(dc, dp, cfun) else setNames(rep(NA_real_, 9),
        c("A", "W", "J1", "J2", "U1", "U2", "Aplus", "Aminus", "crosses"))
      fitted_variance <- if (model == "known_isotropic") 1 else if (fitted_ok) fit$psill + fit$nugget else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(cfg, scheme, diagnostic_model = model,
        as.list(di), exact, matched_size_gap = matched, training_size_effect = training_effect,
        true_J2, bound_lower = lower, bound_upper = upper, proxy = fitted_variance * di[["J2"]],
        fitted_ok = available, fitted_range = if (fitted_ok) fit$range else NA_real_,
        fitted_signal = if (fitted_ok) fit$signal_prop else NA_real_,
        fitted_variance, fit_error = if (is.null(fit$error)) "" else fit$error,
        fit_warnings = paste(unique(fit_warnings), collapse = " | "),
        split = split_fingerprint(folds), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

started <- proc.time()[[3]]
batch_dir <- file.path("revision/results", paste0(stem, "-batches"))
dir.create(batch_dir, showWarnings = FALSE)
batches <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 48L))
for (b in seq_along(batches)) {
  path <- file.path(batch_dir, sprintf("batch-%03d.rds", b))
  if (file.exists(path)) next
  ans <- parallel::mclapply(batches[[b]], function(i) tryCatch(run_one(i),
    error = function(e) list(error = conditionMessage(e), config = grid[i, ])),
    mc.cores = workers, mc.set.seed = FALSE)
  saveRDS(ans, path)
  cat("Completed batch", b, "of", length(batches), "at", format(Sys.time(), tz = "UTC"), "\n")
  flush.console()
}
records <- unlist(lapply(list.files(batch_dir, full.names = TRUE, pattern = "[.]rds$"), readRDS), recursive = FALSE)
errors <- Filter(function(x) !is.data.frame(x), records)
out <- do.call(rbind, Filter(is.data.frame, records))
save_run(out, stem, grid)
bundle <- readRDS(file.path("revision/results", paste0(stem, ".rds")))
bundle$metadata$executed_source_md5 <- source_snapshot
bundle$metadata$package_source_md5 <- package_snapshot
bundle$metadata$dependence_implementation <- dependence_implementation
bundle$metadata$package_version <- as.character(utils::packageVersion("spLeakage"))
bundle$metadata$estimation_scope <- "One field at training locations; exact risk labels refer to fixed geometry-only schemes, not risk conditional on the fitted field."
saveRDS(bundle, file.path("revision/results", paste0(stem, ".rds")))
saveRDS(errors, file.path("revision/results", paste0(stem, "-errors.rds")))
cat("Completed", nrow(out), "rows with", length(errors), "failed configurations in", proc.time()[[3]] - started, "seconds.\n")
if (length(errors)) print(errors)
if (length(errors)) stop("Sensitivity study contains failed configurations; inspect before summarizing.")
