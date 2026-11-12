# Post-hoc fold-count sensitivity; declared before outputs in
# revision/audit/fold-count-sensitivity-protocol.json. Primary studies are unchanged.
source("revision/R/local_libraries.R")
source("revision/R/common.R")
stopifnot(requireNamespace("spLeakage", quietly = TRUE), requireNamespace("CAST", quietly = TRUE))
pilot <- "--pilot" %in% commandArgs(trailingOnly = TRUE)
grid <- expand.grid(range = c(.05, .2), signal = c(.3, .7, .95),
  design = c("random", "clustered"), target = c("domain", "shifted"),
  replicate = seq_len(if (pilot) 1L else 80L), stringsAsFactors = FALSE)
grid$n <- if (pilot) 80L else 160L; grid$id <- seq_len(nrow(grid))
grid$seed <- (if (pilot) 1200000L else 1210000L) + grid$id
grid$response_seed <- grid$seed + 10000000L
stem <- if (pilot) "nested-k10-pilot" else "nested-estimation-k10"
workers <- if (pilot) 1L else as.integer(Sys.getenv("SLURM_NTASKS", "2"))
write.csv(grid, file.path("revision/results", paste0(stem, "-configuration.csv")), row.names = FALSE)
source_snapshot <- tools::md5sum(c("revision/exploratory/run_nested_k10.R", "revision/R/common.R", "revision/R/local_libraries.R", "revision/audit/fold-count-sensitivity-protocol.json"))
package_snapshot <- tools::md5sum(list.files("package/R", full.names = TRUE, pattern = "[.]R$"))
fitter_body <- paste(deparse(body(get(".fit_exponential", asNamespace("spLeakage")))), collapse = "\n")

run_one <- function(i) {
  cfg <- grid[i, ]; set.seed(cfg$seed)
  xy <- sampling_locations(cfg$n, cfg$design)
  target <- target_locations(cfg$target, if (pilot) 8L else 15L)
  candidates <- make_candidates(xy, k = 10L)
  points <- function(x) sf::st_as_sf(as.data.frame(x), coords = 1:2, crs = 3857)
  km <- suppressMessages(CAST::knndm(points(xy), predpoints = points(target), k = 10L))
  candidates$knndm <- Map(function(te, tr) list(test = te, train = tr), km$indx_test, km$indx_train)
  for (folds in candidates) {
    stopifnot(identical(sort(unlist(lapply(folds, `[[`, "test"))), seq_len(cfg$n)))
    if (any(vapply(folds, function(f) length(f$train) < 2L, logical(1)))) stop("Infeasible fold")
  }
  dd <- pair_distance(xy); dq <- pair_distance(target, xy)
  Sigma <- cfg$signal * exp(-dd / cfg$range) + diag(1 - cfg$signal, cfg$n)
  diag(Sigma) <- 1
  Bq <- cfg$signal * exp(-dq / cfg$range)
  # Isolate response simulation from RNG resets inside fold-construction code.
  set.seed(cfg$response_seed)
  y <- as.numeric(t(chol(Sigma)) %*% rnorm(cfg$n))
  true_target_mean <- kriging_prediction(y, Sigma, Bq)
  true_target_variance <- 1 - variance_reduction(Sigma, Bq)
  stopifnot(min(true_target_variance) >= -1e-8)
  true_target_variance <- pmax(true_target_variance, 0)
  fit_records <- list()
  fit_on <- function(tr, name) {
    warnings <- character()
    fitted <- withCallingHandlers(spLeakage::estimate_dependence(
      data.frame(x = xy[tr, 1], y = xy[tr, 2], response = y[tr]), "response", coords = c("x", "y")),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    stopifnot(is.finite(fitted$range), fitted$range > 0, fitted$nugget >= 0,
              fitted$psill >= 0, fitted$nugget + fitted$psill > 0)
    fit_records[[length(fit_records) + 1L]] <<- data.frame(id = cfg$id, name, n = length(tr),
      range = fitted$range, nugget = fitted$nugget, psill = fitted$psill,
      signal = fitted$signal_prop, warnings = paste(unique(warnings), collapse = " | "))
    fitted
  }
  predict_fitted <- function(tr, cross_distances, fit) {
    S <- fit$psill * exp(-dd[tr, tr, drop = FALSE] / fit$range) + diag(fit$nugget, length(tr))
    diag(S) <- fit$psill + fit$nugget
    b <- fit$psill * exp(-cross_distances / fit$range)
    kriging_prediction(y[tr], S, b)
  }
  full_fit <- fit_on(seq_len(cfg$n), "full")
  fitted_target_prediction <- predict_fitted(seq_len(cfg$n), dq, full_fit)
  target_risks <- c(known = mean(true_target_variance),
    estimated = mean((fitted_target_prediction - true_target_mean)^2 + true_target_variance))
  stopifnot(all(is.finite(target_risks)), target_risks[["estimated"]] >= target_risks[["known"]] - 1e-10)
  dp <- apply(dq, 1, min)
  rows <- list()
  for (scheme in names(candidates)) {
    folds <- candidates[[scheme]]
    known_prediction <- estimated_prediction <- dc <- numeric(cfg$n)
    expected_known_cv <- numeric(cfg$n)
    for (fi in seq_along(folds)) {
      f <- folds[[fi]]
      S <- Sigma[f$train, f$train, drop = FALSE]
      b <- Sigma[f$test, f$train, drop = FALSE]
      known_prediction[f$test] <- kriging_prediction(y[f$train], S, b)
      expected_known_cv[f$test] <- 1 - variance_reduction(S, b)
      fold_fit <- fit_on(f$train, paste(scheme, fi, sep = ":"))
      cross_d <- dd[f$test, f$train, drop = FALSE]
      estimated_prediction[f$test] <- predict_fitted(f$train, cross_d, fold_fit)
      dc[f$test] <- apply(cross_d, 1, min)
    }
    known_diagnostic <- distance_summary(dc, dp, function(h) cfg$signal * exp(-h / cfg$range))
    fitted_diagnostic <- distance_summary(dc, dp, full_fit$rho)
    stopifnot(all(is.finite(known_prediction)), all(is.finite(estimated_prediction)))
    cv_scores <- c(known = mean((y - known_prediction)^2), estimated = mean((y - estimated_prediction)^2))
    for (learner in names(target_risks)) {
      rows[[length(rows) + 1L]] <- data.frame(cfg, learner, scheme,
        conditional_target_risk = target_risks[[learner]], cv_score = cv_scores[[learner]],
        score_error = target_risks[[learner]] - cv_scores[[learner]],
        known_expected_cv = mean(expected_known_cv),
        known_exact_gap = target_risks[["known"]] - mean(expected_known_cv),
        A = known_diagnostic[["A"]], W = known_diagnostic[["W"]],
        J1_known = known_diagnostic[["J1"]], J2_known = known_diagnostic[["J2"]],
        U1_known = known_diagnostic[["U1"]], U2_known = known_diagnostic[["U2"]],
        J1_fitted = fitted_diagnostic[["J1"]], J2_fitted = fitted_diagnostic[["J2"]],
        U1_fitted = fitted_diagnostic[["U1"]], U2_fitted = fitted_diagnostic[["U2"]],
        fitted_range = full_fit$range, fitted_signal = full_fit$signal_prop,
        fitted_variance = full_fit$psill + full_fit$nugget,
        min_train = min(vapply(folds, function(f) length(f$train), integer(1))),
        split = split_fingerprint(folds), stringsAsFactors = FALSE)
    }
  }
  list(results = do.call(rbind, rows), fits = do.call(rbind, fit_records))
}

started <- proc.time()[[3]]
batch_dir <- file.path("revision/results", paste0(stem, "-batches"))
dir.create(batch_dir, showWarnings = FALSE)
batches <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 24L))
for (b in seq_along(batches)) {
  path <- file.path(batch_dir, sprintf("batch-%03d.rds", b))
  if (file.exists(path)) next
  out <- parallel::mclapply(batches[[b]], function(i) tryCatch(run_one(i),
    error = function(e) list(error = conditionMessage(e), configuration = grid[i, ])),
    mc.cores = workers, mc.set.seed = FALSE)
  saveRDS(out, path)
  cat("Completed batch", b, "of", length(batches), "at", format(Sys.time(), tz = "UTC"), "\n")
  flush.console()
}
all <- unlist(lapply(list.files(batch_dir, full.names = TRUE, pattern = "[.]rds$"), readRDS), recursive = FALSE)
errors <- Filter(function(x) !is.null(x$error), all)
good <- Filter(function(x) is.null(x$error), all)
results <- do.call(rbind, lapply(good, `[[`, "results"))
fits <- do.call(rbind, lapply(good, `[[`, "fits"))
save_run(results, stem, grid)
bundle <- readRDS(file.path("revision/results", paste0(stem, ".rds")))
bundle$fits <- fits
bundle$metadata$executed_source_md5 <- source_snapshot
bundle$metadata$package_source_md5 <- package_snapshot
bundle$metadata$fitter_body <- fitter_body
bundle$metadata$risk_reference <- "Exact conditional target risk given full training responses under the true Gaussian model; no target outcomes used."
saveRDS(bundle, file.path("revision/results", paste0(stem, ".rds")))
saveRDS(errors, file.path("revision/results", paste0(stem, "-errors.rds")))
cat("Completed", nrow(results), "evaluations from", length(good), "configurations;", length(errors),
    "failed configurations;", nrow(fits), "fits; elapsed seconds", proc.time()[[3]] - started, "\n")
if (length(errors)) { print(errors); stop("Nested study contains failed configurations; inspect before analysis.") }
