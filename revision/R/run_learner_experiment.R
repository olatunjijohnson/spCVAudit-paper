source("revision/R/local_libraries.R")
source("revision/R/common.R")
stopifnot(requireNamespace("CAST", quietly = TRUE), requireNamespace("ranger", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
pilot <- "--pilot" %in% args
grid <- expand.grid(range = c(0.05, 0.2), design = c("random", "clustered"),
                     target = c("domain", "shifted"),
                     replicate = seq_len(if (pilot) 1L else 60L), stringsAsFactors = FALSE)
grid$id <- seq_len(nrow(grid))
grid$seed <- (if (pilot) 510000L else 610000L) + grid$id
grid$n <- if (pilot) 100L else 200L
stem <- if (pilot) "learners-pilot" else "learners"
stem <- paste0(stem, Sys.getenv("REVISION_RUN_SUFFIX", ""))
workers <- if (pilot) 1L else as.integer(Sys.getenv("SLURM_NTASKS", "2"))
stopifnot(is.finite(workers), workers >= 1L)
source_snapshot <- tools::md5sum(c("revision/R/common.R", "revision/R/local_libraries.R",
                                  "revision/R/run_learner_experiment.R"))
write.csv(grid, file.path("revision/results", paste0(stem, "-configuration.csv")), row.names = FALSE)

run_one <- function(i) {
  cfg <- grid[i, ]
  set.seed(cfg$seed)
  xy <- sampling_locations(cfg$n, cfg$design)
  target <- target_locations(cfg$target, if (pilot) 12L else 20L)
  allxy <- rbind(xy, target); nall <- nrow(allxy)
  R <- correlation_matrix(allxy, range = cfg$range)
  L <- chol(R + diag(1e-10, nall))
  fields <- t(L) %*% matrix(rnorm(nall * 2), nall, 2)
  X <- fields[, 1]; Z <- fields[, 2]
  epsilon <- rnorm(cfg$n, sd = sqrt(0.1))
  candidates <- make_candidates(xy, k = 5, buffer = 0.05)
  as_points <- function(z) sf::st_as_sf(as.data.frame(z), coords = 1:2, crs = 3857)
  km <- suppressMessages(CAST::knndm(as_points(xy), predpoints = as_points(target), k = 5))
  candidates$knndm <- Map(function(te, tr) list(test = te, train = tr), km$indx_test, km$indx_train)
  D <- pair_distance(allxy)
  obs <- seq_len(cfg$n); tar <- seq.int(cfg$n + 1, nall)
  target_distance <- apply(D[tar, obs, drop = FALSE], 1, min)
  cv_dist <- lapply(candidates, function(folds) {
    if (any(vapply(folds, function(f) length(f$train) < 2, logical(1)))) stop("Infeasible fold")
    d <- rep(NA_real_, cfg$n)
    for (f in folds) d[f$test] <- apply(D[f$test, f$train, drop = FALSE], 1, min)
    d
  })
  rows <- list()
  for (alpha in c(0, 0.3, 0.6, 0.9)) {
    residual <- max(0, 0.9 - alpha)
    true_surface <- sqrt(alpha) * X + sqrt(residual) * Z
    y <- true_surface[obs] + epsilon
    y_residual <- y - sqrt(alpha) * X[obs]
    learner_names <- c("known_mean", "linear_X", "forest_X", "idw_coordinates", "oracle_kriging")
    predict_one <- function(learner, tr, te, fit_seed) {
      switch(learner,
        known_mean = sqrt(alpha) * X[te],
        linear_X = {
          coefficients <- lm.fit(cbind(1, X[tr]), y[tr])$coefficients
          as.numeric(cbind(1, X[te]) %*% coefficients)
        },
        forest_X = {
          fit <- ranger::ranger(Y ~ X, data = data.frame(Y = y[tr], X = X[tr]),
            num.trees = if (pilot) 100L else 300L, min.node.size = 5,
            seed = fit_seed, num.threads = 1, write.forest = TRUE)
          as.numeric(predict(fit, data.frame(X = X[te]), num.threads = 1)$predictions)
        },
        idw_coordinates = {
          weights <- 1 / pmax(D[te, tr, drop = FALSE], 1e-12)^2
          as.numeric(weights %*% y[tr] / rowSums(weights))
        },
        oracle_kriging = {
          S <- residual * R[tr, tr, drop = FALSE] + diag(0.1, length(tr))
          sqrt(alpha) * X[te] + kriging_prediction(y_residual[tr], S, residual * R[te, tr, drop = FALSE])
        }, stop("Unknown learner"))
    }
    for (li in seq_along(learner_names)) {
      learner <- learner_names[li]
      target_prediction <- predict_one(learner, obs, tar, cfg$seed + li * 1000L)
      target_loss <- mean((target_prediction - true_surface[tar])^2) + 0.1
      for (si in seq_along(candidates)) {
        scheme <- names(candidates)[si]; folds <- candidates[[si]]
        cv_prediction <- rep(NA_real_, cfg$n)
        for (fi in seq_along(folds)) {
          f <- folds[[fi]]
          cv_prediction[f$test] <- predict_one(learner, f$train, f$test,
                                                cfg$seed + li * 1000L + si * 100L + fi)
        }
        stopifnot(all(is.finite(cv_prediction)), all(is.finite(target_prediction)))
        cv_loss <- mean((y - cv_prediction)^2)
        marginal <- distance_summary(cv_dist[[scheme]], target_distance,
                                    function(h) 0.9 * exp(-h / cfg$range))
        conditional <- distance_summary(cv_dist[[scheme]], target_distance,
                              function(h) residual / (residual + 0.1) * exp(-h / cfg$range))
        rows[[length(rows) + 1L]] <- data.frame(cfg, alpha = alpha, learner = learner,
          scheme = scheme, cv_loss = cv_loss, target_loss = target_loss,
          gap = target_loss - cv_loss, A = marginal[["A"]], W = marginal[["W"]],
          J1_marginal = marginal[["J1"]], J2_marginal = marginal[["J2"]],
          U1_marginal = marginal[["U1"]], U2_marginal = marginal[["U2"]],
          J1_conditional = conditional[["J1"]], J2_conditional = conditional[["J2"]],
          U1_conditional = conditional[["U1"]], U2_conditional = conditional[["U2"]],
          one_neighbour_conditional_gap = (1 - alpha) * conditional[["J2"]],
          crosses = marginal[["crosses"]],
          mean_train = mean(vapply(folds, function(f) length(f$train), integer(1))),
          split = split_fingerprint(folds), stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

started <- proc.time()[[3]]
batch_dir <- file.path("revision/results", paste0(stem, "-batches"))
dir.create(batch_dir, showWarnings = FALSE)
batches <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 16L))
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
save_run(out, stem, config = grid)
bundle <- readRDS(file.path("revision/results", paste0(stem, ".rds")))
bundle$metadata$executed_source_md5 <- source_snapshot
bundle$metadata$target_evaluation <- "finite target latent squared loss plus independent noise variance 0.1"
saveRDS(bundle, file.path("revision/results", paste0(stem, ".rds")))
saveRDS(errors, file.path("revision/results", paste0(stem, "-errors.rds")))
cat("Completed", nrow(out), "learner/scheme/alpha evaluations from", nrow(grid),
    "field configurations with", length(errors), "failures in", proc.time()[[3]] - started, "seconds.\n")
if (length(errors)) print(errors)
