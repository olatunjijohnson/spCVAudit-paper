# A quick environment check before an exact replay of the fitted-covariance studies.
# These three configurations exposed numerical-platform sensitivity during the
# public-release audit. They are regression checks, not exclusions or new settings.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
source("revision/R/local_libraries.R")
source("revision/R/common.R")
reference <- readRDS(file.path(args[1], "revision/results/sensitivity-rng.rds"))
rows <- lapply(c(1L, 1077L, 1428L), function(id) {
  cfg <- reference$metadata$config[reference$metadata$config$id == id, ]
  stopifnot(nrow(cfg) == 1L)
  set.seed(cfg$seed)
  xy <- sampling_locations(cfg$n, cfg$design)
  S <- cfg$signal * correlation_matrix(xy, range = cfg$range,
      axis_ratio = cfg$axis_ratio, angle = pi / 6) + diag(1 - cfg$signal, cfg$n)
  set.seed(cfg$response_seed)
  response <- as.numeric(t(chol(S)) %*% rnorm(cfg$n))
  fit <- suppressWarnings(spLeakage::estimate_dependence(
    data.frame(x = xy[, 1], y = xy[, 2], z = response), "z", coords = c("x", "y")))
  saved <- subset(reference$results, id == cfg$id & diagnostic_model == "estimated_isotropic")[1, ]
  expected <- c(range = saved$fitted_range, signal = saved$fitted_signal,
                variance = saved$fitted_variance)
  actual <- c(range = fit$range, signal = fit$signal_prop, variance = fit$psill + fit$nugget)
  data.frame(id = id, parameter = names(expected), expected = unname(expected),
             actual = unname(actual), absolute_difference = unname(abs(expected - actual)))
})
result <- do.call(rbind, rows)
dir.create("revision/audit", recursive = TRUE, showWarnings = FALSE)
write.csv(result, "revision/audit/reference-environment-probe.csv", row.names = FALSE)
print(result)
if (any(!is.finite(result$absolute_difference)) || any(result$absolute_difference > 1e-10))
  stop("This numerical environment does not reproduce the reference covariance fits. ",
       "No study has been run. See REPRODUCING.md and revision/audit/reference-environment-probe.csv.",
       call. = FALSE)
cat("Reference fitting checks passed; full computation can proceed.\n")
