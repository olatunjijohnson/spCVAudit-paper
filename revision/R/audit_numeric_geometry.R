# Check the retained study domains against direct coordinate differences.
# The Gram-matrix distance identity can lose relative precision for extremely
# close points; this audit measures its absolute error in the actual experiments.
source("revision/R/local_libraries.R")
source("revision/R/common.R")
pkgload::load_all("package", quiet=TRUE)
direct_distance <- function(x, y = x) {
  total <- matrix(0, nrow(x), nrow(y))
  for (j in seq_len(ncol(x))) total <- total + outer(x[, j], y[, j], "-")^2
  sqrt(total)
}
compare <- function(x, q) {
  d <- pair_distance(x); truth <- direct_distance(x)
  dq <- pair_distance(q, x); tq <- direct_distance(q, x)
  self <- max(abs(diag(d)))
  diagonal_error <- max(abs(diag(d) - diag(truth)))
  max_error <- max(abs(d - truth), abs(dq - tq))
  diag(d) <- diag(truth) <- Inf
  ni <- max.col(-d, ties.method = "first")
  qi <- max.col(-dq, ties.method = "first")
  # Regular target grids have exact nearest-neighbour ties. Compare the
  # distance attained by the selected neighbour, not only its tied index.
  training_excess <- truth[cbind(seq_len(nrow(d)), ni)] - apply(truth, 1, min)
  target_excess <- tq[cbind(seq_len(nrow(dq)), qi)] - apply(tq, 1, min)
  c(max_absolute_distance_error = max_error, maximum_self_distance = self,
    diagonal_error = diagonal_error,
    minimum_distinct_training_distance = min(truth),
    training_nearest_index_changes = sum(max.col(-d, ties.method = "first") != max.col(-truth, ties.method = "first")),
    target_nearest_index_changes = sum(qi != max.col(-tq, ties.method = "first")),
    maximum_training_nearest_distance_excess = max(training_excess),
    maximum_target_nearest_distance_excess = max(target_excess))
}
workers <- as.integer(Sys.getenv("SLURM_NTASKS", "2"))
records <- list()
for (study in c("gp-reference", "learners", "sensitivity-rng", "nested-estimation-rng")) {
  grid <- read.csv(file.path("revision/results", paste0(study, "-configuration.csv")))
  out <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
    cfg <- grid[i, ]; set.seed(cfg$seed)
    x <- sampling_locations(cfg$n, cfg$design)
    q <- target_locations(cfg$target, if (study == "learners") 20L else 15L)
    if ("axis_ratio" %in% names(cfg) && cfg$axis_ratio != 1) {
      angle <- pi / 6
      rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2)
      x <- x %*% rotation; q <- q %*% rotation
      x[, 2] <- x[, 2] / cfg$axis_ratio; q[, 2] <- q[, 2] / cfg$axis_ratio
    }
    # The learner GP is sampled over training and target jointly.
    values <- if (study == "learners") compare(rbind(x, q), q) else compare(x, q)
    data.frame(study, id = cfg$id, as.list(values))
  }, mc.cores = workers, mc.set.seed = FALSE)
  records[[study]] <- do.call(rbind, out)
}
app <- readRDS("revision/results/application-jura.rds")
details <- Filter(function(d) d$results$repeat_id[1] == 1L, app$details)
detail <- Filter(function(d) d$results$target[1] == "published", details)[[1]]
records$application <- data.frame(study = "application-jura", id = 1L,
  as.list(compare(detail$geometry$xy, detail$geometry$target_xy)))
out <- do.call(rbind, records)
write.csv(out, "revision/results/numeric-geometry-audit.csv", row.names = FALSE)
summary <- aggregate(out[, 3:ncol(out)], out["study"], max)
summary$minimum_distinct_training_distance <- vapply(summary$study,
  function(study) min(out$minimum_distinct_training_distance[out$study == study]), 0.0)
write.csv(summary, "revision/results/numeric-geometry-audit-summary.csv", row.names = FALSE)
print(summary)
stopifnot(max(out$max_absolute_distance_error) < 1e-7,
          max(out$maximum_self_distance) < 1e-12,
          max(out$maximum_training_nearest_distance_excess) < 1e-10,
          max(out$maximum_target_nearest_distance_excess) < 1e-10)

# Verify the public API against saved Jura geometry and both covariance choices.
rows <- list()
for (d in details) for (model in c("ordinary", "trend")) {
  dep <- d$full_models[[model]]$dependence; geo <- d$geometry
  for (scheme in names(geo$folds)) {
    api <- spLeakage::detect_leakage(geo$xy, geo$folds[[scheme]],
      spLeakage::prediction_target(grid = geo$target_xy), dependence = dep)
    ref <- distance_summary(geo$dc[[scheme]], geo$dp, dep$rho)
    mapping <- c(A = "A", W = "W", SLI_rho = "J1", SLI_rho2 = "J2", W_rho = "U1")
    for (name in names(mapping)) rows[[length(rows) + 1L]] <- data.frame(
      target = d$results$target[1], dependence = model, scheme, quantity = name,
      package_value = api[[name]], independent_value = ref[[mapping[[name]]]],
      absolute_difference = abs(api[[name]] - ref[[mapping[[name]]]]))
  }
}
parity <- do.call(rbind, rows)
write.csv(parity, "revision/results/application-api-parity.csv", row.names = FALSE)
stopifnot(max(parity$absolute_difference) < 1e-7)
cat("Public application API maximum discrepancy:", max(parity$absolute_difference), "\n")
