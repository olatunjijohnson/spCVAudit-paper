source("revision/R/local_libraries.R")
source("revision/R/common.R")
stopifnot(requireNamespace("CAST", quietly = TRUE), requireNamespace("sf", quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
pilot <- "--pilot" %in% args
grid <- read.csv(if (pilot) "revision/results/gp-pilot-configuration.csv" else
                   "revision/results/gp-reference-configuration.csv", stringsAsFactors = FALSE)
stem <- if (pilot) "gp-knndm-pilot" else "gp-knndm"
stem <- paste0(stem, Sys.getenv("REVISION_RUN_SUFFIX", ""))
workers <- if (pilot) 1L else as.integer(Sys.getenv("SLURM_NTASKS", "2"))
stopifnot(is.finite(workers), workers >= 1L)
started <- proc.time()[[3]]
source_snapshot <- tools::md5sum(c("revision/R/common.R", "revision/R/local_libraries.R", "revision/R/run_gp_knndm.R"))

run_one <- function(i) {
  cfg <- grid[i, ]
  set.seed(cfg$seed)
  xy <- sampling_locations(cfg$n, cfg$design)
  target <- target_locations(cfg$target, side = if (pilot) 11L else 15L)
  # Consume exactly the same candidate-generation random draws as the paired run.
  ignored <- make_candidates(xy, k = 5, buffer = 0.05)
  as_points <- function(z) sf::st_as_sf(as.data.frame(z), coords = 1:2, crs = 3857)
  warnings <- character()
  val <- tryCatch(withCallingHandlers({
    fit <- CAST::knndm(as_points(xy), predpoints = as_points(target), k = 5)
    folds <- Map(function(te, tr) list(test = te, train = tr), fit$indx_test, fit$indx_train)
    score <- gaussian_scheme(xy, target, folds, cfg$range, cfg$signal, cfg$family)
    # Verify that our distribution calculation agrees with the package statistic.
    stopifnot(abs(score[["W"]] - fit$W) < 1e-7)
    data.frame(cfg, scheme = "knndm", as.list(score), CAST_W = fit$W,
                q_clusters = fit$q, split = split_fingerprint(folds), stringsAsFactors = FALSE)
  }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
  error = function(e) list(error = conditionMessage(e), config = cfg))
  if (is.data.frame(val)) val$warnings <- paste(unique(warnings), collapse = " | ")
  val
}

# Save each batch so an interruption cannot destroy all completed evaluations.
batches <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 128))
batch_dir <- file.path("revision/results", paste0(stem, "-batches"))
dir.create(batch_dir, showWarnings = FALSE)
for (b in seq_along(batches)) {
  destination <- file.path(batch_dir, sprintf("batch-%03d.rds", b))
  if (file.exists(destination)) next
  result <- parallel::mclapply(batches[[b]], run_one, mc.cores = workers,
                                mc.set.seed = FALSE)
  saveRDS(result, destination)
  cat("Completed batch", b, "of", length(batches), "at", format(Sys.time(), tz = "UTC"), "\n")
  flush.console()
}
records <- unlist(lapply(list.files(batch_dir, full.names = TRUE, pattern = "[.]rds$"), readRDS), recursive = FALSE)
good <- Filter(is.data.frame, records)
errors <- Filter(function(x) !is.data.frame(x), records)
out <- do.call(rbind, good)
save_run(out, stem, config = grid)
bundle <- readRDS(file.path("revision/results", paste0(stem, ".rds")))
bundle$metadata$executed_source_md5 <- source_snapshot
bundle$metadata$CAST_version <- as.character(utils::packageVersion("CAST"))
saveRDS(bundle, file.path("revision/results", paste0(stem, ".rds")))
saveRDS(errors, file.path("revision/results", paste0(stem, "-errors.rds")))
cat("Completed", nrow(out), "evaluations with", length(errors), "failures in",
    proc.time()[[3]] - started, "seconds.\n")
if (length(errors)) print(errors)
