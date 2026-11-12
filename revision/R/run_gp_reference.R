source("revision/R/common.R")
args <- commandArgs(trailingOnly = TRUE)
pilot <- "--pilot" %in% args
reps <- if (pilot) 1L else 40L
grid <- expand.grid(n = if (pilot) c(80L, 180L) else c(100L, 250L),
                    range = c(0.05, 0.2), signal = if (pilot) c(0.3, 0.9) else c(0, 0.3, 0.7, 0.95),
                    design = c("random", "clustered"), target = c("domain", "shifted"),
                    family = c("exponential", "matern32"), replicate = seq_len(reps),
                    stringsAsFactors = FALSE)
grid$id <- seq_len(nrow(grid))
grid$seed <- (if (pilot) 210000L else 310000L) + grid$id
stem <- if (pilot) "gp-pilot" else "gp-reference"
stem <- paste0(stem, Sys.getenv("REVISION_RUN_SUFFIX", ""))
write.csv(grid, file.path("revision/results", paste0(stem, "-configuration.csv")), row.names = FALSE)

run_one <- function(i) {
  cfg <- grid[i, ]
  set.seed(cfg$seed)
  xy <- sampling_locations(cfg$n, cfg$design)
  target <- target_locations(cfg$target, side = if (pilot) 11L else 15L)
  candidates <- make_candidates(xy, k = 5, buffer = 0.05)
  records <- lapply(names(candidates), function(scheme) {
    f <- candidates[[scheme]]
    started <- proc.time()[[3]]
    val <- tryCatch(gaussian_scheme(xy, target, f, cfg$range, cfg$signal, cfg$family),
                    error = function(e) e)
    if (inherits(val, "error")) return(list(error = conditionMessage(val), scheme = scheme, config = cfg))
    data.frame(cfg, scheme = scheme, as.list(val), seconds = proc.time()[[3]] - started,
               split = split_fingerprint(f), stringsAsFactors = FALSE, check.names = FALSE)
  })
  records
}

started <- proc.time()[[3]]
workers <- if (pilot) 1L else as.integer(Sys.getenv("SLURM_NTASKS", "2"))
stopifnot(is.finite(workers), workers >= 1L)
ans <- parallel::mclapply(seq_len(nrow(grid)), run_one, mc.cores = workers, mc.set.seed = FALSE)
records <- unlist(ans, recursive = FALSE)
errors <- Filter(function(z) !is.data.frame(z), records)
good <- Filter(is.data.frame, records)
out <- do.call(rbind, good)
save_run(out, stem, config = grid)
saveRDS(errors, file.path("revision/results", paste0(stem, "-errors.rds")))
cat("Configurations:", nrow(grid), "; successful scheme evaluations:", nrow(out),
    "; failures:", length(errors), "; elapsed seconds:", proc.time()[[3]] - started, "\n")
if (length(errors)) print(errors)
print(aggregate(cbind(exact, J1, J2, W) ~ design + target + scheme, out, mean), row.names = FALSE)
