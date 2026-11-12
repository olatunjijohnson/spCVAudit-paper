args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
stems <- if (length(args) > 2L) args[-c(1, 2)] else
  c("gp-reference", "gp-knndm", "learners", "sensitivity-rng", "application-jura", "nested-estimation-rng")
rows <- list()
for (stem in stems) {
  original <- readRDS(file.path(args[1], "revision/results", paste0(stem, ".rds")))
  replay <- readRDS(file.path(args[2], "revision/results", paste0(stem, ".rds")))
  a <- original$results
  b <- replay$results
  key <- intersect(c("id", "scheme", "alpha", "learner", "diagnostic_model", "dependence", "model"), names(a))
  ordering <- function(z) do.call(order, z[key])
  a <- a[ordering(a), ]; b <- b[ordering(b), ]
  stopifnot(nrow(a) == nrow(b), setequal(names(a), names(b)))
  for (column in setdiff(names(a), c("seconds", "fit_warnings"))) {
    x <- a[[column]]; y <- b[[column]]
    stopifnot(identical(is.na(x), is.na(y)))
    if (is.numeric(x)) {
      ok <- isTRUE(all.equal(x, y, tolerance = 1e-10, scale = 1, check.attributes = FALSE))
      error <- if (all(is.na(x))) 0 else max(abs(x - y), na.rm = TRUE)
    } else {
      ok <- identical(as.character(x), as.character(y)); error <- NA_real_
    }
    rows[[length(rows) + 1L]] <- data.frame(study = stem, column, n = length(x), equal = ok, max_absolute_difference = error)
  }
  err <- readRDS(file.path(args[2], "revision/results", paste0(stem, "-errors.rds")))
  stopifnot(length(err) == 0L)
  if (!is.null(original$fits)) {
    rows[[length(rows) + 1L]] <- data.frame(study = stem, column = "all_fit_records", n = nrow(original$fits),
      equal = isTRUE(all.equal(original$fits, replay$fits, tolerance = 1e-10, scale = 1, check.attributes = FALSE)),
      max_absolute_difference = NA_real_)
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(args[2], "revision/results/full-reproduction-comparison.csv"), row.names = FALSE)
print(aggregate(equal ~ study, out, all))
print(aggregate(max_absolute_difference ~ study, out, max, na.rm = TRUE))
stopifnot(all(out$equal))
