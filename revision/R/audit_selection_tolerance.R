# Post-hoc numerical sensitivity: the primary tolerance remains 1e-12.
# This audit does not alter diagnostics, response draws or primary selections.
z <- readRDS("revision/results/nested-estimation-rng.rds")$results
rows <- list()
for (learner in unique(z$learner)) {
  groups <- split(z[z$learner == learner, ], z$id[z$learner == learner])
  cells <- vapply(groups, function(d) paste(d$range[1], d$signal[1], d$design[1], d$target[1]), "")
  for (tolerance in c(0, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8)) {
    values <- vapply(groups, function(d) {
      j <- abs(d$J2_fitted); w <- abs(d$W)
      selected_j <- j - min(j) <= tolerance
      selected_w <- w - min(w) <= tolerance
      c(difference = mean(abs(d$score_error)[selected_j]) - mean(abs(d$score_error)[selected_w]),
        J2_tied = sum(selected_j) > 1L, W_tied = sum(selected_w) > 1L)
    }, numeric(3))
    cell_values <- split(values["difference", ], cells)
    rows[[length(rows) + 1L]] <- data.frame(learner, tolerance,
      mean_difference = mean(values["difference", ]),
      mcse = sqrt(sum(vapply(cell_values, function(x) length(x) * var(x), 0.0))) / ncol(values),
      J2_tie_fraction = mean(values["J2_tied", ]), W_tie_fraction = mean(values["W_tied", ]))
  }
}
out <- do.call(rbind, rows)
primary <- read.csv("revision/results/summaries/nested-paired-comparisons.csv")
primary <- subset(primary, criterion == "J2_fitted" & baseline == "W" & stratum == "all")
check <- subset(out, tolerance == 1e-12)
primary <- primary[match(check$learner, primary$learner), ]
stopifnot(max(abs(check$mean_difference - primary$mean_absolute_error_difference)) < 1e-12,
          max(abs(check$mcse - primary$mcse)) < 1e-12)
write.csv(out, "revision/audit/selection-tolerance.csv", row.names = FALSE)
print(out, row.names = FALSE)
cat("Numerical sensitivity only; primary tolerance and estimand remain unchanged.\n")
