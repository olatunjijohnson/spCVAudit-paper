source("revision/R/common.R")
stem <- paste0("nested-estimation", Sys.getenv("REVISION_NESTED_SUFFIX", "-rng"))
bundle <- readRDS(file.path("revision/results", paste0(stem, ".rds")))
z <- bundle$results
expected <- read.csv(file.path("revision/results", paste0(stem, "-configuration.csv")))
errors <- readRDS(file.path("revision/results", paste0(stem, "-errors.rds")))
stopifnot(nrow(z) == 8L * length(unique(z$id)), all(table(z$id) == 8L),
  !anyDuplicated(z[c("id", "learner", "scheme")]),
  max(abs(z$score_error - z$conditional_target_risk + z$cv_score)) < 1e-8,
  all(is.finite(z$score_error)), all(z$cv_score >= 0), all(z$conditional_target_risk >= 0))
dir.create("revision/results/summaries", showWarnings = FALSE)
factor_columns <- c("range", "signal", "design", "target")
cell <- interaction(expected[factor_columns], drop = TRUE)
groups <- function(d, columns) split(seq_len(nrow(d)), interaction(d[columns], drop = TRUE, lex.order = TRUE))
fixed_mcse <- function(x, ids) {
  g <- split(x, cell[match(ids, expected$id)], drop = TRUE)
  if (any(lengths(g) < 2)) return(NA_real_)
  sqrt(sum(vapply(g, function(y) length(y) * var(y), 0.0))) / length(x)
}
success <- expected
success$completed <- success$id %in% z$id
failure <- aggregate(completed ~ range + signal + design + target, success,
                     function(x) c(attempted = length(x), successful = sum(x)))
write.csv(failure, "revision/results/summaries/nested-completion.csv", row.names = FALSE)
if (length(errors)) warning("Failure counts are nonzero; summaries apply to successful configurations and need explicit qualification.")

criteria <- c("A", "W", "J1_known", "J2_known", "U1_known", "U2_known",
              "J1_fitted", "J2_fitted", "U1_fitted", "U2_fitted", "fixed_random", "fixed_knndm", "oracle")
selection <- lapply(groups(z, c("id", "learner")), function(ii) {
  d <- z[ii, ]
  stopifnot(diff(range(d$conditional_target_risk)) < 1e-10)
  scores <- vapply(criteria, function(cr) switch(cr,
    fixed_random = as.numeric(d$scheme != "random"), fixed_knndm = as.numeric(d$scheme != "knndm"),
    oracle = abs(d$score_error), abs(d[[cr]])), numeric(4))
  selected <- sweep(scores, 2, apply(scores, 2, min), "-") <= 1e-12
  weights <- sweep(selected, 2, colSums(selected), "/")
  data.frame(d[rep(1, length(criteria)), c(names(expected), "learner"), drop = FALSE],
    criterion = criteria, score_error = colSums(weights * d$score_error),
    absolute_error = colSums(weights * abs(d$score_error)),
    squared_error = colSums(weights * d$score_error^2),
    known_reference_gap = colSums(weights * d$known_exact_gap),
    target_risk = d$conditional_target_risk[1], n_tied = colSums(selected),
    selected = vapply(seq_along(criteria), function(j) paste(sort(d$scheme[selected[, j]]), collapse = ";"), ""),
    row.names = NULL)
})
selection <- do.call(rbind, selection)
write.csv(selection, "revision/results/summaries/nested-selection-by-field.csv", row.names = FALSE)
controls <- subset(selection, learner == "known" & criterion %in%
  c("A", "W", "J1_known", "J2_known", "U1_known", "U2_known", "fixed_random", "fixed_knndm"))
calibration <- lapply(split(controls, controls$criterion), function(d) {
  delta <- d$score_error - d$known_reference_gap
  se <- fixed_mcse(delta, d$id)
  data.frame(criterion = d$criterion[1], n = nrow(d),
    mean_observed_minus_analytic_gap = mean(delta), mcse = se, z_score = mean(delta) / se)
})
calibration <- do.call(rbind, calibration)
write.csv(calibration, "revision/results/summaries/nested-known-risk-calibration.csv", row.names = FALSE)
if (any(abs(calibration$z_score) > 4.5)) warning("An analytic-risk calibration check exceeds 4.5 Monte Carlo standard errors; inspect it.")
summary <- list()
for (columns in list(c("learner", "criterion"), c("learner", "criterion", "design", "target"),
                     c("learner", "criterion", factor_columns))) {
  for (ii in groups(selection, columns)) {
    d <- selection[ii, ]; rmse <- sqrt(mean(d$squared_error))
    summary[[length(summary) + 1L]] <- data.frame(
      stratum = paste(paste(columns, d[1, columns], sep = "="), collapse = ";"),
      learner = d$learner[1], criterion = d$criterion[1], n = nrow(d),
      mean_signed_error = mean(d$score_error), signed_mcse = fixed_mcse(d$score_error, d$id),
      mean_absolute_error = mean(d$absolute_error), absolute_mcse = fixed_mcse(d$absolute_error, d$id),
      rmse, rmse_mcse = fixed_mcse(d$squared_error, d$id) / (2 * rmse),
      mean_target_risk = mean(d$target_risk), tie_fraction = mean(d$n_tied > 1))
  }
}
summary <- do.call(rbind, summary)
write.csv(summary, "revision/results/summaries/nested-selection.csv", row.names = FALSE)
paired <- list()
for (learner in unique(selection$learner)) for (baseline in c("W", "fixed_random", "fixed_knndm", "J2_known")) {
  b <- selection[selection$learner == learner & selection$criterion == baseline, ]
  for (criterion in setdiff(criteria, c(baseline, "oracle"))) {
    d <- selection[selection$learner == learner & selection$criterion == criterion, ]
    stopifnot(!anyDuplicated(b$id), !anyDuplicated(d$id), setequal(b$id, d$id))
    difference <- d$absolute_error - b$absolute_error[match(d$id, b$id)]
    for (columns in list(character(), c("design", "target"), factor_columns)) {
      rows <- if (length(columns)) groups(d, columns) else list(all = seq_len(nrow(d)))
      for (ii in rows) {
        se <- fixed_mcse(difference[ii], d$id[ii]); estimate <- mean(difference[ii])
        paired[[length(paired) + 1L]] <- data.frame(learner, criterion, baseline,
          stratum = if (length(columns)) paste(paste(columns, d[ii[1], columns], sep = "="), collapse = ";") else "all",
          n = length(ii), mean_absolute_error_difference = estimate, mcse = se,
          lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)
      }
    }
  }
}
paired <- do.call(rbind, paired)
write.csv(paired, "revision/results/summaries/nested-paired-comparisons.csv", row.names = FALSE)
# Post-hoc description of dependence on the fixed factor-grid weighting. Every
# omitted cell is reported; this is not a confirmatory data-dependent exclusion.
a <- subset(selection, learner == "estimated" & criterion == "J2_fitted")
b <- subset(selection, learner == "estimated" & criterion == "W")
delta <- a$absolute_error - b$absolute_error[match(a$id, b$id)]
influence <- lapply(groups(a, factor_columns), function(ii) {
  keep <- setdiff(seq_len(nrow(a)), ii)
  data.frame(omitted_cell = paste(paste(factor_columns, a[ii[1], factor_columns], sep = "="), collapse = ";"),
    omitted_n = length(ii), omitted_mean_difference = mean(delta[ii]),
    remaining_n = length(keep), remaining_mean_difference = mean(delta[keep]),
    remaining_mcse = fixed_mcse(delta[keep], a$id[keep]),
    status = "Post-hoc descriptive grid sensitivity; no selection-adjusted interval")
})
write.csv(do.call(rbind, influence), "revision/results/summaries/nested-cell-influence.csv", row.names = FALSE)
by_scheme <- lapply(groups(z, c("learner", "scheme", factor_columns)), function(ii) {
  d <- z[ii, ]
  data.frame(d[1, c("learner", "scheme", factor_columns)], n = nrow(d),
    mean_error = mean(d$score_error), mcse = fixed_mcse(d$score_error, d$id),
    mean_absolute_error = mean(abs(d$score_error)),
    conditional_risk = mean(d$conditional_target_risk), mean_cv = mean(d$cv_score),
    mean_known_exact_gap = mean(d$known_exact_gap))
})
write.csv(do.call(rbind, by_scheme), "revision/results/summaries/nested-by-cell.csv", row.names = FALSE)
saveRDS(list(selection = selection, summary = summary, comparisons = paired,
             input_md5 = tools::md5sum(file.path("revision/results", paste0(stem, ".rds")))),
        "revision/results/nested-analysis.rds")
cat("Analysed", length(unique(z$id)), "of", nrow(expected), "configurations.\n")
print(summary[summary$stratum == paste0("learner=", summary$learner, ";criterion=", summary$criterion),
              c("learner", "criterion", "mean_signed_error", "mean_absolute_error", "absolute_mcse", "rmse", "tie_fraction")], row.names = FALSE)
print(subset(paired, stratum == "all" & criterion == "J2_fitted" & baseline == "W"))
