source("revision/R/common.R")
z <- readRDS("revision/results/learners.rds")$results
expected <- read.csv("revision/results/learners-configuration.csv")
stopifnot(nrow(z) == 80L * nrow(expected), all(table(z$id) == 80L),
          !anyDuplicated(z[c("id", "alpha", "learner", "scheme")]),
          all(is.finite(z$gap)), max(abs(z$gap - z$target_loss + z$cv_loss)) < 1e-9)
geometry_groups <- split(seq_len(nrow(z)), interaction(z$id, z$scheme, drop = TRUE))
stopifnot(all(vapply(geometry_groups, function(ii) diff(range(z$J1_marginal[ii])) < 1e-12, logical(1))),
          all(abs(z$J1_conditional[z$alpha == .9]) < 1e-12),
          all(abs(z$target_loss[z$alpha == .9 & z$learner == "known_mean"] - .1) < 1e-8),
          all(abs(z$target_loss[z$alpha == .9 & z$learner == "oracle_kriging"] - .1) < 1e-8))
dir.create("revision/results/summaries", showWarnings = FALSE)
group_rows <- function(d, columns) split(seq_len(nrow(d)), interaction(d[columns], drop = TRUE, lex.order = TRUE))
factors <- c("range", "design", "target", "alpha", "learner", "scheme")
summary <- lapply(group_rows(z, factors), function(ii) {
  d <- z[ii, ]; n <- nrow(d); se <- sd(d$gap) / sqrt(n)
  data.frame(d[1, factors, drop = FALSE], n, mean_gap = mean(d$gap), mcse_gap = se,
    lower = mean(d$gap) - 1.96 * se, upper = mean(d$gap) + 1.96 * se,
    mean_cv = mean(d$cv_loss), mean_target = mean(d$target_loss),
    mean_J1_marginal = mean(d$J1_marginal), mean_J2_marginal = mean(d$J2_marginal),
    mean_J1_conditional = mean(d$J1_conditional), mean_J2_conditional = mean(d$J2_conditional),
    mean_abs_realized_gap = mean(abs(d$gap)),
    mean_conditional_proxy = mean(d$one_neighbour_conditional_gap))
})
summary <- do.call(rbind, summary)
write.csv(summary, "revision/results/summaries/learners-by-cell.csv", row.names = FALSE)

selection <- list()
criteria <- c("A", "W", "J1_marginal", "J2_marginal", "U1_marginal", "U2_marginal",
              "J1_conditional", "J2_conditional", "U1_conditional", "U2_conditional", "fixed_random", "fixed_knndm", "oracle")
for (ii in group_rows(z, c("id", "alpha", "learner"))) {
  d <- z[ii, ]; stopifnot(nrow(d) == 4L, diff(range(d$target_loss)) < 1e-10)
  for (criterion in criteria) {
    score <- switch(criterion, fixed_random = as.numeric(d$scheme != "random"),
      fixed_knndm = as.numeric(d$scheme != "knndm"), oracle = abs(d$gap), abs(d[[criterion]]))
    tied <- which(score <= min(score) + 1e-12)
    selection[[length(selection) + 1L]] <- data.frame(d[1, c(names(expected), "alpha", "learner"), drop = FALSE],
      criterion, mean_abs_realized_gap = mean(abs(d$gap[tied])), signed_gap = mean(d$gap[tied]),
      oracle_abs_gap = min(abs(d$gap)), n_tied = length(tied), selected = paste(sort(d$scheme[tied]), collapse = ";"))
  }
}
selection <- do.call(rbind, selection)
write.csv(selection, "revision/results/summaries/learners-selection-by-field.csv", row.names = FALSE)
ss <- lapply(group_rows(selection, c("range", "design", "target", "alpha", "learner", "criterion")), function(ii) {
  d <- selection[ii, ]; n <- nrow(d)
  data.frame(d[1, c("range", "design", "target", "alpha", "learner", "criterion"), drop = FALSE], n,
    mean_abs_realized_gap = mean(d$mean_abs_realized_gap),
    mcse = sd(d$mean_abs_realized_gap) / sqrt(n), mean_signed_gap = mean(d$signed_gap),
    mean_regret = mean(d$mean_abs_realized_gap - d$oracle_abs_gap))
})
write.csv(do.call(rbind, ss), "revision/results/summaries/learners-selection.csv", row.names = FALSE)

# Changes across alpha are paired within the same geometry and underlying fields.
contrasts <- list()
for (ii in group_rows(z, c("id", "learner", "scheme"))) {
  d <- z[ii, ]; a <- d[d$alpha == 0, ]; b <- d[d$alpha == .9, ]
  contrasts[[length(contrasts) + 1L]] <- data.frame(a[, c("id", "range", "design", "target", "learner", "scheme")],
    change_gap = b$gap - a$gap, change_marginal_J1 = b$J1_marginal - a$J1_marginal,
    change_conditional_J1 = b$J1_conditional - a$J1_conditional)
}
contrasts <- do.call(rbind, contrasts)
cc <- lapply(group_rows(contrasts, c("range", "design", "target", "learner", "scheme")), function(ii) {
  d <- contrasts[ii, ]
  data.frame(d[1, c("range", "design", "target", "learner", "scheme")],
    n = nrow(d), mean_change_gap = mean(d$change_gap), mcse = sd(d$change_gap) / sqrt(nrow(d)),
    marginal_index_change = max(abs(d$change_marginal_J1)))
})
write.csv(do.call(rbind, cc), "revision/results/summaries/learners-alpha-contrast.csv", row.names = FALSE)
cat("Validated", nrow(z), "evaluations from", nrow(expected), "independent fields.\n")
cat("Marginal geometry diagnostics are invariant across alpha and learners.\n")
print(subset(summary, range == .2 & design == "clustered" & target == "domain" & scheme == "random"))
