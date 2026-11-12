# Summaries are paired by independently generated configuration, never by fold.
source("revision/R/common.R")
reference <- readRDS("revision/results/gp-reference.rds")$results
knndm <- readRDS("revision/results/gp-knndm.rds")$results
expected <- read.csv("revision/results/gp-reference-configuration.csv")
stopifnot(nrow(reference) == 3L * nrow(expected), nrow(knndm) == nrow(expected),
          !anyDuplicated(knndm$id), setequal(knndm$id, expected$id))
shared <- intersect(names(reference), names(knndm))
z <- rbind(reference[, shared], knndm[, shared])
stopifnot(all(table(z$id) == 4L), !anyDuplicated(paste(z$id, z$scheme)),
          all(is.finite(z$exact)), all(z$risk_cv >= -1e-8), all(z$risk_target >= -1e-8),
          max(abs(z$exact - (z$risk_target - z$risk_cv))) < 1e-8,
          max(abs(z$J2) - z$U2) < 1e-8,
          max(abs(z$A) - z$W) < 1e-8,
          max(abs(z$exact[z$signal == 0])) < 1e-8)
for (column in names(expected)) {
  stopifnot(isTRUE(all.equal(z[[column]], expected[[column]][match(z$id, expected$id)],
                            check.attributes = FALSE)))
}
dir.create("revision/results/summaries", showWarnings = FALSE)
mean_se <- function(x) c(mean = mean(x), mcse = if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_)
cell <- interaction(expected[c("n", "range", "signal", "design", "target", "family")], drop = TRUE)
fixed_grid_mcse <- function(x, ids) {
  groups <- split(x, cell[match(ids, expected$id)], drop = TRUE)
  # The factor grid is fixed. Only variation within independently replicated
  # factor cells contributes Monte Carlo uncertainty to an equally weighted mean.
  if (any(lengths(groups) < 2L)) return(NA_real_)
  sqrt(sum(vapply(groups, function(v) length(v) * var(v), 0.0))) / length(x)
}
cluster_ratio <- function(numerator, denominator, ids) {
  totals <- matrix(0, nrow(expected), 2)
  a <- rowsum(cbind(numerator, denominator), ids, reorder = FALSE)
  totals[match(as.integer(rownames(a)), expected$id), ] <- a
  if (sum(totals[, 2]) == 0) return(c(estimate = NA_real_, mcse = NA_real_))
  ratio <- sum(totals[, 1]) / sum(totals[, 2])
  residual <- totals[, 1] - ratio * totals[, 2]
  c(estimate = ratio, mcse = fixed_grid_mcse(residual, expected$id) / mean(totals[, 2]))
}
group_rows <- function(d, columns) {
  if (!length(columns)) return(list(all = seq_len(nrow(d))))
  split(seq_len(nrow(d)), interaction(d[columns], drop = TRUE, lex.order = TRUE))
}
sign_rows <- list()
strata <- list(character(), "scheme", c("design", "target"), c("scheme", "design", "target"),
               c("range", "signal"), "family", "n", "crosses",
               c("n", "range", "signal", "design", "target", "family", "scheme"))
for (columns in strata) for (ii in group_rows(z, columns)) {
  d <- z[ii, ]
  label <- if (length(columns)) paste(paste(columns, d[1, columns], sep = "="), collapse = ";") else "all"
  for (threshold in c(.005, .01, .02, .05)) for (diagnostic in c("A", "J1", "J2")) {
    eligible <- abs(d$exact) > threshold
    called <- abs(d[[diagnostic]]) > 1e-12
    correct <- as.numeric(sign(d[[diagnostic]]) == sign(d$exact))
    # Delta-method ratio SE clusters schemes by configuration and stratifies
    # the sampling variance by the fixed experimental factor cells.
    ms <- cluster_ratio(correct * eligible, as.numeric(eligible), d$id)
    sign_rows[[length(sign_rows) + 1L]] <- data.frame(stratum = label, threshold, diagnostic,
      n_configurations = length(unique(d$id)), eligible_configurations = length(unique(d$id[eligible])),
      eligible_evaluations = sum(eligible), coverage = if (any(eligible)) mean(called[eligible]) else NA_real_,
      accuracy = ms[["estimate"]], mcse = ms[["mcse"]])
  }
}
write.csv(do.call(rbind, sign_rows), "revision/results/summaries/reference-sign.csv", row.names = FALSE)

near_zero <- list()
for (columns in list(character(), "scheme", c("design", "target"), "signal")) {
  for (ii in group_rows(z, columns)) {
    d <- z[ii, ]; label <- if (length(columns)) paste(paste(columns, d[1, columns], sep = "="), collapse = ";") else "all"
    for (diagnostic in c("J1", "J2")) for (threshold in c(.005, .01, .02)) {
      small <- abs(d[[diagnostic]]) <= threshold
      near_zero[[length(near_zero) + 1L]] <- data.frame(stratum = label, diagnostic, threshold,
        n = nrow(d), n_small = sum(small), small_fraction = mean(small),
        large_gap_given_small = if (any(small)) mean(abs(d$exact[small]) > .05) else NA_real_,
        maximum_abs_gap_given_small = if (any(small)) max(abs(d$exact[small])) else NA_real_)
    }
  }
}
write.csv(do.call(rbind, near_zero), "revision/results/summaries/reference-cancellation.csv", row.names = FALSE)

# Descriptive rank correlations use the same numerical resolution as selection
# ties. Average ranks are assigned to sorted groups with full span <= tolerance.
# Exact floating-point ties otherwise differ across numerically equal replays.
rank_with_tolerance <- function(x, tolerance = 1e-12) {
  stopifnot(all(is.finite(x)), is.finite(tolerance), tolerance >= 0)
  order_x <- order(x)
  value <- x[order_x]
  result <- numeric(length(x))
  first <- 1L
  while (first <= length(x)) {
    last <- first
    while (last < length(x) && value[last + 1L] - value[first] <= tolerance) last <- last + 1L
    result[order_x[first:last]] <- (first + last) / 2
    first <- last + 1L
  }
  result
}

selection <- list()
diagnostics <- c("A", "W", "J1", "J2", "U1", "U2", "fixed_random", "fixed_knndm", "oracle")
for (d in split(z, z$id)) {
  oracle <- min(abs(d$exact))
  meta <- d[1, names(expected), drop = FALSE]
  for (criterion in diagnostics) {
    scores <- switch(criterion, fixed_random = as.numeric(d$scheme != "random"),
                     fixed_knndm = as.numeric(d$scheme != "knndm"), oracle = abs(d$exact), abs(d[[criterion]]))
    tied <- which(scores <= min(scores) + 1e-12)
    score_rank <- rank_with_tolerance(scores)
    risk_rank <- rank_with_tolerance(abs(d$exact))
    # Average over tied choices; target error never breaks a diagnostic tie.
    selection[[length(selection) + 1L]] <- data.frame(meta, criterion,
      abs_bias = mean(abs(d$exact[tied])), signed_bias = mean(d$exact[tied]),
      regret = mean(abs(d$exact[tied])) - oracle, oracle_abs_bias = oracle,
      n_tied = length(tied), selected = paste(sort(d$scheme[tied]), collapse = ";"),
      rank_correlation = if (sd(score_rank) > 0 && sd(risk_rank) > 0)
        cor(score_rank, risk_rank) else NA_real_)
  }
}
selection <- do.call(rbind, selection)
write.csv(selection, "revision/results/summaries/reference-selection-by-configuration.csv", row.names = FALSE)
selection_summary <- list()
for (columns in list("criterion", c("criterion", "design", "target"), c("criterion", "range", "signal"),
                    c("criterion", "family"), c("criterion", "n"))) {
  for (ii in group_rows(selection, columns)) {
    d <- selection[ii, ]; ms <- c(mean = mean(d$abs_bias), mcse = fixed_grid_mcse(d$abs_bias, d$id))
    rs <- c(mean = mean(d$regret), mcse = fixed_grid_mcse(d$regret, d$id))
    selection_summary[[length(selection_summary) + 1L]] <- data.frame(
      stratum = paste(paste(columns, d[1, columns], sep = "="), collapse = ";"),
      criterion = d$criterion[1], n = nrow(d), abs_bias = ms[["mean"]], abs_bias_mcse = ms[["mcse"]],
      regret = rs[["mean"]], regret_mcse = rs[["mcse"]],
      mean_rank_correlation = mean(d$rank_correlation, na.rm = TRUE),
      rank_coverage = mean(is.finite(d$rank_correlation)))
  }
}
write.csv(do.call(rbind, selection_summary), "revision/results/summaries/reference-selection.csv", row.names = FALSE)
paired <- list()
for (baseline in c("W", "fixed_random", "fixed_knndm")) {
  b <- selection[selection$criterion == baseline, ]
  for (criterion in setdiff(diagnostics, c(baseline, "oracle"))) {
    d <- selection[selection$criterion == criterion, ]; diff <- d$abs_bias - b$abs_bias[match(d$id, b$id)]
    for (columns in list(character(), c("design", "target"), c("range", "signal"),
                         c("n", "range", "signal", "design", "target", "family"))) {
      for (ii in group_rows(d, columns)) {
        ms <- c(mean = mean(diff[ii]), mcse = fixed_grid_mcse(diff[ii], d$id[ii]))
        paired[[length(paired) + 1L]] <- data.frame(criterion, baseline,
          stratum = if (length(columns)) paste(paste(columns, d[ii[1], columns], sep = "="), collapse = ";") else "all",
          n = length(ii), difference = ms[["mean"]], mcse = ms[["mcse"]],
          lower = ms[["mean"]] - 1.96 * ms[["mcse"]], upper = ms[["mean"]] + 1.96 * ms[["mcse"]])
      }
    }
  }
}
write.csv(do.call(rbind, paired), "revision/results/summaries/reference-paired-comparisons.csv", row.names = FALSE)
proxy <- aggregate(cbind(abs_proxy_error = abs(z$J2 - z$exact),
                         squared_proxy_error = (z$J2 - z$exact)^2,
                         abs_remainder = abs(z$remainder)),
                   z[c("n", "range", "signal", "design", "target", "family", "scheme")], mean)
write.csv(proxy, "revision/results/summaries/reference-proxy-error.csv", row.names = FALSE)
saveRDS(list(results = z, selection = selection, sources = tools::md5sum(c(
  "revision/results/gp-reference.rds", "revision/results/gp-knndm.rds",
  "revision/R/analyse_reference.R"))), "revision/results/reference-analysis.rds")
cat("Validated", nrow(z), "paired scheme evaluations from", nrow(expected), "configurations.\n")
print(subset(do.call(rbind, sign_rows), stratum == "all" & threshold == .01))
print(subset(do.call(rbind, paired), stratum == "all"))
