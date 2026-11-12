source("revision/R/common.R")
suffix <- Sys.getenv("REVISION_ESTIMATION_SUFFIX", "-rng")
stem <- paste0("revision/results/sensitivity", suffix)
z <- readRDS(paste0(stem, ".rds"))$results
expected <- read.csv(paste0(stem, "-configuration.csv"))
stopifnot(nrow(z) == nrow(expected) * 8L, all(table(z$id) == 8L),
          all(z$training_size_effect <= 1e-8),
          max(abs(z$exact - z$matched_size_gap - z$training_size_effect)) < 1e-8,
          all(z$exact >= z$bound_lower - 1e-8), all(z$exact <= z$bound_upper + 1e-8))
dir.create("revision/results/summaries", showWarnings = FALSE)
groups <- function(d, names) split(seq_len(nrow(d)), interaction(d[names], drop = TRUE, lex.order = TRUE))
factor_cell <- interaction(expected[c("range", "signal", "design", "target", "axis_ratio")], drop = TRUE)
ratio_mcse <- function(numerator, denominator, ids) {
  totals <- matrix(0, nrow(expected), 2)
  a <- rowsum(cbind(numerator, denominator), ids, reorder = FALSE)
  totals[match(as.integer(rownames(a)), expected$id), ] <- a
  ratio <- sum(totals[, 1]) / sum(totals[, 2])
  residual <- totals[, 1] - ratio * totals[, 2]
  strata <- split(residual, factor_cell)
  sqrt(sum(vapply(strata, function(x) length(x) * var(x), 0.0))) / sum(totals[, 2])
}
factors <- c("range", "signal", "design", "target", "axis_ratio", "scheme", "diagnostic_model")
summary <- lapply(groups(z, factors), function(ii) {
  d <- z[ii, ]; eligible <- abs(d$exact) > .01; fitted <- d$fitted_ok
  sign_known <- sign(d$true_J2) == sign(d$exact)
  sign_est <- sign(d$J2) == sign(d$exact)
  certified <- d$bound_lower > 0 | d$bound_upper < 0
  data.frame(d[1, factors, drop = FALSE], n = nrow(d), fit_success_fraction = mean(fitted),
    eligible = sum(eligible), accuracy_J2 = if (any(eligible & fitted)) mean(sign_est[eligible & fitted]) else NA_real_,
    accuracy_true_J2 = if (any(eligible)) mean(sign_known[eligible]) else NA_real_,
    mean_abs_proxy_error = mean(abs(d$proxy - d$exact), na.rm = TRUE),
    mean_exact_gap = mean(d$exact), mean_matched_size_gap = mean(d$matched_size_gap),
    mean_training_effect = mean(d$training_size_effect),
    mean_abs_training_effect = mean(abs(d$training_size_effect)),
    sign_changed_after_matching_size = mean(sign(d$exact) != sign(d$matched_size_gap)),
    certificate_coverage = mean(certified),
    certificate_correct = if (any(certified)) mean(ifelse(d$bound_lower[certified] > 0, 1, -1) == sign(d$exact[certified])) else NA_real_)
})
write.csv(do.call(rbind, summary), "revision/results/summaries/sensitivity-by-cell.csv", row.names = FALSE)
fit_rows <- z[z$scheme == "random" & z$diagnostic_model == "estimated_isotropic", ]
write.csv(fit_rows[, c(names(expected), "fitted_ok", "fitted_range", "fitted_signal", "fitted_variance", "fit_error", "fit_warnings")],
          "revision/results/summaries/covariance-fits.csv", row.names = FALSE)
headline <- lapply(groups(z, c("axis_ratio", "diagnostic_model")), function(ii) {
  d <- z[ii, ]; eligible <- abs(d$exact) > .01 & d$fitted_ok
  certified <- d$bound_lower > 0 | d$bound_upper < 0
  data.frame(d[1, c("axis_ratio", "diagnostic_model")], n = nrow(d),
    eligible = sum(eligible), sign_accuracy = mean(sign(d$J2[eligible]) == sign(d$exact[eligible])),
    sign_mcse = ratio_mcse(as.numeric(sign(d$J2) == sign(d$exact)) * eligible, as.numeric(eligible), d$id),
    true_J2_sign_accuracy = mean(sign(d$true_J2[eligible]) == sign(d$exact[eligible])),
    true_J2_sign_mcse = ratio_mcse(as.numeric(sign(d$true_J2) == sign(d$exact)) * eligible, as.numeric(eligible), d$id),
    mean_abs_proxy_error = mean(abs(d$proxy - d$exact), na.rm = TRUE),
    certificate_coverage = mean(certified), mean_abs_training_effect = mean(abs(d$training_size_effect)),
    certificate_mcse = ratio_mcse(as.numeric(certified), rep(1, nrow(d)), d$id),
    fit_success_fraction = mean(d$fitted_ok))
})
headline <- do.call(rbind, headline)
write.csv(headline, "revision/results/summaries/sensitivity-headlines.csv", row.names = FALSE)
small_rows <- list()
for (ii in groups(z, c("axis_ratio", "diagnostic_model"))) {
  d <- z[ii, ]; eligible <- abs(d$exact) > .01 & d$fitted_ok
  for (diagnostic in c("J1", "J2")) for (threshold in c(1e-12, .005, .01, .02)) {
    small <- is.finite(d[[diagnostic]]) & abs(d[[diagnostic]]) <= threshold
    called <- is.finite(d[[diagnostic]]) & !small
    small_rows[[length(small_rows) + 1L]] <- data.frame(d[1, c("axis_ratio", "diagnostic_model")],
      diagnostic, threshold, n = nrow(d), n_small = sum(small),
      large_gap_given_small = if (any(small)) mean(abs(d$exact[small]) > .05) else NA_real_,
      maximum_abs_gap_given_small = if (any(small)) max(abs(d$exact[small])) else NA_real_,
      eligible = sum(eligible), call_coverage = mean(called[eligible]),
      called_sign_accuracy = if (any(eligible & called))
        mean(sign(d[[diagnostic]][eligible & called]) == sign(d$exact[eligible & called])) else NA_real_)
  }
}
write.csv(do.call(rbind, small_rows), "revision/results/summaries/sensitivity-near-zero.csv", row.names = FALSE)
cat("Validated", nrow(expected), "sensitivity configurations.\n")
print(headline)
