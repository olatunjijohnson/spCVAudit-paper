# Complete-factor displays use the same retained results as the main figures.
source("revision/R/local_libraries.R")
library(ggplot2)
theme_set(theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), legend.position = "bottom",
  legend.title = element_blank(), strip.background = element_rect(fill = "#f4f4f4"),
  plot.title = element_text(face = "bold", size = 12),
  plot.subtitle = element_text(size = 10), plot.caption = element_text(hjust = 0, size = 9)))
cols <- c(random = "#0072B2", block = "#E69F00", buffered = "#D55E00", knndm = "#009E73")
scheme_names <- c(random = "Random", block = "Blocks", buffered = "Buffered", knndm = "kNNDM")
parse_strata <- function(d) {
  entries <- strsplit(d$stratum, ";", fixed = TRUE)
  for (name in unique(sub("=.*", "", unlist(entries)))) {
    d[[paste0("factor_", name)]] <- vapply(entries, function(x) {
      hit <- x[startsWith(x, paste0(name, "="))]
      if (length(hit) == 1L) substring(hit, nchar(name) + 2L) else NA_character_
    }, "")
  }
  d
}

reference <- read.csv("revision/results/summaries/reference-paired-comparisons.csv")
reference <- parse_strata(subset(reference, criterion == "J2" & baseline == "W" & grepl("family=", stratum)))
stopifnot(nrow(reference) == 128L, all(reference$n == 40L))
reference$signal <- factor(reference$factor_signal, c("0", "0.3", "0.7", "0.95"))
reference$case <- paste(ifelse(reference$factor_design == "random", "Uniform", "Clustered"),
                      ifelse(reference$factor_target == "domain", "in-domain", "shifted"), sep = ": ")
reference$sample_label <- paste("n =", reference$factor_n)
reference$scale_label <- paste("Scale =", reference$factor_range)
pdf("revision/figures/S1-reference-all.pdf", width = 7.2, height = 7.5, useDingbats = FALSE)
for (family in c("exponential", "matern32")) {
  d <- reference[reference$factor_family == family, ]
  p <- ggplot(d, aes(difference, signal, colour = case, shape = case)) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "#888888") +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = .12,
                  position = position_dodge(.65), orientation = "y") +
    geom_point(position = position_dodge(.65), size = 1.8) +
    facet_grid(sample_label ~ scale_label, scales = "free_x") +
    scale_colour_manual(values = c("#222222", "#777777", "#6B4C9A", "#B37700")) +
    labs(title = paste("Known-covariance selection:", if (family == "exponential") "exponential" else "Matérn 3/2"),
         subtitle = "Every factor cell; 40 independent configurations per point",
         x = "Difference in mean absolute exact bias: |J2| minus W", y = "Signal fraction",
         caption = "Bars: ±1.96 Monte Carlo SE of paired differences. Negative values favour |J2|.\nAxis scales differ by covariance scale; zero-signal differences are exactly zero.") +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  print(p)
}
dev.off()

sensitivity <- readRDS("revision/results/sensitivity-rng.rds")$results
fit <- subset(sensitivity, diagnostic_model == "estimated_isotropic" & scheme == "random")
stopifnot(nrow(fit) == 1440L, all(fit$fitted_ok))
fit$geometry <- ifelse(fit$axis_ratio == 1, "Isotropic", "Anisotropic (axis ratio 0.25)")
fit$scale_label <- paste("Generating scale =", fit$range)
fit$signal_label <- factor(fit$signal, c(.3, .7, .95))
fit$design_label <- ifelse(fit$design == "random", "Uniform", "Clustered")
fit$range_ratio <- fit$fitted_range / fit$range
truth <- unique(fit[c("geometry", "scale_label", "signal_label", "signal")])
pdf("revision/figures/S3-covariance-fits.pdf", width = 7.2, height = 7.5, useDingbats = FALSE)
for (quantity in c("fitted_signal", "range_ratio")) {
  p <- ggplot(fit, aes(signal_label, .data[[quantity]], fill = design_label)) +
    geom_boxplot(width = .65, outlier.size = .65, outlier.alpha = .4) +
    facet_grid(scale_label ~ geometry) + scale_fill_manual(values = c("#bbbbbb", "#ffffff")) +
    labs(x = "Generating signal fraction", fill = NULL,
         title = "Estimated isotropic variogram parameters",
         subtitle = "All 1,440 independent fits; targets pooled within each displayed group")
  if (quantity == "fitted_signal") {
    p <- p + geom_point(data = truth, aes(signal_label, signal), inherit.aes = FALSE,
                       shape = 4, size = 2.5, stroke = 1) +
      labs(y = "Estimated signal fraction",
           caption = "Crosses mark the generating signal fraction. Boxes: median and interquartile range;\nwhiskers: 1.5 × IQR. Numerical convergence does not establish a well-identified covariance.")
  } else {
    p <- p + scale_y_log10() + geom_hline(yintercept = 1, linetype = "dotted") +
      labs(y = "Estimated scale / generating scale (log axis)",
           caption = "All finite fitted scales are retained. Under anisotropy, the denominator is the\ngenerating scale parameter, not a true isotropic range. Large ratios can reflect\nweak identification.")
  }
  print(p)
}
dev.off()

nested <- read.csv("revision/results/summaries/nested-paired-comparisons.csv")
nested <- parse_strata(subset(nested, baseline == "W" & grepl("signal=", stratum) &
  criterion %in% c("J2_fitted", "J2_known", "fixed_knndm")))
stopifnot(nrow(nested) == 144L, all(nested$n == 80L))
nested$case <- paste(ifelse(nested$factor_design == "random", "Uniform", "Clustered"),
                    ifelse(nested$factor_target == "domain", "in-domain", "shifted"), sep = ": ")
nested$setting <- factor(paste(nested$factor_range, nested$factor_signal, sep = " / "),
  rev(c("0.05 / 0.3", "0.05 / 0.7", "0.05 / 0.95", "0.2 / 0.3", "0.2 / 0.7", "0.2 / 0.95")))
nested$choice <- factor(nested$criterion, c("J2_fitted", "J2_known", "fixed_knndm"),
                        c("|J2| fitted", "|J2| known", "Fixed kNNDM"))
pdf("revision/figures/S4-nested-all.pdf", width = 7.2, height = 8, useDingbats = FALSE)
for (learner in c("estimated", "known")) {
  d <- nested[nested$learner == learner, ]
  print(ggplot(d, aes(mean_absolute_error_difference, setting, colour = choice, shape = choice)) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "#888888") +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = .12,
                  position = position_dodge(.6), orientation = "y") +
    geom_point(position = position_dodge(.6), size = 1.8) + facet_wrap(~case, ncol = 2) +
    scale_colour_manual(values = c("#222222", "#888888", cols[["knndm"]])) +
    labs(title = paste(if (learner == "estimated") "Estimated" else "Known", "covariance prediction"),
         subtitle = "Every factor cell; 80 independent datasets per point",
         x = "CV estimation MAE difference from selection by W", y = "Generating scale / signal",
         caption = "Bars: ±1.96 Monte Carlo SE of paired differences. Negative values favour the listed rule.\nDiagnostic weights may be fitted or known, separately from the covariance used for prediction."))
}
dev.off()

status <- system2("python3", "revision/R/make_jura_figures.py --supplement")
stopifnot(status == 0L)
cat("Generated supplementary reference, covariance-fit, nested and Jura displays.\n")
