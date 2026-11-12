# Post-hoc paired comparison specified in fold-count-sensitivity-protocol.json.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
baseline <- normalizePath(args[1], mustWork = TRUE)
old <- readRDS(file.path(baseline, "revision/results/nested-estimation-rng.rds"))
new <- readRDS("revision/results/nested-estimation-k10.rds")
configuration <- read.csv("revision/results/nested-estimation-k10-configuration.csv")
old_configuration <- read.csv(file.path(baseline, "revision/results/nested-estimation-rng-configuration.csv"))
stopifnot(identical(configuration, old_configuration), nrow(configuration) == 1920L,
          length(readRDS("revision/results/nested-estimation-k10-errors.rds")) == 0L,
          nrow(new$results) == 15360L, nrow(new$fits) == 78720L,
          identical(old$metadata$fitter_body, new$metadata$fitter_body))
factors <- c("range", "signal", "design", "target")
cell <- interaction(configuration[factors], drop = TRUE)
fixed_mcse <- function(x, ids) {
  strata <- split(x, cell[match(ids, configuration$id)], drop = TRUE)
  stopifnot(all(lengths(strata) >= 2L))
  sqrt(sum(vapply(strata, function(v) length(v) * var(v), 0.0))) / length(x)
}
groups <- function(d) split(seq_len(nrow(d)), interaction(d[factors], drop = TRUE, lex.order = TRUE))

# The response and full fitted predictor must be identical across fold counts.
key <- function(d) paste(d$id, d$learner, d$scheme, sep = ":")
match_old <- match(key(new$results), key(old$results))
stopifnot(!anyNA(match_old), !anyDuplicated(match_old))
verification <- list()
for (column in c("conditional_target_risk", "fitted_range", "fitted_signal", "fitted_variance")) {
  difference <- max(abs(new$results[[column]] - old$results[[column]][match_old]))
  stopifnot(difference <= 1e-12)
  verification[[length(verification) + 1L]] <- data.frame(quantity = column, maximum_absolute_difference = difference)
}
a <- old$fits[old$fits$name == "full", ]; b <- new$fits[new$fits$name == "full", ]
stopifnot(nrow(a) == 1920L, nrow(b) == 1920L, setequal(a$id, b$id))
for (column in c("range", "nugget", "psill", "signal")) {
  difference <- max(abs(b[[column]] - a[[column]][match(b$id, a$id)]))
  stopifnot(difference <= 1e-12)
  verification[[length(verification) + 1L]] <- data.frame(quantity = paste0("full_fit_", column), maximum_absolute_difference = difference)
}
write.csv(do.call(rbind, verification), "revision/results/fold-count-verification.csv", row.names = FALSE)

five <- read.csv(file.path(baseline, "revision/results/summaries/nested-selection-by-field.csv"))
ten <- read.csv("revision/results/summaries/nested-selection-by-field.csv")
stopifnot(identical(names(five), names(ten)), nrow(five) == 49920L, nrow(ten) == 49920L)
five$fold_count <- 5L; ten$fold_count <- 10L
selection <- rbind(five, ten)
contrasts <- list(); paired_fields <- list()
for (learner in c("estimated", "known")) {
  for (fold_count in c(5L, 10L)) {
    d <- selection[selection$learner == learner & selection$fold_count == fold_count, ]
    control <- d[d$criterion == "W", ]
    for (criterion in setdiff(unique(d$criterion), c("W", "oracle"))) {
      candidate <- d[d$criterion == criterion, ]
      stopifnot(nrow(candidate) == 1920L, !anyDuplicated(candidate$id), setequal(candidate$id, control$id))
      difference <- candidate$absolute_error - control$absolute_error[match(candidate$id, control$id)]
      for (ii in c(list(all = seq_len(nrow(candidate))), groups(candidate))) {
        label <- if (length(ii) == nrow(candidate)) "all" else
          paste(paste(factors, candidate[ii[1], factors], sep = "="), collapse = ";")
        estimate <- mean(difference[ii]); se <- fixed_mcse(difference[ii], candidate$id[ii])
        contrasts[[length(contrasts) + 1L]] <- data.frame(learner, fold_count, criterion,
          baseline = "W", stratum = label, n = length(ii), difference = estimate, mcse = se,
          lower = estimate - 1.96 * se, upper = estimate + 1.96 * se,
          candidate_mae = mean(candidate$absolute_error[ii]),
          W_mae = mean(control$absolute_error[match(candidate$id[ii], control$id)]))
      }
      if (criterion == "J2_fitted") {
        paired_fields[[length(paired_fields) + 1L]] <- data.frame(
          candidate[c("id", "learner", factors)], fold_count, difference)
      }
    }
  }
}
contrasts <- do.call(rbind, contrasts)
write.csv(contrasts, "revision/results/fold-count-comparisons.csv", row.names = FALSE)
paired_fields <- do.call(rbind, paired_fields)
change <- list(); change_summary <- list()
for (learner in c("estimated", "known")) {
  a <- paired_fields[paired_fields$learner == learner & paired_fields$fold_count == 5L, ]
  b <- paired_fields[paired_fields$learner == learner & paired_fields$fold_count == 10L, ]
  d <- b[c("id", "learner", factors)]
  d$delta_five <- a$difference[match(b$id, a$id)]
  d$delta_ten <- b$difference
  d$change <- d$delta_ten - d$delta_five
  change[[length(change) + 1L]] <- d
  for (ii in c(list(all = seq_len(nrow(d))), groups(d))) {
    label <- if (length(ii) == nrow(d)) "all" else
      paste(paste(factors, d[ii[1], factors], sep = "="), collapse = ";")
    estimate <- mean(d$change[ii]); se <- fixed_mcse(d$change[ii], d$id[ii])
    change_summary[[length(change_summary) + 1L]] <- data.frame(learner, stratum = label, n = length(ii),
      delta_five = mean(d$delta_five[ii]), delta_ten = mean(d$delta_ten[ii]), change = estimate,
      mcse = se, lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)
  }
}
write.csv(do.call(rbind, change), "revision/results/fold-count-change-by-field.csv", row.names = FALSE)
write.csv(do.call(rbind, change_summary), "revision/results/fold-count-change-summary.csv", row.names = FALSE)
dir.create("revision/figures", showWarnings = FALSE)
pdf("revision/figures/fold-count-comparison.pdf", width = 8.5, height = 9, useDingbats = FALSE)
for (learner in c("estimated", "known")) {
  d <- contrasts[contrasts$learner == learner & contrasts$criterion == "J2_fitted" & contrasts$stratum != "all", ]
  labels <- unique(d$stratum)
  five <- d[d$fold_count == 5L, ]; ten <- d[d$fold_count == 10L, ]
  five <- five[match(labels, five$stratum), ]; ten <- ten[match(labels, ten$stratum), ]
  labels <- gsub("range=", "scale ", labels, fixed = TRUE)
  labels <- gsub("signal=", "signal ", labels, fixed = TRUE)
  labels <- gsub("design=", "", labels, fixed = TRUE)
  labels <- gsub("target=", "", labels, fixed = TRUE)
  labels <- gsub(";", " / ", labels, fixed = TRUE)
  par(mar = c(8.3, 17, 4, 1), las = 1, cex = .78)
  y <- rev(seq_along(labels)); limits <- range(c(d$lower, d$upper, 0))
  plot(NA, xlim = limits, ylim = c(.5, length(y) + .5), yaxt = "n", xlab = "", ylab = "",
       main = paste(if (learner == "estimated") "Estimated" else "Known", "covariance prediction"))
  abline(v = 0, lty = 3, col = "gray50")
  axis(2, at = y, labels = labels, tick = FALSE)
  segments(five$lower, y + .12, five$upper, y + .12, col = "#444444")
  points(five$difference, y + .12, pch = 16, col = "#444444")
  segments(ten$lower, y - .12, ten$upper, y - .12, col = "#0072B2")
  points(ten$difference, y - .12, pch = 17, col = "#0072B2")
  legend("bottomleft", legend = c("Five folds", "Ten folds"), col = c("#444444", "#0072B2"), pch = c(16, 17), bty = "n")
  mtext("CV estimation MAE difference: fitted |J2| minus W", side = 1, line = 2.4)
  mtext("Post-hoc sensitivity; same 80 fields per cell at both fold counts.", side = 1, line = 5.2, cex = .8)
  mtext("Bars: +/- 1.96 Monte Carlo SE. Negative values favour fitted |J2|.", side = 1, line = 6.5, cex = .8)
}
dev.off()
print(subset(contrasts, stratum == "all" & criterion %in% c("J2_fitted", "fixed_knndm", "fixed_random")), row.names = FALSE)
print(subset(do.call(rbind, change_summary), stratum == "all"), row.names = FALSE)
cat("All 1920 paired fields and 24 cells retained; 78720 ten-fold fits; full-data fits and target risks unchanged.\n")
