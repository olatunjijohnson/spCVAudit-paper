# Figure-only entry point, using retained paired comparison values.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
contrasts <- read.csv(args[1])
pdf(args[2], width = 8.5, height = 9, useDingbats = FALSE)
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
print(sessionInfo())
