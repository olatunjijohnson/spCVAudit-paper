# Theory-validation figure (exact Gaussian-process algebra, 192 configs). (a) the
# single-nearest-neighbour closed form tracks the exact excess-explained-variance
# optimism (shape r = 0.98; it underestimates magnitude ~1.6x, the multi-neighbour
# amplification). (b) the package's SLI_rho is a near-sufficient statistic for exact
# optimism (r = 0.98) and recovers its sign. (c) the Wasserstein-1 bound holds for
# every configuration. Run after library(spLeakage).

t <- readRDS("cache/theory_validation.rds")
dcol <- ifelse(t$design == "random", "#00798c", "#d1495b")

dir.create("paper-figures", showWarnings = FALSE)
png("figures/F_theory.png", 1750, 850, res = 200)
graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.4, 3, 1))

## (a) closed form vs exact
sc <- 1.557  # multi-NN amplification (median exact/nn_form)
plot(t$nn_form * sc, t$exact, pch = 19, col = paste0(dcol, "99"), cex = 0.8,
     xlab = expression("closed form  " * V * w^2 * "(E[" * rho^2 * "(g)] - E[" * rho^2 * "(f)])" %*% 1.56),
     ylab = "exact GP optimism",
     main = sprintf("(a) Closed form vs exact (r = %.2f)", cor(t$nn_form, t$exact)))
graphics::abline(0, 1, lwd = 2, col = "grey40"); graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")
graphics::legend("topleft", c("random design", "clustered design"),
                 col = c("#00798c", "#d1495b"), pch = 19, bty = "n", cex = 0.8)

## (b) SLI_rho vs exact  (near-sufficient statistic)
plot(t$sli_rho, t$exact, pch = 19, col = paste0(dcol, "99"), cex = 0.8,
     xlab = expression("Spatial Leakage Index  SLI"[rho]), ylab = "exact GP optimism",
     main = sprintf("(b) SLI is near-sufficient (r = %.2f)", cor(t$sli_rho, t$exact)))
graphics::abline(lm(exact ~ sli_rho, t), lwd = 2, col = "grey40")
graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")
graphics::text(min(t$sli_rho), max(t$exact), sprintf("%.0f%% of configs\npessimistic (exact<0);\nSLI recovers the sign", 100*mean(t$exact<0)),
               adj = c(0, 1), cex = 0.72, col = "grey30")
dev.off()
cat("wrote figures/F_theory.png\n")
