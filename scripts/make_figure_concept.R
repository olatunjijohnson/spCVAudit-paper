# Conceptual schematic (new Figure 1). (a) Under random CV a test point usually has a
# nearby training neighbour (short reach), but the wall-to-wall map must predict at
# grid cells far from any sample (long reach). (b) The two nearest-neighbour-distance
# distributions therefore differ; the signed gap between them, converted to retained
# spatial correlation, is the Spatial Leakage Index. Run after library(spLeakage).

suppressMessages(library(spLeakage))
set.seed(7)

# a clustered sample in the unit square (a few tight clusters)
K <- 6; per <- 12
ctr <- cbind(runif(K, 0.12, 0.88), runif(K, 0.12, 0.88))
pts <- do.call(rbind, lapply(seq_len(K), function(k)
  cbind(rnorm(per, ctr[k, 1], 0.05), rnorm(per, ctr[k, 2], 0.05))))
pts <- pts[pts[, 1] > 0 & pts[, 1] < 1 & pts[, 2] > 0 & pts[, 2] < 1, ]
n <- nrow(pts)

# deployment grid (wall-to-wall)
gg <- as.matrix(expand.grid(seq(0.03, 0.97, length.out = 22),
                            seq(0.03, 0.97, length.out = 22)))

# NN distances: test->train under random 10-fold; grid->sample at deployment
folds <- sample(rep_len(1:10, n))
gobs <- numeric(n)
for (k in 1:10) { te <- which(folds == k); tr <- which(folds != k)
  for (i in te) gobs[i] <- sqrt(min(colSums((t(pts[tr, , drop = FALSE]) - pts[i, ])^2))) }
fpred <- apply(gg, 1, function(p) sqrt(min(colSums((t(pts) - p)^2))))

dir.create("paper-figures", showWarnings = FALSE)
png("figures/F0_concept.png", 2000, 900, res = 200)
graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3, 1.4), cex.main = 0.95)

## (a) the map with the two reaches
plot(gg, pch = 3, col = "grey80", cex = 0.5, xlab = "", ylab = "", axes = FALSE,
     xlim = c(0, 1), ylim = c(0, 1), main = "(a) CV reach vs deployment reach")
graphics::box(col = "grey60")
graphics::points(pts, pch = 19, col = "grey35", cex = 0.9)
# a test point and its nearest training neighbour (short)
ti <- which.max(gobs > stats::quantile(gobs, 0.2) & gobs < stats::quantile(gobs, 0.5))
tr <- which(folds != folds[ti]); nn <- tr[which.min(colSums((t(pts[tr, ]) - pts[ti, ])^2))]
graphics::points(pts[ti, 1], pts[ti, 2], pch = 21, bg = "#d1495b", col = "black", cex = 1.6)
graphics::arrows(pts[ti, 1], pts[ti, 2], pts[nn, 1], pts[nn, 2], length = 0.06, lwd = 2, col = "#d1495b")
# a deployment grid cell far from any sample (long)
gi <- which.max(fpred); ns <- which.min(colSums((t(pts) - gg[gi, ])^2))
graphics::points(gg[gi, 1], gg[gi, 2], pch = 22, bg = "#00798c", col = "black", cex = 1.6)
graphics::arrows(gg[gi, 1], gg[gi, 2], pts[ns, 1], pts[ns, 2], length = 0.06, lwd = 2, col = "#00798c", lty = 1)
graphics::legend("topleft", c("sample (training)", "deployment grid",
                 "test point -> nearest train", "map cell -> nearest sample"),
       pch = c(19, 3, 21, 22), col = c("grey35", "grey70", "black", "black"),
       pt.bg = c(NA, NA, "#d1495b", "#00798c"), bty = "n", cex = 0.68)

## (b) the two ECDFs and the signed gap
xs <- seq(0, max(fpred) * 1.02, length.out = 400)
Fo <- ecdf(gobs)(xs); Fp <- ecdf(fpred)(xs)
graphics::par(mar = c(4.2, 4.4, 3, 1.4), cex.main = 0.95)
plot(xs, Fo, type = "n", xlab = "nearest-neighbour distance", ylab = "cumulative fraction",
     main = "(b) Leakage = the signed gap between reaches", ylim = c(0, 1))
graphics::polygon(c(xs, rev(xs)), c(Fo, rev(Fp)), col = "#edae4933", border = NA)
graphics::lines(xs, Fo, lwd = 3, col = "#d1495b")   # CV / Gobs = red
graphics::lines(xs, Fp, lwd = 3, col = "#00798c")   # deployment / Gpred = teal
graphics::legend("bottomright", c(expression(hat(G)[obs] ~ "(CV: test->train)"),
                 expression(hat(G)[pred] ~ "(deployment: map->sample)")),
       lwd = 3, col = c("#d1495b", "#00798c"), bty = "n", cex = 0.78)
graphics::text(stats::median(xs) * 1.05, 0.5, "CV reaches\nnearer\n(optimistic\nleakage)",
               cex = 0.72, col = "grey30")
dev.off()
cat("wrote figures/F0_concept.png\n")
