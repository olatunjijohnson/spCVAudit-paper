# Build Figure 5 for the MEE paper: the SDM case study (5 panels).
# ---------------------------------------------------------------------------
# READS cache/case_study_sdm.rds -- it does NOT recompute the audit. The
# previous version re-derived everything independently, which is how a figure and
# its own case study drift apart (the same duplication that let meta_audit.R and
# trend_analysis.R diverge). Run case_study_sdm.R first.
#
# Panels:
#  (a) per-point leakage map            -- WHERE the split leaks
#  (b) CV vs deployment NN-distance ECDFs -- WHY (the mechanism)
#  (c) signed SLI by candidate scheme   -- WHAT to do instead
#  (d) Brier: reported vs honest        -- WHAT it costs, in the ecologist's currency
#  (e) the whole assemblage             -- that it is not one cherry-picked species
#
# Run with: library(spLeakage); source("scripts/make_figure_sdm.R")

suppressMessages(library(spLeakage))
dir.create("figures", showWarnings = FALSE)

RDS <- "cache/case_study_sdm.rds"
if (!file.exists(RDS)) stop("Run cache/case_study_sdm.R first (missing ", RDS, ").")
res <- readRDS(RDS)
lk  <- res$leak; opt <- res$opt; rk <- res$rank; agg <- res$assemblage
FOCAL <- res$focal; REPS <- res$reps

teal <- "#00798c"; red <- "#d1495b"; grey <- "grey80"

draw <- function() {
  op <- par(mar = c(4.2, 4.2, 3, 1), mgp = c(2.4, 0.8, 0), cex.main = 1.12, font.main = 2)
  on.exit(par(op))
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 2, byrow = TRUE))

  ## (a) leakage map
  lp <- lk$leak_point; co <- lk$coords; keep <- !is.na(lp)
  rng <- max(abs(lp[keep]))
  pal <- grDevices::colorRampPalette(c(teal, "grey92", red))(100)
  idx <- as.integer(cut(lp[keep], breaks = seq(-rng, rng, length.out = 101),
                        include.lowest = TRUE))
  plot(co[keep, 1] / 1e3, co[keep, 2] / 1e3, col = pal[idx], pch = 19, cex = 0.8,
       asp = 1, xlab = "easting (km)", ylab = "northing (km)",
       main = "(a) Where the random split leaks")
  # Honest legend: on this survey EVERY site leaks (the transects are tightly
  # clustered, so no point has a negative contribution). A two-category
  # leaking/not-leaking key would advertise a class with no members.
  legend("topright", bty = "n", pch = 19, col = c(red, "grey92"),
         legend = c(sprintf("strongest leak (%+.2f)", max(lp[keep])),
                    sprintf("weakest leak (%+.2f)", min(lp[keep]))), cex = 0.8)
  legend("bottomleft", bty = "n", cex = 0.78, text.font = 3,
         legend = sprintf("all %d sites leak", sum(keep)))

  ## (b) NN-distance ECDFs
  e <- lk$ecdf
  plot(e$r / 1e3, e$Gobs, type = "s", col = red, lwd = 2.5, ylim = c(0, 1),
       xlab = "nearest-neighbour distance (km)", ylab = "ECDF",
       main = "(b) CV vs deployment reach")
  lines(e$r / 1e3, e$Gpred, type = "s", col = teal, lwd = 2.5)
  abline(v = lk$phi / 1e3, lty = 3, col = "grey50")
  legend("bottomright", bty = "n", lwd = 2.5, col = c(red, teal, "grey50"),
         lty = c(1, 1, 3), legend = c("CV (test->train)", "deployment (pred->sample)",
         expression(phi)), cex = 0.8)

  ## (c) signed SLI by scheme
  rr <- rk$ranking[match(c("random", "block", "buffered"), rk$ranking$scheme), ]
  cols <- ifelse(rr$sli_rho > 0.02, red, ifelse(rr$sli_rho < -0.02, teal, grey))
  bp <- barplot(rr$sli_rho, names.arg = rr$scheme, col = cols, border = NA,
                ylim = range(c(rr$sli_rho, 0)) * 1.3,
                ylab = expression(SLI[rho] ~ "(signed)"),
                main = "(c) Which scheme matches deployment")
  abline(h = 0, col = "grey40")
  text(bp, rr$sli_rho, sprintf("%+.2f", rr$sli_rho),
       pos = ifelse(rr$sli_rho >= 0, 3, 1), cex = 0.8, xpd = NA)
  text(bp[2], 0, "best\n(matched)", pos = 1, cex = 0.7, col = teal, xpd = NA)

  ## (d) Brier: reported (random) vs honest (block).
  ## Medians over the REPS random splits, with IQR -- NOT a single draw, so this
  ## panel agrees with (e) by construction. A single draw of awt06 moves between
  ## +32% and +38%; quoting one would contradict the assemblage panel.
  fd <- res$focal_draws
  vals <- c(stats::median(fd$brier_cv), stats::median(fd$brier_block))
  lo   <- c(stats::quantile(fd$brier_cv, .25), stats::quantile(fd$brier_block, .25))
  hi   <- c(stats::quantile(fd$brier_cv, .75), stats::quantile(fd$brier_block, .75))
  bp2 <- barplot(vals, names.arg = c("random CV\n(reported)", "block CV\n(honest)"),
                 col = c(red, teal), border = NA, ylim = c(0, max(hi) * 1.35),
                 ylab = "Brier score (lower = better)",
                 main = sprintf("(d) Reported vs honest (%s, %d splits)", FOCAL, REPS))
  arrows(bp2, lo, bp2, hi, angle = 90, code = 3, length = 0.04, col = "grey30")
  text(bp2, hi, sprintf("%.3f", vals), pos = 3, cex = 0.85, xpd = NA)
  text(mean(bp2), max(hi) * 1.24,
       sprintf("optimism %+.0f%%  [%+.0f, %+.0f]",
               100 * stats::median(fd$optimism_rel),
               100 * stats::quantile(fd$optimism_rel, .25),
               100 * stats::quantile(fd$optimism_rel, .75)), cex = 0.9, font = 2)

  ## (e) the whole assemblage: optimism per species, median over REPS splits
  a <- agg[order(agg$optimism), ]
  cols <- ifelse(a$optimism > 0.02, red, grey)
  bord <- ifelse(a$spid == FOCAL, "black", NA)
  par(mar = c(4.2, 5.2, 3, 1))
  bp3 <- barplot(100 * a$optimism, horiz = TRUE, names.arg = a$spid, las = 1,
                 col = cols, border = bord,
                 xlim = c(min(0, 100 * min(a$opt_lo)) - 3, 100 * max(a$opt_hi) * 1.12),
                 xlab = "optimism of random CV (% of honest Brier)",
                 main = sprintf("(e) The whole assemblage: %d of %d species optimistic (median %+.0f%%, %d random splits each)",
                                sum(a$optimism > 0), nrow(a),
                                100 * stats::median(a$optimism), REPS),
                 cex.names = 0.8, cex.main = 1.0)
  arrows(100 * a$opt_lo, bp3, 100 * a$opt_hi, bp3, angle = 90, code = 3,
         length = 0.02, col = "grey30", lwd = 1.1)
  abline(v = 0, col = "grey40")
  fi <- which(a$spid == FOCAL)
  lab <- sprintf(" %+.0f%%", 100 * a$optimism)
  lab[fi] <- sprintf(" %+.0f%%  <- worked in panels (a)-(d)", 100 * a$optimism[fi])
  text(100 * a$opt_hi, bp3, lab, pos = 4, cex = 0.72,
       font = ifelse(seq_len(nrow(a)) == fi, 3, 1), xpd = NA)
}

for (f in c("figures/F5_sdm.png")) {
  grDevices::png(f, width = 2400, height = 1750, res = 200)
  draw()
  grDevices::dev.off()
  message("wrote ", f)
}
