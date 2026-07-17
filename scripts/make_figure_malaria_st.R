# Figure: the Kenya Pf-malaria data are spatio-temporal, and a purely spatial audit is
# incomplete. READS cache/case_study_kenya_st.rds -- it does not recompute the audit,
# so the figure and the case study cannot drift. Run case_study_kenya_st.R first.
#
# Panels:
#  (a) survey locations coloured by year (31-year pooling)
#  (b) the epoch ladder: w-hat is destroyed by POOLING, not by small n
#  (c) multi-channel audit under naive random CV -- the joint space-time index dwarfs
#      both margins
#  (d) no marginal fix repairs the joint leak
#
# Run with: Rscript scripts/make_figure_malaria_st.R

suppressMessages({library(spLeakage); library(sf)})

RDS <- "cache/case_study_kenya_st.rds"
if (!file.exists(RDS)) stop("Run cache/case_study_kenya_st.R first (missing ", RDS, ").")
res <- readRDS(RDS)
g   <- readRDS("cache/kenya_pr.rds")
d   <- res$draws
med <- function(v) stats::median(v)

teal <- "#00798c"; red <- "#d1495b"; amber <- "#edae49"; grey <- "grey80"

draw <- function() {
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3, 1), mgp = c(2.5, 0.8, 0),
            cex.main = 1.1, font.main = 2)
  on.exit(par(op))

  ## (a) locations by year
  yr <- g$year_start
  pal <- grDevices::colorRampPalette(c(teal, amber, red))(diff(range(yr)) + 1)
  plot(g$longitude, g$latitude, col = pal[yr - min(yr) + 1], pch = 19, cex = 0.45,
       asp = 1, xlab = "longitude", ylab = "latitude",
       main = sprintf("(a) %d surveys, %d sites, %d-%d", nrow(g),
                      res$n_sites, res$years[1], res$years[2]))
  legend("topright", bty = "n", pch = 19, col = c(teal, amber, red), cex = 0.8,
         legend = c(res$years[1], round(mean(res$years)), res$years[2]))

  ## (b) epoch ladder: pooling, not n, destroys the signal
  el <- res$epoch_ladder
  cols <- ifelse(el$w_hat > 0.05, teal, grey)
  bp <- barplot(el$w_hat, names.arg = sprintf("%s\n(n=%d)", el$window, el$n),
                col = cols, border = NA, ylim = c(0, max(el$w_hat) * 1.45),
                ylab = expression("estimated signal proportion  " * hat(w)),
                main = "(b) Pooling destroys the signal, not small n",
                cex.names = 0.7)
  text(bp, el$w_hat, sprintf("%.2f", el$w_hat), pos = 3, cex = 0.8, xpd = NA)
  text(bp, el$w_hat, sprintf("%.0f km", el$range_km), pos = 3, cex = 0.62,
       offset = 1.4, col = "grey35", xpd = NA)
  text(mean(bp[1:3]), max(el$w_hat) * 1.32, "n falls 5x, signal stays 0",
       cex = 0.72, font = 3, col = "grey30", xpd = NA)

  ## (c) four channels under naive random CV
  v <- c(`spatial\nSLI` = med(d$sli_rho), `co-location\nfraction` = med(d$group_frac),
         `space-only\nSLI` = med(d$sli_space), `time-only\nSLI` = med(d$sli_time),
         `JOINT\nspace-time` = med(d$sli_st))
  cols <- c(red, red, grey, amber, red)
  bp2 <- barplot(v, col = cols, border = NA, las = 1, cex.names = 0.62,
                 ylim = c(min(0, min(v)) - 0.05, max(v) * 1.3),
                 ylab = "signed index / fraction",
                 main = "(c) The joint channel dwarfs both margins")
  abline(h = 0, col = "grey40")
  text(bp2, v, sprintf("%+.2f", v), pos = ifelse(v >= 0, 3, 1), cex = 0.78, xpd = NA)
  text(bp2[5], max(v) * 1.18,
       sprintf("%.0fx the larger margin", med(d$sli_st) / med(d$sli_time)),
       cex = 0.72, font = 3, col = "grey30", xpd = NA)

  ## (d) no marginal fix repairs the joint leak
  f <- res$fixes
  bp3 <- barplot(f$joint_st, names.arg = sub("_", "\n", f$scheme),
                 col = c(red, amber, red), border = NA,
                 ylim = c(0, max(f$joint_st) * 1.4), cex.names = 0.75,
                 ylab = "joint space-time SLI",
                 main = "(d) Neither marginal fix repairs it")
  text(bp3, f$joint_st, sprintf("%+.3f", f$joint_st), pos = 3, cex = 0.8, xpd = NA)
  text(bp3, 0, sprintf("co-loc %.0f%%\nahead %.0f%%", 100 * f$co_location,
                       100 * f$lookahead), pos = 3, cex = 0.6, col = "white", xpd = NA)
  text(mean(bp3), max(f$joint_st) * 1.28,
       "each fix clears its own channel, not the joint one",
       cex = 0.7, font = 3, col = "grey30", xpd = NA)
}

for (fp in c("figures/F6_malaria_st.png")) {
  grDevices::png(fp, width = 2200, height = 1800, res = 200)
  draw()
  grDevices::dev.off()
  message("wrote ", fp)
}
