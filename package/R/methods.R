# S3 methods and accessors for leakage_diagnosis objects.

#' Extract the Spatial Leakage Index
#'
#' @param x A `leakage_diagnosis` object.
#' @param type `"rho"` (dependence form, default), `"rho2"` (squared-correlation
#'   contrast for `k = 1`), or `"d"` (distance form).
#' @return A single numeric SLI value.
#' @examples
#' set.seed(1)
#' n <- 150
#' d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
#' g <- as.matrix(expand.grid(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)))
#' lk <- detect_leakage(d, split = sample(rep_len(1:10, n)),
#'                      target = prediction_target(grid = g, type = "grid"),
#'                      response = "z", coords = c("x", "y"))
#' sli(lk)             # dependence form: signed, bounded, the headline index
#' sli(lk, type = "d") # distance form: the bridge to NNDM's W
#' @export
sli <- function(x, type = c("rho", "d", "rho2")) {
  stopifnot(inherits(x, "leakage_diagnosis"))
  type <- match.arg(type)
  switch(type, rho = x$SLI_rho, d = x$SLI_d, rho2 = x$SLI_rho2)
}

# Human-readable verdict from the signed index.
.leak_verdict <- function(sli_rho, delta) {
  dir <- if (sli_rho > 1e-12) "higher mean correlation in CV"
         else if (sli_rho < -1e-12) "lower mean correlation in CV"
         else "zero signed contrast; inspect distributions"
  cross <- if (abs(delta) < 1 - 1e-12 && abs(delta) > 0)
    " (ECDFs cross; inspect plot)" else ""
  paste0(dir, cross)
}

#' @export
print.leakage_diagnosis <- function(x, ...) {
  cat("<leakage_diagnosis>\n")
  cat(sprintf("  target            : %s   |  n = %d, test = %d, folds = %d%s\n",
              x$target_type, x$n, x$n_test, x$n_folds,
              if (isTRUE(x$anisotropic)) "  [anisotropic]" else ""))
  ci_r <- if (!is.null(x$SLI_rho_ci))
    sprintf("  90%% parameter interval [%+.3f, %+.3f]", x$SLI_rho_ci[1], x$SLI_rho_ci[2]) else ""
  ci_d <- if (!is.null(x$SLI_d_ci))
    sprintf("  90%% parameter interval [%+.3f, %+.3f]", x$SLI_d_ci[1], x$SLI_d_ci[2]) else ""
  cat(sprintf("  SLI_rho (signed)  : %+.3f   [%s]%s\n",
              x$SLI_rho, .leak_verdict(x$SLI_rho, x$delta), ci_r))
  cat(sprintf("  SLI_d  (signed)   : %+.3f   (A = %+.4g, phi = %.4g)%s\n",
              x$SLI_d, x$A, x$phi, ci_d))
  cat(sprintf("  retained corr.    : c_obs = %.3f vs c_pred = %.3f\n",
              x$c_obs, x$c_pred))
  cat(sprintf("  W (NNDM) / delta  : %.4g / %+.2f\n", x$W, x$delta))
  if (!is.null(x$SLI_rho2) && is.finite(x$SLI_rho2))
    cat(sprintf("  squared contrast  : %+.3f (one-neighbour model reference)\n", x$SLI_rho2))
  if (isTRUE(x$crosses)) cat("  distance ECDFs cross: signed contrasts can cancel.\n")
  cat("  These are mismatch diagnostics, not estimates of predictive bias.\n")
  invisible(x)
}

#' @export
summary.leakage_diagnosis <- function(object, ...) {
  print(object)
  if (length(object$fold_sli) > 1L) {
    cat("\n  per-fold contributions (mean rho excess):\n")
    fs <- object$fold_sli
    for (k in seq_along(fs)) {
      cat(sprintf("    fold %-3s : %+.3f\n", names(fs)[k], fs[k]))
    }
  }
  invisible(object)
}

#' Plot a leakage diagnosis
#'
#' @param x A `leakage_diagnosis` object.
#' @param which `"ecdf"` (observed vs target NN-distance ECDFs) or `"map"`
#'   (per-point contribution to the correlation contrast, not prediction error).
#' @param ... Passed to the underlying plot call.
#' @return `x`, invisibly.
#' @export
plot.leakage_diagnosis <- function(x, which = c("ecdf", "map"), ...) {
  which <- match.arg(which)
  if (which == "ecdf") {
    e <- x$ecdf
    plot(e$r, e$Gobs, type = "s", col = "#d1495b", lwd = 2,
         xlab = "nearest-neighbour distance", ylab = "ECDF",
         main = "Observed (test->train) vs target (pred->sample)",
         ylim = c(0, 1), ...)
    graphics::lines(e$r, e$Gpred, type = "s", col = "#00798c", lwd = 2)
    graphics::abline(v = x$phi, lty = 3, col = "grey50")
    graphics::legend("bottomright", bty = "n",
                     col = c("#d1495b", "#00798c", "grey50"),
                     lwd = c(2, 2, 1), lty = c(1, 1, 3),
                     legend = c("observed (CV)", "target (deployment)", "phi"))
  } else {
    co <- x$coords; lp <- x$leak_point
    keep <- !is.na(lp)
    rng <- max(abs(lp[keep]))
    if (!is.finite(rng) || rng == 0) rng <- 1
    pal <- grDevices::colorRampPalette(c("#00798c", "grey85", "#d1495b"))(100)
    idx <- as.integer(cut(lp[keep], breaks = seq(-rng, rng, length.out = 101),
                          include.lowest = TRUE))
    # extra right margin for the colour key
    op <- graphics::par(mar = graphics::par("mar") + c(0, 0, 0, 6))
    on.exit(graphics::par(op), add = TRUE)
    plot(co[keep, 1], co[keep, 2], col = pal[idx], pch = 19, cex = 1.1,
         xlab = "x", ylab = "y", asp = 1,
         main = "Validation correlation relative to the target mean", ...)
    graphics::points(co[keep, 1], co[keep, 2], pch = 21, col = "grey40",
                     cex = 1.1, lwd = 0.4)   # thin outline for contrast

    # vertical colour key in the right margin, red = leaking (test easier than
    # deployment), teal = pessimistic (test harder), grey = well matched.
    usr <- graphics::par("usr")
    kx <- usr[2] + (usr[2] - usr[1]) * 0.04          # key bar x position
    kw <- (usr[2] - usr[1]) * 0.03
    ky <- seq(usr[3] + (usr[4] - usr[3]) * 0.18,
              usr[4] - (usr[4] - usr[3]) * 0.18, length.out = 100)
    graphics::rect(kx, ky[-100], kx + kw, ky[-1], col = pal, border = NA, xpd = NA)
    graphics::rect(kx, ky[1], kx + kw, ky[100], border = "grey40", xpd = NA)
    tx <- kx + kw * 1.5
    graphics::text(tx, ky[100], sprintf("%+.2f", rng), adj = 0, cex = 0.72, xpd = NA)
    graphics::text(tx, mean(ky), "0", adj = 0, cex = 0.72, xpd = NA)
    graphics::text(tx, ky[1], sprintf("%+.2f", -rng), adj = 0, cex = 0.72, xpd = NA)
    graphics::text(tx, ky[100], "higher", adj = c(0, 1.8), cex = 0.68,
                   col = "#d1495b", font = 2, xpd = NA)
    graphics::text(tx, ky[1], "lower", adj = c(0, -0.9), cex = 0.68,
                   col = "#00798c", font = 2, xpd = NA)
    graphics::text(kx, usr[4], expression("SLI"[rho] * " contribution"),
                   adj = c(0, -0.6), cex = 0.72, font = 2, xpd = NA)
  }
  invisible(x)
}
