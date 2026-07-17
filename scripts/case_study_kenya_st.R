# Spatio-temporal case study: Pf malaria prevalence mapping, Kenya
# ---------------------------------------------------------------------------
# The data are SPATIO-TEMPORAL (surveys 1985-2016), so a purely spatial audit is
# incomplete. This script runs the geographic, grouped/co-location, temporal
# (lookahead) and joint space-time channels, and shows how the correct validation
# depends on the DECLARED estimand/target (contemporary map vs forecast) -- the
# paper's thesis, made operational on real data.
#
# WHY KENYA AND NOT NIGERIA. The earlier version used Nigeria (n = 66), which was not a
# scientific choice but a data-availability artifact: MAP returns 2,672 Nigeria records,
# of which 2,476 are DHS surveys whose cluster coordinates MAP cannot share, leaving
# only 66 with both coordinates and prevalence. That ceiling is hard without DHS
# credentials. A scan of MAP's public geolocated PfPR found Kenya (1,607 records at
# 1,230 sites) and Tanzania (1,429 at 657) far better served.
#
# NOTE the reason for the swap is NOT that Kenya has a spatial signal where Nigeria did
# not -- pooled over 31 years Kenya's signal proportion is also ~0. It is that n = 1,607
# lets us DEMONSTRATE why, which n = 66 never could. The epoch ladder below shows the
# signal is destroyed by POOLING, not by small n: it stays at ~0 for 1985-2016
# (n=1607), 2005+ (n=973) and 2010+ (n=319), and only reappears in a single survey year
# (2010, n=255: w-hat = 0.35, range 29 km -- a plausible malaria correlation range).
# The Nigeria case could assert this; the Kenya case can prove it.
#
# Design notes:
#  - REPEATED SPLITS. The split under audit is random, so every reported number is a
#    median over REPS independent random 10-fold draws with an interquartile range.
#  - The lookahead fraction is reported but NOT headlined: a random split of any
#    time-stamped data puts nearly every test point after some training point, so the
#    number is close to arithmetic and says little about these data. The joint
#    space-time index is the informative quantity, and it is the one no marginal
#    channel recovers.
#
# Reproducible offline from cache/kenya_pr.rds (built by the getPR call below).
# Run with: library(spLeakage); source("scripts/case_study_kenya_st.R")

suppressMessages({library(spLeakage); library(sf)})
set.seed(1)

REPS <- 20L
CACHE <- "cache/kenya_pr.rds"

if (!file.exists(CACHE)) {
  suppressMessages(library(malariaAtlas))
  d <- suppressWarnings(malariaAtlas::getPR(country = "Kenya", species = "Pf"))
  ok <- !is.na(d$longitude) & !is.na(d$latitude) & !is.na(d$pr) & !is.na(d$year_start)
  saveRDS(as.data.frame(d)[ok, c("longitude", "latitude", "pr", "examined",
                                 "positive", "year_start")], CACHE)
}
g <- readRDS(CACHE)
sfd <- sf::st_as_sf(g, coords = c("longitude", "latitude"), crs = 4326)
xy  <- sf::st_coordinates(sfd)
dat <- data.frame(x = xy[, 1], y = xy[, 2], pr = g$pr, year = g$year_start)

message(sprintf("Kenya Pf prevalence: n=%d at %d sites, %d-%d",
                nrow(g), nrow(unique(g[, 1:2])), min(g$year_start), max(g$year_start)))

## ---- deployment target: a contemporary wall-to-wall prevalence surface ------
gr <- as.matrix(expand.grid(
  x = seq(min(xy[, 1]), max(xy[, 1]), length.out = 45),
  y = seq(min(xy[, 2]), max(xy[, 2]), length.out = 45)))
tgt <- prediction_target(grid = gr, type = "grid")

dep <- estimate_dependence(sfd, response = "pr")

res <- list(n = nrow(g), n_sites = nrow(unique(g[, 1:2])), reps = REPS,
            years = range(g$year_start), dependence = dep,
            nn_index = spLeakage:::.nn_index(xy, list(geographic = TRUE)))

## ---- one audit of one random split, across every channel --------------------
audit_once <- function(seed) {
  set.seed(seed)
  folds <- sample(rep_len(1:10, nrow(sfd)))
  lk <- detect_leakage(sfd, split = folds, target = tgt, response = "pr",
                       dependence = dep)
  gl <- detect_group_leakage(dat, folds, coords = c("x", "y"))
  tl <- detect_temporal_leakage(dat, folds, time = "year")
  st <- detect_st_leakage(dat, split = folds, time = "year", coords = c("x", "y"))
  data.frame(rep = seed, sli_rho = lk$SLI_rho, c_obs = lk$c_obs, c_pred = lk$c_pred,
             group_frac = gl$frac_leaked, lookahead = tl$lookahead_frac,
             sli_space = st$sli_space, sli_time = st$sli_time, sli_st = st$sli_st)
}
draws <- do.call(rbind, lapply(seq_len(REPS), audit_once))
res$draws <- draws

med <- function(v) c(m = stats::median(v), lo = stats::quantile(v, .25),
                     hi = stats::quantile(v, .75))

## ---- is the lost signal about n, or about pooling? --------------------------
# The epoch ladder. If w-hat ~ 0 were a small-n problem it would recover as n grows;
# instead it is flat at ~0 across n = 1607/973/319 and only returns within a single
# survey year. Pooling epochs of differing method and endemicity is what destroys it.
res$epoch_ladder <- do.call(rbind, lapply(
  list(list("1985-2016", rep(TRUE, nrow(g))), list("2005+", g$year_start >= 2005),
       list("2010+", g$year_start >= 2010), list("2010 only", g$year_start == 2010)),
  function(e) {
    d <- tryCatch(estimate_dependence(sfd[e[[2]], ], response = "pr"),
                  error = function(err) NULL)
    if (is.null(d)) return(NULL)
    data.frame(window = e[[1]], n = sum(e[[2]]), w_hat = d$signal_prop,
               range_km = d$practical_range / 1000)
  }))

## ---- the fix follows from the declaration ----------------------------------
# Each candidate scheme against every channel. The point of the joint space-time
# channel is here: NEITHER marginal fix repairs it. Forward-chaining removes lookahead
# entirely yet leaves the joint index high, because co-located repeat visits still sit
# close in both space and time across the time blocks; group k-fold removes co-location
# entirely yet leaves the joint index essentially untouched, because it does nothing
# about time. A leak that is only visible jointly needs a scheme that is aware of both.
set.seed(1)
schemes <- list(random = sample(rep_len(1:10, nrow(dat))),
                temporal_kfold = temporal_kfold(dat, k = 10, time = "year"),
                group_kfold = group_kfold(dat, k = 10, coords = c("x", "y")))
res$fixes <- do.call(rbind, lapply(names(schemes), function(nm) {
  f <- schemes[[nm]]
  data.frame(scheme = nm,
             joint_st = detect_st_leakage(dat, split = f, time = "year",
                                          coords = c("x", "y"))$sli_st,
             co_location = detect_group_leakage(dat, f, coords = c("x", "y"))$frac_leaked,
             lookahead = detect_temporal_leakage(dat, f, time = "year")$lookahead_frac)
}))

# Declaration B: the target is a CONTEMPORARY map -> the pooling is itself the problem;
# restrict to a single survey year and the audit becomes a real, purely spatial one.
res$n_recent <- sum(g$year_start == 2010)
sub <- sfd[g$year_start == 2010, ]
set.seed(1); fr <- sample(rep_len(1:10, nrow(sub)))
res$recent_dep  <- estimate_dependence(sub, response = "pr")
res$recent_leak <- detect_leakage(sub, split = fr, target = tgt, response = "pr")

## ---- report ----------------------------------------------------------------
cat(sprintf("\n=== Kenya Pf prevalence: n=%d at %d sites, %d-%d (%d random splits) ===\n",
            res$n, res$n_sites, res$years[1], res$years[2], REPS))
cat(sprintf("  pooled signal proportion w-hat : %.3f  (pure nugget; see epoch ladder)\n",
            dep$signal_prop))
cat(sprintf("  pooled practical range         : %.0f km\n", dep$practical_range / 1000))
cat(sprintf("  Clark-Evans nn_index           : %.3f\n", res$nn_index))

cat("\n--- epoch ladder: is w-hat ~ 0 about n, or about POOLING? ---\n")
print(transform(res$epoch_ladder, w_hat = round(w_hat, 3),
                range_km = round(range_km)), row.names = FALSE)
cat("  ^ flat at ~0 as n falls 1607 -> 973 -> 319, then the signal RETURNS at n=255\n")
cat("    within a single survey year. Not a sample-size problem: pooling epochs of\n")
cat("    differing method and endemicity is what destroys the spatial signal.\n")

cat("\n--- naive random 10-fold, four channels (median [IQR] over splits) ---\n")
for (v in c("sli_rho", "group_frac", "lookahead", "sli_space", "sli_time", "sli_st")) {
  s <- med(draws[[v]])
  cat(sprintf("  %-11s %+7.3f  [%+.3f, %+.3f]\n", v, s[1], s[2], s[3]))
}
cat("\n  ^ the JOINT space-time index dwarfs both margins: a random split routinely\n")
cat("    places a training record close in BOTH space and time (the same site a year\n")
cat("    apart), while a forecast must predict a future time across all locations.\n")
cat("    Neither marginal channel reveals it. The lookahead fraction is near-arithmetic\n")
cat("    for any random split of time-stamped data and is reported, not headlined.\n")

cat("\n--- NEITHER marginal fix repairs the joint space-time leak ---\n")
print(transform(res$fixes, joint_st = sprintf("%+.3f", joint_st),
                co_location = sprintf("%.1f%%", 100 * co_location),
                lookahead = sprintf("%.1f%%", 100 * lookahead)), row.names = FALSE)
cat("  ^ temporal_kfold() drives lookahead to 0% yet leaves the joint index high (repeat\n")
cat("    visits still sit close in BOTH space and time across the time blocks);\n")
cat("    group_kfold() drives co-location to 0% and barely moves it (it ignores time).\n")
cat("    A leak visible only jointly needs a scheme aware of both -- which is the case\n")
cat("    for measuring the joint channel rather than assuming a marginal fix suffices.\n")

cat("\n--- the fix follows from the DECLARED target ---\n")
cat(sprintf("  contemporary map -> restrict to a single survey year (2010, n=%d):\n", res$n_recent))
cat(sprintf("      w-hat %.3f (signal recovered), range %.0f km, SLI_rho %+.3f\n",
            res$recent_dep$signal_prop, res$recent_dep$practical_range / 1000,
            res$recent_leak$SLI_rho))
cat("      -- now a real spatial audit, on genuine autocorrelation, to correct with block CV.\n")

saveRDS(res, "cache/case_study_kenya_st.rds")
invisible(res)
