# Ecological case study: species distribution models (presence/absence)
# ---------------------------------------------------------------------------
# Audits the spatial validation of SDMs for the Australian Wet Tropics bird
# assemblage from the Elith et al. (2006) / Elith et al. (2020) NCEAS benchmark,
# distributed in `disdat`: 340 structured survey sites, 20 real bird species,
# 13 environmental predictors, GDA94 / MGA zone 55.
#
# Replaces the previous blockCV-based script, which used a SIMULATED species over
# continental Australia and audited it with a coordinates-only IDW interpolator
# that never saw the covariates. Here the species are real and the learner is a
# random forest on the actual predictors -- i.e. an SDM an ecologist would fit.
#
# Deployment target: the 10,000 disdat background points, which carry covariates
# and span a WIDER extent than the survey sites -- a genuine wall-to-wall map.
#
# Design decisions (all load-bearing for the paper's claims):
#  1. REPEATED SPLITS. The split being audited is random, so a single draw is a
#     single realisation: awt06 moved between +32% and +38% across two seeds.
#     Every reported number is a median over REPS independent random 10-fold
#     draws, with an interquartile range.
#  2. REAL LEARNER. Random forest on the covariates, not the package default
#     `spLeakage:::.idw_predictor()`, which is coordinates-only and is documented as being
#     chosen because it "clearly exhibits leakage". Using it would let the
#     learner, not the data, drive the result.
#  3. SHARED GEOMETRY, SCOPED CLAIM. All 20 species share the same 340 sites, so
#     they are NOT 20 independent datasets: g_obs and f_pred are identical across
#     species and only the fitted variogram differs. The assemblage result is a
#     scoped claim about ONE survey, never a field-wide census.
#  4. FOCAL SPECIES = a worked example, NOT an "average" species. The full
#     assemblage is reported alongside, so the reader sees the whole spread and
#     nothing rests on the choice. awt06 is used because it is the best-powered
#     species (187 presences, balanced prevalence) and exercises every channel.
#  5. sf INPUT with a declared CRS, so the case study passes the package's own
#     hygiene check rather than tripping it in a printed scorecard.
#
# Reproducible; requires Suggests: disdat, ranger, sf.
# Run with: library(spLeakage); source("scripts/case_study_sdm.R")

suppressMessages({library(spLeakage); library(disdat); library(sf)})
set.seed(1)

REGION <- "AWT"
GROUP  <- "bird"
FOCAL  <- "awt06"
N_GRID <- 4000L
REPS   <- 20L
NCORE  <- max(1L, min(12L, parallel::detectCores() - 1L))
covn   <- disPredictors(REGION)

## ---- data: real presence/absence survey + matched environment --------------
pa <- disPa(REGION, group = GROUP)
ev <- disEnv(REGION, group = GROUP)
stopifnot(identical(pa$siteid, ev$siteid))
d0 <- cbind(pa[, c("siteid", "x", "y")], ev[, covn])
d0 <- sf::st_as_sf(d0, coords = c("x", "y"), crs = sf::st_crs(disCRS(REGION)))

spp  <- setdiff(names(pa), c("siteid", "x", "y", "group"))
prev <- vapply(spp, function(s) mean(pa[[s]] > 0), numeric(1))
spp_use <- spp[prev >= 0.15 & prev <= 0.60]   # drop near-absent / near-ubiquitous

message(sprintf("%s/%s: %d sites, %d species (%d usable); CRS: %s",
                REGION, GROUP, nrow(pa), length(spp), length(spp_use),
                sf::st_crs(d0)$input))

## ---- deployment target: wall-to-wall map over the background points ---------
bg  <- disBg(REGION)
bg  <- bg[stats::complete.cases(bg[, covn]), ]
gi  <- sample(nrow(bg), min(N_GRID, nrow(bg)))
tgt <- prediction_target(grid = as.matrix(bg[gi, c("x", "y")]), type = "grid")

## ---- the learner an ecologist would actually fit ---------------------------
rf_predictor <- function(covariates, response, ntree = 300L) {
  function(train, test) {
    tr <- sf::st_drop_geometry(train)[, c(response, covariates), drop = FALSE]
    tr[[response]] <- factor(tr[[response]], levels = c(0, 1))
    m <- ranger::ranger(
      stats::as.formula(paste(response, "~", paste(covariates, collapse = " + "))),
      data = tr, probability = TRUE, num.trees = ntree, seed = 1, verbose = FALSE)
    stats::predict(m, data = sf::st_drop_geometry(test), verbose = FALSE)$predictions[, "1"]
  }
}
rf <- rf_predictor(covn, "occ")

## ---- one audit of one random split, for one species ------------------------
audit_once <- function(s, seed) {
  set.seed(seed)
  d <- d0; d$occ <- as.integer(pa[[s]] > 0)
  folds <- sample(rep_len(1:10, nrow(d)))
  lk <- detect_leakage(d, split = folds, target = tgt, response = "occ")
  op <- estimate_optimism(d, folds, response = "occ", metric = "brier",
                          control = "block", predict_fun = rf)
  data.frame(spid = s, rep = seed, sli_rho = lk$SLI_rho,
             brier_cv = op$E_cv, brier_block = op$E_control,
             optimism_rel = op$optimism_rel)
}

## ---- (1) assemblage-wide, over REPS random splits per species --------------
message(sprintf("auditing %d species x %d random splits on %d cores...",
                length(spp_use), REPS, NCORE))
draws <- do.call(rbind, parallel::mclapply(spp_use, function(s)
  do.call(rbind, lapply(seq_len(REPS), function(i) audit_once(s, i))),
  mc.cores = NCORE))

agg <- do.call(rbind, lapply(split(draws, draws$spid), function(z) data.frame(
  spid = z$spid[1], prev = mean(as.integer(pa[[z$spid[1]]] > 0)),
  n_pres = sum(pa[[z$spid[1]]] > 0),
  sli_rho = stats::median(z$sli_rho),
  brier_cv = stats::median(z$brier_cv), brier_block = stats::median(z$brier_block),
  optimism = stats::median(z$optimism_rel),
  opt_lo = stats::quantile(z$optimism_rel, .25),
  opt_hi = stats::quantile(z$optimism_rel, .75),
  frac_optimistic = mean(z$optimism_rel > 0))))
agg <- agg[order(-agg$optimism), ]

## ---- (2) the focal species, worked through every channel -------------------
d <- d0; d$occ <- as.integer(pa[[FOCAL]] > 0)
set.seed(1); folds <- sample(rep_len(1:10, nrow(d)))   # one draw, for the maps

res <- list(assemblage = agg, draws = draws, focal = FOCAL, reps = REPS,
            focal_draws = draws[draws$spid == FOCAL, ])

res$leak  <- detect_leakage(d, split = folds, target = tgt, response = "occ", n_boot = 200)
res$opt   <- estimate_optimism(d, folds, response = "occ", metric = "brier",
                               control = "block", predict_fun = rf)
res$feat  <- detect_feature_leakage(d, folds, covariates = covn,
                                    newdata = bg[sample(nrow(bg), 3000), covn])
res$grp   <- tryCatch(detect_group_leakage(d, folds), error = function(e) NULL)
# Refit-free emulator: OUT OF SCOPE for this audit, and it now says so rather than
# guessing. Every learner in the emulator's calibration grid is a function of
# coordinates alone, so it has never seen a covariate-driven SDM. Unguarded it answers
# anyway and is wrong: it returned +10.1% (90% CI +3.1, +17.1) reporting itself in-AOA,
# against a refit route giving +36.9% with covariates and +56.3% coordinates-only.
# Declaring the covariates makes the refusal structural rather than a matter of luck
# with the AOA distance -- the same move as the `time` guard.
res$emul <- tryCatch(
  predict_optimism(d, split = folds, target = tgt, response = "occ", model = "rf",
                   covariates = covn),
  error = function(e) paste("emulator declined (by construction):",
                            conditionMessage(e)))
res$rec   <- recommend_validation(d, estimand = "prediction", design = "clustered",
                                  target = "grid")
res$rank  <- rank_cv_schemes(d, tgt, response = "occ")
res$audit <- audit_workflow(d, split = folds, target = tgt, response = "occ")

## ---- report ----------------------------------------------------------------
cat(sprintf("\n=== AWT bird assemblage: %d species, RF on covariates, median over %d random splits ===\n",
            nrow(agg), REPS))
print(data.frame(spid = agg$spid, prev = round(agg$prev, 2), n_pres = agg$n_pres,
                 sli_rho = round(agg$sli_rho, 3),
                 brier_cv = round(agg$brier_cv, 4),
                 brier_block = round(agg$brier_block, 4),
                 optimism = sprintf("%+.0f%% [%+.0f,%+.0f]", 100 * agg$optimism,
                                    100 * agg$opt_lo, 100 * agg$opt_hi),
                 always_opt = sprintf("%.0f%%", 100 * agg$frac_optimistic)),
      row.names = FALSE)
cat(sprintf("\noptimistic in %d / %d species (median split); median %+.1f%%, range [%+.0f%%, %+.0f%%]\n",
            sum(agg$optimism > 0), nrow(agg), 100 * stats::median(agg$optimism),
            100 * min(agg$optimism), 100 * max(agg$optimism)))
cat(sprintf("cor(SLI_rho, optimism) across the assemblage = %+.2f\n",
            stats::cor(agg$sli_rho, agg$optimism)))
cat("  ^ shared geometry: g_obs/f_pred identical across species, so SLI varies only\n",
    "  via the fitted variogram while optimism also scales with Vw^2. The index audits\n",
    "  a split for a given response; it does not rank responses. (Theory, eq. closed form.)\n")

fd <- res$focal_draws
cat(sprintf("\n=== focal species: %s (worked example; %d splits) ===\n", FOCAL, REPS))
cat(sprintf("  SLI_rho  : median %+.3f  IQR [%+.3f, %+.3f]\n", stats::median(fd$sli_rho),
            stats::quantile(fd$sli_rho, .25), stats::quantile(fd$sli_rho, .75)))
cat(sprintf("  optimism : median %+.1f%%  IQR [%+.1f%%, %+.1f%%]  (optimistic in %.0f%% of splits)\n",
            100 * stats::median(fd$optimism_rel), 100 * stats::quantile(fd$optimism_rel, .25),
            100 * stats::quantile(fd$optimism_rel, .75), 100 * mean(fd$optimism_rel > 0)))
cat("\n")
print(res$leak); cat("\n"); print(res$opt); cat("\n"); print(res$feat); cat("\n")
if (!is.null(res$grp))  { print(res$grp);  cat("\n") }
if (!is.null(res$emul)) { print(res$emul); cat("\n") }
print(res$rank); cat("\n"); print(res$audit)

saveRDS(res, "cache/case_study_sdm.rds")
invisible(res)
