# Listing 1 from the paper, verbatim and runnable.
# Audits the SDM validation of the case study (one focal species, one split); the
# case study itself (scripts/case_study_sdm.R) averages 20 splits over all species.
# Requires: spLeakage, disdat, ranger, sf. Produces SLI_rho = +0.17 and ~ +37% optimism.

library(spLeakage); library(disdat)
set.seed(1)
covs <- disPredictors("AWT")

# the data an ecologist already has (Australian Wet Tropics birds)
pa <- disPa("AWT", group = "bird")            # presence/absence surveys
ev <- disEnv("AWT", group = "bird")           # environment at those sites
d  <- sf::st_as_sf(cbind(pa[, c("siteid", "x", "y")], ev[, covs]),
                   coords = c("x", "y"), crs = sf::st_crs(disCRS("AWT")))
d$occ <- as.integer(pa$awt06 > 0)             # one focal species
folds <- sample(rep_len(1:10, nrow(d)))       # the naive random split
bg  <- disBg("AWT")                           # background points = where the map is used
tgt <- prediction_target(grid = as.matrix(bg[, c("x", "y")]), type = "grid")

lk <- detect_leakage(d, split = folds, target = tgt,      # signed leakage index
                     response = "occ", n_boot = 200)
lk; plot(lk)                                              # score, CI, ECDFs, leak map

# a random forest on the predictors -- the ecologist's own model
rf <- function(train, test) {
  tr <- sf::st_drop_geometry(train); tr$occ <- factor(tr$occ, levels = c(0, 1))
  m  <- ranger::ranger(occ ~ ., data = tr[, c("occ", covs)], probability = TRUE)
  predict(m, sf::st_drop_geometry(test))$predictions[, "1"]
}
estimate_optimism(d, folds, response = "occ", metric = "brier",   # % inflation (refit)
                  control = "block", predict_fun = rf)            #   -> ~ +37%

audit_workflow(d, split = folds, target = tgt,            # all channels at once
               response = "occ")                          # (pass time= if time-stamped)
recommend_validation(d, estimand = "prediction",          # what to do instead
                     design = "clustered", target = "grid")
