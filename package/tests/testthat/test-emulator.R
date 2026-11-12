test_that("the shipped emulator loads and is well-formed", {
  em <- spLeakage:::.get_emulator()
  expect_s3_class(em, "optimism_emulator")
  expect_true(em$n_train > 200)
  expect_true(all(c("sli_rho", "range_rel", "signal", "nn_index") %in% em$features))
  expect_setequal(unique(em$model), c("idw", "rf", "gam"))
  expect_setequal(unique(em$response_type), c("gaussian", "poisson", "binomial"))
})

# A central, in-distribution feature vector for a given (model, response) cell:
# the column-wise median of that cell's (unstandardised) training features.
central_feats <- function(em, model = "idw", rt = "gaussian") {
  ix <- which(em$model == model & em$response_type == rt)
  X <- t(t(em$Xs[ix, , drop = FALSE]) * em$scale + em$center)
  stats::setNames(apply(X, 2L, stats::median), em$features)
}

test_that("an in-distribution feature vector gets an in-AOA finite estimate", {
  skip_if_not_installed("xgboost")
  em <- spLeakage:::.get_emulator()
  res <- spLeakage:::.emu_predict_one(em, central_feats(em), "idw", "gaussian")
  expect_true(res$in_aoa)
  expect_true(is.finite(res$optimism_rel))
  expect_true(res$optimism_rel > -1 && res$optimism_rel < 1)
})

test_that("higher SLI predicts higher optimism (the C1 -> C2 relationship)", {
  skip_if_not_installed("xgboost")
  em <- spLeakage:::.get_emulator()
  raw <- t(em$Xs) * em$scale + em$center
  # Compare real training rows, NOT a synthetic query that holds every other feature
  # at its median while pushing sli_rho to a tail. The features are not independent
  # -- clustered designs carry high sli_rho, so cor(sli_rho, nn_index) ~= -0.5 -- and
  # such a marginal query lands off the manifold, where the AOA guard correctly
  # refuses to answer. Varying sli_rho along the data is the honest version of this
  # test, and the guard refusing the off-manifold version is a feature, not a bug.
  o <- order(raw["sli_rho", ])
  lo <- raw[.EMU_FEATURES, o[round(0.1 * length(o))]]
  hi <- raw[.EMU_FEATURES, o[round(0.9 * length(o))]]
  r_lo <- spLeakage:::.emu_predict_one(em, lo, "idw", "gaussian")
  r_hi <- spLeakage:::.emu_predict_one(em, hi, "idw", "gaussian")
  expect_true(r_lo$in_aoa && r_hi$in_aoa)
  expect_gt(r_hi$optimism_rel, r_lo$optimism_rel)
})

test_that("emulator refuses an extreme out-of-AOA query", {
  skip_if_not_installed("xgboost")
  em <- spLeakage:::.get_emulator()
  x <- central_feats(em); x["range_rel"] <- 50; x["sli_rho"] <- 0.9
  res <- spLeakage:::.emu_predict_one(em, x, "idw", "gaussian")
  expect_false(res$in_aoa)
})

test_that("predict_optimism runs end-to-end and returns a valid object", {
  skip_if_not_installed("xgboost")
  d <- sim_gp(n = 300L, range = 0.2, signal = 0.6, design = "clustered")
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  # May fall outside the AOA (different generator than the training study); the
  # AOA-refusal warning is expected behaviour, so allow it.
  pr <- suppressWarnings(predict_optimism(d, random_folds(nrow(d)), tgt,
                         response = "z", coords = c("x", "y"), model = "idw"))
  expect_s3_class(pr, "optimism_prediction")
  expect_equal(pr$response_type, "gaussian")
  expect_true(is.logical(pr$in_aoa))
  if (pr$in_aoa) expect_true(is.finite(pr$optimism_rel))
})

# --- scope guards (added with the space-time guard) ---------------------------

test_that("predict_optimism refuses space-time queries by construction", {
  skip_if_not_installed("xgboost")
  set.seed(1)
  d <- data.frame(x = runif(80), y = runif(80), year = sample(2000:2010, 80, TRUE))
  d$z <- rnorm(80)
  tgt <- prediction_target(grid = as.matrix(expand.grid(x = seq(0, 1, .2),
                                                        y = seq(0, 1, .2))),
                           type = "grid")
  folds <- sample(rep_len(1:5, nrow(d)))
  # The AOA cannot see a time dimension, so the guard must be structural.
  expect_error(
    predict_optimism(d, split = folds, target = tgt, response = "z",
                     coords = c("x", "y"), time = "year"),
    "purely spatial"
  )
  # Same structural reason: the calibration grid's learners are coordinate-only, so
  # "the response is covariate-driven" is invisible to the AOA too.
  expect_error(
    predict_optimism(d, split = folds, target = tgt, response = "z",
                     coords = c("x", "y"), covariates = c("a", "b")),
    "coordinates alone"
  )
})

test_that(".emulator_features conforms to the .EMU_FEATURES contract", {
  set.seed(1)
  d <- data.frame(x = runif(60), y = runif(60)); d$z <- rnorm(60)
  tgt <- prediction_target(grid = as.matrix(expand.grid(x = seq(0, 1, .25),
                                                        y = seq(0, 1, .25))),
                           type = "grid")
  f <- .emulator_features(d, sample(rep_len(1:5, 60)), tgt, "z", c("x", "y"))
  # sli_d is deliberately absent: unstable (phi -> 0) and not in the training table.
  expect_identical(names(f), .EMU_FEATURES)
  expect_false("sli_d" %in% names(f))
})
