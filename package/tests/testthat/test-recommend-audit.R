test_that("probability sampling does not certify a predictive CV scheme", {
  rec <- recommend_validation(estimand = "prediction", design = "probability",
                              target = "grid")
  expect_true(is.na(rec$spatial_cv_appropriate))
  expect_match(rec$rationale, "smaller training")
})

test_that("clustering alone does not certify spatial CV", {
  rec <- recommend_validation(estimand = "prediction", design = "clustered",
                              target = "grid")
  expect_true(is.na(rec$spatial_cv_appropriate))
  expect_match(rec$avoid, "proof")
})

test_that("population estimand + probability sample -> design-based, not spatial", {
  rec <- recommend_validation(estimand = "population", design = "probability",
                              target = "grid")
  expect_false(rec$spatial_cv_appropriate)
  expect_match(paste(rec$recommended, collapse = " "), "[Dd]esign-based")
})

test_that("unknown design remains unknown", {
  rec <- recommend_validation(estimand = "prediction", design = "unknown",
                              target = "grid")
  expect_equal(rec$assumed_design, "unknown")
  expect_true(is.na(rec$spatial_cv_appropriate))
})

test_that("clustering flag is computed from data but does not set the design", {
  d <- sim_spatial_clustered(n = 400L)
  rec <- recommend_validation(d, estimand = "prediction", design = "unknown",
                              target = "grid", coords = c("x", "y"))
  expect_true(rec$clustering$clustered)       # geometry IS clustered
  expect_equal(rec$design, "unknown")         # but design stays unknown
})

test_that("audit_workflow and report_leakage assemble a scorecard", {
  d <- sim_spatial_clustered(n = 400L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  au <- audit_workflow(d, split = random_folds(nrow(d)), target = tgt,
                       response = "z", coords = c("x", "y"))
  expect_s3_class(au, "workflow_audit")
  expect_true(au$flags$leakage_grade %in% c("A", "B", "C", "D", "F"))

  lk <- detect_leakage(d, random_folds(nrow(d)), tgt, response = "z", coords = c("x", "y"))
  opt <- estimate_optimism(d, random_folds(nrow(d)), response = "z", coords = c("x", "y"))
  rec <- recommend_validation(estimand = "prediction", design = "clustered", target = "grid")
  rep <- report_leakage(lk, optimism = opt, recommendation = rec)
  expect_s3_class(rep, "leakage_report")
  expect_output(print(rep), "spLeakage report")
})

test_that("audit_workflow adds temporal channels when a time column is supplied", {
  d <- sim_spatial_clustered(n = 300L)
  d$year <- sample(2000:2018, nrow(d), replace = TRUE)   # space-time data
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  fold <- random_folds(nrow(d))
  au_sp <- audit_workflow(d, split = fold, target = tgt, response = "z",
                          coords = c("x", "y"))
  expect_true(is.na(au_sp$flags$lookahead_frac))          # no time -> not computed
  expect_null(au_sp$temporal)

  au_st <- audit_workflow(d, split = fold, target = tgt, response = "z",
                          coords = c("x", "y"), time = "year")
  expect_s3_class(au_st$temporal, "temporal_leakage")
  expect_s3_class(au_st$st, "st_leakage")
  expect_true(is.finite(au_st$flags$lookahead_frac))
  # a random split of time-stamped data trains on the future for most test points
  expect_gt(au_st$flags$lookahead_frac, 0.5)
  expect_output(print(au_st), "temporal lookahead")
})
