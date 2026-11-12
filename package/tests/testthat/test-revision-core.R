test_that("crossing distances distinguish signed cancellation from transport", {
  ar <- spLeakage:::.sli_areas(c(0, 2), 1)
  expect_equal(ar$A, 0)
  expect_equal(ar$W, 1)
  expect_true(ar$crosses)
  expect_false(spLeakage:::.sli_areas(c(0, 1), c(2, 3))$crosses)
  expect_false(spLeakage:::.sli_areas(c(0, 1), c(0, 1))$crosses)
})

test_that("nonphysical variogram fits are refitted instead of clipped", {
  h <- seq(.05, .6, length.out = 15)
  vg <- data.frame(dist = h, gamma = -.1 + 1 - exp(-h / .2), np = 40)
  fitted <- spLeakage:::.fit_exponential(vg)
  co <- fitted$coef
  fitted_loss <- sum(vg$np * (vg$gamma - co$nugget - co$psill * (1 - exp(-h / co$range)))^2)
  clipped_loss <- sum(vg$np * (vg$gamma - (1 - exp(-h / .2)))^2)
  expect_gte(co$nugget, 0)
  expect_gte(co$psill, 0)
  expect_gt(co$range, 0)
  expect_lt(fitted_loss, clipped_loss / 2)
  expect_null(fitted$vcov)
})

test_that("a constant response produces zero fitted signal without a fabricated covariance", {
  vg <- data.frame(dist = seq(.05, .6, length.out = 15), gamma = 0, np = 40)
  fitted <- spLeakage:::.fit_exponential(vg)
  expect_equal(fitted$coef$nugget + fitted$coef$psill, 0)
  expect_gt(fitted$coef$range, 0)
  expect_null(fitted$vcov)
})

test_that("bounded fit roundoff cannot produce a negative nugget", {
  # Synthetic nested-CV pilot: the bounded optimizer returned about -1e-18.
  vg <- data.frame(
    dist = c(.0174166149201063, .044482346267583, .0724806242276635,
      .102209213265044, .131067302260208, .161505552556981, .188739273472322,
      .218720854717328, .247872245454855, .276013318698493, .304983047134764,
      .334498562271149, .365647502587037, .394962954944214, .421207796404485),
    gamma = c(.0718540819172004, .140432722074197, .288194101055274,
      .363569650800861, .532038688563101, .87964126980158, 1.29405343401255,
      1.21526250010067, 1.28581986039573, 1.65346264369948, 1.64357533204122,
      1.35462202838517, 2.16886647849344, 1.89285627569227, 1.40666906480299),
    np = c(45, 108, 96, 95, 78, 77, 75, 84, 86, 61, 74, 62, 52, 57, 67))
  fit <- spLeakage:::.fit_exponential(vg)
  expect_gte(fit$coef$nugget, 0)
  expect_gte(fit$coef$psill, 0)
  expect_gt(fit$coef$range, 0)
  expect_equal(fit$coef$nugget, 0, tolerance = 1e-12)
  expect_lte(fit$coef$psill / (fit$coef$psill + fit$coef$nugget), 1)
})

test_that("constrained dependence fitting respects coordinate and response units", {
  h <- seq(.05, .6, length.out = 15)
  vg <- data.frame(dist = h, gamma = -.1 + 1 - exp(-h / .2), np = 40)
  small <- spLeakage:::.fit_exponential(vg)$coef
  vg$dist <- vg$dist * 1000
  vg$gamma <- vg$gamma * 100
  large <- spLeakage:::.fit_exponential(vg)$coef
  expect_equal(large$range / 1000, small$range, tolerance = 1e-6)
  expect_equal(large$nugget / 100, small$nugget, tolerance = 1e-6)
  expect_equal(large$psill / 100, small$psill, tolerance = 1e-6)
})

test_that("geographic nearest neighbours remain geographic without s2", {
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))
  query <- matrix(c(-179, 89), nrow = 1)
  ref <- rbind(c(179, 89), c(-179, 85))
  d <- spLeakage:::.nn_dist(query, ref, list(geographic = TRUE, crs = sf::st_crs(4326)))
  expect_lt(d, 10000)
})

test_that("geographic buffer caps use metres and preserve the requested buffer", {
  grid <- expand.grid(x = seq(-3, -1, length.out = 12), y = seq(52, 53, length.out = 10))
  points <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)
  xy <- as.matrix(grid)
  info <- list(geographic = TRUE, crs = sf::st_crs(points))
  extent <- spLeakage:::.domain_diameter(xy, info)
  expect_gt(extent, 100000)
  expect_lt(extent, 200000)
  expect_equal(extent, max(as.numeric(sf::st_distance(points))))
  dep <- list(rho = function(h) exp(-h / 5000), practical_range = 15000)
  set.seed(42)
  base <- spLeakage:::.parse_split(spatial_block_cv(xy, 3), nrow(xy))
  reference <- lapply(base, function(f) {
    dd <- sf::st_distance(points[f$train, ], points[f$test, ])
    nearest <- apply(matrix(as.numeric(dd), length(f$train)), 1, min)
    list(test = f$test, train = f$train[nearest > 15000])
  })
  expect_lt(sum(lengths(lapply(reference, `[[`, "train"))), sum(lengths(lapply(base, `[[`, "train"))))
  target <- prediction_target(grid = points)
  expected <- detect_leakage(points, reference, target, dependence = dep)
  ranked <- rank_cv_schemes(points, target, k = 3, schemes = "buffered", dependence = dep, seed = 42)
  expect_equal(ranked$ranking$W, expected$W, tolerance = 1e-10)
})

test_that("default geographic IDW respects proximity across the date line", {
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))
  predict <- spLeakage:::.idw_predictor(c("x", "y"), "z",
    info = list(geographic = TRUE, crs = sf::st_crs(4326)))
  train <- data.frame(x = c(179, -179), y = c(89, 85), z = c(1, 0))
  expect_gt(predict(train, data.frame(x = -179, y = 89)), .999)
  planar <- spLeakage:::.idw_predictor(c("x", "y"), "z")
  expect_equal(planar(data.frame(x = c(0, 2), y = 0, z = c(0, 1)),
                      data.frame(x = 1, y = 0)), .5)
})

test_that("pair subsampling preserves pair selection without allocating all pairs", {
  xy <- cbind(seq_len(10), 0)
  set.seed(719)
  expected <- utils::combn(10, 2)[, sample.int(45, 12), drop = FALSE]
  set.seed(719)
  actual <- spLeakage:::.pair_geometry(xy, list(geographic = FALSE), max_pairs = 12)
  expect_identical(rbind(actual$i, actual$j), expected)
  # More than five billion possible pairs: only the requested 500 are created.
  large <- spLeakage:::.pair_geometry(cbind(seq_len(100001), 0),
                                      list(geographic = FALSE), max_pairs = 500)
  expect_length(large$h, 500)
  expect_true(all(large$i < large$j))
  expect_equal(large$h, large$j - large$i)
  expect_equal(nrow(unique(cbind(large$i, large$j))), 500)
  single <- spLeakage:::.pair_geometry(xy, list(geographic = FALSE), max_pairs = 1)
  expect_length(single$h, 1)
})

test_that("incompatible declared coordinate systems cannot silently mix", {
  d <- sf::st_as_sf(data.frame(x = c(-2, -1, 0, 1), y = 53), coords = c("x", "y"), crs = 4326)
  target <- prediction_target(grid = sf::st_transform(d, 3857))
  dep <- list(rho = function(h) exp(-h / 10000), practical_range = 30000)
  expect_error(detect_leakage(d, c(1, 1, 2, 2), target, dependence = dep), "CRS differ")
  expect_error(prediction_target(grid = matrix(numeric(), ncol = 2)), "at least one")
})

test_that("first and second correlation moments can have opposite signs", {
  cv <- rep(0.5, 5)
  target <- c(rep(0.8, 3), rep(0.01, 2))
  expect_equal(mean(cv) - mean(target), 0.016)
  expect_equal(mean(cv^2) - mean(target^2), -0.13404)
  expect_true(spLeakage:::.sli_areas(-log(cv), -log(target))$crosses)
})

test_that("independent co-located noise does not create perfect correlation", {
  d <- sim_spatial(n = 120, noise = 0.5)
  dep <- estimate_dependence(d, "z", coords = c("x", "y"))
  expect_equal(dep$rho(0), dep$psill / (dep$psill + dep$nugget))
  expect_equal(dep$rho(0), dep$rho(1e-12), tolerance = 1e-8)
})

test_that("pure observation noise gives no correlation contrast at repeats", {
  d <- data.frame(x = c(0, 0, 1, 1), y = 0)
  dep <- list(rho = function(h) rep(0, length(h)), practical_range = 1)
  target <- prediction_target(grid = matrix(c(0.5, 0), nrow = 1))
  lk <- detect_leakage(d, c(1, 2, 1, 2), target,
                       dependence = dep, coords = c("x", "y"))
  expect_equal(lk$SLI_rho, 0)
  expect_equal(lk$SLI_rho2, 0)
  expect_gt(lk$W, 0)
  expect_output(print(lk), "not estimates of predictive bias")
})

test_that("the squared contrast uses consistent correlation scaling", {
  d <- data.frame(x = c(0, 1, 2, 3), y = 0)
  dep <- list(rho = function(h) 0.4 * exp(-h), practical_range = 3)
  target <- prediction_target(grid = matrix(c(0.5, 0), nrow = 1))
  lk <- detect_leakage(d, c(1, 2, 1, 2), target,
                       dependence = dep, coords = c("x", "y"))
  expect_equal(sli(lk, "rho2"), 0.4^2 * (exp(-2) - exp(-1)))
})

test_that("unavailable parameter covariance does not fabricate an interval", {
  d <- data.frame(x = c(0, 1, 2, 3), y = 0)
  dep <- list(rho = function(h) 0.4 * exp(-h), practical_range = 3,
              coef = list(nugget = 0.6, psill = 0.4, range = 1), vcov = NULL)
  target <- prediction_target(grid = matrix(c(0.5, 0), nrow = 1))
  expect_warning(lk <- detect_leakage(d, c(1, 2, 1, 2), target,
    dependence = dep, coords = c("x", "y"), n_boot = 20), "no interval returned")
  expect_null(lk$SLI_rho_ci)
  expect_null(lk$boot)
})

test_that("invalid splits cannot silently reuse observations", {
  d <- sim_spatial(30)
  dep <- list(rho = function(h) exp(-h), practical_range = 3)
  target <- prediction_target(grid = grid_target(3))
  audit <- function(f) detect_leakage(d, f, target, dependence = dep, coords = c("x", "y"))
  expect_error(audit(list(test = 1:10, train = 5:30)), "must not overlap")
  expect_error(audit(list(1:10, 1:10)), "more than once")
  expect_error(audit(list(test = c(1, 31))), "valid training and test")
})

test_that("target candidate generation does not reset the caller random stream", {
  d <- sim_spatial(80)
  xy <- spLeakage:::.extract_coords(d, c("x", "y"))
  info <- list(crs = NA, geographic = FALSE)
  dep <- list(rho = function(h) exp(-h), practical_range = 3)
  target <- prediction_target(grid = grid_target(4))
  run <- function() spLeakage:::.matched_control(d, target, "z", c("x", "y"),
                                                 4, xy, info, dep)
  set.seed(88)
  run(); after_first <- .Random.seed
  run(); after_second <- .Random.seed
  expect_false(identical(after_first, after_second))
  set.seed(88)
  run()
  expect_identical(.Random.seed, after_first)
})
