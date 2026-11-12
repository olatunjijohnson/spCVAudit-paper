source("revision/R/common.R")

rows <- list()
add <- function(name, observed, expected, tolerance = 1e-10) {
  stopifnot(is.finite(observed), abs(observed - expected) <= tolerance)
  rows[[length(rows) + 1L]] <<- data.frame(check = name, observed = observed,
                                           expected = expected, tolerance = tolerance)
}

# A one-observation conditional variance can be evaluated directly.
v <- 1.7; w <- 0.6; h <- 0.2; a <- 0.3
c0 <- w * exp(-h / a)
ev <- variance_reduction(matrix(v, 1), matrix(v * c0, 1))
add("one-neighbour explained variance", ev, v * c0^2)

# Distinct co-located replicates, with independent measurement errors.
sigma2 <- 0.7; tau2 <- 0.3
for (k in c(1, 2, 5, 10)) {
  S <- matrix(sigma2, k, k) + diag(tau2, k)
  e <- variance_reduction(S, matrix(sigma2, 1, k))
  add(paste("co-located replicates", k), e, k * sigma2^2 / (k * sigma2 + tau2))
}

# A crossing distribution reverses the sign after squaring.
cv <- rep(0.5, 5); target <- c(rep(0.8, 3), rep(0.01, 2))
di <- distance_summary(-log(cv), -log(target), function(h) exp(-h))
add("crossing first moment", di[["J1"]], 0.016)
add("crossing second moment", di[["J2"]], -0.13404)
add("crossing flag", di[["crosses"]], 1)

# Equal mean distances are not equal distributions.
di <- distance_summary(c(0, 2), 1, function(h) exp(-h))
add("signed area cancellation", di[["A"]], 0)
add("unsigned transport despite cancellation", di[["W"]], 1)

# Actual nested training sets: the held-out point lies between two observations,
# but deployment lies beyond the sample. Both nearest distances are one.
xy <- cbind(c(-1, 0, 1), 0); q <- matrix(c(2, 0), 1)
R <- correlation_matrix(xy, range = 1)
S <- 0.8 * R + diag(0.2, 3)
ec <- variance_reduction(S[c(1, 3), c(1, 3)], S[2, c(1, 3), drop = FALSE])
eq <- variance_reduction(S, 0.8 * correlation_matrix(q, xy, range = 1))
di <- distance_summary(1, 1, function(h) 0.8 * exp(-h))
add("identical nearest distances first moment", di[["J1"]], 0)
add("identical nearest distances second moment", di[["J2"]], 0)
stopifnot(ec - eq > 0.01)
counterexample <- data.frame(J1 = di[["J1"]], J2 = di[["J2"]],
                             full_optimism = ec - eq, cv_ev = ec, target_ev = eq)
write.csv(counterexample, "revision/results/neighbourhood-counterexample.csv", row.names = FALSE)

# Exhaustive leave-one-out CV on the same three observations gives a stronger
# example: every observation is tested once, as in the main estimand.
loo <- lapply(seq_len(nrow(xy)), function(i) list(test = i, train = setdiff(seq_len(nrow(xy)), i)))
loo_result <- gaussian_scheme(xy, q, loo, range = 1, signal = .8)
for (name in c("A", "W", "J1", "J2", "U1", "U2")) add(paste("full LOO equal-distance", name), loo_result[[name]], 0)
stopifnot(loo_result[["exact"]] > .02)
write.csv(as.data.frame(as.list(loo_result)), "revision/results/loo-neighbourhood-counterexample.csv", row.names = FALSE)

# Ordering and transport inequalities over many reproducible distributions.
set.seed(11401)
for (i in seq_len(100)) {
  cv <- rexp(20); target <- cv + runif(20)
  w <- runif(1); a <- runif(1, 0.05, 2)
  z <- distance_summary(cv, target, function(h) w * exp(-h / a))
  stopifnot(z[["J1"]] >= -1e-12, z[["J2"]] >= -1e-12,
            abs(z[["A"]] - z[["W"]]) < 1e-10)
  cv <- rexp(17); target <- rexp(23)
  z <- distance_summary(cv, target, function(h) w * exp(-h / a))
  stopifnot(abs(z[["J2"]]) <= 2 * w^2 / a * z[["W"]] + 1e-12,
            abs(z[["J2"]]) <= z[["U2"]] + 1e-12,
            abs(z[["J1"]]) <= z[["U1"]] + 1e-12,
            abs(z[["A"]]) <= z[["W"]] + 1e-12)
}
add("ordered and crossing-distribution inequality trials", 200, 200)

# Cross-check the independent reference implementation against the package.
pkgload::load_all("package", quiet = TRUE)
set.seed(11703)
d <- sampling_locations(40, "random"); q <- target_locations("domain", 7)
folds <- make_candidates(d)$random
dep <- list(rho = function(h) 0.6 * exp(-h / 0.2), practical_range = 0.6)
lk <- detect_leakage(d, folds, prediction_target(grid = q), dependence = dep)
ref <- gaussian_scheme(d, q, folds, range = 0.2, signal = 0.6)
for (pair in list(c("A", "A"), c("W", "W"), c("SLI_rho", "J1"),
                  c("SLI_rho2", "J2"), c("W_rho", "U1"))) {
  add(paste("package/reference", pair[1]), lk[[pair[1]]], ref[[pair[2]]], 1e-7)
}
out <- do.call(rbind, rows)
save_run(out, "theory-checks")
print(out, row.names = FALSE)
cat("All algebraic checks passed. Original single-holdout counterexample:\n")
print(counterexample, row.names = FALSE)
cat("Exhaustive leave-one-out counterexample:\n")
print(loo_result)
