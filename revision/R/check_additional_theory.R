source("revision/R/common.R")
set.seed(728110)
rows <- list()
for (trial in seq_len(200)) {
  n <- sample(2:15, 1); xy <- matrix(runif(2 * n), n, 2)
  qxy <- matrix(runif(8), 4, 2)
  signal <- runif(1, .05, .95); nugget <- 1 - signal; range <- runif(1, .03, .5)
  S <- signal * correlation_matrix(xy, range = range) + diag(nugget, n)
  # The identities below assume marginal variance exactly one. BLAS rounding
  # in the squared-distance formula can leave tiny positive self-distances;
  # enforce the stated variance rather than relaxing the identity tolerance.
  diag(S) <- 1
  B <- signal * correlation_matrix(qxy, xy, range = range)
  full <- variance_reduction(S, B)
  for (i in seq_len(nrow(qxy))) {
    nn <- which.min(pair_distance(qxy[i, , drop = FALSE], xy))
    rest <- setdiff(seq_len(n), nn)
    s <- S[rest, nn]; b1 <- B[i, nn]
    residual_b <- B[i, rest] - s * b1
    conditional <- S[rest, rest, drop = FALSE] - tcrossprod(s)
    extra <- as.numeric(crossprod(residual_b, solve(conditional, residual_b)))
    bound <- min(sum(residual_b^2) / nugget, signal - b1^2)
    stopifnot(abs(full[i] - b1^2 - extra) < 1e-8, extra >= -1e-9,
              extra <= bound + 1e-8,
              min(eigen(conditional, symmetric = TRUE, only.values = TRUE)$values) >= nugget - 1e-9)
    rows[[length(rows) + 1L]] <- data.frame(trial, point = i, full = full[i],
      one_neighbour = b1^2, extra, bound, identity_error = abs(full[i] - b1^2 - extra))
  }
  cv_cor <- runif(20); target_cor <- runif(30)
  J1 <- mean(cv_cor) - mean(target_cor)
  J2 <- mean(cv_cor^2) - mean(target_cor^2)
  population_variance <- function(x) mean((x - mean(x))^2)
  decomposed <- (mean(cv_cor) + mean(target_cor)) * J1 +
    population_variance(cv_cor) - population_variance(target_cor)
  stopifnot(abs(J2 - decomposed) < 1e-12)
}
out <- do.call(rbind, rows)
write.csv(out, "revision/results/additional-theory-checks.csv", row.names = FALSE)
cat("Verified 800 Schur-complement remainders/bounds and 200 moment decompositions.\n")
cat("Maximum numerical identity error:", max(out$identity_error), "\n")
