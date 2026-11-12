# Reconstruct the exact RNG path before Gaussian draws in the two fitted studies.
# Distinct fingerprints establish absence of exact stream resets/reuse, not a
# mathematical proof that a pseudorandom generator supplies independent draws.
source("revision/R/local_libraries.R")
source("revision/R/common.R")
stopifnot(requireNamespace("CAST", quietly = TRUE), requireNamespace("digest", quietly = TRUE))
historical <- "--historical" %in% commandArgs(trailingOnly = TRUE)
suffix <- if (historical) "-final" else "-rng"
studies <- paste0(c("sensitivity", "nested-estimation"), suffix)
records <- list(); moments <- list()
workers <- as.integer(Sys.getenv("SLURM_NTASKS", "2"))
for (study in studies) {
  grid <- read.csv(file.path("revision/results", paste0(study, "-configuration.csv")))
  out <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
    cfg <- grid[i, ]; set.seed(cfg$seed)
    xy <- sampling_locations(cfg$n, cfg$design)
    target <- target_locations(cfg$target, 15L)
    ignored <- make_candidates(xy)
    points <- function(x) sf::st_as_sf(as.data.frame(x), coords = 1:2, crs = 3857)
    km <- suppressMessages(CAST::knndm(points(xy), predpoints = points(target), k = 5L))
    post_cast_state <- digest::digest(.Random.seed, algo = "sha256")
    if (!historical) set.seed(cfg$response_seed)
    state <- digest::digest(.Random.seed, algo = "sha256")
    innovations <- rnorm(cfg$n)
    list(record = data.frame(study, id = cfg$id, seed = cfg$seed, n = cfg$n,
      post_cast_state_sha256 = post_cast_state, state_sha256 = state,
      innovation_sha256 = digest::digest(innovations, algo = "sha256")),
      innovations = innovations)
  }, mc.cores = workers, mc.set.seed = FALSE)
  d <- do.call(rbind, lapply(out, `[[`, "record"))
  write.csv(d, file.path("revision/results", paste0("random-stream-", study, ".csv")), row.names = FALSE)
  if (!historical) stopifnot(!anyDuplicated(d$state_sha256), !anyDuplicated(d$innovation_sha256), !anyDuplicated(d$seed))
  m <- do.call(rbind, lapply(out, `[[`, "innovations"))
  # Consecutive configurations are independent draws in the declared protocol.
  lag_cor <- cor(as.vector(m[-nrow(m), ]), as.vector(m[-1, ]))
  moments[[study]] <- data.frame(study, configurations = nrow(m), draws = length(m),
    unique_post_cast_states = length(unique(d$post_cast_state_sha256)),
    unique_response_states = length(unique(d$state_sha256)),
    unique_innovation_vectors = length(unique(d$innovation_sha256)),
    mean = mean(m), variance = var(as.vector(m)), adjacent_stream_correlation = lag_cor)
  records[[study]] <- d
}
d <- do.call(rbind, records)
if (!historical) stopifnot(!anyDuplicated(d$state_sha256), !anyDuplicated(d$innovation_sha256))
tag <- if (historical) "-historical" else ""
write.csv(d, paste0("revision/results/random-stream-audit", tag, ".csv"), row.names = FALSE)
write.csv(do.call(rbind, moments), paste0("revision/results/random-stream-moments", tag, ".csv"), row.names = FALSE)
cat(if (historical) "Historical duplicate-stream diagnosis:" else "Verified distinct response streams:", nrow(d), "configurations.\n")
print(do.call(rbind, moments))

# Independently simulate from a conditional target distribution to check the
# risk calculation for a fixed predictor that differs from the conditional mean.
set.seed(1234501L)
x <- sampling_locations(20, "random"); q <- target_locations("domain", 4)
allx <- rbind(x, q); n <- nrow(x); nq <- nrow(q)
Sall <- .7 * exp(-pair_distance(allx) / .2) + diag(.3, n + nq)
Strain <- Sall[seq_len(n), seq_len(n)]
B <- Sall[n + seq_len(nq), seq_len(n)]
Sqq <- Sall[n + seq_len(nq), n + seq_len(nq)]
y <- as.numeric(t(chol(Strain)) %*% rnorm(n))
conditional_mean <- as.numeric(B %*% solve(Strain, y))
conditional_covariance <- Sqq - B %*% solve(Strain, t(B))
prediction <- rep(mean(y), nq)
formula <- mean((prediction - conditional_mean)^2 + diag(conditional_covariance))
draws <- sweep(t(chol(conditional_covariance)) %*% matrix(rnorm(nq * 20000L), nq), 1, conditional_mean, "+")
loss <- colMeans((draws - prediction)^2)
se <- sd(loss) / sqrt(length(loss))
check <- data.frame(analytic_risk = formula, simulated_mean = mean(loss), mcse = se,
                     z_score = (mean(loss) - formula) / se)
write.csv(check, "revision/results/conditional-risk-check.csv", row.names = FALSE)
print(check)
stopifnot(abs(check$z_score) < 4.5)
