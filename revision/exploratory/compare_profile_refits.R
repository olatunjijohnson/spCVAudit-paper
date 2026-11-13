# Compare retained fitted studies with the profile-search refits declared in
# revision/audit/profile-refit-protocol.json. Runs the unchanged analysis scripts
# on both versions and restores the retained summaries afterwards.
# Usage: Rscript revision/exploratory/compare_profile_refits.R
out_dir <- "revision/results/profile-refit"
dir.create(out_dir, showWarnings = FALSE)
summaries <- "revision/results/summaries"
backup <- tempfile("summaries-"); dir.create(backup)
retained_files <- c(list.files(summaries, full.names = TRUE), "revision/results/nested-analysis.rds")
retained_files <- retained_files[file.exists(retained_files)]
stopifnot(file.copy(retained_files, backup))
on.exit({
  file.copy(file.path(backup, basename(retained_files)), dirname(retained_files), overwrite = TRUE)
  unlink(backup, recursive = TRUE)
}, add = TRUE)

runs <- data.frame(
  label = c("nested-rng", "nested-rng-profile", "nested-k10", "nested-k10-profile",
            "sensitivity-rng", "sensitivity-rng-profile"),
  script = rep(c("revision/R/analyse_nested_estimation.R", "revision/R/analyse_sensitivity.R"), c(4, 2)),
  variable = rep(c("REVISION_NESTED_SUFFIX", "REVISION_ESTIMATION_SUFFIX"), c(4, 2)),
  suffix = c("-rng", "-rng-profile", "-k10", "-k10-profile", "-rng", "-rng-profile"),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]
  status <- system2("Rscript", c("--vanilla", r$script), env = paste0(r$variable, "=", r$suffix),
                    stdout = file.path(out_dir, paste0(r$label, ".log")), stderr = file.path(out_dir, paste0(r$label, ".log")))
  if (status != 0) stop(r$label, " analysis failed; see ", file.path(out_dir, paste0(r$label, ".log")))
  target <- file.path(out_dir, r$label); dir.create(target, showWarnings = FALSE)
  prefix <- if (startsWith(r$label, "nested")) "^nested-" else "^(sensitivity-|covariance-fits)"
  stopifnot(file.copy(list.files(summaries, prefix, full.names = TRUE), target, overwrite = TRUE))
}

# Headline selection results, defined as in revision/R/make_manuscript_numbers.R.
cell <- "range=0.05;signal=0.95;design=random;target=shifted"
nested_headline <- function(label) {
  d <- file.path(out_dir, label)
  s <- read.csv(file.path(d, "nested-selection.csv"))
  p <- read.csv(file.path(d, "nested-paired-comparisons.csv"))
  inf <- read.csv(file.path(d, "nested-cell-influence.csv"))
  key <- p$learner == "estimated" & p$criterion == "J2_fitted" & p$baseline == "W"
  all <- p[key & p$stratum == "all", ]
  base <- s[s$stratum == "learner=estimated;criterion=W", "mean_absolute_error"]
  data.frame(study = label, mae_W = base,
    mae_J2 = s[s$stratum == "learner=estimated;criterion=J2_fitted", "mean_absolute_error"],
    delta = all$mean_absolute_error_difference, mcse = all$mcse, lower = all$lower, upper = all$upper,
    reduction_percent = -100 * all$mean_absolute_error_difference / base,
    strong_cell_delta = p[key & p$stratum == cell, "mean_absolute_error_difference"],
    other_cells_delta = inf[inf$omitted_cell == cell, "remaining_mean_difference"])
}
headline <- do.call(rbind, lapply(grep("^nested", runs$label, value = TRUE), nested_headline))
write.csv(headline, file.path(out_dir, "nested-headlines.csv"), row.names = FALSE)
sens <- do.call(rbind, lapply(c("sensitivity-rng", "sensitivity-rng-profile"), function(label)
  cbind(study = label, read.csv(file.path(out_dir, label, "sensitivity-headlines.csv")))))
write.csv(sens, file.path(out_dir, "sensitivity-headlines.csv"), row.names = FALSE)

# Fit-level changes, matched on configuration and fit name.
fit_change <- function(retained, profile, keys) {
  m <- merge(retained, profile, by = keys, suffixes = c("_retained", "_profile"))
  stopifnot(nrow(m) == nrow(retained), nrow(m) == nrow(profile))
  changed <- (abs(m$range_profile - m$range_retained) > 1e-10 * pmax(1, m$range_retained) |
    abs(m$signal_profile - m$signal_retained) > 1e-10) %in% TRUE
  list(rows = m, summary = data.frame(fits = nrow(m), changed = sum(changed),
    changed_fraction = mean(changed),
    median_abs_log_range_change = if (any(changed)) median(abs(log(m$range_profile[changed] / m$range_retained[changed]))) else 0,
    max_abs_signal_change = max(abs(m$signal_profile - m$signal_retained), na.rm = TRUE)))
}
fits <- list()
for (s in c("rng", "k10")) {
  a <- readRDS(sprintf("revision/results/nested-estimation-%s.rds", s))$fits
  b <- readRDS(sprintf("revision/results/nested-estimation-%s-profile.rds", s))$fits
  fc <- fit_change(a, b, c("id", "name"))
  fits[[s]] <- cbind(study = paste0("nested-", s), fc$summary)
}
sens_fit <- function(stem) {
  z <- readRDS(stem)$results
  z <- z[z$scheme == "random" & z$diagnostic_model == "estimated_isotropic", ]
  data.frame(id = z$id, name = "full", range = z$fitted_range, signal = z$fitted_signal)
}
fc <- fit_change(sens_fit("revision/results/sensitivity-rng.rds"),
                 sens_fit("revision/results/sensitivity-rng-profile.rds"), c("id", "name"))
fits$sensitivity <- cbind(study = "sensitivity-rng", fc$summary)
fits <- do.call(rbind, fits)
write.csv(fits, file.path(out_dir, "fit-changes.csv"), row.names = FALSE)

cat("\nVariogram fits changed by the profile search:\n"); print(fits, row.names = FALSE)
cat("\nFitted |J2| versus W, estimated covariance:\n"); print(headline, row.names = FALSE, digits = 4)
cat("\nSensitivity headlines:\n"); print(sens[sens$diagnostic_model == "estimated_isotropic", ], row.names = FALSE, digits = 4)
