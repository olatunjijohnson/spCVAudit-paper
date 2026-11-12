# Post-hoc descriptive analysis declared in EXPERIMENTS.md on 29 October 2026.
# It does not change the primary MAE comparison or discard any configuration.
selection <- read.csv("revision/results/summaries/nested-selection-by-field.csv")
a <- subset(selection, criterion == "J2_fitted")
b <- subset(selection, criterion == "W")
key <- function(d) paste(d$id, d$learner, sep = ":")
b <- b[match(key(a), key(b)), ]
stopifnot(identical(key(a), key(b)), nrow(a) == 3840L)
factors <- c("range", "signal", "design", "target")
cell_key <- function(d) apply(d[factors], 1, function(x) paste(x, collapse = ":"))
field <- a[c("id", factors, "learner", "target_risk")]
field$cell <- cell_key(field)
field$J2_selected <- a$selected; field$W_selected <- b$selected
field$J2_ties <- a$n_tied; field$W_ties <- b$n_tied
field$J2_absolute_error <- a$absolute_error; field$W_absolute_error <- b$absolute_error
field$difference <- a$absolute_error - b$absolute_error
field$J2_signed_error <- a$score_error; field$W_signed_error <- b$score_error
stopifnot(all(is.finite(field$difference)), max(abs(a$target_risk - b$target_risk)) < 1e-12)

summarise_difference <- function(d, label) {
  cells <- split(d$difference, d$cell)
  stopifnot(length(unique(lengths(cells))) == 1L)
  mcse <- sqrt(sum(vapply(cells, function(x) length(x) * var(x), 0.0))) / nrow(d)
  # Removing one replicate leaves all factor cells equally weighted. This is
  # different from simply dropping a row from the pooled mean.
  cell_mean <- vapply(cells, mean, 0.0)
  cell_n <- lengths(cells)
  omitted <- mean(d$difference) - (d$difference - cell_mean[d$cell]) /
    (length(cells) * (cell_n[d$cell] - 1L))
  data.frame(learner = d$learner[1], stratum = label, n = nrow(d),
    mean_difference = mean(d$difference), mcse,
    J2_lower_loss = sum(d$difference < -1e-12),
    equal_loss = sum(abs(d$difference) <= 1e-12),
    W_lower_loss = sum(d$difference > 1e-12),
    same_selected_set = sum(d$J2_selected == d$W_selected),
    J2_random_only = sum(d$J2_selected == "random"),
    W_random_only = sum(d$W_selected == "random"),
    minimum = min(d$difference), q05 = unname(quantile(d$difference, .05)),
    q25 = unname(quantile(d$difference, .25)), median = median(d$difference),
    q75 = unname(quantile(d$difference, .75)), q95 = unname(quantile(d$difference, .95)),
    maximum = max(d$difference),
    single_replicate_omission_min = min(omitted), single_replicate_omission_max = max(omitted))
}

schemes <- c("random", "block", "buffered", "knndm")
choice_weights <- function(selected) {
  t(vapply(strsplit(selected, ";", fixed = TRUE), function(s) {
    stopifnot(all(s %in% schemes), !anyDuplicated(s))
    as.numeric(schemes %in% s) / length(s)
  }, numeric(length(schemes))))
}
summaries <- frequencies <- list()
for (learner in unique(field$learner)) {
  d <- field[field$learner == learner, ]
  groups <- c(list(all = seq_len(nrow(d))), split(seq_len(nrow(d)), d$cell))
  for (label in names(groups)) {
    z <- d[groups[[label]], ]
    summaries[[length(summaries) + 1L]] <- summarise_difference(z, label)
    for (criterion in c("J2", "W")) {
      selected <- z[[paste0(criterion, "_selected")]]
      weight <- choice_weights(selected)
      stopifnot(max(abs(rowSums(weight) - 1)) < 1e-12)
      frequencies[[length(frequencies) + 1L]] <- data.frame(learner, stratum = label,
        criterion, scheme = schemes, n = nrow(z), selection_probability = colMeans(weight),
        effective_selections = colSums(weight), unique_selections = vapply(schemes,
          function(s) sum(selected == s), integer(1)))
    }
  }
}
summaries <- do.call(rbind, summaries)
frequencies <- do.call(rbind, frequencies)

# The independent descriptive calculation must reproduce the primary paired
# comparison for both predictors, including its stratified Monte Carlo SE.
primary <- read.csv("revision/results/summaries/nested-paired-comparisons.csv")
primary <- subset(primary, criterion == "J2_fitted" & baseline == "W" & stratum == "all")
check <- summaries[summaries$stratum == "all", ]
primary <- primary[match(check$learner, primary$learner), ]
stopifnot(max(abs(check$mean_difference - primary$mean_absolute_error_difference)) < 1e-12,
          max(abs(check$mcse - primary$mcse)) < 1e-12)

raw <- readRDS("revision/results/nested-estimation-rng.rds")$results
scores <- aggregate(raw[c("conditional_target_risk", "cv_score", "score_error", "known_expected_cv")],
  raw[c(factors, "learner", "scheme")], mean)
scores$absolute_error <- aggregate(abs(raw$score_error), raw[c(factors, "learner", "scheme")], mean)$x
write.csv(field, "revision/results/summaries/nested-behaviour-by-field.csv", row.names = FALSE)
write.csv(summaries, "revision/results/summaries/nested-behaviour-by-cell.csv", row.names = FALSE)
write.csv(frequencies, "revision/results/summaries/nested-choice-frequencies.csv", row.names = FALSE)
write.csv(scores, "revision/results/summaries/nested-candidate-scores.csv", row.names = FALSE)
print(check, row.names = FALSE)
cat("Post-hoc descriptions retain all 24 factor cells and all 1,920 configurations per predictor.\n",
    "Losses and choice probabilities average uniformly over numerical diagnostic ties.\n", sep = "")
