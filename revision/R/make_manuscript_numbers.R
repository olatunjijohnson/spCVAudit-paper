# Every generated number carries a row-level source in the companion CSV.
dir.create("revision/paper-outputs", showWarnings = FALSE)
records <- list()
add <- function(name, value, digits, source, selector, percent = FALSE) {
  stopifnot(length(value) == 1L, is.finite(value))
  formatted <- formatC(value * if (percent) 100 else 1, format = "f", digits = digits)
  records[[name]] <<- data.frame(macro = name, value, formatted, source, selector, percent)
}
read <- function(stem) read.csv(file.path("revision/results/summaries", paste0(stem, ".csv")))
sig <- read("reference-sign")
for (d in c("A", "J1", "J2")) {
  row <- subset(sig, stratum == "all" & threshold == .01 & diagnostic == d)
  label <- c(A = "A", J1 = "JOne", J2 = "JTwo")[[d]]
  add(paste0("Sign", label), row$accuracy, 2, "reference-sign.csv", paste(d, "all, threshold .01"), TRUE)
  add(paste0("SignSE", label), row$mcse, 2, "reference-sign.csv", paste(d, "all, threshold .01 MCSE"), TRUE)
}
sel <- read("reference-selection")
allsel <- sel[sel$stratum == paste0("criterion=", sel$criterion), ]
for (criterion in c("A", "W", "J1", "J2", "U1", "U2", "fixed_random", "fixed_knndm", "oracle")) {
  row <- allsel[allsel$criterion == criterion, ]
  name <- c(A = "A", W = "W", J1 = "JOne", J2 = "JTwo", U1 = "UOne", U2 = "UTwo", fixed_random = "Random", fixed_knndm = "Knn", oracle = "Oracle")[[criterion]]
  add(paste0("Bias", name), row$abs_bias, 5, "reference-selection.csv", row$stratum)
}
paired <- read("reference-paired-comparisons")
r <- subset(paired, criterion == "J2" & baseline == "W" & stratum == "all")
add("BiasDelta", r$difference, 5, "reference-paired-comparisons.csv", "J2 vs W, all")
add("BiasDeltaSE", r$mcse, 5, "reference-paired-comparisons.csv", "J2 vs W, all MCSE")
add("BiasReduction", -r$difference / allsel$abs_bias[allsel$criterion == "W"], 1,
    "reference-selection.csv + reference-paired-comparisons.csv", "(W - J2)/W, all", TRUE)
sens <- read("sensitivity-headlines")
fits <- read("covariance-fits")
for (a in c(1, .25)) {
  d <- fits[fits$axis_ratio == a, ]
  for (boundary in c("Low", "High")) {
    count <- if (boundary == "Low") sum(d$fitted_signal < 1e-7) else sum(d$fitted_signal > 1 - 1e-7)
    add(paste0("Boundary", if (a == 1) "Iso" else "Aniso", boundary), count, 0,
        "covariance-fits.csv", paste("axis ratio", a, if (boundary == "Low") "signal < 1e-7" else "signal > 1 - 1e-7"))
  }
}
small <- read("sensitivity-near-zero")
for (a in c(1, .25)) {
  r <- subset(small, axis_ratio == a & diagnostic_model == "estimated_isotropic" & diagnostic == "J2" & threshold == .005)
  add(paste0("SmallGap", if (a == 1) "Iso" else "Aniso"), r$large_gap_given_small, 1,
      "sensitivity-near-zero.csv", paste(a, "estimated isotropic J2 threshold .005, large |O| > .05"), TRUE)
  add(paste0("SmallCount", if (a == 1) "Iso" else "Aniso"), r$n_small, 0,
      "sensitivity-near-zero.csv", paste(a, "estimated isotropic J2 threshold .005, n_small"))
}
for (a in c(.25, 1)) for (m in c("estimated_isotropic", "known_isotropic")) {
  r <- sens[sens$axis_ratio == a & sens$diagnostic_model == m, ]
  suffix <- paste0(if (a == 1) "Iso" else "Aniso", if (m == "estimated_isotropic") "Est" else "Nominal")
  add(paste0("Sens", suffix), r$sign_accuracy, 1, "sensitivity-headlines.csv", paste(a, m), TRUE)
  if (m == "known_isotropic") {
    add(paste0("Sens", if (a == 1) "Iso" else "Aniso", "True"), r$true_J2_sign_accuracy, 1, "sensitivity-headlines.csv", paste(a, "true J2"), TRUE)
    add(paste0("Cert", if (a == 1) "Iso" else "Aniso"), r$certificate_coverage, 1, "sensitivity-headlines.csv", paste(a, "certificate coverage"), TRUE)
    add(paste0("Size", if (a == 1) "Iso" else "Aniso"), r$mean_abs_training_effect, 3, "sensitivity-headlines.csv", paste(a, "mean absolute training-size effect"))
  }
}
counter <- read.csv("revision/results/loo-neighbourhood-counterexample.csv")
add("NeighbourCounterexample", counter$exact, 5, "loo-neighbourhood-counterexample.csv", "exact")
learn <- read("learners-by-cell")
for (learner in c("forest_X", "idw_coordinates", "oracle_kriging", "linear_X")) for (a in c(0, .9)) {
  r <- subset(learn, range == .2 & design == "clustered" & target == "domain" & scheme == "random")
  r <- r[r$learner == learner & r$alpha == a, ]
  name <- paste0(c(forest_X = "Forest", idw_coordinates = "IDW", oracle_kriging = "Kriging", linear_X = "Linear")[[learner]], if (a == 0) "Zero" else "Nine")
  add(name, r$mean_gap, 3, "learners-by-cell.csv", paste(learner, a, "range .2, clustered, domain, random"))
  add(paste0(name, "SE"), r$mcse_gap, 3, "learners-by-cell.csv", paste(learner, a, "range .2, clustered, domain, random MCSE"))
}
r <- subset(learn, range == .2 & design == "clustered" & target == "domain" & scheme == "random" & learner == "forest_X" & alpha == 0)
add("ConstantMarginal", r$mean_J2_marginal, 3, "learners-by-cell.csv", "range .2, clustered, domain, random marginal J2")
app <- read("application-jura-candidates")
for (target in c("published", "uniform", "north")) for (scheme in c("random", "block", "buffered", "knndm")) {
  r <- app[app$target == target & app$model == "ordinary" & app$scheme == scheme, ]
  stem <- paste0("Jura", c(published = "Published", uniform = "Uniform", north = "North")[[target]],
                 c(random = "Random", block = "Block", buffered = "Buffer", knndm = "Knn")[[scheme]])
  for (v in c("cv_mse_mean", "J2_mean", "W_mean", "absolute_error_mean")) {
    name <- paste0(stem, c(cv_mse_mean = "Score", J2_mean = "JTwo", W_mean = "W", absolute_error_mean = "MAE")[[v]])
    add(name, r[[v]], if (v == "W_mean") 4 else 3, "application-jura-candidates.csv", paste(target, "ordinary", scheme, v))
  }
  if (scheme == "random") add(paste0("Jura", c(published="Published",uniform="Uniform",north="North")[[target]], "Target"),
    r$target_mse_mean, 3, "application-jura-candidates.csv", paste(target, "ordinary target_mse_mean"))
}
nested <- read("nested-selection")
np <- read("nested-paired-comparisons")
for (criterion in c("W", "J2_fitted", "fixed_knndm", "fixed_random")) {
  r <- nested[nested$stratum == paste0("learner=estimated;criterion=", criterion), ]
  name <- c(W = "W", J2_fitted = "JTwo", fixed_knndm = "Knn", fixed_random = "Random")[[criterion]]
  add(paste0("NestedMAE", name), r$mean_absolute_error, 5, "nested-selection.csv", r$stratum)
  add(paste0("NestedBias", name), r$mean_signed_error, 5, "nested-selection.csv", paste(r$stratum, "mean signed error"))
}
r <- subset(np, learner == "estimated" & criterion == "J2_fitted" & baseline == "W" & stratum == "all")
for (v in c("mean_absolute_error_difference", "mcse", "lower", "upper")) {
  name <- c(mean_absolute_error_difference = "Delta", mcse = "SE", lower = "Lower", upper = "Upper")[[v]]
  add(paste0("Nested", name), r[[v]], 5, "nested-paired-comparisons.csv", paste("estimated, J2_fitted vs W, all", v))
}
base <- nested[nested$stratum == "learner=estimated;criterion=W", ]
add("NestedReduction", -r$mean_absolute_error_difference / base$mean_absolute_error, 1,
    "nested-selection.csv + nested-paired-comparisons.csv", "estimated: (W - J2_fitted)/W", TRUE)
r <- subset(np, learner == "estimated" & criterion == "J2_fitted" & baseline == "W" &
  stratum == "design=random;target=shifted")
add("NestedShiftDelta", r$mean_absolute_error_difference, 5, "nested-paired-comparisons.csv",
    "estimated, J2_fitted vs W, design=random;target=shifted")
add("NestedShiftSE", r$mcse, 5, "nested-paired-comparisons.csv",
    "estimated, J2_fitted vs W, design=random;target=shifted MCSE")
cell <- "range=0.05;signal=0.95;design=random;target=shifted"
r <- np[np$learner == "estimated" & np$criterion == "J2_fitted" & np$baseline == "W" & np$stratum == cell, ]
add("NestedStrongCellDelta", r$mean_absolute_error_difference, 5, "nested-paired-comparisons.csv", paste("estimated, J2_fitted vs W,", cell))
influence <- read("nested-cell-influence")
r <- influence[influence$omitted_cell == cell, ]
add("NestedOtherCellsDelta", r$remaining_mean_difference, 5, "nested-cell-influence.csv", paste("omitted cell", cell))
behaviour <- read("nested-behaviour-by-cell")
r <- subset(behaviour, learner == "estimated" & stratum == "all")
for (v in c("same_selected_set", "J2_lower_loss", "W_lower_loss", "single_replicate_omission_min", "single_replicate_omission_max")) {
  name <- c(same_selected_set = "SameChoices", J2_lower_loss = "Wins", W_lower_loss = "Losses",
            single_replicate_omission_min = "OmissionMin", single_replicate_omission_max = "OmissionMax")[[v]]
  add(paste0("Nested", name), r[[v]], if (grepl("omission", v)) 5 else 0,
      "nested-behaviour-by-cell.csv", paste("estimated, all,", v))
}
r <- subset(behaviour, learner == "estimated" & stratum == "0.05:0.95:random:shifted")
add("NestedStrongRandomChoices", r$J2_random_only, 0, "nested-behaviour-by-cell.csv",
    "estimated, scale .05, signal .95, uniform, shifted; J2 selects random alone")
add("NestedStrongWRandomChoices", r$W_random_only, 0, "nested-behaviour-by-cell.csv",
    "estimated, scale .05, signal .95, uniform, shifted; W selects random alone")
scores <- read("nested-candidate-scores")
r <- subset(scores, learner == "estimated" & range == .05 & signal == .95 & design == "random" & target == "shifted")
add("NestedStrongTargetRisk", r$conditional_target_risk[r$scheme == "random"], 3,
    "nested-candidate-scores.csv", "estimated, scale .05, signal .95, uniform, shifted; mean conditional target risk")
for (scheme in c("random", "knndm")) add(paste0("NestedStrongCV", if (scheme == "random") "Random" else "Knn"),
  r$cv_score[r$scheme == scheme], 3, "nested-candidate-scores.csv",
  paste("estimated, scale .05, signal .95, uniform, shifted; mean CV score", scheme))
out <- do.call(rbind, records)
write.csv(out, "revision/results/manuscript-number-provenance.csv", row.names = FALSE)
writeLines(c("% Generated by revision/R/make_manuscript_numbers.R; do not edit.",
             sprintf("\\newcommand{\\%s}{%s}", out$macro, out$formatted)), "revision/paper-outputs/numbers.tex")

# Table 2 compares selection rules and benchmarks within sampling/target groups.
lines <- c("% Generated from the paired reference analysis.",
  "\\begin{tabular}{llrrrrr}", "\\toprule",
  "Sampling & Target & \\shortstack{Selection\\\\by $W$} & \\shortstack{Selection\\\\by $|J_2|$} & \\shortstack{Fixed\\\\kNNDM} & \\shortstack{Best of\\\\four} & \\shortstack{Difference\\\\(MCSE)} \\\\", "\\midrule")
for (design in c("random", "clustered")) for (target in c("domain", "shifted")) {
  suffix <- paste0("design=", design, ";target=", target)
  values <- vapply(c("W", "J2", "fixed_knndm", "oracle"), function(cr) {
    sel$abs_bias[sel$stratum == paste0("criterion=", cr, ";", suffix)]
  }, 0.0)
  r <- paired[paired$stratum == suffix & paired$criterion == "J2" & paired$baseline == "W", ]
  lines <- c(lines, sprintf("%s & %s & %.4f & %.4f & %.4f & %.4f & %.5f (%.5f) \\\\",
    if (design == "random") "Uniform" else "Clustered", if (target == "domain") "In-domain" else "Shifted",
    values[1], values[2], values[3], values[4], r$difference, r$mcse))
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), "revision/paper-outputs/table-reference.tex")
cat("Generated", nrow(out), "traceable manuscript values and one table.\n")
