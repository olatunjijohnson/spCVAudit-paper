source("revision/R/local_libraries.R")
source("revision/R/common.R")
library(ggplot2)
library(patchwork)
dir.create("revision/figures", showWarnings = FALSE)
cols <- c(random = "#0072B2", block = "#E69F00", buffered = "#D55E00", knndm = "#009E73", target = "#202020")
scheme_names <- c(random = "Random", block = "Blocks", buffered = "Buffered", knndm = "kNNDM", target = "Target")
theme_set(theme_bw(base_size = 10.5, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ededed"),
        plot.title = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "#f4f4f4"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(face = "bold", size = 13)))
scheme_scale <- function() scale_colour_manual(values = cols, breaks = names(cols), labels = scheme_names, drop = FALSE)
save_plot <- function(p, stem, width = 7.2, height = 6) {
  ggsave(file.path("revision/figures", paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path("revision/figures", paste0(stem, ".png")), p, width = width, height = height,
         dpi = 220, bg = "white", device = grDevices::png, type = "cairo")
}
ecdf_frame <- function(d, label) {
  knots <- sort(unique(c(0, d)))
  data.frame(distance = knots, probability = ecdf(d)(knots), scheme = label)
}
extend_curves <- function(d, xmax = max(d$distance)) {
  rbind(d, data.frame(distance = xmax, probability = 1, scheme = unique(d$scheme)))
}

# Figure 1 uses an illustrative seed, unrelated to the evaluated configurations.
set.seed(20261029)
xy <- sampling_locations(100, "clustered")
q <- target_locations("domain", 11)
folds <- make_candidates(xy)
dd <- pair_distance(xy)
cv_d <- lapply(folds, function(ff) {
  v <- numeric(nrow(xy))
  for (f in ff) v[f$test] <- apply(dd[f$test, f$train, drop = FALSE], 1, min)
  v
})
qd <- apply(pair_distance(q, xy), 1, min)
map_example <- rbind(data.frame(x = q[, 1], y = q[, 2], type = "Target grid"),
                     data.frame(x = xy[, 1], y = xy[, 2], type = "Observations"))
p1 <- ggplot(map_example, aes(x, y, shape = type, colour = type)) + geom_point(size = 1.8) +
  scale_shape_manual(values = c(Observations = 16, "Target grid" = 3)) +
  scale_colour_manual(values = c(Observations = "#202020", "Target grid" = "#888888")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "Horizontal coordinate", y = "Vertical coordinate", title = "Sample and target")
curves <- do.call(rbind, c(lapply(c("random", "block"), function(s) ecdf_frame(cv_d[[s]], s)), list(ecdf_frame(qd, "target"))))
p2 <- ggplot(extend_curves(curves), aes(distance, probability, colour = scheme, linetype = scheme)) +
  geom_step(linewidth = .8) + scheme_scale() +
  scale_linetype_manual(values = c(random = "solid", block = "longdash", target = "solid"), labels = scheme_names) +
  labs(x = "Nearest training distance", y = "Cumulative probability", title = "Distance distributions") +
  guides(linetype = "none", colour = guide_legend(nrow = 1,
    override.aes = list(linetype = c("solid", "longdash", "solid")))) +
  theme(legend.key.width = grid::unit(1.2, "cm"))
h <- seq(0, .8, length.out = 300)
weights <- rbind(data.frame(h, weight = exp(-h / .15) / .15, curve = "First contrast"),
                 data.frame(h, weight = 2 * exp(-2 * h / .15) / .15, curve = "Squared contrast"))
p3 <- ggplot(weights, aes(h, weight, linetype = curve)) + geom_line(linewidth = .8, colour = "#6B4C9A") +
  scale_linetype_manual(values = c("First contrast" = "longdash", "Squared contrast" = "solid"),
    labels = c("First contrast" = "J1 (dashed)", "Squared contrast" = "J2 (solid)")) +
  theme(legend.key.width = grid::unit(1.2, "cm")) +
  labs(x = "Distance", y = "Weight on CDF difference", title = "Dependence weighting") +
  annotate("text", x = .47, y = 10.7, label = "Exponential correlation\nscale = 0.15, signal = 1", size = 3)
cross <- rbind(ecdf_frame(rep(-log(.5), 5), "random"),
               ecdf_frame(c(rep(-log(.8), 3), rep(-log(.01), 2)), "target"))
p4 <- ggplot(extend_curves(cross, 5), aes(distance, probability, colour = scheme, linetype = scheme)) + geom_step(linewidth = .8) +
  scale_colour_manual(values = c(random = "#777777", target = "#202020"), breaks = c("random", "target"), labels = c("Validation (dashed)", "Target (solid)")) +
  scale_linetype_manual(values = c(random = "longdash", target = "solid"), breaks = c("random", "target"), labels = c("Validation (dashed)", "Target (solid)")) +
  theme(legend.key.width = grid::unit(1.2, "cm"), legend.text = element_text(size = 9)) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 1)) +
  annotate("label", x = 2.8, y = .29, label = "J1 = +0.016\nJ2 = -0.134", size = 3.2, fill = "white", linewidth = 0) +
  labs(x = "Distance (unit exponential scale)", y = "Cumulative probability", title = "A crossing counterexample")
save_plot(((p1 + p2) + plot_layout(widths = c(1, 1))) / (p3 + p4) + plot_annotation(tag_levels = "A"), "F1-diagnostic", 7.2, 6.4)

# Figure 2: main factor grid and separately generated sensitivity extension.
sig <- read.csv("revision/results/summaries/reference-sign.csv")
sig <- subset(sig, threshold == .01 & stratum %in% c("all", "scheme=random", "scheme=block", "scheme=buffered", "scheme=knndm"))
sig$group <- factor(sig$stratum, levels = rev(c("all", "scheme=random", "scheme=block", "scheme=buffered", "scheme=knndm")),
                    labels = rev(c("All schemes", "Random", "Blocks", "Buffered", "kNNDM")))
sig$lower <- pmax(0, sig$accuracy - 1.96 * sig$mcse)
sig$upper <- pmin(1, sig$accuracy + 1.96 * sig$mcse)
diagcols <- c(A = "#888888", J1 = "#6B4C9A", J2 = "#222222")
dodge <- position_dodge(width = .55)
p1 <- ggplot(sig, aes(100 * accuracy, group, colour = diagnostic, shape = diagnostic)) +
  geom_errorbar(aes(xmin = 100 * lower, xmax = 100 * upper), width = .12, position = dodge, orientation = "y") +
  geom_point(position = dodge, size = 2) + scale_colour_manual(values = diagcols) +
  scale_x_continuous(limits = c(75, 101), breaks = c(80, 90, 100)) +
  labs(x = "Sign agreement (%)", y = NULL, title = "Sign agreement")
sel <- read.csv("revision/results/summaries/reference-selection.csv")
sel <- sel[sel$stratum == paste0("criterion=", sel$criterion) & sel$criterion != "fixed_random", ]
sel$choice <- factor(sel$criterion, levels = rev(c("A", "W", "J1", "J2", "U1", "U2", "fixed_knndm", "oracle")),
  labels = rev(c("|A|", "W", "|J1|", "|J2|", "U1", "U2", "Fixed kNNDM", "Best of four")))
p2 <- ggplot(sel, aes(abs_bias, choice)) +
  geom_errorbar(aes(xmin = abs_bias - 1.96 * abs_bias_mcse, xmax = abs_bias + 1.96 * abs_bias_mcse), width = .15, orientation = "y") +
  geom_point(size = 2) + scale_x_continuous(limits = c(.032, .046), breaks = c(.035, .04, .045)) +
  labs(x = "Mean absolute expected bias", y = "Selection rule or benchmark", title = "Choosing a CV scheme")
sens <- read.csv("revision/results/summaries/sensitivity-headlines.csv")
sens$setting <- ifelse(sens$axis_ratio == 1, "Isotropic field", "Anisotropic field")
sens$weight <- ifelse(sens$diagnostic_model == "known_isotropic", "Nominal isotropic", "Estimated isotropic")
true <- sens[sens$diagnostic_model == "known_isotropic", ]
true$weight <- "True covariance metric"; true$sign_accuracy <- true$true_J2_sign_accuracy
true$sign_mcse <- true$true_J2_sign_mcse
sens <- rbind(sens, true)
sens <- sens[!(sens$axis_ratio == 1 & sens$weight == "Nominal isotropic"), ]
sens$case <- paste(ifelse(sens$axis_ratio == 1, "Isotropic", "Anisotropic"),
  c("True covariance metric" = "true", "Nominal isotropic" = "nominal", "Estimated isotropic" = "fitted")[sens$weight], sep = ": ")
sens$case <- factor(sens$case, rev(c("Isotropic: true", "Isotropic: fitted", "Anisotropic: true", "Anisotropic: nominal", "Anisotropic: fitted")))
p3 <- ggplot(sens, aes(100 * sign_accuracy, case)) + geom_point(size = 2.3) +
  geom_errorbar(aes(xmin = 100 * pmax(0, sign_accuracy - 1.96 * sign_mcse),
                    xmax = 100 * pmin(1, sign_accuracy + 1.96 * sign_mcse)), orientation = "y", width = .15) +
  scale_x_continuous(limits = c(85, 102), breaks = c(85, 90, 95, 100)) +
  labs(x = "Sign agreement for J2 (%)", y = NULL, title = "Estimation and anisotropy")
np <- read.csv("revision/results/summaries/nested-paired-comparisons.csv")
np <- subset(np, learner == "estimated" & baseline == "W" & stratum == "all" &
  criterion %in% c("J1_fitted", "J2_fitted", "U1_fitted", "U2_fitted", "J2_known", "fixed_knndm"))
np$criterion <- factor(np$criterion, rev(c("J1_fitted", "J2_fitted", "U1_fitted", "U2_fitted", "J2_known", "fixed_knndm")),
  rev(c("|J1| (fitted)", "|J2| (fitted)", "U1 (fitted)", "U2 (fitted)", "|J2| (known)", "Fixed kNNDM")))
p4 <- ggplot(np, aes(mean_absolute_error_difference, criterion)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "#888888") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = .15, orientation = "y") + geom_point(size = 2) +
  labs(x = "MAE difference from W", y = NULL, title = "CV estimation error")
save_plot((p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A"), "F2-reference", 7.2, 6.4)

# Figure 3: displayed slice is explicit; all factor cells are in S2.
learn <- read.csv("revision/results/summaries/learners-by-cell.csv")
learner_names <- c(known_mean = "Known mean from X", linear_X = "Linear regression on X", forest_X = "Random forest on X",
                  idw_coordinates = "Inverse-distance weighting", oracle_kriging = "Known residual kriging")
learn$learner_label <- factor(learn$learner, names(learner_names), learner_names)
learner_plot <- function(d, ncol = 2) {
  ggplot(d, aes(alpha, mean_gap, colour = scheme, group = scheme)) +
    geom_hline(yintercept = 0, colour = "#aaaaaa", linetype = "dotted") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .025, linewidth = .4,
                  position = position_dodge(width = .035)) +
    geom_line(linewidth = .6) + geom_point(size = 1.7) + scheme_scale() +
    facet_wrap(~learner_label, ncol = ncol, drop = TRUE) +
    scale_x_continuous(breaks = c(0, .3, .6, .9)) +
    labs(x = expression("Variance explained by the predictor (" * alpha * ")"),
         y = "Mean target loss minus CV loss")
}
d <- subset(learn, range == .2 & design == "clustered" & target == "domain" & learner != "known_mean" & scheme %in% c("random", "knndm"))
p <- learner_plot(d) + guides(colour = guide_legend(override.aes = list(linewidth = 1)))
save_plot(p, "F3-learners", 7.2, 5.2)

# Figure 4: descriptive variation within a single survey, not repeated populations.
status <- system2("python3", "revision/R/make_jura_figures.py --main")
stopifnot(status == 0L)

# Supplementary pages preserve the full learner factorial grid.
pdf("revision/figures/S2-learners-all.pdf", width = 7.2, height = 8, family = "sans", useDingbats = FALSE)
for (r in c(.05, .2)) for (design in c("random", "clustered")) for (target in c("domain", "shifted")) {
  d <- learn[learn$range == r & learn$design == design & learn$target == target, ]
  design_label <- if (design == "random") "Uniform" else "Clustered"
  target_label <- if (target == "domain") "in-domain" else "shifted"
  print(learner_plot(d) + labs(title = paste0("Scale ", r, " | ", design_label, " sampling | ", target_label, " target")))
}
dev.off()
cat("Generated four main figures and full-factor learner supplement.\n")
