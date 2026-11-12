# Reproduce every figure and table in the manuscript.
# ---------------------------------------------------------------------------
#   source("run_all.R")
#
# REBUILD = FALSE (default): reuse cache/, rebuild only the figures. ~3 minutes, and
#   reproduces every published number exactly.
# REBUILD = TRUE: re-run the simulations and case studies end to end. ~4 hours on 11
#   cores. Expect Monte Carlo noise in the last reported digit of simulation-based
#   results; the exact-algebra results (Table 3, Table 4) are noise-free and will match.
#
# Run from the repository root.

REBUILD <- FALSE

stopifnot(dir.exists("scripts"), dir.exists("cache"))
dir.create("figures", showWarnings = FALSE)

if (!requireNamespace("spLeakage", quietly = TRUE)) {
  stop("Install the package under test first:\n",
       '  remotes::install_github("olatunjijohnson/spLeakage")', call. = FALSE)
}

needed <- c("disdat", "malariaAtlas", "ranger", "mgcv", "fields", "xgboost", "sf")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages used by the paper's analyses: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

run <- function(f) {
  message("\n=== ", f, " ===")
  t0 <- Sys.time()
  source(file.path("scripts", f), echo = FALSE)
  message(sprintf("--- %s done in %.1f min", f,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---- expensive: only with REBUILD ------------------------------------------
if (REBUILD) {
  run("validate_theory.R")        # -> cache/theory_validation.rds   (Fig 2, Table 3)
  run("benchmark_nndm.R")         # -> cache/benchmark_nndm.rds      (Table 4)
  run("simulate_optimism.R")      # -> cache/optimism_sim_table.rds  (Table 1, Fig 3-4)
  run("validate_independent.R")   # -> cache/independent_validation.rds  (S3.4)
  run("case_study_sdm.R")         # -> cache/case_study_sdm.rds      (Fig 5, Table 6)
  run("case_study_kenya_st.R")    # -> cache/case_study_kenya_st.rds (Table 5, Fig 6, Table 6)
} else {
  have <- file.exists(file.path("cache", c(
    "theory_validation.rds", "benchmark_nndm.rds", "optimism_sim_table.rds",
    "independent_validation.rds", "case_study_sdm.rds", "case_study_kenya_st.rds")))
  if (!all(have)) stop("cache/ is incomplete; set REBUILD <- TRUE", call. = FALSE)
  message("Using cached results (set REBUILD <- TRUE to recompute).")
}

## ---- figures: always rebuilt, always FROM the cache -------------------------
# Every figure script reads the saved .rds rather than recomputing its audit, so a
# figure and the number quoted beside it in the text cannot drift apart.
run("make_figure_concept.R")      # Fig 1
run("make_figure_theory.R")       # Fig 2
run("make_figures.R")             # Fig 3, Fig 4
run("make_figure_sdm.R")          # Fig 5
run("make_figure_malaria_st.R")   # Fig 6

## ---- report the headline numbers, so a reader can check them against the paper ----
message("\n================ published numbers, recomputed from cache ================")
t <- readRDS("cache/theory_validation.rds")
message(sprintf("Table 3  closed form vs exact optimism   r = %.3f   (paper: 0.975)",
                cor(t$nn_form, t$exact)))
message(sprintf("Table 3  SLI_rho vs exact optimism       r = %.3f   (paper: 0.976)",
                cor(t$sli_rho, t$exact)))
message(sprintf("Table 3  Wasserstein bound holds         %.0f%%      (paper: 100%%)",
                100 * mean(abs(t$nn_form) <= t$bound + 1e-9)))
b <- readRDS("cache/benchmark_nndm.rds")
message(sprintf("Table 4  signed SLI vs signed optimism   r = %+.3f  (paper: +0.977)",
                cor(b$sli_rho, b$exact)))
message(sprintf("Table 4  unsigned W vs signed optimism   r = %+.3f  (paper: -0.090)",
                cor(b$W, b$exact)))
s <- readRDS("cache/case_study_sdm.rds")
message(sprintf("Fig 5e   species optimistic              %d/%d, median %+.0f%%  (paper: 14/15, +27%%)",
                sum(s$assemblage$optimism > 0), nrow(s$assemblage),
                100 * median(s$assemblage$optimism)))
k <- readRDS("cache/case_study_kenya_st.rds")
message(sprintf("Fig 6c   joint space-time SLI            %+.3f    (paper: +0.502)",
                median(k$draws$sli_st)))
message(sprintf("Fig 6c   SLI_rho == co-located fraction  %+.3f == %+.3f  (paper: identical)",
                median(k$draws$sli_rho), median(k$draws$group_frac)))
message("=========================================================================")
