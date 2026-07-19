# Reproduction code — *Auditing spatial information leakage in ecological prediction*

Everything needed to reproduce every figure and table in the manuscript
(Johnson, submitted to *Methods in Ecology and Evolution*), which describes the
[**spLeakage**](https://github.com/olatunjijohnson/spLeakage) R package.

This repository is the **paper**; the package lives in its own repository. Nothing here
is required to *use* `spLeakage` — install that from its own repo and read its vignettes.
What is here is the code that produced the published numbers.

---

## Quick start

```r
# 1. the package under test
remotes::install_github("olatunjijohnson/spLeakage")

# 2. dependencies used only by the paper's analyses
install.packages(c("disdat", "malariaAtlas", "ranger", "mgcv", "fields",
                   "xgboost", "sf", "geoR", "sp"))

# 3. reproduce everything from cached results (~3 minutes)
source("run_all.R")

# ...or rebuild the expensive results from scratch (~4 hours, 11 cores)
source("run_all.R")   # after setting REBUILD <- TRUE at the top
```

By default `run_all.R` **reuses `cache/`** and rebuilds only the figures. That is the
fast path and it reproduces every published number exactly. Set `REBUILD <- TRUE` to
re-run the simulations and case studies end to end; expect ~4 hours, and expect Monte
Carlo noise in the last reported digit of the simulation-based results.

---

## What produces what

Every float in the manuscript, and the script that makes it:

| Manuscript float | Script | Cached result | Cost to rebuild |
|---|---|---|---|
| **Fig. 1** — the leakage concept | `scripts/make_figure_concept.R` | — (self-contained) | seconds |
| **Fig. 2**, **Table 3** — theory vs exact optimism | `scripts/validate_theory.R` → `scripts/make_figure_theory.R` | `theory_validation.rds` | ~2 min |
| **Table 1** — the realisation generator | `scripts/simulate_optimism.R` (describes it) | `optimism_sim_table.rds` | ~3 h |
| **Table 2** — design × target recommendation | *static* — no computation; follows the cited literature | — | — |
| **Table 4** — signed SLI vs NNDM's unsigned *W* | `scripts/benchmark_nndm.R` | `benchmark_nndm.rds` | ~3 min |
| **Fig. 3**, **Fig. 4** — emulator calibration; design effect | `scripts/simulate_optimism.R` → `scripts/make_figures.R` | `optimism_sim_table.rds` | ~3 h |
| **Fig. 5** — the SDM case study | `scripts/case_study_sdm.R` → `scripts/make_figure_sdm.R` | `case_study_sdm.rds` | ~5 min |
| **Table 5**, **Fig. 6** — the malaria case study | `scripts/case_study_kenya_st.R` → `scripts/make_figure_malaria_st.R` | `case_study_kenya_st.rds`, `kenya_pr.rds` | ~10 min |
| **Table 6** — the two case studies compared | both case-study scripts | both `.rds` | — |
| **Listing 1** — the paper's worked code, verbatim | `scripts/paper_listing.R` | (uses `disdat`) | seconds |
| §3.4 — real fields vs independent observations | `scripts/validate_independent.R` | `independent_validation.rds` | ~15 min |

Figure scripts **read the cached `.rds`** rather than recomputing an audit, so a figure
and the number quoted beside it in the text cannot drift apart.

`make_figures.R` also emits `figures/F2_sli_vs_optimism.png`, which is **not** used in
the manuscript — it is a diagnostic of the index-versus-optimism relationship kept for
inspection. Every other figure in `figures/` appears in the paper.

---

## Data

All of it is public; none is redistributed here except a small cache.

| Data | Source | Note |
|---|---|---|
| Australian Wet Tropics birds (340 sites, 20 species, 13 predictors) | `disdat` — the Elith et al. (2006, 2020) NCEAS benchmark | fetched from the package at run time |
| Kenya *P. falciparum* prevalence (1,607 records, 1,230 sites, 1985–2016) | `malariaAtlas::getPR()` | cached in `cache/kenya_pr.rds` so the analysis runs offline; delete it to re-fetch |
| meuse, Paraná, ca20 | `sp`, `geoR` | fetched at run time |

**On country choice.** The malaria case uses Kenya, not Nigeria, and this is worth
stating because it looks arbitrary and is not. MAP returns 2,672 Nigeria records but
2,476 are DHS surveys whose cluster coordinates MAP cannot share, leaving **66** with
both coordinates and prevalence — a hard ceiling without DHS credentials. Kenya has
1,607. The reason this matters scientifically is that n = 66 cannot distinguish "the
spatial signal is missing because the data are thin" from "the spatial signal is missing
because pooling 31 years of surveys destroys it". At n = 1,607 the epoch ladder settles
it (see `case_study_kenya_st.R`): it is the pooling.

---

## Reproducibility notes

- **Seeds** are fixed in every script. The case studies additionally report medians over
  20 independent random splits, because the split being audited *is* random: a single
  draw of one bird species moved between +32% and +38%.
- **`sessionInfo()`** from the run that produced the manuscript is in
  `sessionInfo.txt`.
- **The simulation is the only slow step.** If you only want to check the case studies
  or the theory, run those scripts directly; they do not depend on it.
- **Exact vs realisation designs.** Two simulation designs appear in the paper and they
  do different jobs. `validate_theory.R` and `benchmark_nndm.R` use the *exact* design:
  optimism is computed straight from the covariance algebra, with no field realised and
  no learner involved, so it carries no Monte Carlo noise. `simulate_optimism.R` uses
  the *realisation* design: it simulates fields, samples them, fits a learner, and
  measures the optimism actually realised. The first says what optimism *is*; the second
  says what a practitioner can *recover* of it.

## Licence

MIT, matching the package.
