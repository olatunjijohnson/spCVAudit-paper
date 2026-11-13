# Reproducing the current paper

Run commands from the repository root. Linux is the tested platform; parallel
computation uses R's `mclapply`. Use one worker on platforms without fork support.
No script downloads data or changes the installed user-wide `spLeakage` package.

## Environment

The simulations used R 4.5.0, CAST 1.1.1 and ranger 0.16.0. The current Jura and
rendering work used R 4.6.1 with those same CAST and ranger versions. Individual
RDS bundles retain their original `sessionInfo()`. `environment/runtime-lock.csv`
records the complete historical R 4.5.0 dependency closure; source hashes and
non-CRAN provenance accompany it. These records predate the Jura replacement.

The automated retained-result check uses Ubuntu 24.04 and R 4.5.0. This also
avoids compiling the archived ranger version with newer default C++ toolchains
that are incompatible with that release. For the original dependency closure,
use the source-pinned environment builder below.

The pipeline requires R packages `sf`, `CAST`, `ranger`, `gstat`, `pkgload`,
`jsonlite`, `ggplot2`, `patchwork` and `digest`, including their dependencies.
The runner checks the CAST and ranger versions and installs the vendored package
into `revision/csf/library` before analysis. It never installs the moving GitHub
package head. `local_libraries.R` can use already installed libraries from earlier
R minor versions as a fallback; use the locked environment for a controlled rerun.

With **R 4.5.0** on PATH, rebuild the historical source-pinned library using:

```bash
python3 environment/build_runtime.py environment/runtime-lock.csv build/r45-runtime \
  --source-cache build/r-sources \
  --provenance environment/runtime-source-provenance.csv \
  --source-hashes environment/runtime-source-hashes.csv
export R_LIBS_USER="$PWD/build/r45-runtime/library"
```

The builder downloads versioned source archives, verifies their recorded hashes,
and records installation logs. It requires native build tools and the system
libraries listed in `environment/runtime-source-hashes.csv` (including the spatial
stack used by `sf` and `terra`). The library lock is not an operating-system image.
An existing source cache can be supplied with `--cache-only`. Download and install
times depend on the platform; failed builds leave logs and can use `--resume`.

### Numerical environment for exact replay

The original simulation environment was Rocky Linux 9.7, R 4.5.0, with R's
bundled BLAS and LAPACK 3.12.1. Matching R and package versions alone is not
sufficient to reproduce every nonlinear covariance fit. The public-release
audit found signal-fraction differences exceeding 0.01 in three of the 1,440
diagnostic-study configurations on Ubuntu 24.04; the original environment
reproduced their saved fits. A local
R 4.5.0 source build with bundled numerical libraries also differed. Tiny changes
in the simulated inputs and in optimisation arithmetic can affect this fitter.
The new spCVAudit fitter is a separate implementation and is not substituted here.

With `compute --compare-to`, a quick check of configurations 1, 1077 and 1428 now
runs before the studies. It stops if the reference fits differ by more than
1e-10 and saves `revision/audit/reference-environment-probe.csv`. This avoids a
long run in an environment that cannot reproduce the archived fits. It does not
exclude configurations or change the fitting algorithm. Rendering retained
results is independently checked on Ubuntu and remains the recommended first step.

Python 3.12 was used locally. Set up a project environment with:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r environment/requirements.txt
```

These versions record the working Python environment. The pipeline uses NumPy,
pandas, SciPy, Matplotlib and scikit-learn. Matplotlib uses a writable project
cache and a noninteractive backend. PDF compilation additionally requires
`latexmk`, pdfLaTeX, BibTeX and the TeX packages declared in the manuscript.

## Rebuild from retained results

```bash
python3 reproduce.py render
```

This is the recommended first check. It regenerates the paper's summaries and
figures, verifies the numerical identities and compares the generated tables and
132 formatted manuscript values against the saved references. It then reproduces
the supplementary analyses and checks the Jura implementation against `gstat` and
the public package API. The target-response exclusion check recomputes one full
Jura configuration with changed target responses, which must not affect fitting,
CV scores or scheme choices.

Outputs are `revision/manuscript/main.pdf`, `supplement.pdf`, the files in
`revision/figures`, and tables under `revision/results/summaries` and
`revision/results/supplementary`. Logs and machine-readable check records are
under `revision/logs` and `revision/audit`. Saved scientific RDS inputs are not
changed. Generated PDFs may differ in metadata across environments; numerical
comparisons use an absolute tolerance of 1e-10, while formatted TeX tables and
manuscript values must match exactly.

## Recompute every current study

```bash
python3 reproduce.py compute --output ../spLeakage-paper-rerun --workers 8 --compare-to .
```

The output directory must not exist. The runner copies the code and Jura data,
installs the included package, and runs:

| Paper analysis | Computation | Main saved result |
| --- | --- | --- |
| Section 3.2, known covariance | `run_gp_reference.R`, `run_gp_knndm.R` | `gp-reference.rds`, `gp-knndm.rds` |
| Section 3.3, estimated diagnostics | `run_sensitivity.R` | `sensitivity-rng.rds` |
| Section 3.4, covariance refitting | `run_nested_estimation.R` | `nested-estimation-rng.rds` |
| Section 3.5, use of explanatory information | `run_learner_experiment.R` | `learners.rds` |
| Jura application | `run_jura_application.R` | `application-jura.rds` |
| Supplementary fold-count comparison | `run_nested_k10.R` | `nested-estimation-k10.rds` |
| Jura fitting sensitivity | `run_jura_application.R --original-fit` | `application-jura-original.csv` |

It then runs the same rendering and verification as above. `--compare-to .`
compares regenerated primary scientific records against the retained inputs;
the output checks additionally compare derived summaries, supplementary fold
comparisons and Jura resampling with their retained references. All seeds and
setting grids are specified in the computation scripts and configuration CSVs.
The two `-rng` studies reset the response seed separately from fold construction.

The Jura application uses the profile-search refinement implemented in
`jura_profile_estimator.R`. The earlier simulations use the original package
covariance fitter. That distinction is intentional and is discussed in the paper;
updating the simulations to a different fitter would be a new experiment.

## Individual checks and supplementary calculations

After a primary render, these commands may be run separately:

```bash
python3 scripts/fold_count.py
Rscript --vanilla revision/R/audit_selection_tolerance.R
python3 scripts/jura_sensitivity.py
python3 scripts/check_jura_fitting.py
Rscript --vanilla scripts/check_jura_implementation.R
Rscript --vanilla revision/R/application_walkthrough.R
python3 scripts/check_outputs.py
```

The fold-count script uses an isolated temporary directory, so ten-fold summaries
cannot replace the primary five-fold results. Jura resampling keeps the model and
CV scores fixed, reconstructs the original coordinate clusters, and reproduces
2,000 resamples per fixed target and mean model. These are conditional sensitivity
summaries, not independent datasets or confidence intervals for population risk.

The fold-count script also regenerates `revision/exploratory/fold-count-report.pdf`
and `fold-count-comparison.pdf`, both cited in the supplementary material.
The tie audit writes `revision/audit/selection-tolerance.csv`; its 12 rows cover
six numerical tolerances for each covariance choice and are checked against the
retained audit. These outputs are included in the root reproduction command.

The original-fitting Jura CSV is a retained input for rendering only. Full
`compute` regenerates it with the package fitter before comparing it with the
profile-search fit. Independent SciPy optimisation checks the full-model fitting
objective for both mean models at the original target split.

## Scope of verification

For a small computation check using the existing pilot grids, with no retained
simulation results as inputs, run:

```bash
python3 scripts/smoke_simulations.py --output build/simulation-smoke
```

The output directory must be new. This checks all six simulation entry points,
including kNNDM and ten-fold refitting. It does not replace a full experiment.

`provenance/imported-files.json` records the imported research files before
repository integration. Original simulation records were validated in earlier
full reruns described in the supplementary material. Those historical checks
are not a claim that the newly assembled repository has been rerun in full.
The dated assembly checks are in `docs/VERIFICATION.md`; subsequent public-release
checks are reported in `docs/VALIDATION.md`.
Selected historical validation receipts are retained under
`provenance/historical-checks`; their original paths refer to the research archive.

Old analyses are preserved only in Git history and the
`before-public-cleanup-20261112` tag; they are not dependencies. `revision/EXPERIMENTS.md`
preserves the historical study protocol, including superseded plans; the current
paper and the table above define the active analyses.
