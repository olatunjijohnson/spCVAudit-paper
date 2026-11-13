# Nearest-neighbour distances and prediction error in spatial cross-validation

Public research code and manuscript for Olatunji Johnson's paper. This research
draft has not been submitted or independently reviewed.

- [Paper](revision/manuscript/main.pdf) and [supplementary material](revision/manuscript/supplement.pdf).
- [Reproduction instructions](REPRODUCING.md), including dependencies and full recomputation.
- [spCVAudit R package](https://github.com/olatunjijohnson/spCVAudit) and
  [website with tutorials](https://olatunjijohnson.github.io/spCVAudit/).
- [Relationship between the packages and the paper](docs/PACKAGE-COVERAGE.md).
- [Jura data and attribution](revision/data/JURA-SOURCE.txt).
- [Verification of the public repository](docs/VALIDATION.md).

The repository contains the analysis code, seeds, retained results, figure scripts,
LaTeX source and the fixed `spLeakage` source used by the reproduction scripts.
The package source is included because development versions can change.
The new **spCVAudit** package provides a focused public interface for the paper's
diagnostics and Gaussian calculations. The fixed source under `package/` remains
**spLeakage 0.1.0.9004**. The original simulation runs recorded version 0.1.0.9003;
that [source archive](provenance/original-package/) is also preserved, and its
variogram fitting function is identical. Use the pinned reproduction environment
for the paper; use spCVAudit for a new analysis.

After installing the dependencies described in [REPRODUCING.md](REPRODUCING.md), run:

```bash
python3 reproduce.py render
```

This rebuilds all main and supplementary figures, generated manuscript numbers,
tables and PDFs from retained results. It also verifies the theory, package/API
agreement, Jura predictions and target-response exclusion, target-cluster
resampling, covariance-fitting sensitivity and the paired fold-count comparison.
Regenerated numerical outputs are compared with the retained paper results.

To rerun all current studies and supplementary analyses in a new directory:

```bash
python3 reproduce.py compute --output ../spLeakage-paper-rerun --workers 8 --compare-to .
```

Full computation is substantially slower than rendering; the simulation studies
were originally run on the University of Manchester's CSF. The `compute` command
requires a new output directory and does not reuse saved simulation checkpoints.
Use `--no-pdf` if LaTeX is unavailable; figures and numerical checks still run.

## Contents

| Location | Contents |
| --- | --- |
| `revision/R/` | Current theory, simulation, application and figure scripts |
| `revision/exploratory/` | Ten-fold sensitivity study and paired comparison |
| `revision/data/` | Jura measurements, original split and source attribution |
| `revision/results/` | Scientific result bundles, configurations and failure records |
| `revision/manuscript/` | Current paper and supplementary material |
| `revision/figures/` | Figures 1–4 and supplementary Figures S1–S5 |
| `package/` | Fixed package source, version 0.1.0.9004 |
| `environment/` | R version lock, source hashes, builder and Python dependencies |
| `scripts/` | Supplementary reproduction and verification |
| `provenance/` | Import checksums and independent retained output references |

Current analyses use Jura. Superseded applications have been removed from the
current checkout; the `before-public-cleanup-20261112` tag preserves the earlier
repository state for historical reference. They are not inputs to reproduction.

The imported-file manifest records the original workspace bytes. Subsequent
integration edits are identified by Git history. RDS metadata also records the
original executed source checksums, package version, seeds and session information.

Project code follows the MIT licence used by `spLeakage`; third-party dependencies
and the Jura data retain their upstream terms and attribution. See
[DATA-SOURCES.md](DATA-SOURCES.md). This repository is public. Cite a release or
commit when reproducing the draft; a publication DOI has not yet been assigned.
