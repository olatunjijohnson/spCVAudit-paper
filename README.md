# Nearest-neighbour distances and prediction error in spatial cross-validation

Research code and manuscript for Olatunji Johnson's paper, current draft dated
12 November 2026. This draft has not been submitted or independently reviewed.

- [Paper](revision/manuscript/main.pdf) and [supplementary material](revision/manuscript/supplement.pdf).
- [Reproduction instructions](REPRODUCING.md), including dependencies and full recomputation.
- [Which paper methods are in the spLeakage package?](docs/PACKAGE-COVERAGE.md).
- [Submission audit and outstanding items](docs/submission-audit.pdf).
- [Jura data and attribution](revision/data/JURA-SOURCE.txt).

The repository contains the analysis code, seeds, retained results, figure scripts,
LaTeX source and the exact `spLeakage` working-tree snapshot used for this paper.
The package source is included because the development version installed from
[the package repository](https://github.com/olatunjijohnson/spLeakage) can change.
The package and paper repository serve different purposes: several theoretical
calculations and study-specific selection rules are implemented in the research code.

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
| `legacy/2026-submission/` | Earlier repository contents, preserved for history |

The legacy directory is excluded from every current reproduction command. Its
earlier application, title and submission notes do not describe the current paper.
Current analyses use Jura; the earlier bird analysis and the exploratory PM2.5
comparison are not inputs to the current paper.

The imported-file manifest records the original workspace bytes. Subsequent
integration edits are identified by Git history. RDS metadata also records the
original executed source checksums, package version, seeds and session information.

Project code follows the MIT licence used by `spLeakage`; third-party dependencies
and the Jura data retain their upstream terms and attribution. See
[DATA-SOURCES.md](DATA-SOURCES.md). The repository is private during preparation;
a permanent public archive and DOI have not yet been assigned.
