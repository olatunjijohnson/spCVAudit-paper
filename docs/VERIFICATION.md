# Initial repository verification

This is the dated assembly record. See [the current validation record](VALIDATION.md)
for subsequent public-release checks and recomputation.

Checks run locally on 12 November 2026 with R 4.6.1, CAST 1.1.1, ranger 0.16.0
and the Python versions in `environment/requirements.txt`.

## Fresh simulation computations

`python3 scripts/smoke_simulations.py --output build/simulation-smoke` completed
all six existing pilot grids, with no saved simulation inputs and no failed
configurations. It produced 1,664 candidate/model evaluations:

| Study | Evaluations |
| --- | ---: |
| Gaussian reference | 192 |
| Paired kNNDM reference | 64 |
| Learner experiment | 640 |
| Covariance-diagnostic sensitivity | 384 |
| Five-fold covariance refitting | 192 |
| Ten-fold covariance refitting | 192 |

The machine-readable receipt is `revision/audit/simulation-smoke.json`.

## Full render from retained scientific results

`python3 reproduce.py render` completed successfully, including every supplementary
stage. It rebuilt the 21-page main paper, 18-page supplementary material, Figures
1–4 and Figures S1–S5. The final LaTeX logs had no undefined references or citations.

- All 132 formatted manuscript values and both generated tables matched exactly.
- All 36 retained summary, fold-count, resampling and manuscript-output comparisons
  passed, using 1e-10 absolute tolerance for numerical CSV entries.
- The theory checks included 800 additional-neighbour decompositions and bounds,
  plus 200 moment decompositions; the maximum identity error was 1.84e-16.
- All 1,920 paired fields in the fold-count comparison passed the checks for
  identical full-data fits and target risks across five and ten folds.
- Jura verification reconstructed all 960 evaluations and 5,040 covariance-fit
  records from 120 configurations. CV and target MSE reconstruction agreed to
  1.06e-15; the largest distance reconstruction difference was 2.27e-13 km.
- Ordinary and universal kriging predictions agreed with `gstat` within 7.35e-14.
  Changing only target responses left fitting, CV scores and scheme choices unchanged.
- A fresh original-fitter Jura configuration matched its retained result within
  5.00e-16. Independent SciPy optimisation checks passed for both full mean models.
- The public package API reproduced A, W, J1, J2 and U1 for all four primary Jura
  CV schemes, including the externally supplied kNNDM folds.

Receipts and text copies of the execution logs are in
[`provenance/verification-20261112`](../provenance/verification-20261112).
These are checks of the assembled repository; historical validation receipts are
kept separately under `provenance/historical-checks`.

The subsequent submission audit restored the two cited fold-count PDF reports
and connected the tie-tolerance audit to the root workflow. Regeneration of those
outputs passed; all **37** output comparisons now pass with zero numerical
differences. The manuscript clarifies that the anisotropic reference uses the
most correlated neighbour, and that target response values (rather than target
locations) are excluded from fitting and diagnostic selection. The numerical
results and formatted manuscript values are unchanged.

## Scope

These small computations test the current scripts and dependency integration.
They do not constitute a new full-scale simulation rerun. The retained full
results and selected historical validation records are included separately.
The source-pinned R 4.5.0 environment was not rebuilt during repository assembly;
the builder and its historical validation receipts are provided for that purpose.
