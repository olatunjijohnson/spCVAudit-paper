# Public repository verification

Checked on 12–13 November 2026 when preparing the public repository.

All current studies were computed afresh from the public source in a new
workspace, without retained simulation inputs or cached batches. The six primary
studies pass all **179** raw-result comparisons; adding the fresh ten-fold study
gives **211** passing comparisons. All **29** columns of the original-fitting
Jura sensitivity also agree. The original numerical environment reproduces the
five- and ten-fold covariance-fit records exactly.

The fresh outputs were transferred with verified hashes and used to rebuild the
paper locally. All **37** derived-output comparisons pass. The rebuilt main and
supplementary PDFs have exactly the same page text as the release PDFs. The
[fresh-replay records](../provenance/fresh-replay-20261113/) describe the separate
computation and rendering stages, including recovery from a missing plotting
dependency on the cluster.

- The complete `python3 reproduce.py render` workflow passed, including main and
  supplementary figures, manuscript tables, theory checks, Jura prediction
  checks, target-response exclusion, fitting sensitivity, and the paired
  fold-count comparison.
- All **37** retained numerical-output and formatted-table comparisons passed.
- The updated main paper and supplementary material compile successfully with
  their public repository links. No unresolved references or overfull lines
  remain in these two documents.
- The new spCVAudit API agrees with **984** Gaussian-reference and retained Jura
  cases. The maximum absolute difference is **1.72e-14**. The comparison uses
  the same supplied covariance, including the separately labelled residual
  covariance sensitivity; it does not replace the paper's fitting procedure.
- The public package website passed **28** desktop/mobile browser checks;
  **656** internal documentation links resolved.

Earlier machine-readable checks are under `provenance/public-release-20261112/`.
The GitHub Actions workflow independently regenerates and verifies the retained
results on Linux. Its live status is available in the repository's Actions tab.

The hosted run [43342785310](https://github.com/olatunjijohnson/spCVAudit-paper/actions/runs/43342785310)
passed all 37 comparisons. The theory check now enforces its stipulated unit
variance exactly, and the target-response exclusion check compares two fresh
runs in the same environment. Fresh-fit differences from saved Jura results are
reported separately. They are required to pass in exact compute mode; retained
result rendering instead requires the fixed-fit and paired invariance checks.

Exact recomputation of nonlinear covariance fits also depends on the numerical
environment. Three regression configurations reproduce exactly under the
original Rocky Linux/R environment and differ under the tested Ubuntu runtime.
The complete Ubuntu attempt failed 19 of the 179 comparisons in the fitted
studies; its outputs were not used to replace the paper results. The
[paired probe records](../provenance/public-release-20261112/) and the
[full comparison record](../provenance/fresh-replay-20261113/ubuntu-primary-comparison.csv)
are preserved.
`compute --compare-to` now checks these cases before running the studies; see
[REPRODUCING.md](../REPRODUCING.md). The retained results and fitting algorithm
are preserved.

The paper is an author research draft. Software verification does not constitute
independent scientific peer review or establish robustness to every covariance
fitting procedure. The limitations stated in the manuscript remain applicable.
