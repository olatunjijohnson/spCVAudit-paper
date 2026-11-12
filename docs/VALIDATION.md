# Public repository verification

Checked on 12 November 2026 when preparing the public repository.

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

Machine-readable results are under `provenance/public-release-20261112/`.
The GitHub Actions workflow independently regenerates and verifies the retained
results on Linux. Its live status is available in the repository's Actions tab.

The paper is an author research draft. Software verification does not constitute
independent scientific peer review or establish robustness to every covariance
fitting procedure. The limitations stated in the manuscript remain applicable.
