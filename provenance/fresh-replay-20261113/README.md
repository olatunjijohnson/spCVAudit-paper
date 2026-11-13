# Fresh public-repository replay, 13 November 2026

All six primary studies, the ten-fold sensitivity study and the original-fitting
Jura sensitivity were computed afresh. The computation used public repository
commit `19e3e505f1574eab31ea5f367129f3352c854ce3` in a new directory, without
copying retained results or cached simulation batches. All 54 recorded scientific
source and data files agree with the release source. Validation and local
rendering use commit `6536e32e44a7106b6fa4ad827ec5f8634cb622e7`; intervening changes
concern validation and documentation, not the scientific fitting procedures.

## Computation and comparisons

- The six primary studies pass all **179** raw-result/fit-record comparisons.
- Including the fresh ten-fold study gives **211** passing comparisons. Its
  numeric records agree exactly, including all 78,720 covariance fits.
- All **29** retained columns of the original-fitting Jura sensitivity agree
  across 960 evaluations. The largest difference is 7.23e-13.
- The absolute comparison tolerance is 1e-10. Timing and warning text are excluded
  from scientific numeric comparisons. No configurations were removed.
- The transfer record hashes all 37 collected files and confirms that local
  copies match their freshly computed remote counterparts. It records the files
  before local rendering and comparison overwrite generated outputs.

The reference environment is Rocky Linux 9.7, R 4.5.0, with bundled BLAS and
LAPACK 3.12.1. Compiler and runtime details are in
[`reference-numerical-runtime.txt`](../../environment/reference-numerical-runtime.txt).
The original cluster archive was used only for installed dependencies; source
and outputs were isolated in a new workspace.

The primary computation completed successfully. Its first rendering attempt
stopped at figure generation because NumPy was absent on the cluster. The
supplementary computations then ran in a separate successful batch job. Fresh
result files were transferred, hash-verified and used for local rendering with
the documented Python and TeX dependencies. This is a staged recovery, rather
than a claim that the initial command exited successfully.

An earlier full computation on Ubuntu completed the six studies but failed 19
of the 179 strict comparisons, all in the studies using the original nonlinear
covariance fitter. `ubuntu-primary-comparison.csv` preserves that result.
Matching R and package versions alone did not resolve those differences; the
original numerical environment did. The preflight check in `compute --compare-to`
now detects this incompatibility before long computations start.

The completed local rebuild passes all **37** derived-output comparisons,
including exact formatted manuscript numbers and tables. Fresh Jura replay,
independent kriging prediction and target-response exclusion checks pass. The
rebuilt main and supplementary PDFs have exactly the same page text as the
release PDFs, with no unresolved references or overfull lines.

## Records

- `verification.json`: completion record and scope of the fresh replay.
- `fresh-derived-comparisons.json`: all 37 checks using freshly computed inputs.
- `fresh-jura-fit-replay.csv` and `fresh-jura-implementation.csv`: fresh Jura checks.
- `pdf-text-comparison.json`: the rebuilt and release PDF text comparisons.
- `primary-comparison.csv`: all six primary raw-result comparisons.
- `primary-and-ten-fold-comparison.csv`: the seven-study comparison.
- `original-jura-comparison.json`: complete original-fitting Jura CSV comparison.
- `fit-record-comparison.csv`: per-field checks of all retained covariance-fit
  records; five- and ten-fold fits agree exactly and the largest Jura fit-field
  difference is 2.62e-11.
- `scientific-source-parity.json`: unchanged scientific source and data hashes.
- `transfer.json`: file hashes from the fresh reference run.
- `github-retained-render.json`: the separate successful hosted render check,
  including all 37 numerical-output and formatted-text comparisons.
- `github-jura-fit-replay.csv`: fresh Jura checks from that hosted environment.
- `ubuntu-primary-comparison.csv`: the unsuccessful alternate-environment replay.

The scientific result bundles committed to the repository remain the retained
paper inputs. Fresh computation records validate them; they do not replace them
or change the manuscript's reported results.
