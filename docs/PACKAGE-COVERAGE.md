# Relationship between the paper, spCVAudit and spLeakage

## spCVAudit: the public companion package

[spCVAudit](https://github.com/olatunjijohnson/spCVAudit) provides a focused API
for new analyses. [Its website](https://olatunjijohnson.github.io/spCVAudit/)
contains the function reference and four worked tutorials.

| Paper method | spCVAudit interface |
| --- | --- |
| Target locations and arbitrary nonnegative target weights | `prediction_target()` |
| Explicit fold training/test sets; random, block, buffered and optional CAST kNNDM constructors | `cv_folds()`, `spatial_folds()` |
| A, W, J1, J2, U1, U2 and distribution crossings | `distance_contrasts()`, `audit_cv()` |
| Selection by absolute J2 or another common diagnostic; equal probabilities for tied candidates | `compare_cv_schemes()` |
| Specified exponential and Matern-3/2 covariance | `spatial_covariance()` |
| Estimated exponential covariance with constrained profile search | `estimate_covariance()` |
| Exact known-covariance MSE, full-neighbour decomposition and positive-noise bound | `gaussian_cv()` |
| Conditional target MSE under a specified Gaussian generating distribution | `gaussian_target_mse()` |

The new estimator profiles an explicit bounded weighted variogram objective.
It is not substituted into the retained paper computations: its bin-edge
handling and allowance for exactly zero spatial signal differ from the frozen
implementation. Estimated parameters do not make a reference calculation an
exact estimate of real-world bias. The package deliberately excludes the old
experimental leakage and error-correction interfaces.

## spLeakage: the fixed reproduction snapshot

Audit of the included source snapshot, version **0.1.0.9004**, on 12 November 2026.
The package implements the central nearest-neighbour diagnostics, but not every
calculation or study-specific rule in the paper. The repository contains the
additional research implementations needed to reproduce the results.

| Paper idea | Package support | Research implementation or limitation |
| --- | --- | --- |
| Training and target nearest-neighbour distances; signed distance contrast and Wasserstein discrepancy | `prediction_target()`, `detect_leakage()`; fields `A`, `W`, `gobs`, `fpred` | `common.R::distance_summary()` independently calculates these quantities |
| Correlation contrast J1 | `detect_leakage()$SLI_rho`, or `sli(x, "rho")` | Paper name J1 corresponds to the legacy API name SLI_rho |
| Squared-correlation contrast J2 | `detect_leakage()$SLI_rho2`, or `sli(x, "rho2")`, for `k = 1` | This is the paper's nearest-neighbour contrast; the experimental `k > 1` calculation is not its full-neighbour theory |
| Unsigned correlation discrepancy U1 | `detect_leakage()$W_rho` for `k = 1` | U2 is implemented by `distance_summary()` in the research code, but is not a returned public package diagnostic |
| Random, spatial-block and buffered folds | Public fold constructors and `rank_cv_schemes()` | Paper-specific settings and buffer rules are in `common.R` and the study scripts; package defaults can differ |
| kNNDM | `detect_leakage()` accepts externally constructed fold lists | The paper constructs these with `CAST::knndm()`; the package ranking helper does not construct kNNDM |
| Select among four schemes by smallest absolute J2 | **Not supported by `rank_cv_schemes()`** | Research analysis scripts evaluate all four schemes; the package helper accepts `W`, `W_rho` and `abs_sli` (absolute J1) |
| Equal weighting of evaluation losses for tied schemes | Not the ranking helper's behaviour | Explicitly implemented in the paper's analysis scripts; it averages tied losses, not their signed errors |
| General target distribution Q or unequal target weights | Target-coordinate API uses equal weights | The empirical studies also use equal weights; arbitrary weighted-Q calculations from the theory are not exposed as an API |
| Exact one-neighbour MSE relation, full-neighbour decomposition and bounds | No dedicated public functions | `common.R`, `check_theory.R`, `check_additional_theory.R` |
| Exact conditional target MSE with refitted covariance | No public oracle-risk function | `run_nested_estimation.R`; this uses the simulation's known generating distribution for evaluation only |
| Covariance estimation | `estimate_dependence()` implements the original exponential variogram fitter | The Jura profile-search refinement is in `jura_profile_estimator.R`, not the package API |
| Ordinary/universal kriging application and learner experiments | Research workflows rather than one public package entry point | `jura_common.R`, `jura_profile_fit.R`, `run_learner_experiment.R` |

`application_walkthrough.R` verifies that the package and paper code return the
same A, W, J1, J2 and U1 for the Jura folds. The full reproduction additionally
checks the theoretical identities and compares application predictions with
`gstat`. Legacy `estimate_optimism()` and `deleak_estimate()` are not implementations
of the paper's exact full-neighbour MSE decomposition or a general bias correction.

The snapshot is preserved byte-for-byte in its scientific source. New features
are implemented and tested separately in spCVAudit. Reproduction does not install
or depend on the latest development version of either package.
