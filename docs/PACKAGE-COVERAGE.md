# Relationship between the paper and spLeakage

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

No missing package features were added as part of preparing this repository.
Priorities for a later package update would be a four-scheme selection API with
absolute-J2 ranking and explicit tie handling, followed by U2 and an explicit
choice of covariance fitter. Any such update should be validated separately
against this fixed reproduction snapshot.
