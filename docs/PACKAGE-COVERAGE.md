# Relationship between the paper and spCVAudit

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

## The frozen research implementation

The paper's results were computed with the research package preserved unchanged
in `package/`, together with the study scripts in `revision/R/`.
`application_walkthrough.R` verifies that this implementation and the paper code
return the same A, W, J1, J2 and U1 for the Jura folds. The full reproduction
additionally checks the theoretical identities and compares application
predictions with `gstat`. The frozen implementation is not maintained; new
analyses should use spCVAudit. Reproduction does not install or depend on the
latest development version of either package.
