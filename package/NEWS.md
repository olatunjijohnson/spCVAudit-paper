# spLeakage 0.1.0.9004 (geographic units and bounded pair allocation)

* Automatic buffer limits and target-matched buffer candidates now use metres
  for geographic coordinates, instead of comparing metre distances with a
  degree-based coordinate extent. Geographic extent is the largest observed
  geodesic separation; planar inputs retain the bounding-box diagonal.
* The default inverse-distance predictor in legacy score-comparison functions
  now uses geodesic distances for geographic inputs, including across the date
  line. Coordinate k-means is still a planar construction and is documented as
  such; geographic data should be projected when metric blocks are required.
* Variogram pair subsampling allocates only the requested pairs. It no longer
  constructs all quadratic pairs before reducing them to `max_pairs`. The
  sampled pair order remains consistent with the previous algorithm for the
  same seed. Regression checks include a five-billion-pair case using 500 pairs.

# spLeakage 0.1.0.9003 (bounded-fit roundoff)

* Enforce exact feasibility after bounded variogram optimisation when a parameter
  lies outside its bound only by floating-point roundoff. Larger violations are
  errors. A nested-CV pilot exposed a nugget of approximately -1e-18; its synthetic
  variogram is retained as a regression test. This does not clip unconstrained
  nonphysical fits or supply an uncertainty covariance for constrained fits.

# spLeakage 0.1.0.9002 (interpretation and documentation)

* Nonphysical unconstrained variogram fits are now refitted under the existing
  nonnegative parameter bounds, instead of clipping negative parameters and
  retaining an incompatible parameter covariance. Constant responses return zero
  fitted signal. Failed constrained optimization is reported explicitly.
  Estimated-dependence sensitivity results and the application were regenerated.
  The constrained fit is scaled to avoid coordinate-unit-dependent convergence.
* Geographic nearest neighbours are selected using geographic distances even
  when `sf`'s s2 backend is disabled. Audits reject incompatible declared CRSs.
  The target constructor rejects empty/nonfinite coordinates and documents its
  equal-location weighting and rectangular interpolation target.
* `recommend_validation()` no longer declares random CV unbiased because a
  sample is probabilistic, declares spatial CV appropriate because locations
  are clustered, or silently substitutes a clustered design for an unknown one.
  Prediction recommendations list candidate assessments and return indeterminate
  spatial-CV suitability. Fixed-map assessment requires a separate validation
  sample and the assumptions of its survey estimator.
* Workflow reports no longer print arbitrary letter grades as validity judgments.
  Historical grade fields remain in returned objects for compatibility.
* The validation-location heuristic now identifies its weights and precision
  calculations as model-dependent; spatial spread does not establish known
  inclusion probabilities or a design-unbiased estimator.
* README and tutorials explain the revised mismatch interpretation and the
  restrictions of legacy score-comparison and experimental functions.

# spLeakage 0.1.0.9001 (development revision)

The core interpretation is now validation–target mismatch. A positive signed
correlation contrast is not a universal estimate of predictive optimism.

* `detect_leakage()` now reports `SLI_rho2`, the squared-correlation contrast for
  one-neighbour prediction, `W_rho`, and whether distance ECDFs cross. The existing
  `SLI_rho` remains the first-moment contrast; its definition is not silently
  replaced by the squared form.
* `estimate_dependence()` treats its nugget as independent observation noise.
  Distinct co-located observations therefore have correlation equal to the signal
  fraction, rather than one. Models with shared observation errors need a
  different dependence specification. Historical results involving co-location
  must be recomputed.
* Missing parameter covariance no longer triggers an arbitrary 20% perturbation
  interval. Parameter propagation conditions on the fitted family and geometry;
  it does not quantify uncertainty in predictive bias. Propagation with the
  heuristic `k > 1` option is rejected.
* Split checks reject training/test overlap, invalid row indices, empty folds and
  repeated assessment of a row. Audit repeated CV one repeat at a time.
* Geographic k-neighbour calculations now use geographic distances for every k.
* `rank_cv_schemes()` defaults to unsigned Wasserstein mismatch. Its `criterion`
  argument retains signed-contrast ranking for explicit comparisons. Ranking
  preserves the caller's random-number state when using its explicit seed.
* Removed an internal global random-seed reset in target comparison selection.
  Historical repeated-sampling studies using that route require regeneration.
* Target comparison selection now minimizes unsigned distance mismatch. Legacy
  `estimate_optimism()` and `deleak_estimate()` print CV-score contrasts and
  descriptive resampling intervals; they no longer describe a comparison score
  or a rescaled published score as established deployment accuracy. Historical
  result-field names remain for compatibility.

The emulator and multi-channel extensions are outside the focused methodological
paper. Their presence in this package is not evidence of general accuracy
correction. The redevelopment record and fresh experiments are in `../revision/`.
