# spLeakage

Audit spatial validation against a declared prediction target.

`spLeakage` compares the distances from validation observations to their fold's
training data with the distances from target locations to the full training
sample. It reports signed correlation contrasts, unsigned distribution mismatch
and whether the distance distributions cross.

A positive signed contrast means greater average nearest-neighbour correlation
under validation, given the specified dependence model. It does not by itself
establish predictive optimism. A value near zero can arise through cancellation
or omit information supplied by other neighbours and covariates.

The focused workflow is:

1. Declare the equally weighted locations that define the prediction target.
2. Inspect validation and target distance distributions.
3. Compare the original correlation contrast (`SLI_rho`) with the squared
   contrast (`SLI_rho2`), unsigned mismatch (`W`) and distribution crossings.
4. Compare feasible validation schemes and state their training sizes.
5. Evaluate predictive error using a reference appropriate to the actual learner
   and target, where such a reference is available.

The squared contrast has an exact risk interpretation for a specified
one-neighbour linear predictor with known mean and covariance. Applying it to
full kriging or another learning algorithm requires further assumptions or
empirical evidence. The development record and fresh experiments are in
[`../revision/`](../revision/README.md).

## Example

```r
library(spLeakage)

target <- prediction_target(grid = prediction_locations, type = "grid")
audit <- detect_leakage(data, split = folds, target = target,
                        response = "y", coords = c("x", "y"))
audit
plot(audit)

# Compare candidate distance mismatches. The default criterion is W.
rank_cv_schemes(data, target, response = "y", coords = c("x", "y"))
```

Use an appropriate projected coordinate system or supply an `sf` object with a
correct CRS. The target should describe the intended predictions, including its
weights. The current constructor weights supplied locations equally; nonuniform
target integration requires separate calculations. Parameter propagation, when available, conditions on the chosen
covariance family and geometry; it is not an interval for predictive bias.

## Development status

The current development version is recorded in [DESCRIPTION](DESCRIPTION).
See [NEWS.md](NEWS.md) for behavior changes. Historical results that depended on the former seed reset, zero-distance
correlation override or arbitrary parameter-uncertainty fallback require
recomputation. A previous package-check claim is not evidence about this revision;
current checks are recorded with the redevelopment outputs.

`estimate_optimism()` and `deleak_estimate()` retain their historical names and
fields but compare CV scores. A blocked or target-matched score is not established
deployment truth. Their resampling intervals and rescaled reported scores must
not be presented as corrected map accuracy. The emulator, multi-channel
extensions and validation-location heuristic remain experimental and are outside
the focused methodological paper.

## Installation and documentation

From the parent project workspace, install the working source with:

```sh
R CMD INSTALL package
```

The core package imports `sf`. Optional functions require the packages listed in
`Suggests`; install only those needed for the intended analysis.

- `vignette("spLeakage")`: auditing validation–target mismatch.
- `vignette("quantify-correct")`: comparing scores and interpreting legacy APIs.
- `vignette("channels")`: additional descriptive checks and their limits.

## License

MIT © Olatunji Johnson
