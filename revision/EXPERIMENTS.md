# Experiments specified before the main runs

27 October 2026. Pilots check implementation and computational cost. Pilot results will be stored separately from the main runs. Changes to the design after inspecting results will be recorded here, with their reason.

## Questions

1. When does the dependence-weighted signed contrast describe the direction of the known model-based error gap, and when does it fail?
2. Does it add useful information beyond the mean distance gap, squared-correlation contrast and unsigned distribution mismatch?
3. Can the same spatial diagnostic have different predictive implications for different learners and different amounts of covariate-explained variation?
4. Does selecting a CV scheme by a small signed contrast improve target-risk estimation, or does cancellation defeat the rule?

## A. Algebra and counterexamples

Executable assertions will verify the one-neighbour identity, the co-located replicate formula, the stochastic-order sign result in numerical examples, the signed-area/mean identity and the correct transport bound. Deliberately crossing distributions will produce opposite first- and second-moment signs. Non-coincident neighbourhoods will demonstrate different full prediction variances at the same nearest-neighbour distance.

These examples establish existence and algebraic correctness. They do not estimate how common failures are in ecological applications.

## B. Exact Gaussian reference experiment

Generate independent spatial sampling configurations. Calculate expected MSE directly from the covariance matrices, avoiding simulation noise in this reference comparison. The replication unit is an independently generated configuration. It is not a test point or a fold.

Core factors: random and clustered locations; short and longer dependence ranges; zero, moderate and strong spatial signal; two training sample sizes; a target covering the sampled domain and a shifted target; exponential and Matérn covariance. Geometric anisotropy is a separate sensitivity experiment. The exact grid and replicate count will be saved before the main execution after a timing pilot.

Each configuration uses the same candidate schemes: random folds, spatial blocks, buffered blocks and, if the official implementation can be made available, kNNDM. All schemes predict the same target from the full training sample. Record training sizes and any infeasible folds. Never silently omit a failed fold and call the resulting test sample comparable.

Compare signed J₁, J₂ and mean-distance gap A; use W and transformed unsigned distances U₁ and U₂ as mismatch diagnostics and scheme-selection criteria, not as signed predictors. Also retain the unweighted latent-correlation contrasts to expose the effect of inconsistent signal scaling. The known Gaussian one-neighbour risk proxy is VJ₂; V is one in the normalized reference experiment.

Primary outcomes: sign agreement outside |O| ≤ 0.01V, false reassurance when a small signed contrast coexists with |O| > 0.05V, absolute error of the actual one-neighbour risk proxy, and the excess absolute risk bias of the selected scheme relative to the best available candidate. Report sensitivity to the negligible-effect threshold. Pooled correlation is secondary. Present results separately by design, range, signal, sample size, target and covariance family.

The one-neighbour transport inequality must hold to numerical tolerance. Full-kriging bound violations are reported as limitations of extrapolating that inequality. Test whether CDF ordering predicts one-neighbour sign agreement, keeping it distinct from full-kriging agreement.

## C. Covariates, residual dependence and learning algorithms

Generate spatial covariate X and residual field Z independently, with the same covariance family. Let

Y = sqrt(alpha) X + sqrt(1 − tau² − alpha) Z + epsilon,

with independent noise variance tau² and alpha varying from zero to 1 − tau². Marginal spatial covariance is then held constant across alpha, while residual spatial covariance given X changes. Within a replicate, retain the same geometry and underlying standardized fields across alpha to allow paired comparisons. Replicates use independently drawn fields and observation errors.

Assess a known mean-function predictor, fitted linear regression on X, a random forest using X, a spatial learner, and oracle residual kriging as a reference. Predict the same target from each full fitted model. The target loss integrates squared error against the known latent response surface and adds tau², so target observation noise need not be simulated merely to estimate mean loss. The finite target approximates the stated spatial distribution.

Compare marginal dependence diagnostics with known conditional dependence diagnostics. Estimated dependence is a separate analysis, using only training observations; any residual estimation must state how the mean fit was obtained. The experiment does not assume that a residual variogram automatically validates CV.

Analyze realized CV-minus-target loss gaps across independent field replicates. Distinguish these from the exact expected gaps in experiment B. A scheme-selection method may use locations, target covariates and the declared dependence model, but never target response values. An oracle selection using target loss is only a reference for scoring regret.

## D. Continuous-response application (updated 12 November 2026)

Use natural-log cadmium in the Swiss Jura, with the published 259 training and 100 reserved validation sites. Compare ordinary kriging and a coordinate-trend sensitivity, the four five-fold CV schemes, and selection by W and |J2|. Retain 40 repetitions per target and all three targets (published, random and northern). The pre-analysis design and post-hoc numerical-search protocol are preserved in [the original application protocol](../provenance/application-protocols-20261112/application-comparison-protocol.txt) and [the covariance-search protocol](../provenance/application-protocols-20261112/covariance-search-protocol.txt).

The manuscript uses the checked profile-search results, reproduced through `R/run_jura_application.R` in `results/application-jura.rds`. The observed target responses enter evaluation only. These repeated analyses of one soil survey do not estimate expected bias across independent fields or area-wide map error. The original exploratory comparison and both covariance-fitting variants are preserved. The dated application specification below is superseded and retained as research history.

## Reproducibility and reporting rules

- Save configuration tables, individual seeds, split fingerprints, software versions and source hashes.
- Keep independent field/configuration replicates distinct from repeated partitions of the same field.
- Report all planned configurations, failures and exclusions. Do not select a favorable covariance range after viewing results.
- Use paired within-replicate comparisons. Intervals across independent fields quantify Monte Carlo uncertainty; they do not represent general ecological uncertainty.
- Keep development pilots and subsequent sensitivity analyses identifiable.
- Rebuild manuscript numbers and figures from machine-readable results.
- If dependence weighting does not outperform simpler diagnostics, report that result and narrow the contribution accordingly.

## Execution decisions after the timing pilot

The 64-configuration Gaussian pilot completed 192 scheme evaluations in approximately 2.7 seconds, without infeasible folds. It is retained as `gp-pilot` and excluded from main summaries. The main Gaussian grid is fixed at n = 100 or 250; exponential range parameter 0.05 or 0.2 (the same scale parameter convention, not matched practical ranges, for Matérn 3/2); signal fraction 0, 0.3, 0.7 or 0.95; random or clustered design; domain or shifted target; exponential or Matérn 3/2 covariance; and 40 independent configurations per cell. This gives 5,120 configurations. Each target uses 225 midpoint grid locations. Candidate folds use k = 5 and a buffer of 0.05 in unit-domain coordinates.

The initial main run uses random, block and buffered schemes. CAST is not installed, and the execution environment disables shell network access, so an attempted CRAN installation could not access the repository index. This is an environment limitation, not evidence that CAST is unavailable for this R version. Official kNNDM comparison remains a tracked task; no substitute algorithm will be labelled as an evaluation of CAST.

**Dependency resolution before the kNNDM main run:** CAST 1.1.1 was found in the existing R 4.5 library and loads under the current interpreter. The official default hierarchical kNNDM implementation was tested on the same 64 pilot configurations: all succeeded in approximately 46 seconds, and every returned W value agreed with the independent calculation. The main paired kNNDM run therefore uses all 5,120 pre-specified configurations, k = 5 and the other CAST defaults. Sampling locations and preceding candidate-construction random draws are reproduced from each configuration's seed. This fourth scheme is generated without using any target response. Results remain separated from the pilot.

**Learner experiment timing and specification:** an eight-field pilot with n = 100, 144 target points and 100-tree forests completed 640 learner/alpha/scheme evaluations in approximately 18 seconds. Checks confirmed that the marginal diagnostic is identical across alpha and learners within each field/scheme; the conditional diagnostic vanishes at alpha = 0.9; the known-mean CV score is invariant to the scheme; and oracle prediction reaches the 0.1 noise floor when residual spatial variance is zero. The main run uses n = 200, 400 target points, 300-tree forests and 60 independent fields in each of eight range/design/target cells (480 fields). Each field evaluates alpha = 0, 0.3, 0.6 and 0.9, five learners and four schemes. This gives 38,400 evaluations, grouped into 480 independent field/configuration units. The known-mean learner deliberately uses no observed response: it is a control showing that residual dependence alone does not establish proximity-driven error bias.

## Continuation sensitivity specification — 28 October 2026

This extension is specified after a preliminary aggregate inspection of the three-scheme main Gaussian results and before running the sensitivity study. It therefore remains a separately labelled analysis. Use 1,440 independent configurations: n = 160, exponential range 0.05 or 0.2, signal fraction 0.3/0.7/0.95, random or clustered sampling, domain or shifted target, axis ratio 1 or 0.25 at angle pi/6, and 30 replicates per cell. Seeds start at 810001. The pilot uses distinct seeds starting at 710001, n = 80 and one replicate per cell.

The four existing validation schemes remain unchanged. For every scheme, decompose the full-model target-minus-CV risk into (a) the target risk using each fold's identical training subset, weighted by its test count, minus its CV risk and (b) the improvement from training on the full sample. This controls training size without treating the full-model deployment estimand as an error.

Compare the isotropic diagnostic with true range and signal against the package's isotropic variogram fit to one training field. Under anisotropy the isotropic diagnostic is deliberately misspecified; the actual one-neighbour covariance contrast is retained as a reference. Target responses are never supplied to covariance fitting. The model-averaged exact risk labels apply to the fixed geometry-only schemes. Selecting a scheme after fitting a covariance to responses does not make these labels the conditional risk of that response-adaptive selection procedure. Report fitting errors, boundary/extreme parameters and all configurations.

Compute the Schur-complement sign bound from known covariance and report both correctness and the fraction of cases in which it makes a sign determination. Its use of all-neighbour covariance must be explicit; it is not a nearest-neighbour-only diagnostic.

For selection comparisons, average loss over tied diagnostic minima. Target loss never breaks a tie. Monte Carlo standard errors pair schemes within each configuration and estimate simulation variance within the fixed factor cells. The headline sign statistic uses the predefined |O| > 0.01 zone, with thresholds 0.005, 0.02 and 0.05 retained as sensitivity analyses.

## Historical application specification — 28 October 2026 (superseded)

The primary illustration retains AWT species awt06, which was selected before this redevelopment because its 187 presences among 340 survey sites provide a relatively balanced response. Use all 13 supplied environmental predictors and a probability random forest with 300 trees, mtry = 3, minimum terminal node size = 10, fixed seeds and one thread per fitted forest; no tuning against any target outcome. Brier loss is the mean squared probability error. Background points are unlabelled, never absences.

For the map illustration, declare the 10,000 supplied background points with equal weights as the finite target and use all 340 survey sites for training. Compare five-fold random, k-means spatial blocks, blocks buffered by a fixed 5 km, and official kNNDM. Repeat fold generation 20 times for awt06. Report all 20 species for the first shared fold-generation seed as supporting descriptions, without treating species as independent spatial datasets or selecting species by a large score discrepancy.

Add two explicitly limited survey-holdout tasks for awt06: 85 survey sites selected uniformly without replacement, and the 85 northernmost sites selected by northing (site ID breaks ties). The remaining 255 sites form the training sample. Target labels are excluded from fold construction, dependence fitting and scheme selection. Repeat the uniform partition 20 times; for the fixed northern holdout repeat only CV fold generation 20 times while holding the final full-training forest seed fixed. These are repeated partitions or fits within one survey, not 20 independent fields or a probability sample of the background map. Report medians and quartiles as descriptive variation rather than confidence intervals for ecological generalization.

Fit the package exponential variogram to the training response and, separately, to random-CV residuals from the specified forest. Both are descriptive dependence specifications for this binary-response application. The Gaussian one-neighbour risk identity is not claimed for its random-forest Brier loss. This application can show how the audit changes with target and dependence specification, and compare scores on held-out survey sites; it cannot establish map-wide optimism or validate an accuracy correction.

### Variogram implementation repair and rerun

Inspection of the initial sensitivity output exposed an implementation problem: a converged unconstrained exponential fit with negative nugget or partial sill was accepted and then clipped to zero, without refitting or discarding its parameter covariance. A deterministic regression case confirms that this is not constrained weighted least squares. Version 0.1.0.9002 now rejects nonphysical unconstrained fits and uses the pre-existing bounded optimizer instead. Constant responses return zero fitted signal; optimizer failures are reported. This is a correction of fitting behavior, not a choice of covariance model based on diagnostic performance.

The first sensitivity and application outputs and the exact 0.1.0.9001 study package are preserved. The affected studies are repeated with the same predeclared configurations, seeds, models and scoring rules under suffix `-bounded`; final estimated-dependence results will use those runs. Known-covariance reference and learner experiments do not call this fitter and remain unchanged.

The first bounded refit exposed optimizer convergence failures in package checks, particularly when coordinate magnitudes changed. The bounded objective now normalizes distances and semivariances, uses its analytic gradient and a second constrained optimizer on failure. Bounds and statistical objective are unchanged. A unit-transformation regression check and the full local suite pass. The final retained fitted-dependence runs use suffix `-scaled`; both earlier stages remain available. The full known-covariance and learner studies are unaffected by this fitting correction.

## Nested covariance estimation and response-adaptive selection extension

Specified 28 October 2026, approximately 23:59 UTC, after completing the known-covariance and fitted-diagnostic studies, but before generating this extension's pilot or main outputs. This is an additional study, not an independent preregistration of the whole redevelopment. Its motivation is a remaining inferential gap: fixed-candidate model-average risk labels do not assess CV after response-adaptive diagnostic selection.

**Question.** Does selecting a geometry-defined candidate by a fitted dependence contrast improve the actual CV estimate of a fitted model's conditional target risk when covariance is estimated independently within each training fold?

**Design fixed before execution.** Exponential stationary Gaussian response with known zero mean, marginal variance one, scale 0.05 or 0.2, signal fraction 0.3, 0.7 or 0.95, uniform or six-cluster sampling, in-domain or shifted finite target. These 24 cells each have 80 independent training response/location realisations (1,920 configurations), n=160 and 15-by-15 equally weighted target midpoints. Main seeds are 1210000 + id. Pilot seeds are 1200000 + id with one configuration per cell and n=80, target 8-by-8. The pilot checks runtime and implementation; the main factor grid, metrics and replicates are fixed now.

**Candidates and predictions.** The same four five-fold geometry-only candidates are random, k-means blocks, 0.05-buffered blocks and official CAST kNNDM. Compare (a) known-covariance kriging and (b) estimated-covariance kriging, both with known zero mean. The full model fits the variogram using all training observations, while every CV model fits its own variogram using only that fold's training observations. No full-data covariance estimate is reused for a held-out prediction. The package's final scaled fitter is fixed; no method tuning follows from these results.

**Risk reference.** Under the true generating covariance, calculate the exact conditional target mean m_Q and variances v_Q given the realised full training response. For either final predictor p_Q, its conditional target risk is the equally weighted average of (p_Q - m_Q)^2 + v_Q. This integrates target responses analytically and includes independent target noise. It does not require simulated target labels. Observed CV mean squared error uses the held-out training responses. The score error is conditional target risk minus observed CV error; its average across independent training realisations measures model-average optimism, and its absolute/squared values assess the CV estimator against the conditional risk of the fitted predictor.

**Selection and controls.** Compare |A|, W, |J1|, |J2|, U1, U2 using a full-training-data fitted variogram, with known-covariance versions as diagnostic controls, fixed random and fixed kNNDM, plus a clearly labelled unattainable oracle using the smallest absolute realised score error. Average exact numerical ties within 1e-12. Target conditional-risk calculations cannot affect candidate construction, diagnostic selection or tuning. Full-data response use for diagnostic selection is deliberate and its selection effects are included in the actual CV score error.

**Reporting.** Report mean signed error, mean absolute error and root mean squared error, paired within a configuration, with Monte Carlo standard errors stratified by fixed factor cell. Do not call mean absolute realised error 'absolute expected bias'. Include all factor-cell results, all failures, fitting warnings and success denominators. Report known-covariance prediction and estimated-covariance prediction separately. A configuration with a failed fit is retained in the error file; failures are not replaced by a different predictor. Do not present this correctly specified exponential experiment as validation for forests, binary data, unknown mean functions, anisotropy or preferential sampling. These results may weaken or strengthen the practical selection claim and will be reported either way.

## Response-stream repair — 29 October 2026, 01:02 UTC

A post-run RNG audit found only 760 unique innovation vectors among the 1,440 sensitivity configurations. A trace establishes the cause: CAST 1.1.1's kNNDM calls `CreateSpacetimeFolds()`, whose default seed is drawn from a small set and then passed to `set.seed()`. Responses simulated after that call can reuse innovations across configurations. The old fitted sensitivity and nested-estimation runs therefore do not support the independent-replicate Monte Carlo uncertainties assumed in their summaries. Their numerical replay does not validate independence. All old results and scripts are retained for audit, including the failed diagnostic logs and a minimal seed-reset trace.

Before generating corrected results, fix the response seed for each configuration at `geometry_seed + 10000000`, separately from fold construction. Retain the exact factor grids, geometry seeds, candidate folds, covariance fitter, predictors and evaluation metrics. Both studies, including pilots, are rerun under suffix `-rng`; the sensitivity and nested `-final` results are superseded. No replacement seed is chosen based on outcomes. Verify the unique innovation fingerprints, reconstruct the same fold fingerprints, and repeat the conditional-risk and independent replay checks. Report the corrected results regardless of direction.

The exact Gaussian reference has no simulated response draws. The learner study generates all Gaussian fields and observation noise before calling CAST and sets explicit forest seeds. The application has observed responses and explicitly seeded forests. Those studies are unaffected by this response-stream defect. Public package source is unchanged by this study-specific repair.

## Documentation audit — 29 October 2026, 02:20 UTC

The implemented and reported sign convention throughout the retained studies is **target loss minus CV loss**. The earlier planning sentence in Section C says the reverse; it is preserved above as history and corrected by this note. No scores or scientific conclusions are changed. The theory now states finite first distance moments, nonempty training sets disjoint from their tests, and nonsingular covariance matrices explicitly.

## Descriptive selection behaviour — 29 October 2026, 02:37 UTC

After inspecting the corrected nested-study mean differences and the all-cell omission analysis, examine selection frequencies, paired loss-difference distributions and single-replicate influence. An exploratory inspection already showed that the strongest cell often selects random folds with fitted J2 and separated folds with W. This analysis is post hoc. Report every factor cell for both learners, include numerical tie averaging, and keep the predeclared primary MAE comparison unchanged. Single-replicate omission summaries retain equal factor-cell weights; they describe influence and are not confirmatory exclusions. No new responses, tuned criteria or selected additional replications will be generated.

A subsequent numerical check varies the absolute tie tolerance from exact equality through 1e-8 for fitted J2 and W. This addresses possible rank changes when a fitted signal fraction makes all weighted contrasts tiny. The primary 1e-12 rule is unchanged. All six tolerances and both predictors are reported as a post-hoc audit, without choosing a replacement based on performance.

## Interpretation and metadata clarification — 29 October 2026, 03:28 UTC

The known-covariance setting permits exact target-risk calculation. Its criterion historically labelled `oracle` is therefore described as an exact-risk benchmark, rather than inherently unattainable under those assumptions. The fitted study's corresponding truth-based benchmark uses conditional risk withheld from the diagnostic rules. This changes neither selection nor any scientific value. The abstract now explicitly includes the response-variance multiplier in the one-neighbour identity. The learner's known mean conditions on the covariate alone, distinct from conditioning on all training responses.

Repeated uniform survey holdouts change the training presence count. The application summary now reports its median and range rather than copying the first repetition's count. All non-count summary columns and all 119 manuscript numeric macros remain unchanged. This is a descriptive metadata correction, not a change to the experiment or scoring rule.


## Descriptive rank reproducibility repair — 29 October 2026, 05:08 UTC

The isolated replay agrees with every scientific input column to at most 1.70e-13, but its summary comparison exposed unstable exact-tie Spearman correlations when effectively identical candidate risks differ in their last floating-point bits. These rank correlations are supplementary CSV diagnostics and do not enter candidate selection, primary errors, Monte Carlo intervals or manuscript numbers. Before regenerating them, fix an absolute tie tolerance of 1e-12 on both ranking variables, matching the tolerance already used for diagnostic selection. Sort values; each tie group contains values no more than that tolerance above its smallest member, and receives their average rank. This avoids chaining a sequence of small adjacent gaps into a group with a larger total span. The original script and two affected summary tables are archived. All other analyses remain fixed, and both original and isolated raw results will be checked under this rule.
