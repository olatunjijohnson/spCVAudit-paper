"""Reproduce Jura covariance-optimisation checks and the fitting sensitivity."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
from jura_sensitivity import selection_weights

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "revision/results/supplementary/jura"
EXPORT = ROOT / "build/jura-export"
old = pd.read_csv(ROOT / "revision/results/application-jura-original.csv")
new = pd.read_csv(EXPORT / "results.csv")
key = ["dataset", "target", "repeat_id", "model", "scheme"]
old = old.sort_values(key).reset_index(drop=True)
new = new.sort_values(key).reset_index(drop=True)
assert old[key].equals(new[key])
for column in ["A", "W", "buffer_km", "min_train", "mean_train", "n_train", "n_target", "split_seed", "fold_seed"]:
    np.testing.assert_allclose(old[column], new[column], rtol=0, atol=1e-12)
fits = pd.read_csv(EXPORT / "fits.csv")
assert (fits.profile_objective <= fits.original_objective + 1e-12).all()
details = json.loads((EXPORT / "geometry.json").read_text())
geometry = next(g for g in details if g["target"] == "published" and g["repeat_id"] == 1)
checks = []
for model in ("ordinary", "trend"):
    v = geometry["variograms"][model]
    empirical = pd.DataFrame(v["empirical"])
    h = empirical.dist.to_numpy() / empirical.dist.max()
    y = empirical.gamma.to_numpy() / empirical.gamma.max()
    w = empirical.np.to_numpy() / empirical.np.sum()

    def objective(p):
        return float(np.sum(w * (p[0] + p[1] * (-np.expm1(-h / np.exp(p[2]))) - y) ** 2))

    search = differential_evolution(objective, [(0, 1), (1e-8, 2), (np.log(1e-6), np.log(5))],
                                    seed=1600000 + len(checks), popsize=16, maxiter=500,
                                    tol=1e-9, polish=True)
    actual = objective([v["nugget"] / empirical.gamma.max(), v["psill"] / empirical.gamma.max(),
                        np.log(v["range_km"] / empirical.dist.max())])
    assert actual <= search.fun + 1e-8
    checks.append(dict(model=model, R_profile_objective=actual,
                       independent_scipy_objective=search.fun, difference=actual-search.fun))
pd.DataFrame(checks).to_csv(OUT / "independent-optimiser-checks.csv", index=False)
comparisons = []
groups = ["target", "model", "repeat_id"]
for values, frame in new.groupby(groups):
    reference = old
    for col, value in zip(groups, values):
        reference = reference[reference[col] == value]
    for rule in ("W", "J2", "fixed_random", "fixed_knndm"):
        comparisons.append(dict(zip(groups, values)) | dict(
            rule=rule, original_error=float(selection_weights(reference, rule) @ reference.absolute_error),
            profile_error=float(selection_weights(frame, rule) @ frame.absolute_error)))
comparison = pd.DataFrame(comparisons)
summary = comparison.groupby(["target", "model", "rule"])[["original_error", "profile_error"]].mean()
summary["difference"] = summary.profile_error - summary.original_error
summary.to_csv(OUT / "fitting-sensitivity.csv")
record = dict(covariance_fits=len(fits), all_profile_objectives_no_worse=True,
              improved_fits=int(fits.improved.sum()), independent_checks=len(checks),
              largest_absolute_mae_change=float(summary.difference.abs().max()))
(OUT / "fitting-verification.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
