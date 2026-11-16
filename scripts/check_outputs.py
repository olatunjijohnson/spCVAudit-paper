"""Compare regenerated numerical outputs with the author's retained paper results."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("compare", ROOT / "revision/compare_summaries.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
checks = []
for reference, actual in [
    (ROOT / "provenance/expected-summaries", ROOT / "revision/results/summaries"),
    (ROOT / "provenance/expected-fold-count", ROOT / "revision/results/supplementary/fold-count")]:
    for file in sorted(reference.glob("*.csv")):
        target = actual / file.name
        result = module.compare_csv(file, target, 1e-10) if target.exists() else {"equal": False, "reason": "Missing output"}
        checks.append({"file": str(target.relative_to(ROOT)), **result})
for name in ("numbers.tex", "table-reference.tex", "table-jura.tex"):
    checks.append({"file": name, "equal": (ROOT / "provenance" / name).read_bytes() ==
                   (ROOT / "revision/paper-outputs" / name).read_bytes(), "comparison": "exact formatted text"})
for name in ("target-spatial-blocks.csv", "target-block-sensitivity.csv"):
    checks.append({"file": name, **module.compare_csv(ROOT / "provenance" / name,
                   ROOT / "revision/results/supplementary/jura" / name, 1e-10)})
checks.append({"file": "selection-tolerance.csv", **module.compare_csv(
    ROOT / "provenance/selection-tolerance.csv", ROOT / "revision/audit/selection-tolerance.csv", 1e-10)})
record = {"all_equal": all(x["equal"] for x in checks), "absolute_tolerance": 1e-10, "checks": checks}
output = ROOT / "revision/audit/current-reproduction-check.json"
output.write_text(json.dumps(record, indent=2) + "\n")
if not record["all_equal"]:
    raise RuntimeError("Reproduction comparison failed; inspect " + str(output))
print("Verified", len(checks), "retained numerical outputs and formatted manuscript tables.")
