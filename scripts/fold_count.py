"""Analyse ten-fold results separately, preserving the primary five-fold summaries."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "revision/results/supplementary/fold-count"
OUT.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix="spleakage-fold-count-") as temp:
    work = Path(temp)
    for folder in ("revision/R", "revision/exploratory"):
        shutil.copytree(ROOT / folder, work / folder)
    result = work / "revision/results"
    (result / "summaries").mkdir(parents=True)
    for suffix in (".rds", "-configuration.csv", "-errors.rds"):
        name = "nested-estimation-k10" + suffix
        shutil.copy2(ROOT / "revision/results" / name, result / name)
    env = os.environ.copy()
    env["REVISION_NESTED_SUFFIX"] = "-k10"
    env["R_LIBS_USER"] = os.pathsep.join(filter(None, [str(ROOT / "revision/csf/library"),
                                                    env.get("R_LIBS_USER")]))
    for command in (["Rscript", "--vanilla", "revision/R/analyse_nested_estimation.R"],
                    ["Rscript", "--vanilla", "revision/exploratory/compare_fold_counts.R", str(ROOT)]):
        subprocess.run(command, cwd=work, env=env, check=True)
    for path in result.glob("fold-count-*.csv"):
        shutil.copy2(path, OUT / path.name)
    shutil.copytree(result / "summaries", OUT / "ten-fold-summaries", dirs_exist_ok=True)
print("Wrote paired five-versus-ten-fold comparisons:", OUT)
