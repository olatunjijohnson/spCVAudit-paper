"""Analyse ten-fold results separately, preserving the primary five-fold summaries."""
import os
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--no-pdf", action="store_true", help="Omit LaTeX compilation of the report")
args = parser.parse_args()
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
subprocess.run(["Rscript", "--vanilla", "revision/exploratory/plot_fold_counts.R",
                str(OUT / "fold-count-comparisons.csv"),
                "revision/exploratory/fold-count-comparison.pdf"], cwd=ROOT, check=True)
subprocess.run([sys.executable, "revision/exploratory/make_fold_count_report.py", str(OUT),
                "revision/exploratory/fold-count-report.tex"], cwd=ROOT, check=True)
if not args.no_pdf:
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
                    "fold-count-report.tex"], cwd=ROOT / "revision/exploratory", check=True)
print("Wrote paired five-versus-ten-fold comparisons:", OUT)
