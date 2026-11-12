"""Run the existing small simulation grids in an isolated directory, without saved results."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, help="New directory; otherwise use a temporary directory")
args = parser.parse_args()
temporary = tempfile.TemporaryDirectory(prefix="spleakage-smoke-") if args.output is None else None
work = Path(temporary.name) if temporary else args.output.expanduser().resolve()
if not temporary:
    work.mkdir(parents=True, exist_ok=False)
env = os.environ.copy()
libraries = subprocess.check_output(["Rscript", "--vanilla", "-e",
    'source("revision/R/local_libraries.R");cat(.libPaths(),sep="\\n")'], cwd=ROOT, text=True).strip().splitlines()
env.update(R_LIBS_USER=os.pathsep.join(libraries), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", SLURM_NTASKS="1")
record = {"started_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "stages": [], "status": "running"}
receipt = ROOT / "revision/audit/simulation-smoke.json"
try:
    for folder in ("package", "revision/R", "revision/exploratory"):
        shutil.copytree(ROOT / folder, work / folder, ignore=shutil.ignore_patterns("__pycache__"))
    for folder in ("revision/results", "revision/audit", "logs"):
        (work / folder).mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "revision/audit/fold-count-sensitivity-protocol.json", work / "revision/audit")
    for script in ("R/run_gp_reference.R", "R/run_gp_knndm.R", "R/run_learner_experiment.R",
                   "R/run_sensitivity.R", "R/run_nested_estimation.R", "exploratory/run_nested_k10.R"):
        name = Path(script).stem
        print("Smoke computation:", name, flush=True)
        with (work / "logs" / (name + ".log")).open("w") as log:
            subprocess.run(["Rscript", "--vanilla", "revision/" + script, "--pilot"],
                           cwd=work, env={**env, "REVISION_RUN_SUFFIX": ""},
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        record["stages"].append({"script": script, "exit_code": 0})
        receipt.write_text(json.dumps(record, indent=2) + "\n")
    verification = r'''
files <- list.files('revision/results',pattern='-errors[.]rds$',full.names=TRUE)
stopifnot(length(files)==6L,all(vapply(files,function(p) length(readRDS(p))==0L,logical(1))))
bundles <- c('gp-pilot','gp-knndm-pilot','learners-pilot','sensitivity-pilot','nested-estimation-pilot','nested-k10-pilot')
counts <- vapply(bundles,function(s) nrow(readRDS(paste0('revision/results/',s,'.rds'))$results),integer(1))
stopifnot(all(counts>0L))
jsonlite::write_json(as.list(counts),'revision/audit/smoke-counts.json',auto_unbox=TRUE)
'''
    subprocess.run(["Rscript", "--vanilla", "-e", verification], cwd=work, env=env, check=True)
    record.update(status="complete", evaluation_counts=json.loads((work / "revision/audit/smoke-counts.json").read_text()))
except Exception as error:
    record.update(status="failed", error=str(error), work_directory=str(work))
    raise
finally:
    record["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    receipt.write_text(json.dumps(record, indent=2) + "\n")
    if temporary:
        temporary.cleanup()
print("All six simulation pilot computations completed with no failed configurations.")
