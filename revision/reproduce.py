#!/usr/bin/env python3
"""Rebuild the research outputs, or compute all studies in a new directory.

The formatted values and tables cited in the paper are written to
revision/paper-outputs; the manuscript text is not part of this repository.

Run from any directory. No dependencies are downloaded. Computation requires
the installed R dependencies described in REPRODUCING.md. Existing batch files
are never copied into a compute run.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


STUDIES = (
    ("run_gp_reference.R", "gp-reference", ""),
    ("run_gp_knndm.R", "gp-knndm", ""),
    ("run_learner_experiment.R", "learners", ""),
    ("run_sensitivity.R", "sensitivity-rng", "-rng"),
    ("run_jura_application.R", "application-jura", ""),
    ("run_nested_estimation.R", "nested-estimation-rng", "-rng"),
)
ANALYSES = (
    "check_theory.R", "check_additional_theory.R", "analyse_reference.R",
    "analyse_learners.R", "analyse_sensitivity.R", "analyse_jura_application.R",
    "analyse_nested_estimation.R", "analyse_selection_behaviour.R", "make_figures.R",
    "make_supplement_figures.R", "make_manuscript_numbers.R",
)


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def checksum(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inventory(root, paths):
    return {str(p.relative_to(root)): {"sha256": checksum(p), "bytes": p.stat().st_size}
            for p in sorted(set(paths)) if p.is_file()}


def source_files(root):
    files = [root / "revision/reproduce.py", root / "revision/EXPERIMENTS.md"]
    for folder in ("revision/R", "revision/data", "package/R", "package/man", "package/tests"):
        files.extend((root / folder).rglob("*"))
    for name in ("DESCRIPTION", "NAMESPACE", "LICENSE", "NEWS.md"):
        files.append(root / "package" / name)
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("render", "compute"))
    parser.add_argument("--output", type=Path, help="Required new directory for compute")
    parser.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_NTASKS", "1")))
    parser.add_argument("--compare-to", type=Path, help="Compare regenerated scientific results with this project root")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.mode == "compute" and args.output is None:
        parser.error("compute requires --output pointing to a new directory")
    if args.mode == "render" and (args.output is not None or args.compare_to is not None):
        parser.error("--output and --compare-to apply only to compute")
    if not shutil.which("Rscript") or not shutil.which("R"):
        parser.error("R and Rscript must be on PATH")

    source = Path(__file__).resolve().parents[1]
    root = source
    env = os.environ.copy()
    for name in tuple(env):
        if name.startswith("REVISION_"):
            del env[name]
    env.update(SLURM_NTASKS=str(args.workers), OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    # Resolve installed dependency libraries before changing working directory.
    libraries = subprocess.check_output(
        ["Rscript", "--vanilla", "-e",
         'source("revision/R/local_libraries.R"); cat(.libPaths(), sep="\\n")'],
        cwd=source, env=env, text=True).strip().splitlines()

    if args.mode == "compute":
        root = args.output.expanduser().resolve()
        if root.exists():
            parser.error("compute output already exists; use a new directory to exclude stale checkpoints")
        root.mkdir(parents=True)
        shutil.copytree(source / "package", root / "package",
                        ignore=shutil.ignore_patterns(".git", ".DS_Store", "*.Rcheck"))
        for folder in ("R", "data", "exploratory"):
            shutil.copytree(source / "revision" / folder, root / "revision" / folder)
        for folder in ("scripts", "environment", "provenance", "docs"):
            shutil.copytree(source / folder, root / folder)
        for name in ("reproduce.py", "README.md", "REPRODUCING.md", "DATA-SOURCES.md", "CITATION.cff", "LICENSE"):
            shutil.copy2(source / name, root / name)
        (root / "revision/audit").mkdir(parents=True)
        shutil.copy2(source / "revision/audit/fold-count-sensitivity-protocol.json",
                     root / "revision/audit/fold-count-sensitivity-protocol.json")
        for name in ("reproduce.py", "EXPERIMENTS.md", "REPRODUCING.md"):
            shutil.copy2(source / "revision" / name, root / "revision" / name)
        library = root / "revision/csf/library"
        library.mkdir(parents=True)
        env["R_LIBS_USER"] = os.pathsep.join([str(library)] + libraries)

    library = root / "revision/csf/library"
    library.mkdir(parents=True, exist_ok=True)
    env["R_LIBS_USER"] = os.pathsep.join([str(library)] + libraries)
    env["MPLCONFIGDIR"] = str(root / "build/matplotlib")
    for folder in ("results", "results/summaries", "figures", "logs", "audit"):
        (root / "revision" / folder).mkdir(parents=True, exist_ok=True)
    tag = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = root / "revision/audit" / ("reproduction-" + tag + ".json")
    record = {"started_utc": utc(), "mode": args.mode, "workers": args.workers,
              "source_root": str(source), "output_root": str(root),
              "command": sys.argv, "inputs": inventory(root, source_files(root)), "stages": []}
    if args.mode == "render":
        inputs = []
        for _, stem, _ in STUDIES:
            for suffix in (".rds", "-configuration.csv", "-errors.rds"):
                if stem == "gp-knndm" and suffix == "-configuration.csv":
                    continue  # This paired run uses the Gaussian reference grid.
                path = root / "revision/results" / (stem + suffix)
                if not path.is_file():
                    parser.error("Missing retained input: " + str(path))
                inputs.append(path)
        record["inputs"].update(inventory(root, inputs))

    def save_record():
        manifest_path.write_text(json.dumps(record, indent=2) + "\n")

    def stage(name, command, cwd=None, extra_env=None):
        stage_env = env.copy()
        stage_env.update(extra_env or {})
        log = root / "revision/logs" / (tag + "-" + name + ".log")
        item = {"name": name, "command": command, "started_utc": utc(),
                "cwd": str(cwd or root), "environment_overrides": extra_env or {},
                "log": str(log.relative_to(root))}
        record["stages"].append(item)
        save_record()
        print(utc(), name, "->", log, flush=True)
        with log.open("w") as handle:
            result = subprocess.run(command, cwd=cwd or root, env=stage_env,
                                    stdout=handle, stderr=subprocess.STDOUT)
        item.update(finished_utc=utc(), exit_code=result.returncode, log_sha256=checksum(log))
        save_record()
        if result.returncode:
            print("\n".join(log.read_text(errors="replace").splitlines()[-20:]),
                  file=sys.stderr, flush=True)
            raise RuntimeError(f"{name} failed; inspect {log}")

    try:
        stage("environment", ["Rscript", "--vanilla", "revision/R/check_environment.R"])
        stage("dependencies", ["Rscript", "--vanilla", "-e",
              'source("revision/R/local_libraries.R"); '
              'stopifnot(packageVersion("CAST")=="1.1.1", packageVersion("ranger")=="0.16.0"); '
              'for(p in c("sf","ggplot2","patchwork","digest","gstat","pkgload","jsonlite")) stopifnot(requireNamespace(p,quietly=TRUE))'])
        stage("install", ["R", "CMD", "INSTALL", "--library=revision/csf/library", "package"])
        if args.mode == "compute":
            if args.compare_to is not None:
                stage("reference-environment", ["Rscript", "--vanilla",
                      "revision/R/check_reference_fitter.R", str(args.compare_to.expanduser().resolve())])
            for script, stem, suffix in STUDIES:
                stage(stem, ["Rscript", "--vanilla", "revision/R/" + script],
                      extra_env={"REVISION_RUN_SUFFIX": suffix})
            stage("streams", ["Rscript", "--vanilla", "revision/R/audit_random_streams.R"])
            stage("geometry", ["Rscript", "--vanilla", "revision/R/audit_numeric_geometry.R"])
            if args.compare_to is not None:
                stage("comparison", ["Rscript", "--vanilla", "revision/R/compare_reproduction.R",
                      str(args.compare_to.expanduser().resolve()), str(root)] + [s[1] for s in STUDIES])
        for script in ANALYSES:
            stage(script.removesuffix(".R"), ["Rscript", "--vanilla", "revision/R/" + script])
        outputs = []
        for _, stem, _ in STUDIES:
            outputs.extend((root / "revision/results").glob(stem + ".*"))
            outputs.extend((root / "revision/results").glob(stem + "-configuration.csv"))
            outputs.extend((root / "revision/results").glob(stem + "-errors.rds"))
        for pattern in ("*-analysis.rds", "*theory-checks.*", "*neighbourhood-counterexample.csv",
                        "manuscript-number-provenance.csv"):
            outputs.extend((root / "revision/results").glob(pattern))
        if args.mode == "compute":
            for pattern in ("random-stream*-rng.csv", "random-stream-audit.csv", "random-stream-moments.csv",
                            "numeric-geometry-audit*.csv", "application-api-parity.csv",
                            "conditional-risk-check.csv", "full-reproduction-comparison.csv"):
                outputs.extend((root / "revision/results").glob(pattern))
        outputs.extend((root / "revision/results/summaries").glob("*.*"))
        outputs.extend((root / "revision/figures").glob("*.*"))
        outputs.extend((root / "revision/paper-outputs").glob("*.tex"))
        record.update(finished_utc=utc(), status="complete", outputs=inventory(root, outputs))
    except Exception as error:
        record.update(finished_utc=utc(), status="failed", error=str(error))
        save_record()
        raise
    save_record()
    print("Completed. Provenance:", manifest_path, flush=True)


if __name__ == "__main__":
    main()
