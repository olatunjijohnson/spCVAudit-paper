#!/usr/bin/env python3
"""Reproduce the paper's results, including its supplementary analyses."""
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("render", "compute"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--no-pdf", action="store_true", help="Skip the exploratory fold-count report PDF")
    parser.add_argument("--compare-to", type=Path)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    root = args.output.expanduser().resolve() if args.mode == "compute" and args.output else source
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               SLURM_NTASKS=str(args.workers), MPLCONFIGDIR=str(root / "build/matplotlib"))
    primary_args = [args.mode, "--workers", str(args.workers)]
    if args.output:
        primary_args += ["--output", str(args.output.expanduser().resolve())]
    if args.compare_to:
        primary_args += ["--compare-to", str(args.compare_to.expanduser().resolve())]
    subprocess.run([sys.executable, str(source / "revision/reproduce.py"), *primary_args],
                   cwd=source, env=env, check=True)
    commands = []
    if args.mode == "compute":
        commands += [["Rscript", "--vanilla", "revision/exploratory/run_nested_k10.R"],
                     ["Rscript", "--vanilla", "revision/R/run_jura_application.R", "--original-fit"]]
    fold_command = [sys.executable, "scripts/fold_count.py"]
    if args.no_pdf:
        fold_command.append("--no-pdf")
    implementation_command = ["Rscript", "--vanilla", "scripts/check_jura_implementation.R"]
    if args.mode == "compute" and args.compare_to is not None:
        implementation_command.append("--strict-reference")
    commands += [fold_command,
                 ["Rscript", "--vanilla", "revision/R/audit_selection_tolerance.R"],
                 [sys.executable, "scripts/jura_sensitivity.py"],
                 [sys.executable, "scripts/check_jura_fitting.py"],
                 implementation_command,
                 ["Rscript", "--vanilla", "revision/R/application_walkthrough.R"],
                 [sys.executable, "scripts/check_outputs.py"]]
    for command in commands:
        print("Running", " ".join(command), flush=True)
        subprocess.run(command, cwd=root, env=env, check=True)
    print("Completed primary and supplementary reproduction:", root, flush=True)


if __name__ == "__main__":
    main()
