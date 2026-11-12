#!/usr/bin/env python3
"""Generate the post-hoc fold-count report and its numeric provenance."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("results", type=Path)
parser.add_argument("output", type=Path, help="Output TeX file")
args = parser.parse_args()
comparisons = list(csv.DictReader((args.results / "fold-count-comparisons.csv").open()))
changes = list(csv.DictReader((args.results / "fold-count-change-summary.csv").open()))
verification = list(csv.DictReader((args.results / "fold-count-verification.csv").open()))
assert len(verification) == 8 and all(float(x["maximum_absolute_difference"]) == 0 for x in verification)


def comparison(learner, folds, criterion="J2_fitted"):
    rows = [r for r in comparisons if r["learner"] == learner and r["fold_count"] == str(folds)
            and r["criterion"] == criterion and r["stratum"] == "all"]
    assert len(rows) == 1 and rows[0]["n"] == "1920"
    return rows[0]


def change(learner):
    rows = [r for r in changes if r["learner"] == learner and r["stratum"] == "all"]
    assert len(rows) == 1 and rows[0]["n"] == "1920"
    return rows[0]


five, ten, delta = comparison("estimated", 5), comparison("estimated", 10), change("estimated")
number = lambda x: f"{float(x):.5f}"
percent = lambda row: f"{-100 * float(row['difference']) / float(row['W_mae']):.1f}"
assert float(ten["lower"]) < 0 < float(ten["upper"])
assert float(delta["lower"]) < 0 < float(delta["upper"])
table = []
benchmark_table = []
values = {}
for learner in ("estimated", "known"):
    for folds in (5, 10):
        r = comparison(learner, folds)
        values[f"{learner}_{folds}"] = r
        table.append(f"{learner.capitalize()} & {folds} & " + " & ".join(
            "$" + number(r[k]) + "$" for k in ("W_mae", "candidate_mae", "difference", "mcse")) + r" \\")
        random = comparison(learner, folds, "fixed_random")
        knndm = comparison(learner, folds, "fixed_knndm")
        benchmark_table.append(f"{learner.capitalize()} & {folds} & ${number(random['candidate_mae'])}$ & ${number(knndm['candidate_mae'])}$" + r" \\")
        values[f"{learner}_{folds}_fixed_random"] = random
        values[f"{learner}_{folds}_fixed_knndm"] = knndm
    values[f"{learner}_paired_change"] = change(learner)

template = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=24mm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern,amsmath,booktabs}
\usepackage[colorlinks=true,allcolors=blue]{hyperref}
\hypersetup{pdftitle={spLeakage: ten-fold sensitivity check},pdfauthor={Prepared for Olatunji Johnson},pdfsubject={Paired fold-count comparison; research draft}}
\title{Ten-fold sensitivity check}
\author{spLeakage: author review addendum}
\date{1 November 2026}
\begin{document}
\maketitle
\section*{What changes with ten folds?}
For prediction using estimated covariance, choosing a CV scheme by fitted $|J_2|$ gave a ten-fold MAE difference of \( @TENDELTA@ \) relative to selection by $W$ (MCSE \( @TENSE@ \)). This corresponds to a @TENPERCENT@\% reduction in MAE. The interval formed as the difference plus or minus 1.96 MCSE is \([@TENLOW@,@TENHIGH@]\), which includes zero. The primary five-fold result is a @FIVEPERCENT@\% reduction. Both comparisons use the same simulated fields.

\begin{center}\small
\begin{tabular}{llrrrr}\toprule
Covariance & Folds & $W$ MAE & Fitted $|J_2|$ MAE & Difference & MCSE\\\midrule
@TABLE@
\bottomrule\end{tabular}
\end{center}
The first column identifies the covariance used by the prediction model. In every row, selection by $|J_2|$ uses correlation estimated from the full training dataset. Both rules choose among the same four candidates at the stated fold count. MAE is the mean absolute difference between the observed CV score and the final prediction model's true conditional target risk.

To assess whether fold count changes the relative benefit, we compare the ten-fold and five-fold MAE differences within each field. With estimated-covariance prediction, the paired change is \( @CHANGE@ \) (MCSE \( @CHANGESE@ \); interval \([@CHANGELOW@,@CHANGEHIGH@]\)). This leaves the change in relative benefit unresolved. Comparing whether two separate intervals include zero would not answer that question.

For comparison, always using one fixed scheme gives:
\begin{center}\small
\begin{tabular}{llrr}\toprule
Covariance & Folds & Random-CV MAE & kNNDM MAE\\\midrule
@BENCHMARKS@
\bottomrule\end{tabular}
\end{center}
Changing fold count changes both the training subsets and their spatial arrangement. These results do not isolate a training-size effect or establish an optimal fold count. The accompanying two-page figure shows all 24 simulation settings for prediction with estimated and known covariance. The result tables retain all 13 selection rules.

\newpage
\section*{Design and interpretation}
This check was added after the primary findings and further reading of the fold-count discussion in \href{https://doi.org/10.5194/gmd-17-5897-2024}{Linnenbrink et al. (2024)}. Its protocol was recorded before its pilot and main results were generated. It supplements the primary studies without changing their settings or results.

The check reuses 1,920 datasets: 80 fields in each of 24 combinations of two covariance scales, three signal fractions, two sampling patterns and two targets. Sampling is uniform or clustered; the target lies inside the sampled domain or is shifted. Each dataset has 160 training and 225 target locations. Balanced random folds, $k$-means blocks, buffered blocks and CAST kNNDM all use ten folds. The buffer, covariance model, fitting procedure, seeds and scoring rules are unchanged. Ties within $10^{-12}$ are evaluated by averaging candidate losses, without using target error to choose a candidate.

Because the fields are reused, this is a paired sensitivity analysis rather than an independent replication. MCSEs use paired differences and variation between repetitions within each setting. They describe simulation precision for this fixed collection of settings. They do not describe uncertainty across a new ecological population or account for all decisions made after the primary results.

The 24-configuration pilot completed 984 covariance fits. The main run completed 78,720 fits with no failed configuration. Eight recorded comparisons confirm identical full-data covariance parameters and conditional target risks between fold counts. Thus the final prediction model and its target risk stay fixed while the validation procedure changes. A separate Python calculation reproduces the paired mean changes and MCSEs.

\section*{Files and reproduction}
\begingroup\raggedright
The protocol is \nolinkurl{revision/audit/fold-count-sensitivity-protocol.json}; the changes needed for the ten-fold run are recorded in \nolinkurl{revision/audit/fold-count-source-changes.patch}. Raw results, fit records, summaries and execution logs are under \nolinkurl{revision/csf/fold-count-23933572/}. That job used the isolated R library rebuilt from saved sources and spLeakage 0.1.0.9004.

To rerun the analysis, use \nolinkurl{revision/csf/fold-count-sensitivity.slurm} from a research-delivery root with a completed runtime and research archive. It creates a new directory, runs the pilot and main study, and compares them with the saved five-fold results. This report is generated from the comparison CSVs by \nolinkurl{revision/exploratory/make_fold_count_report.py}. Its adjacent JSON records the numerical values and source hashes. A separate plotting entry point regenerates the figure from the saved comparison table.

The original run and checks are documented in \nolinkurl{revision/audit/delivery-fold-count-proof.json}. Later presentation changes are recorded separately and do not alter the scientific results. Independent human scientific review remains outstanding.
\par\endgroup
\end{document}
"""

replacements = {"TENDELTA": number(ten["difference"]), "TENSE": number(ten["mcse"]),
    "TENPERCENT": percent(ten), "TENLOW": number(ten["lower"]), "TENHIGH": number(ten["upper"]),
    "FIVEPERCENT": percent(five), "CHANGE": number(delta["change"]), "CHANGESE": number(delta["mcse"]),
    "CHANGELOW": number(delta["lower"]), "CHANGEHIGH": number(delta["upper"]),
    "TABLE": "\n".join(table), "BENCHMARKS": "\n".join(benchmark_table)}
for name, value in replacements.items():
    template = template.replace("@" + name + "@", value)
assert "@" not in template
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(template)
sources = [args.results / name for name in ("fold-count-comparisons.csv", "fold-count-change-summary.csv", "fold-count-verification.csv")]
args.output.with_suffix(".json").write_text(json.dumps({"values": values,
    "source_hashes": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2) + "\n")
print("Generated the report and numeric provenance from retained comparison tables.")
