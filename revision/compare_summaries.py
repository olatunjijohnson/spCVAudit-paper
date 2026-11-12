#!/usr/bin/env python3
"""Compare regenerated scientific summaries, allowing only numerical roundoff.

Raw results are checked separately by compare_reproduction.R. CSV structure,
text, ordering and missing values must match; finite numeric values have a
declared absolute tolerance. Formatted manuscript numbers must match exactly.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def compare_csv(reference, actual, tolerance):
    with reference.open(newline="") as left, actual.open(newline="") as right:
        a, b = list(csv.reader(left)), list(csv.reader(right))
    if not a or not b or a[0] != b[0] or len(a) != len(b):
        return {"equal": False, "reason": "Different header or row count"}
    maximum = 0.0
    failures = []
    for i, (x, y) in enumerate(zip(a[1:], b[1:]), 2):
        if len(x) != len(a[0]) or len(y) != len(a[0]):
            return {"equal": False, "reason": "Malformed row width", "row": i}
        for j, (old, new) in enumerate(zip(x, y)):
            if old == new:
                continue
            try:
                u, v = float(old), float(new)
            except ValueError:
                equal = False
            else:
                difference = abs(u - v)
                equal = math.isfinite(u) and math.isfinite(v) and difference <= tolerance
                if math.isfinite(difference):
                    maximum = max(maximum, difference)
            if not equal and len(failures) < 10:
                failures.append({"row": i, "column": a[0][j], "reference": old, "actual": new})
    return {"equal": not failures, "rows": len(a)-1,
            "maximum_absolute_difference": maximum, "first_discrepancies": failures}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("actual", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance < 0:
        parser.error("Tolerance must be finite and nonnegative")
    reference, actual = args.reference.resolve(), args.actual.resolve()
    expected = sorted((reference / "revision/results/summaries").glob("*.csv"))
    if not expected:
        raise ValueError("No reference summary tables")
    expected_names = {p.name for p in expected}
    actual_names = {p.name for p in (actual / "revision/results/summaries").glob("*.csv")}
    if expected_names != actual_names:
        raise ValueError("Summary file sets differ")
    expected.append(reference / "revision/results/manuscript-number-provenance.csv")
    rows = []
    for source in expected:
        relative = source.relative_to(reference)
        result = compare_csv(source, actual / relative, args.tolerance)
        rows.append({"file": str(relative), **result})
    for relative in ("revision/manuscript/numbers.tex", "revision/manuscript/table-reference.tex"):
        rows.append({"file": relative, "equal": (reference / relative).read_bytes() == (actual / relative).read_bytes(),
                     "comparison": "exact bytes of formatted manuscript output"})
    record = {"comparison": "Finite numeric CSV values within absolute tolerance; other fields and TeX exact",
              "absolute_tolerance": args.tolerance, "reference": str(reference), "actual": str(actual),
              "all_equal": all(row["equal"] for row in rows), "files": rows}
    output = actual / "revision/audit/summary-reproduction-comparison.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n")
    if not record["all_equal"]:
        raise RuntimeError("Scientific summary discrepancies; inspect " + str(output))
    maximum = max(row.get("maximum_absolute_difference", 0) for row in rows)
    print("All", len(rows), "summary/provenance/manuscript artifacts match; maximum numeric difference", maximum)


if __name__ == "__main__":
    main()
