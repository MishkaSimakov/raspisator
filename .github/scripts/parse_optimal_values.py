#!/usr/bin/env python3
"""Generate lp_optimal_values.csv from the Netlib LP `00readme.txt`.

This is a one-shot generator, run locally and checked into the repo (the
dataset and its readme are git-ignored, so we commit the derived reference
table instead). The CI report does NOT call this script -- it reads the
committed CSV directly.

Two sources of optimal values inside the readme:

  1. The PROBLEM SUMMARY TABLE -- the baseline `Optimal Value` column.
  2. A later correction table listing CPLEX(Sparc) / MINOS(MIPS) values that
     differ from the baseline. Per the user's choice, when a CPLEX value is
     present we prefer it as the reference (`source=cplex`).

Usage:
    parse_optimal_values.py [README_PATH] [OUTPUT_CSV]
"""

import csv
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_README = REPO_ROOT / "resources" / "lp_problems" / "00readme.txt"
DEFAULT_OUTPUT = REPO_ROOT / ".github" / "data" / "lp_optimal_values.csv"

# Scientific-notation floats only (e.g. -4.6475314286E+02). Plain integers such
# as the row/column counts have no exponent and are therefore never matched.
SCI = re.compile(r"[-+]?\d*\.?\d+[eE][-+]?\d+")

# In the CPLEX/MINOS correction table the CPLEX column begins around character
# 13 and the MINOS column around character 34. A value whose match starts left
# of this threshold belongs to the CPLEX column; one starting to its right is
# MINOS-only (and is ignored, keeping the summary-table value).
CPLEX_COLUMN_MAX_START = 30


def parse_summary(lines):
    """Baseline optimal values from the PROBLEM SUMMARY TABLE."""
    values = {}
    in_table = False
    for line in lines:
        # Anchor on the column-header row, not the prose mentions of the table
        # name that appear earlier in the file.
        if not in_table:
            if "Optimal Value" in line:
                in_table = True
            continue
        if "BOUND-TYPE TABLE" in line:
            break

        tokens = line.split()
        if len(tokens) < 2:
            continue
        name = tokens[0]
        # Skip the column header row ("Name Rows Cols ...").
        if name.lower() == "name":
            continue
        matches = SCI.findall(line)
        if not matches:
            continue  # e.g. STANDGUB -> "(see NOTES)", no numeric optimum.
        values[name] = matches[-1]
    return values


def apply_cplex_overrides(lines, values):
    """Override with CPLEX(Sparc) values from the correction table."""
    in_table = False
    for line in lines:
        if "CPLEX(Sparc)" in line and "MINOS(MIPS)" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if "above CPLEX and MINOS" in line:
            break

        tokens = line.split()
        if not tokens:
            continue
        name = tokens[0]
        if name not in values:
            continue
        # First scientific number starting within the CPLEX column wins.
        for m in SCI.finditer(line):
            if m.start() <= CPLEX_COLUMN_MAX_START:
                values[name] = (m.group(), "cplex")
                break
    return values


def main():
    readme = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_README
    output = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUTPUT

    lines = readme.read_text().splitlines()

    # name -> optimal string; mark source as we go.
    summary = parse_summary(lines)
    values = {name: (opt, "summary") for name, opt in summary.items()}
    apply_cplex_overrides(lines, values)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "optimal", "source"])
        for name in sorted(values):
            opt, source = values[name]
            writer.writerow([name, opt, source])

    n_cplex = sum(1 for v in values.values() if v[1] == "cplex")
    print(f"Wrote {len(values)} reference values "
          f"({n_cplex} from CPLEX column) to {output}")


if __name__ == "__main__":
    main()
