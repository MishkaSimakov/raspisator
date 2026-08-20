#!/usr/bin/env python3
"""Classify simplex benchmark results and render a Markdown summary.

Usage: benchmark_report.py [variant]   # variant: primal (default) | dual

Joins the benchmark output (`log/benchmark_simplex_<variant>.csv`) with the
committed reference optima (`.github/data/lp_optimal_values.csv`) and writes a
results table to `$GITHUB_STEP_SUMMARY` (and stdout for local runs).

Classification per problem:
  SOLVED        status OPTIMAL and |obj - ref| / max(1, |ref|) < TOL
  WRONG         status OPTIMAL but objective disagrees with the reference
  NO_REFERENCE  no reference optimum available (e.g. STANDGUB)
  SKIPPED       status SKIPPED -- problem deliberately excluded from the run
                (does not count as a failure)
  <status>      any other non-OPTIMAL status, passed through as a failure bucket
                (INFEASIBLE / UNBOUNDED / ITERATIONS_LIMIT / PHASE1_ERROR /
                 EXCEPTION)

The benchmark reports phase 1 and phase 2 iterations separately; both are shown
per problem.

Timing is reported but never used to gate or classify -- shared CI runners are
too noisy. Iteration count is the deterministic performance proxy.
"""

import csv
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_CSV = REPO_ROOT / ".github" / "data" / "lp_optimal_values.csv"
VARIANTS = ("primal", "dual")

# Relative tolerance ~ 7 significant digits (see plan / user decision).
TOL = 1e-7

EMOJI = {
    "SOLVED": "✅",
    "WRONG": "❌",
    "NO_REFERENCE": "➖",
    "SKIPPED": "⏭️",
}
FAIL_EMOJI = "❌"

# Categories that are not solver failures (excluded from the "Failed" bucket).
NON_FAILURE = ("SOLVED", "WRONG", "NO_REFERENCE", "SKIPPED")


def load_references(path):
    refs = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            refs[row["name"]] = float(row["optimal"])
    return refs


def rel_error(obj, ref):
    return abs(obj - ref) / max(1.0, abs(ref))


def classify(status, objective, ref):
    """Return (category, rel_err_or_None)."""
    if status != "OPTIMAL":
        # SKIPPED and the failure statuses pass through as their own category.
        return status, None
    if ref is None:
        return "NO_REFERENCE", None
    err = rel_error(objective, ref)
    return ("SOLVED" if err < TOL else "WRONG"), err


def fmt_num(x):
    return f"{x:.10g}"


def fmt_err(err):
    return "—" if err is None else f"{err:.2e}"


def main():
    variant = sys.argv[1] if len(sys.argv) > 1 else "primal"
    if variant not in VARIANTS:
        sys.exit(f"Unknown variant {variant!r}, expected one of {VARIANTS}.")

    benchmark_csv = REPO_ROOT / "log" / f"benchmark_simplex_{variant}.csv"
    if not benchmark_csv.exists():
        sys.exit(f"Benchmark output not found: {benchmark_csv}")

    refs = load_references(REFERENCE_CSV)

    rows = []
    with benchmark_csv.open(newline="") as f:
        for r in csv.DictReader(f):
            name = r["name"]
            status = r["status"]
            objective = float(r["objective"])
            ref = refs.get(name)
            category, err = classify(status, objective, ref)
            phase1_iters = int(r["phase1_iterations"])
            phase2_iters = int(r["phase2_iterations"])
            rows.append({
                "name": name,
                "status": status,
                "category": category,
                "objective": objective,
                "ref": ref,
                "err": err,
                "phase1_iters": phase1_iters,
                "phase2_iters": phase2_iters,
                "total_iters": phase1_iters + phase2_iters,
                "time_s": int(r["time"]) / 1e9,
            })

    rows.sort(key=lambda r: r["name"])

    total = len(rows)
    solved = sum(1 for r in rows if r["category"] == "SOLVED")
    wrong = sum(1 for r in rows if r["category"] == "WRONG")
    no_ref = sum(1 for r in rows if r["category"] == "NO_REFERENCE")
    skipped = sum(1 for r in rows if r["category"] == "SKIPPED")
    failed = total - solved - wrong - no_ref - skipped

    # Skipped problems were never attempted, so report against the run total.
    attempted = total - skipped

    # Per-status breakdown of the failure buckets, for the summary line.
    fail_statuses = {}
    for r in rows:
        if r["category"] not in NON_FAILURE:
            fail_statuses[r["category"]] = fail_statuses.get(r["category"], 0) + 1
    fail_detail = ", ".join(f"{k} {v}" for k, v in sorted(fail_statuses.items()))

    out = []
    out.append(f"## Simplex {variant} benchmark\n")
    out.append(
        f"**Solved {solved} / {attempted}** · "
        f"Wrong {wrong} · Failed {failed} · No-ref {no_ref} · Skipped {skipped}\n"
    )
    if fail_detail:
        out.append(f"<sub>Failures: {fail_detail}</sub>\n")
    out.append("")
    out.append("| | Problem | Result | Status | Objective | Reference | Rel. err | Iters P1 | Iters P2 | Iters Σ | Time (s) |")
    out.append("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        mark = EMOJI.get(r["category"], FAIL_EMOJI)
        ref_str = fmt_num(r["ref"]) if r["ref"] is not None else "—"
        obj_str = fmt_num(r["objective"]) if r["status"] == "OPTIMAL" else "—"
        if r["category"] == "SKIPPED":
            p1_str = p2_str = sum_str = time_str = "—"
        else:
            p1_str = str(r["phase1_iters"])
            p2_str = str(r["phase2_iters"])
            sum_str = str(r["total_iters"])
            time_str = f"{r['time_s']:.3f}"
        out.append(
            f"| {mark} | {r['name']} | {r['category']} | {r['status']} | "
            f"{obj_str} | {ref_str} | {fmt_err(r['err'])} | "
            f"{p1_str} | {p2_str} | {sum_str} | {time_str} |"
        )
    report = "\n".join(out) + "\n"

    print(report)
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        with open(summary_path, "a") as f:
            f.write(report)


if __name__ == "__main__":
    main()
