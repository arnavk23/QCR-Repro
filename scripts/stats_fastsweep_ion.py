"""One-sample t-test of the fast-sweep exact-engine results against the
published baseline (111 total, 43 RXX), matching report/draft_paper.tex's
existing reporting format (t, p, Cohen's d, 95% CI) for the cost-aware
objective's total-length and two-qubit-count comparisons, plus the
length-minimizing objective for the table.

Usage:
    python scripts/stats_fastsweep_ion.py results/comparison_exactfast_fast/comparison_ion_trap.csv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PAPER_TOTAL = 111
PAPER_RXX = 43


def load(csv_path: Path, method: str) -> tuple[np.ndarray, np.ndarray, dict]:
    ends, rxx = [], []
    per_type = {"RX": [], "RY": [], "RZ": [], "RXX": []}
    with csv_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["method"] != method:
                continue
            ends.append(int(row["end"]))
            rxx.append(int(row["RXX"]))
            for k in per_type:
                per_type[k].append(int(row[k]))
    arrs = {k: np.array(v, dtype=float) for k, v in per_type.items()}
    return np.array(ends, dtype=float), np.array(rxx, dtype=float), arrs


def report(name: str, values: np.ndarray, paper_value: float) -> None:
    n = len(values)
    mean = values.mean()
    sd = values.std(ddof=1)
    t, p = stats.ttest_1samp(values, paper_value)
    d = (mean - paper_value) / sd
    se = sd / np.sqrt(n)
    ci_lo, ci_hi = stats.t.interval(0.95, df=n - 1, loc=mean, scale=se)
    pct = (paper_value - mean) / paper_value * 100
    print(f"  {name}: mean {mean:.1f} (sd {sd:.1f}), n={n}, vs paper {paper_value:.0f} "
          f"({pct:+.1f}%)")
    print(f"    t={t:.1f}  p={p:.2e}  Cohen's d={d:.2f}  95% CI [{ci_lo:.1f}, {ci_hi:.1f}]")


def main() -> None:
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/comparison_exactfast_fast/comparison_ion_trap.csv")
    for method in ("exact_len", "exact_cost"):
        ends, rxx, counts_mean = load(csv_path, method)
        if len(ends) == 0:
            print(f"-- {method}: no rows --")
            continue
        print(f"-- {method} (n={len(ends)}) --")
        per_type_str = "  ".join(
            f"{k} {counts_mean[k].mean():.1f} ({counts_mean[k].std(ddof=1):.1f})" for k in ("RX", "RY", "RZ", "RXX")
        )
        print(f"  {per_type_str}  total {ends.mean():.1f} ({ends.std(ddof=1):.1f})")
        report("total length", ends, PAPER_TOTAL)
        report("RXX (two-qubit)", rxx, PAPER_RXX)
        print()


if __name__ == "__main__":
    main()
