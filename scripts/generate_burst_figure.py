"""Figure for the matlab_demo-faithful burst reducer finding: NISQ total
gate count, current default (burst_frac=1, i.e. reduce_matlab_burst for the
whole budget then one structural+sweep cleanup pass) vs the previous
default (exhaustive sweep, burst_frac=0), full official protocol (4 qubits,
length 300, 60s budget, max_block_len=10, rz_pass + hybrid on), n=100,
identical seeds and identical hybrid setting for both -- NEW_CSV is the
current results/comparison/comparison_nisq.csv (burst_frac=1, the new
default); OLD_CSV is the pre-burst dataset pulled from git history
(git show HEAD:results/comparison/comparison_nisq.csv, before this pass's
commit), i.e. the exact same 100 seeds and hybrid=True setting with
burst_frac=0 (the burst block is skipped entirely when burst_frac=0, so
this is byte-identical to the code path that produced the paper's
originally-reported 160.9/160.5).

Usage:
    PYTHONPATH=src python scripts/generate_burst_figure.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
NEW_CSV = ROOT / "results" / "comparison" / "comparison_nisq.csv"
OLD_CSV = Path(
    r"C:\Users\kapoo\AppData\Local\Temp\claude\C--Users-kapoo-Downloads-QCR-Repro"
    r"\717f9095-d1f3-46ac-b818-eda5fcf1f0f5\scratchpad\nisq_research\old_comparison_nisq.csv"
)

SWEEP_C = "#616161"   # current default (pure exhaustive sweep, burst_frac=0)
BURST_C = "#2e7d32"   # matlab-faithful burst (burst_frac=1)
PAPER_C = "#c62828"   # published baseline reference line

PUBLISHED_BASELINE = 107


def _load(path: Path, method: str) -> np.ndarray:
    rows = list(csv.DictReader(path.open()))
    by_seed = {int(r["seed"]): int(r["end"]) for r in rows if r["method"] == method}
    return by_seed


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for label, method in [("length-minimizing", "numeric_len"), ("cost-aware", "numeric_cost")]:
        sweep_by_seed = _load(OLD_CSV, method)
        burst_by_seed = _load(NEW_CSV, method)
        seeds = sorted(set(sweep_by_seed) & set(burst_by_seed))
        sweep = np.array([sweep_by_seed[s] for s in seeds], dtype=float)
        burst = np.array([burst_by_seed[s] for s in seeds], dtype=float)

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

        ax = axes[0]
        bp = ax.boxplot(
            [sweep, burst],
            tick_labels=["exhaustive sweep\n(current default)", "matlab-faithful burst\n(burst_frac=1.0)"],
            patch_artist=True, widths=0.55, showmeans=True,
        )
        for patch, color in zip(bp["boxes"], [SWEEP_C, BURST_C]):
            patch.set_facecolor(color)
            patch.set_alpha(0.25)
            patch.set_edgecolor(color)
        for med, color in zip(bp["medians"], [SWEEP_C, BURST_C]):
            med.set_color(color)
            med.set_linewidth(2)
        ax.axhline(PUBLISHED_BASELINE, color=PAPER_C, linestyle="--", linewidth=1.2,
                   label=f"published baseline ({PUBLISHED_BASELINE})")
        ax.set_ylabel("NISQ total gate count")
        ax.set_title(f"Per-circuit distribution (n={len(seeds)}, same seeds)")
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(axis="y", alpha=0.3)

        ax = axes[1]
        rng = np.random.default_rng(0)
        jitter = rng.uniform(-0.06, 0.06, size=len(sweep))
        ax.scatter(np.zeros_like(sweep) + jitter, sweep, s=14, color=SWEEP_C, alpha=0.6)
        ax.scatter(np.ones_like(burst) + jitter, burst, s=14, color=BURST_C, alpha=0.6)
        for i in range(len(sweep)):
            color = BURST_C if burst[i] < sweep[i] else (SWEEP_C if burst[i] > sweep[i] else "#bdbdbd")
            ax.plot([jitter[i], 1 + jitter[i]], [sweep[i], burst[i]], color=color, alpha=0.25, linewidth=0.9)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["sweep", "burst"])
        ax.set_xlim(-0.3, 1.3)
        ax.set_ylabel("NISQ total gate count")
        n_better = int(np.sum(burst < sweep))
        n_worse = int(np.sum(burst > sweep))
        ax.set_title(f"Paired per-circuit: burst wins {n_better}/{len(sweep)}, loses {n_worse}/{len(sweep)}")
        ax.grid(axis="y", alpha=0.3)

        delta = burst.mean() - sweep.mean()
        fig.suptitle(
            f"NISQ ({label}): matlab_demo-faithful burst vs. exhaustive sweep, full official protocol, n={len(seeds)}\n"
            f"mean {sweep.mean():.1f} → {burst.mean():.1f} gates ({delta:+.1f}, {100*delta/sweep.mean():+.1f}%)",
            fontsize=12, fontweight="bold",
        )
        fig.tight_layout(rect=(0, 0, 1, 0.90))
        suffix = method.replace("numeric_", "")
        out = OUT_DIR / f"figure_matlab_burst_nisq_{suffix}.png"
        fig.savefig(out, dpi=150)
        print(f"saved {out}  (n={len(seeds)}, sweep={sweep.mean():.1f}, burst={burst.mean():.1f}, delta={delta:+.1f})")


if __name__ == "__main__":
    main()
