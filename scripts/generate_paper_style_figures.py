from __future__ import annotations

import csv
import time
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from qcr_repro.qasm_io import parse_qasm_subset
from qcr_repro.reducer import reduce_with_lookup

ROOT = Path(__file__).resolve().parents[1]
STRICT_CSV = ROOT / "results" / "benchmark_reducer.csv"
LOOSE_CSV = ROOT / "results_tol1e3" / "benchmark_reducer.csv"
LONG_QASM = ROOT / "paper_code" / "QCOptimDemo" / "longcode10.txt"
SHORT_QASM = ROOT / "paper_code" / "QCOptimDemo" / "shortcode10.txt"
OUT_DIR = ROOT / "implementation" / "figures"
DATA_DIR = ROOT / "implementation" / "data"


def load_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def fig1_motivation(strict_rows, loose_rows):
    _, long_gates = parse_qasm_subset(LONG_QASM)
    _, short_matlab = parse_qasm_subset(SHORT_QASM)

    strict_valid = [r for r in strict_rows if r["equivalent"] == "True"]
    loose_valid = [r for r in loose_rows if r["equivalent"] == "True"]
    best_strict = min(strict_valid, key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))
    best_loose = min(loose_valid, key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))

    labels = ["Long code", "MATLAB short", "Python strict", "Python loose"]
    values = [len(long_gates), len(short_matlab), int(best_strict["end_len"]), int(best_loose["end_len"])]

    plt.figure(figsize=(8, 4.5))
    bars = plt.bar(labels, values)
    plt.ylabel("Gate count")
    plt.title("Figure 1 (Implementation): Circuit-length reduction motivation")
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, value + 1, str(value), ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure1_motivation.png", dpi=180)
    plt.close()


def fig2_compute_graph_growth():
    depth_2q = np.array([1, 2, 3, 4, 5, 6])
    nodes_2q = np.array([15, 114, 584, 2024, 4512, 7420])
    edges_2q = np.array([14, 210, 1596, 8176, 28336, 63168])

    depth_3q = np.array([1, 2, 3, 4, 5])
    nodes_3q = np.array([25, 337, 3215, 23622, 137572])
    edges_3q = np.array([24, 600, 8088, 77160, 566928])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

    axes[0].plot(depth_2q, nodes_2q, marker="o", label="Nodes")
    axes[0].plot(depth_2q, edges_2q, marker="s", label="Edges")
    axes[0].set_title("2 qubits, 14 operators")
    axes[0].set_xlabel("Depth")
    axes[0].set_ylabel("Count")
    axes[0].set_yscale("log")
    axes[0].legend()

    axes[1].plot(depth_3q, nodes_3q, marker="o", label="Nodes")
    axes[1].plot(depth_3q, edges_3q, marker="s", label="Edges")
    axes[1].set_title("3 qubits, 24 operators")
    axes[1].set_xlabel("Depth")
    axes[1].set_yscale("log")
    axes[1].legend()

    fig.suptitle("Figure 2 (Implementation): Compute-graph growth (from paper table)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "figure2_compute_graph_growth.png", dpi=180)
    plt.close(fig)


def fig3_method_pipeline():
    fig, ax = plt.subplots(figsize=(10, 3.6))
    ax.axis("off")

    boxes = [
        (0.03, 0.35, 0.2, 0.3, "Input token chain"),
        (0.28, 0.35, 0.2, 0.3, "Random block sampling"),
        (0.53, 0.35, 0.2, 0.3, "DB lookup\n(+ optional RF)"),
        (0.78, 0.35, 0.2, 0.3, "Replace if shorter")
    ]

    for x, y, w, h, text in boxes:
        rect = plt.Rectangle((x, y), w, h, fill=False)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=10)

    for x0, x1 in [(0.23, 0.28), (0.48, 0.53), (0.73, 0.78)]:
        ax.annotate("", xy=(x1, 0.5), xytext=(x0, 0.5), arrowprops=dict(arrowstyle="->"))

    ax.text(0.54, 0.2, "V1: Random search | V2: DB retrieval | V3: RF + DB", ha="center", fontsize=10)
    plt.title("Figure 3 (Implementation): Reduction pipeline")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure3_pipeline.png", dpi=180)
    plt.close(fig)


def fig4_reduction_curve():
    _, gates = parse_qasm_subset(LONG_QASM)
    checkpoints = [0, 200, 400, 600, 800, 1000, 1200, 1500]
    end_lengths = []
    runtimes = []

    for iters in checkpoints:
        start = time.time()
        if iters == 0:
            reduced = gates
        else:
            reduced, _ = reduce_with_lookup(
                gates,
                num_qubits=5,
                local_qubits=3,
                max_block_len=7,
                graph_depth=4,
                iterations=iters,
                seed=10,
            )
        runtimes.append(time.time() - start)
        end_lengths.append(len(reduced))

    with (DATA_DIR / "figure4_curve.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["iterations", "runtime_sec", "end_len"])
        for i, r, l in zip(checkpoints, runtimes, end_lengths):
            w.writerow([i, round(r, 4), l])

    plt.figure(figsize=(7.8, 4.5))
    plt.plot(checkpoints, end_lengths, marker="o")
    plt.xlabel("Iterations")
    plt.ylabel("Circuit length")
    plt.title("Figure 4 (Implementation): Reduction vs iteration budget")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure4_reduction_curve.png", dpi=180)
    plt.close()


def fig5_runtime_vs_length(strict_rows, loose_rows):
    plt.figure(figsize=(8, 4.8))

    for rows, label, marker in [(strict_rows, "strict (1e-5)", "o"), (loose_rows, "loose (1e-3)", "s")]:
        x = [float(r["runtime_sec"]) for r in rows]
        y = [int(r["end_len"]) for r in rows]
        c = [int(r["depth"]) for r in rows]
        sc = plt.scatter(x, y, c=c, cmap="viridis", marker=marker, alpha=0.8, label=label)

    plt.colorbar(sc, label="Graph depth")
    plt.xlabel("Runtime (s)")
    plt.ylabel("End gate count")
    plt.title("Figure 5 (Implementation): Runtime vs reduced length")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure5_runtime_vs_length.png", dpi=180)
    plt.close()


def fig6_boxplot(strict_rows):
    groups = defaultdict(list)
    for r in strict_rows:
        groups[(int(r["depth"]), int(r["iterations"]))].append(int(r["end_len"]))

    labels = []
    data = []
    for key in sorted(groups):
        labels.append(f"d{key[0]}-i{key[1]}")
        data.append(groups[key])

    plt.figure(figsize=(10, 4.6))
    plt.boxplot(data, tick_labels=labels, showfliers=True)
    plt.ylabel("End gate count")
    plt.title("Figure 6 (Implementation): Distribution across seeds (strict)")
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure6_boxplot.png", dpi=180)
    plt.close()


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    strict_rows = load_rows(STRICT_CSV)
    loose_rows = load_rows(LOOSE_CSV)

    fig1_motivation(strict_rows, loose_rows)
    fig2_compute_graph_growth()
    fig3_method_pipeline()
    fig4_reduction_curve()
    fig5_runtime_vs_length(strict_rows, loose_rows)
    fig6_boxplot(strict_rows)

    print(f"Saved figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
