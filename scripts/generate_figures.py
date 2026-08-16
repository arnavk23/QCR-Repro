"""Generate the full figure set (paper protocol + head-to-head comparison).

Figures 1-6 summarize the paper-protocol sweep; figures 7-9 plot the comparison against the paper's Tables 6/7 numbers.

Usage:
    PYTHONPATH=src python scripts/generate_figures.py"""

from __future__ import annotations

import csv
import statistics
import time
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from qcr_repro.circuits import count_gates, random_circuit
from qcr_repro.qasm import parse_qasm_subset
from qcr_repro.reducer import reduce_with_lookup

ROOT = Path(__file__).resolve().parents[1]
STRICT_CSV = ROOT / "results" / "demo_sweep" / "strict" / "benchmark_reducer.csv"
LOOSE_CSV = ROOT / "results" / "demo_sweep" / "loose" / "benchmark_reducer.csv"
DEEP_CSV = ROOT / "results" / "demo_sweep" / "deep" / "benchmark_reducer.csv"
LONG_QASM = ROOT / "matlab_demo" / "QCOptimDemo" / "longcode10.txt"
SHORT_QASM = ROOT / "matlab_demo" / "QCOptimDemo" / "shortcode10.txt"
COMPARISON_ION_CSV = ROOT / "results" / "comparison" / "comparison_ion_trap.csv"
COMPARISON_NISQ_CSV = ROOT / "results" / "comparison" / "comparison_nisq.csv"
OUT_DIR = ROOT / "figures"
DATA_DIR = ROOT / "report" / "data"

# Paper's "Ours" means (Tables 6/7, 4 qubits, length 300) used in the
# comparison figures.
PAPER_OURS = {
    "ion_trap": {"total": 111.0, "two_qubit": 43.0},
    "nisq": {"total": 107.0, "two_qubit": 43.0},
}

# --------------------------------------------------------------------------- #
# shared design system
# --------------------------------------------------------------------------- #

# Muted, print-friendly palette.
PAPER_C = "#616161"      # paper reference
INPUT_C = "#90a4ae"      # input circuits
OURS_C = "#2e7d32"       # our reducers (length objective)
DEEP_C = "#1b5e20"       # our reducers, larger iteration budget / deeper database
COST_C = "#00838f"       # our reducers (cost-aware objective)
BASE_C = "#c62828"       # baseline compilers (qiskit/BQSKit)
GRID_C = "#d7d7d7"

FONT = "DejaVu Sans"


def _set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": FONT,
            "font.size": 10.5,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 11,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.5,
            "axes.edgecolor": "#444444",
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "grid.color": GRID_C,
            "grid.linestyle": "--",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.7,
            "axes.axisbelow": True,
            "figure.dpi": 200,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
        }
    )


def _finish(fig, out_name: str) -> None:
    """Apply final layout polish and save."""
    fig.tight_layout()
    fig.savefig(OUT_DIR / out_name)
    plt.close(fig)


def _bar_labels(ax, bars, fmt: str = "{:.0f}", dy: float = 0.0, fontsize: int = 9.5,
                color: str = "#333333", va: str = "bottom") -> None:
    """Annotate each bar with its value, offset above the bar top."""
    for bar in bars:
        v = bar.get_height()
        ax.annotate(fmt.format(v), (bar.get_x() + bar.get_width() / 2, v),
                    xytext=(0, 3 + dy), textcoords="offset points",
                    ha="center", va=va, fontsize=fontsize, color=color)


def load_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _csv_ok(path: Path) -> bool:
    """True if the data file exists and is non-empty (figures warn+skip otherwise)."""
    return path.exists() and path.stat().st_size > 0


# --------------------------------------------------------------------------- #
# Figures 1-6: paper-protocol runs
# --------------------------------------------------------------------------- #


def fig1_motivation(strict_rows, loose_rows, deep_rows=()):
    """Endpoints of the demo-circuit runs, with % reduction labels."""
    _, long_gates = parse_qasm_subset(LONG_QASM)
    _, short_matlab = parse_qasm_subset(SHORT_QASM)

    strict_valid = [r for r in strict_rows if r["equivalent"] == "True"]
    loose_valid = [r for r in loose_rows if r["equivalent"] == "True"]
    best_strict = min(strict_valid, key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))
    best_loose = min(loose_valid, key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))

    long_n = len(long_gates)
    labels = ["Python loose (1e-3)", "Python strict (1e-5)", "MATLAB short",
              "Long code (input)"]
    values = [int(best_loose["end_len"]), int(best_strict["end_len"]), len(short_matlab), long_n]
    colors = [COST_C, OURS_C, "#7a5c3e", INPUT_C]

    deep_valid = [r for r in deep_rows if r["equivalent"] == "True"]
    if deep_valid:
        best_deep = min(deep_valid, key=lambda r: int(r["end_len"]))
        labels.insert(2, "Python V2, deep (100k iters, depth 5)")
        values.insert(2, int(best_deep["end_len"]))
        colors.insert(2, DEEP_C)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    bars = ax.barh(labels, values, color=colors, edgecolor="#333333", linewidth=0.6, height=0.62)
    ax.invert_yaxis()
    ax.set_xlabel("Gate count")
    ax.set_title("Figure 1: Circuit-length reduction on the demo circuit")
    ax.grid(axis="x", alpha=0.4)
    ax.grid(axis="y", visible=False)

    for bar, v in zip(bars, values):
        pct = 100.0 * (long_n - v) / long_n if long_n > v else 0.0
        txt = f"{v}" + (f"   ({pct:.1f}% shorter)" if pct > 0 else "   (input)")
        ax.text(bar.get_width() + 3, bar.get_y() + bar.get_height() / 2, txt,
                va="center", ha="left", fontsize=10)

    ax.set_xlim(0, long_n * 1.22)
    _finish(fig, "figure1_motivation.png")


def fig2_database_growth():
    """Compute-graph growth with depth (values from the paper, Table 1)."""
    depth_2q = np.array([1, 2, 3, 4, 5, 6])
    nodes_2q = np.array([15, 114, 584, 2024, 4512, 7420])
    edges_2q = np.array([14, 210, 1596, 8176, 28336, 63168])

    depth_3q = np.array([1, 2, 3, 4, 5])
    nodes_3q = np.array([25, 337, 3215, 23622, 137572])
    edges_3q = np.array([24, 600, 8088, 77160, 566928])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), sharey=True)

    for ax, dep, nodes, edges, title in [
        (axes[0], depth_2q, nodes_2q, edges_2q, "2 qubits, 14 operators"),
        (axes[1], depth_3q, nodes_3q, edges_3q, "3 qubits, 24 operators"),
    ]:
        ax.plot(dep, nodes, marker="o", markersize=5, linewidth=1.8, color=OURS_C,
                label="Nodes (distinct unitaries)")
        ax.plot(dep, edges, marker="s", markersize=5, linewidth=1.8, color=BASE_C,
                label="Edges (adjacent pairs)")
        ax.set_title(title)
        ax.set_xlabel("Graph depth")
        ax.set_yscale("log")
        ax.legend(loc="upper left")
        for x, y in zip(dep, nodes):
            ax.annotate(f"{y:,}", (x, y), xytext=(0, 5), textcoords="offset points",
                        ha="center", fontsize=7.5, color="#333333")

    axes[0].set_ylabel("Count (log scale)")
    fig.suptitle("Figure 2: Compute-graph growth with depth\n(values from the paper, Table 1)",
                 fontsize=12, fontweight="bold")
    _finish(fig, "figure2_database_growth.png")


def _input_composition(gateset: str, n_circuits: int, weights) -> dict[str, float]:
    """Mean per-type gate counts of the random input circuits (same seeds as the
    comparison benchmark, so input statistics match the runs)."""
    counts: dict[str, list] = defaultdict(list)
    for seed in range(n_circuits):
        gates, _ = random_circuit(4, 300, gateset, seed=seed, weights=weights)
        c = count_gates(gates)
        for name, n in c.items():
            counts[name].append(n)
    return {name: statistics.mean(vals) for name, vals in counts.items()}


def fig3_gate_composition():
    """Input vs reduced gate composition (mean counts over the comparison runs)."""
    if not (_csv_ok(COMPARISON_ION_CSV) and _csv_ok(COMPARISON_NISQ_CSV)):
        print("WARNING: figure3 skipped (results/comparison CSVs missing)")
        return
    ion_rows = [r for r in load_rows(COMPARISON_ION_CSV) if r["end"] != "-1"]
    nisq_rows = [r for r in load_rows(COMPARISON_NISQ_CSV) if r["end"] != "-1"]
    n_ion = len({r["seed"] for r in ion_rows}) or 24
    n_nisq = len({r["seed"] for r in nisq_rows}) or 12

    def means(rows, method, types):
        rs = [r for r in rows if r["method"] == method]
        return [statistics.mean(int(r[t]) for r in rs) for t in types]

    ion_types = ["RX", "RY", "RZ", "RXX"]
    nisq_types = ["RX", "RZ", "CZ"]
    ion_in = _input_composition("ion_trap", n_ion, None)
    nisq_in = _input_composition("nisq", n_nisq, {"RX": 1.0, "RZ": 1.0, "CZ": 2.0})

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

    for ax, types, inp, rows, reduced_method, title in [
        (axes[0], ion_types, ion_in, ion_rows, "exact_len", "Ion trap (RX/RY/RZ/RXX)"),
        (axes[1], nisq_types, nisq_in, nisq_rows, "numeric_len", "NISQ (RX/RZ/CZ)"),
    ]:
        x = np.arange(len(types))
        w = 0.34
        in_vals = [inp[t] for t in types]
        red_vals = means(rows, reduced_method, types)
        b1 = ax.bar(x - w / 2, in_vals, w, label="Input", color=INPUT_C,
                    edgecolor="#333333", linewidth=0.5)
        b2 = ax.bar(x + w / 2, red_vals, w, label=f"Reduced ({reduced_method})",
                    color=OURS_C, edgecolor="#333333", linewidth=0.5)
        _bar_labels(ax, b1, fmt="{:.0f}", fontsize=8.5)
        _bar_labels(ax, b2, fmt="{:.0f}", fontsize=8.5)
        ax.set_xticks(x, types)
        ax.set_ylabel("Mean gate count (per circuit)")
        ax.set_title(title)
        ax.legend()
        # percentage reduction annotation per type (above the taller bar's label)
        for xi, (iv, rv) in zip(x, zip(in_vals, red_vals)):
            if iv <= 0:
                continue
            ax.annotate(f"{-100.0 * (rv - iv) / iv:.0f}%", (xi, max(iv, rv)),
                        xytext=(0, 16), textcoords="offset points", ha="center",
                        fontsize=8, color="#555555")

    fig.suptitle("Figure 3: Gate composition of input vs reduced circuits\n"
                 "(mean over the comparison benchmark runs)", fontsize=12, fontweight="bold")
    _finish(fig, "figure3_pipeline.png")


def fig4_reduction_curve():
    """Convergence of end length with the iteration budget (demo circuit)."""
    curve_csv = DATA_DIR / "figure4_curve.csv"
    if curve_csv.exists():
        with curve_csv.open("r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))
        checkpoints = [int(r["iterations"]) for r in rows]
        end_lengths = [int(r["end_len"]) for r in rows]
    else:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
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

        with curve_csv.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["iterations", "runtime_sec", "end_len"])
            for i, r, l in zip(checkpoints, runtimes, end_lengths):
                w.writerow([i, round(r, 4), l])

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.plot(checkpoints, end_lengths, marker="o", markersize=5, linewidth=2,
            color=OURS_C)
    ax.fill_between(checkpoints, end_lengths, end_lengths[0], color=OURS_C,
                    alpha=0.08)
    ax.set_xlabel("Iteration budget")
    ax.set_ylabel("Circuit length")
    ax.set_title("Figure 4: Reduction vs iteration budget (demo circuit)")
    ax.set_ylim(0, max(end_lengths) * 1.08)
    ax.annotate(f"input: {end_lengths[0]}", (checkpoints[0], end_lengths[0]),
                xytext=(8, -16), textcoords="offset points", fontsize=9.5)
    ax.annotate(f"reduced: {end_lengths[-1]}", (checkpoints[-1], end_lengths[-1]),
                xytext=(-30, 10), textcoords="offset points", fontsize=9.5,
                color=OURS_C, fontweight="bold")
    _finish(fig, "figure4_reduction_curve.png")


def fig5_runtime_vs_length(strict_rows, loose_rows):
    """Runtime vs achieved end length; color = tolerance, marker = depth."""
    fig, ax = plt.subplots(figsize=(8.2, 4.8))

    style = [
        (strict_rows, "strict (1e-5)", OURS_C),
        (loose_rows, "loose (1e-3)", BASE_C),
    ]
    depth_markers = {3: "o", 4: "s"}

    for rows, label, color in style:
        for depth in sorted(depth_markers):
            sub = [r for r in rows if int(r["depth"]) == depth]
            if not sub:
                continue
            x = [float(r["runtime_sec"]) for r in sub]
            y = [int(r["end_len"]) for r in sub]
            ax.scatter(x, y, marker=depth_markers[depth], s=46, alpha=0.85,
                       color=color, edgecolor="#333333", linewidth=0.4,
                       label=f"{label}, depth {depth}")

    ax.set_xlabel("Runtime (s)")
    ax.set_ylabel("End gate count")
    ax.set_title("Figure 5: Runtime vs reduced length (color = tolerance)")
    ax.legend(loc="upper right", framealpha=0.95)
    _finish(fig, "figure5_runtime_vs_length.png")


def fig6_boxplot(strict_rows):
    """End-length distribution per (depth, iterations) group, strict setting."""
    groups = defaultdict(list)
    for r in strict_rows:
        groups[(int(r["depth"]), int(r["iterations"]))].append(int(r["end_len"]))

    labels = []
    data = []
    for key in sorted(groups):
        labels.append(f"depth {key[0]}\niters {key[1]}")
        data.append(groups[key])

    fig, ax = plt.subplots(figsize=(10, 4.6))
    bp = ax.boxplot(data, tick_labels=labels, showfliers=True, patch_artist=True,
                    medianprops=dict(color="#c62828", linewidth=1.6))
    for patch in bp["boxes"]:
        patch.set(facecolor="#dcebe2", edgecolor=OURS_C, linewidth=1.1)
    for i, d in enumerate(data, start=1):
        ax.annotate(f"n={len(d)}", (i, max(d)), xytext=(0, 8),
                    textcoords="offset points", ha="center", fontsize=8.5,
                    color="#555555")
    ax.set_ylabel("End gate count")
    ax.set_title("Figure 6: End-length distribution across seeds (strict setting)")
    _finish(fig, "figure6_boxplot.png")


# --------------------------------------------------------------------------- #
# Figures 7-9: head-to-head comparison with the paper
# --------------------------------------------------------------------------- #


def _comparison_stats(csv_path: Path) -> dict[str, dict]:
    """Per-method (mean, std) of total and two-qubit counts from a comparison CSV."""
    rows = [r for r in load_rows(csv_path) if r["end"] != "-1" and int(r["end"]) >= 0]
    stats: dict[str, dict] = {}
    for method in sorted({r["method"] for r in rows}):
        rs = [r for r in rows if r["method"] == method]
        ends = [int(r["end"]) for r in rs]
        twqs = [int(r["twq"]) for r in rs]
        stats[method] = {
            "n": len(rs),
            "total": statistics.mean(ends),
            "total_std": statistics.pstdev(ends),
            "twq": statistics.mean(twqs),
            "twq_std": statistics.pstdev(twqs),
        }
    return stats


def _comparison_bar(fig_title, ylabel, labels, methods, stats, paper_total, out_name):
    """Shared comparison bar chart: paper reference + our methods + baselines."""
    values = [paper_total]
    errs = [0.0]
    colors = [PAPER_C]
    for m in methods:
        st = stats[m]
        values.append(st["total"])
        errs.append(st["total_std"])
        colors.append(OURS_C if m in ("exact_len", "numeric_len")
                      else COST_C if m in ("exact_cost", "numeric_cost")
                      else BASE_C)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    bars = ax.bar(range(len(labels)), values, yerr=errs, capsize=3.5,
                  color=colors, edgecolor="#333333", linewidth=0.6, width=0.66)
    ax.axhline(paper_total, color=PAPER_C, ls="--", lw=1.1, alpha=0.9)
    ax.text(len(labels) - 0.45, paper_total, f" paper mean {paper_total:.0f}",
            va="center", fontsize=8.5, color=PAPER_C)
    _bar_labels(ax, bars, fmt="{:.0f}", fontsize=9.5)

    ax.set_xticks(range(len(labels)), labels, fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title(fig_title)
    ax.set_ylim(0, max(values) * 1.22)
    _finish(fig, out_name)


def fig7_comparison_ion():
    """Total gate counts, ion-trap pool: paper vs our reducers vs qiskit."""
    if not _csv_ok(COMPARISON_ION_CSV):
        print(f"WARNING: figure7 skipped ({COMPARISON_ION_CSV} missing)")
        return
    stats = _comparison_stats(COMPARISON_ION_CSV)
    labels = ["Paper\n'Ours'", "exact_len\n(ours)", "exact_cost\n(ours)",
              "qiskit\nL1", "qiskit\nL2", "qiskit\nL3"]
    methods = ["exact_len", "exact_cost", "qiskit_l1", "qiskit_l2", "qiskit_l3"]
    missing = [m for m in methods if m not in stats]
    if missing:
        print(f"WARNING: figure7 skipped, missing methods in {COMPARISON_ION_CSV}: {missing}")
        return
    _comparison_bar(
        "Figure 7: Ion trap (4 qubits) - comparison with the paper",
        "Mean final gate count", labels, methods, stats,
        PAPER_OURS["ion_trap"]["total"], "figure7_comparison_ion_trap.png",
    )


def fig8_comparison_nisq():
    """Total gate counts, NISQ pool: paper vs our reducers vs qiskit."""
    if not _csv_ok(COMPARISON_NISQ_CSV):
        print(f"WARNING: figure8 skipped ({COMPARISON_NISQ_CSV} missing)")
        return
    stats = _comparison_stats(COMPARISON_NISQ_CSV)
    labels = ["Paper\n'Ours'", "numeric_len\n(ours)", "numeric_cost\n(ours)",
              "qiskit\nL1", "qiskit\nL2", "qiskit\nL3"]
    methods = ["numeric_len", "numeric_cost", "qiskit_l1", "qiskit_l2", "qiskit_l3"]
    missing = [m for m in methods if m not in stats]
    if missing:
        print(f"WARNING: figure8 skipped, missing methods in {COMPARISON_NISQ_CSV}: {missing}")
        return
    _comparison_bar(
        "Figure 8: NISQ (4 qubits) - comparison with the paper",
        "Mean final gate count", labels, methods, stats,
        PAPER_OURS["nisq"]["total"], "figure8_comparison_nisq.png",
    )


def fig9_two_qubit_counts():
    """Two-qubit gate counts (RXX/CZ): hardware-cost objective vs paper."""
    if not (_csv_ok(COMPARISON_ION_CSV) and _csv_ok(COMPARISON_NISQ_CSV)):
        print("WARNING: figure9 skipped (results/comparison CSVs missing)")
        return
    ion = _comparison_stats(COMPARISON_ION_CSV)
    nisq = _comparison_stats(COMPARISON_NISQ_CSV)

    ion_labels = ["Paper", "exact_len", "exact_cost"]
    nisq_labels = ["Paper", "numeric_len", "numeric_cost"]
    missing = [m for m in ["exact_len", "exact_cost"] if m not in ion] + \
              [m for m in ["numeric_len", "numeric_cost"] if m not in nisq]
    if missing:
        print(f"WARNING: figure9 skipped, missing methods: {missing}")
        return

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    panels = [
        (axes[0], ion_labels, ["exact_len", "exact_cost"], ion,
         PAPER_OURS["ion_trap"]["two_qubit"], "Ion trap (RXX)", "Mean RXX count"),
        (axes[1], nisq_labels, ["numeric_len", "numeric_cost"], nisq,
         PAPER_OURS["nisq"]["two_qubit"], "NISQ (CZ)", "Mean CZ count"),
    ]
    for ax, labels, methods, stats, paper_twq, title, ylabel in panels:
        values = [paper_twq]
        errs = [0.0]
        colors = [PAPER_C]
        for m in methods:
            values.append(stats[m]["twq"])
            errs.append(stats[m]["twq_std"])
            colors.append(OURS_C if m.endswith("_len") else COST_C)
        bars = ax.bar(labels, values, yerr=errs, capsize=3.5, color=colors,
                      edgecolor="#333333", linewidth=0.6, width=0.62)
        _bar_labels(ax, bars, fmt="{:.1f}", fontsize=9.5)
        ax.axhline(paper_twq, color=PAPER_C, ls="--", lw=1.0, alpha=0.8)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, max(values) * 1.2)

    fig.suptitle("Figure 9: Two-qubit gate counts - hardware-cost objective",
                 fontsize=12, fontweight="bold")
    _finish(fig, "figure9_two_qubit_counts.png")


def main():
    _set_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    strict_rows = load_rows(STRICT_CSV)
    loose_rows = load_rows(LOOSE_CSV)
    deep_rows = load_rows(DEEP_CSV) if _csv_ok(DEEP_CSV) else []

    fig1_motivation(strict_rows, loose_rows, deep_rows)
    fig2_database_growth()
    fig3_gate_composition()
    fig4_reduction_curve()
    fig5_runtime_vs_length(strict_rows, loose_rows)
    fig6_boxplot(strict_rows)
    fig7_comparison_ion()
    fig8_comparison_nisq()
    fig9_two_qubit_counts()

    print(f"Saved figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
