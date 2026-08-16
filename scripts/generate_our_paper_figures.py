"""Our-data analogs of the paper's Figures 1-10 (see report/results_report.md
Appendix B for the originals). Every figure here is computed from this
repository's own runs, data and circuits -- no numbers are copied from the
paper. Figure 10 (IBM hardware measurement) has no analog: this project has
no quantum-hardware access, and is not fabricated.

Usage:
    PYTHONPATH=src python scripts/generate_our_paper_figures.py
"""

from __future__ import annotations

import csv
import statistics
import time
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.database import ComputeGraph, load_or_build_database
from qcr_repro.exact_database import load_or_build_exact
from qcr_repro.exact_reducer import reduce_circuit_exact
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import (
    reduce_circuit,
    reduce_random_sampling,
    reduce_random_sampling_gated,
)
from qcr_repro.token_pool import TokenPool
from qcr_repro.unitary import equivalent_up_to_global_phase

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
COMPARISON_DIR = ROOT / "results" / "comparison"

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}
NISQ_ANGLES = (-1.5707963267948966, -0.7853981633974483, 0.7853981633974483, 1.5707963267948966)

# --------------------------------------------------------------------------- #
# shared style (matches scripts/generate_figures.py)
# --------------------------------------------------------------------------- #

PAPER_C = "#616161"
INPUT_C = "#90a4ae"
OURS_C = "#2e7d32"
COST_C = "#00838f"
BASE_C = "#c62828"
V1_C = "#c62828"
V2_C = "#8e24aa"
V3_C = "#1565c0"
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
    fig.tight_layout()
    fig.savefig(OUT_DIR / out_name)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# minimal circuit-diagram renderer (no qiskit dependency)
# --------------------------------------------------------------------------- #


def _label(gate: GateInstance) -> str:
    if gate.theta is None:
        return gate.name
    import math

    frac = gate.theta / math.pi
    for num, den, sym in ((1, 2, "π/2"), (-1, 2, "-π/2"), (1, 4, "π/4"), (-1, 4, "-π/4")):
        if abs(frac - num / den) < 1e-6:
            return f"{gate.name}({sym})"
    return f"{gate.name}({gate.theta:.2f})"


def draw_circuit(ax, gates: list[GateInstance], num_qubits: int, title: str) -> None:
    """Draw a token chain as a simple ASAP-scheduled circuit diagram."""
    col_of_wire = [0] * num_qubits
    box_w, box_h, gap = 0.8, 0.6, 0.5
    for gate in gates:
        col = max(col_of_wire[q] for q in gate.qubits) + 1
        x = col * (box_w + gap)
        if len(gate.qubits) == 1:
            q = gate.qubits[0]
            y = num_qubits - 1 - q
            rect = mpatches.FancyBboxPatch(
                (x - box_w / 2, y - box_h / 2), box_w, box_h,
                boxstyle="round,pad=0.02,rounding_size=0.05",
                linewidth=0.9, edgecolor="#333333", facecolor="#e8eef7",
            )
            ax.add_patch(rect)
            ax.text(x, y, _label(gate), ha="center", va="center", fontsize=6.5)
        else:
            q0, q1 = gate.qubits
            y0, y1 = num_qubits - 1 - q0, num_qubits - 1 - q1
            ax.plot([x, x], [min(y0, y1), max(y0, y1)], color="#333333", linewidth=1.3, zorder=1)
            for y in (y0, y1):
                rect = mpatches.FancyBboxPatch(
                    (x - box_w / 2, y - box_h / 2), box_w, box_h,
                    boxstyle="round,pad=0.02,rounding_size=0.05",
                    linewidth=0.9, edgecolor="#333333", facecolor="#f7e8e8", zorder=2,
                )
                ax.add_patch(rect)
            ax.text(x, y0, _label(gate), ha="center", va="center", fontsize=6.5, zorder=3)
        for q in gate.qubits:
            col_of_wire[q] = col

    n_cols = max(col_of_wire) + 1 if gates else 1
    width = n_cols * (box_w + gap) + box_w
    for q in range(num_qubits):
        y = num_qubits - 1 - q
        ax.plot([0, width], [y, y], color="#999999", linewidth=0.8, zorder=0)
        ax.text(-box_w, y, f"q{q}", ha="right", va="center", fontsize=9)
    ax.set_xlim(-box_w * 1.6, width)
    ax.set_ylim(-1, num_qubits)
    ax.axis("off")
    ax.set_title(title, fontsize=10.5, fontweight="bold")


# --------------------------------------------------------------------------- #
# Figure 1 (motivation) / Figure 4 (example steps): before/after circuits
# --------------------------------------------------------------------------- #


def fig1_and_fig4_example_reduction(db_exact) -> None:
    """A small example circuit, reduced with our exact ion-trap engine,
    drawn before and after (analog of the paper's Figures 1 and 4)."""
    num_qubits = 3
    gates, _ = random_circuit(num_qubits, 16, "ion_trap", seed=7)
    u0 = circuit_unitary(num_qubits, gates)
    reduced, stats = reduce_circuit_exact(gates, num_qubits, db_exact, budget_s=5.0, seed=0)
    ok = equivalent_up_to_global_phase(u0, circuit_unitary(num_qubits, reduced), atol=1e-9)
    assert ok, "reduction changed the unitary"

    _set_style()
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.5))
    draw_circuit(axes[0], gates, num_qubits,
                 f"Before: {len(gates)} gates (random ion-trap circuit, seed 7)")
    draw_circuit(axes[1], reduced, num_qubits,
                 f"After: {len(reduced)} gates (our exact reducer, {stats.runtime_sec:.2f}s, exact-verified)")
    fig.suptitle("Our Figure 1/4 analog: example circuit before/after reduction", fontsize=11.5, fontweight="bold")
    _finish(fig, "our_figure1_4_example_reduction.png")
    print(f"  fig1/4  {len(gates)} -> {len(reduced)} gates, unitary_ok={ok}")


# --------------------------------------------------------------------------- #
# Figure 2: compute-graph growth, measured on our own ComputeGraph
# --------------------------------------------------------------------------- #


def fig2_compute_graph_growth() -> None:
    _set_style()
    fig, ax = plt.subplots(figsize=(7.5, 5))
    configs = [
        ("ion trap, 2-wire", "ion_trap", 2, range(1, 7), None),
        ("ion trap, 3-wire", "ion_trap", 3, range(1, 5), None),
        ("NISQ, 2-wire", "nisq", 2, range(1, 6), NISQ_ANGLES),
        ("NISQ, 3-wire", "nisq", 3, range(1, 4), NISQ_ANGLES),
    ]
    colors = [OURS_C, COST_C, BASE_C, "#6a1b9a"]
    for (label, gate_set, wires, depths, angles), color in zip(configs, colors):
        pool = TokenPool(num_qubits=wires, gate_set=gate_set, angles=angles)
        nodes = []
        for depth in depths:
            g = ComputeGraph(pool=pool, max_depth=depth)
            nodes.append(g.num_nodes)
            print(f"  fig2  {label} depth={depth}: {g.num_nodes} nodes")
        ax.plot(list(depths), nodes, marker="o", color=color, label=f"{label} ({len(pool.tokens())} tokens)")
    ax.set_yscale("log")
    ax.set_xlabel("compute-graph depth")
    ax.set_ylabel("nodes (log scale)")
    ax.set_title("Our Figure 2 analog: measured compute-graph growth")
    ax.legend()
    _finish(fig, "our_figure2_compute_graph_growth.png")


# --------------------------------------------------------------------------- #
# Figure 3: our pipeline, schematic
# --------------------------------------------------------------------------- #


def fig3_pipeline() -> None:
    _set_style()
    fig, ax = plt.subplots(figsize=(11, 3.2))
    stages = [
        ("Input\ncircuit", INPUT_C),
        ("Pre-pass\n(algebraic/ZX)", "#cfd8dc"),
        ("Cluster +\n1-wire collapse", "#cfd8dc"),
        ("Exhaustive\nsweep (DB)", OURS_C),
        ("dag_compact /\ntransport shuffle", COST_C),
        ("Escape\n(resample)", "#8e24aa"),
        ("Verified\noutput", BASE_C),
    ]
    x = 0.0
    box_w, box_h, gap = 1.7, 1.1, 0.55
    for i, (text, color) in enumerate(stages):
        rect = mpatches.FancyBboxPatch(
            (x, -box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            linewidth=1.0, edgecolor="#333333", facecolor=color, alpha=0.85,
        )
        ax.add_patch(rect)
        ax.text(x + box_w / 2, 0, text, ha="center", va="center", fontsize=9,
                 color="white" if color in (OURS_C, COST_C, BASE_C, "#8e24aa") else "#222222")
        if i < len(stages) - 1:
            ax.annotate("", xy=(x + box_w + gap, 0), xytext=(x + box_w, 0),
                        arrowprops=dict(arrowstyle="->", color="#333333", lw=1.2))
        x += box_w + gap
    ax.annotate("repeat until budget exhausted\nor no window reduces",
                xy=(x * 0.45, -1.3), ha="center", fontsize=8.5, color="#555555")
    ax.set_xlim(-0.3, x)
    ax.set_ylim(-1.8, 1.2)
    ax.axis("off")
    ax.set_title("Our Figure 3 analog: reduce_circuit pipeline (src/reducer.py)")
    _finish(fig, "our_figure3_pipeline.png")


# --------------------------------------------------------------------------- #
# Figure 5: length vs time, our exhaustive sweep vs V2/V3-style loops
# --------------------------------------------------------------------------- #


def fig5_reduction_vs_time(db) -> None:
    num_qubits = 4
    gates, _ = random_circuit(num_qubits, 100, "ion_trap", seed=3)
    budgets = [0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30]

    curves = {"exhaustive sweep": [], "V2 random sampling": [], "V3 RF-gated": []}
    for b in budgets:
        r, _, _ = reduce_circuit(list(gates), num_qubits, db, budget_s=b, seed=0)
        curves["exhaustive sweep"].append(len(r))
        r2, _ = reduce_random_sampling(list(gates), num_qubits, db, budget_sec=b, seed=0)
        curves["V2 random sampling"].append(len(r2))
        r3, _, _ = reduce_random_sampling_gated(list(gates), num_qubits, db, budget_sec=b, seed=0)
        curves["V3 RF-gated"].append(len(r3))
        print(f"  fig5  budget={b}s: sweep={len(r)} V2={len(r2)} V3={len(r3)}")

    _set_style()
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for (label, ys), color in zip(curves.items(), (OURS_C, V2_C, V3_C)):
        ax.plot(budgets, ys, marker="o", color=color, label=label)
    ax.set_xlabel("computation time (s)")
    ax.set_ylabel("circuit length")
    ax.set_title("Our Figure 5 analog: length vs. time (100-gate ion-trap circuit)")
    ax.legend()
    _finish(fig, "our_figure5_reduction_vs_time.png")


# --------------------------------------------------------------------------- #
# Figure 6: computation-time boxplot, our three loop variants
# --------------------------------------------------------------------------- #


def fig6_boxplot_methods(db, num_circuits: int = 12, target: int = 60, budget_s: float = 15.0) -> None:
    num_qubits = 4
    times = {"exhaustive\nsweep": [], "V2 random\nsampling": [], "V3 RF-gated": []}
    for seed in range(num_circuits):
        gates, _ = random_circuit(num_qubits, 100, "ion_trap", seed=seed)

        t0 = time.perf_counter()
        reduce_circuit(list(gates), num_qubits, db, budget_s=budget_s, seed=seed)
        times["exhaustive\nsweep"].append(time.perf_counter() - t0)

        t0 = time.perf_counter()
        reduce_random_sampling(list(gates), num_qubits, db, budget_sec=budget_s, seed=seed)
        times["V2 random\nsampling"].append(time.perf_counter() - t0)

        t0 = time.perf_counter()
        reduce_random_sampling_gated(list(gates), num_qubits, db, budget_sec=budget_s, seed=seed)
        times["V3 RF-gated"].append(time.perf_counter() - t0)
        print(f"  fig6  seed={seed} done")

    _set_style()
    fig, ax = plt.subplots(figsize=(7, 5))
    data = list(times.values())
    bp = ax.boxplot(data, labels=list(times.keys()), patch_artist=True, widths=0.5)
    for patch, color in zip(bp["boxes"], (OURS_C, V2_C, V3_C)):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
    ax.set_ylabel("wall-clock time (s)")
    ax.set_title(f"Our Figure 6 analog: time distribution, n={num_circuits}, budget={budget_s:.0f}s cap")
    _finish(fig, "our_figure6_boxplot_methods.png")


# --------------------------------------------------------------------------- #
# Figure 7: wire reduction, schematic with a real example block
# --------------------------------------------------------------------------- #


def fig7_wire_reduction(db) -> None:
    num_qubits = 5
    full_block = [
        GateInstance("RX", (0,), 1.5707963267948966),
        GateInstance("RXX", (0, 2), 1.5707963267948966),
        GateInstance("RZ", (2,), 1.5707963267948966),
        GateInstance("RZ", (2,), -1.5707963267948966),
    ]
    wires = sorted({q for g in full_block for q in g.qubits})
    forward = {w: i for i, w in enumerate(wires)}
    local_block = [
        GateInstance(g.name, tuple(sorted(forward[q] for q in g.qubits)), g.theta) for g in full_block
    ]
    reduced_local = db.try_reduce(local_block) or local_block
    reverse = {i: w for w, i in forward.items()}
    reduced_full = [
        GateInstance(g.name, tuple(sorted(reverse[q] for q in g.qubits)), g.theta) for g in reduced_local
    ]

    _set_style()
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    draw_circuit(axes[0], full_block, num_qubits, f"Full register ({num_qubits} wires)")
    draw_circuit(axes[1], local_block, len(wires), f"Wire-reduced ({len(wires)} active wires)")
    draw_circuit(axes[2], reduced_full, num_qubits, "Lifted back, reduced")
    fig.suptitle("Our Figure 7 analog: wire reduction -> local lookup -> lift back "
                 "(ReductionDatabase.try_reduce)", fontsize=10.5, fontweight="bold")
    _finish(fig, "our_figure7_wire_reduction.png")
    print(f"  fig7  {len(full_block)} -> {len(reduced_full)} gates on the active block")


# --------------------------------------------------------------------------- #
# Figures 8/9: boxplots from our own committed comparison runs
# --------------------------------------------------------------------------- #


def _boxplot_from_csv(csv_path: Path, methods_order: list[str], labels: list[str],
                       title: str, out_name: str, input_len: int = 300,
                       num_ours: int = 1) -> None:
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        print(f"  skip {out_name}: {csv_path} missing/empty")
        return
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    data = [[input_len]]
    for m in methods_order:
        vals = [int(r["end"]) for r in rows if r["method"] == m and r.get("ok") in ("True", "1")]
        data.append(vals if vals else [0])

    _set_style()
    fig, ax = plt.subplots(figsize=(8.5, 5))
    bp = ax.boxplot(data, labels=["Input"] + labels, patch_artist=True, widths=0.55)
    n = len(data)
    colors = [INPUT_C] + [BASE_C] * (n - 1 - num_ours) + [OURS_C] * num_ours
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
    ax.set_ylabel("# gates")
    ax.set_title(title)
    _finish(fig, out_name)
    print(f"  wrote {out_name} (n={len(rows)} rows)")


def fig8_boxplot_iontrap() -> None:
    _boxplot_from_csv(
        COMPARISON_DIR / "comparison_ion_trap.csv",
        ["qiskit_l1", "qiskit_l2", "qiskit_l3", "exact_len", "exact_cost"],
        ["Q-L1", "Q-L2", "Q-L3", "exact_len", "exact_cost"],
        "Our Figure 8 analog: ion-trap gate counts (committed 100-run data)",
        "our_figure8_boxplot_iontrap.png",
        num_ours=2,
    )


def fig9_boxplot_nisq() -> None:
    _boxplot_from_csv(
        COMPARISON_DIR / "comparison_nisq.csv",
        ["qiskit_l1", "qiskit_l2", "qiskit_l3", "numeric_len", "numeric_cost"],
        ["Q-L1", "Q-L2", "Q-L3", "numeric_len", "numeric_cost"],
        "Our Figure 9 analog: NISQ gate counts (committed 100-run data)",
        "our_figure9_boxplot_nisq.png",
        num_ours=2,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("== figures 8/9: boxplots from committed comparison CSVs ==")
    fig8_boxplot_iontrap()
    fig9_boxplot_nisq()

    print("== figure 2: measured compute-graph growth ==")
    fig2_compute_graph_growth()

    print("== figure 3: pipeline schematic ==")
    fig3_pipeline()

    print("== building ion-trap databases (numeric + exact) ==")
    db = load_or_build_database("ion_trap", ION_DEPTHS, verbose=True)
    db_exact = load_or_build_exact("ion_trap", ION_DEPTHS)

    print("== figures 1/4: example circuit before/after ==")
    fig1_and_fig4_example_reduction(db_exact)

    print("== figure 7: wire reduction on a real block ==")
    fig7_wire_reduction(db)

    print("== figure 5: length vs time (single 100-gate circuit) ==")
    fig5_reduction_vs_time(db)

    print("== figure 6: time boxplot over multiple circuits (slow) ==")
    fig6_boxplot_methods(db)

    print("\nFigure 10 (IBM hardware measurement) has no analog here: this "
          "project has no quantum-hardware access. Not fabricated.")
    print("Done. See figures/our_figure*.png")


if __name__ == "__main__":
    main()
