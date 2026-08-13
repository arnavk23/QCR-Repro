"""Generate the paper's tables and figures (Rosenhahn, Osborne & Hirche, NJP 27 104509) from data accumulated on this setup.

Covers Table 1 (compute-graph growth), Table 2 (V1/V2/V3 timing), Tables 3-5 (single-circuit examples), Tables 6-7 (100-run statistics), Figures 1-10. Figure 10 (IBM hardware) is a placeholder using the paper's reported numbers.

Usage:
    PYTHONPATH=src python scripts/generate_paper_tables.py --single-circuits | --timing | --figures | --tables | --all"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from qcr_repro.circuits import count_gates, random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.exact_reducer import reduce_circuit_exact, verify_exact
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import (
    reduce_random_sampling,
    reduce_random_sampling_gated,
    reduce_with_database,
    reduce_with_lookup,
)
from qcr_repro.unitary import equivalent_up_to_global_phase

ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "paper_tables"
FIG_DIR = ROOT / "figures" / "paper"
COMPARISON_DIR = ROOT / "results" / "comparison"

# --------------------------------------------------------------------------- #
# Reference values read from the paper (Tables 2-7)
# --------------------------------------------------------------------------- #

# Table 2: mean/std computation time (s) of V1/V2/V3 on 100 length-100 circuits.
PAPER_TABLE2 = {
    "V1-RS": {"mean": 199.0, "std": 351.5},
    "V2-DR": {"mean": 55.0, "std": 96.3},
    "V3-RF": {"mean": 38.0, "std": 39.8},
}

# Tables 3-5: single-circuit examples (per-type gate counts).  The paper's
# exact example circuits are not published, so we re-run the experiment
# (same architecture, length, gate set) on our own instances and report both.
PAPER_TABLE3 = {  # ion trap, 4 qubits, length 300 (RX, RY, RZ, RXX)
    "original": [82, 71, 86, 61], "Q-L1": [36, 43, 63, 61], "Q-L2": [41, 48, 58, 60],
    "Q-L3": [41, 48, 58, 60], "B-L2": [59, 10, 65, 68], "B-L3": [54, 9, 72, 59],
    "B-L4": [69, 0, 79, 58], "Ours": [9, 27, 38, 38],
}
PAPER_TABLE4 = {  # NISQ, 6 qubits, length 300 (RX, RZ, CZ)
    "original": [93, 100, 107], "Q-L1": [64, 66, 93], "Q-L2": [63, 39, 66],
    "Q-L3": [63, 39, 66], "B-L2": [82, 100, 94], "B-L3": [84, 111, 68],
    "B-L4": [90, 132, 60], "Ours": [51, 22, 60],
}
PAPER_TABLE5 = {  # NISQ, 15 qubits, length 500 (RX, RZ, CZ)
    "original": [96, 91, 313], "Q-L1": [74, 74, 285], "Q-L2": [74, 36, 203],
    "Q-L3": [74, 36, 203], "B-L2": [87, 112, 306], "B-L3": [126, 164, 269],
    "B-L4": [293, 395, 253], "Ours": [72, 26, 191],
}

# Tables 6/7: 100-run statistics (means, plus std in brackets).
PAPER_TABLE6 = {  # ion trap, 4 qubits, length 300 (RX, RY, RZ, RXX)
    "In": [78, 83, 78, 59], "Q-L1": [32, 46, 59, 59], "Q-L2": [33, 49, 66, 56],
    "Q-L3": [33, 49, 66, 56], "B-L2": [49, 2, 66, 37], "B-L3": [39, 1, 57, 32],
    "B-L4": [40, 0, 58, 28], "Ours": [10, 29, 29, 43],
}
PAPER_TABLE7 = {  # NISQ, 4 qubits, length 300 (RX, RZ, CZ)
    "In": [108, 109, 82], "Q-L1": [59, 68, 69], "Q-L2": [59, 39, 51],
    "Q-L3": [59, 39, 51], "B-L2": [67, 85, 62], "B-L3": [55, 70, 39],
    "B-L4": [56, 75, 37], "Ours": [45, 19, 43],
}

TABLE3_CFG = {"gateset": "ion_trap", "qubits": 4, "length": 300, "types": ["RX", "RY", "RZ", "RXX"]}
TABLE4_CFG = {"gateset": "nisq", "qubits": 6, "length": 300, "types": ["RX", "RZ", "CZ"]}
TABLE5_CFG = {"gateset": "nisq", "qubits": 15, "length": 500, "types": ["RX", "RZ", "CZ"]}

GATE_TYPES = {"ion_trap": ["RX", "RY", "RZ", "RXX"], "nisq": ["RX", "RZ", "CZ"]}

# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #


def load_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _qiskit_transpile(gates, num_qubits: int, gateset: str, level: int):
    """Transpile with qiskit to the target basis at level L.

    Returns (counts, total, ok) or None if unavailable/failed.
    """
    try:
        from qiskit import QuantumCircuit, transpile
        from qiskit.circuit.library import RXGate, RYGate, RZGate, CZGate, RXXGate
        from qiskit.quantum_info import Operator

        qc = QuantumCircuit(num_qubits)
        rev = lambda qubits: [num_qubits - 1 - q for q in qubits]  # little-endian
        for g in gates:
            if g.name == "RX":
                qc.append(RXGate(g.theta), rev(g.qubits))
            elif g.name == "RY":
                qc.append(RYGate(g.theta), rev(g.qubits))
            elif g.name == "RZ":
                qc.append(RZGate(g.theta), rev(g.qubits))
            elif g.name == "RXX":
                qc.append(RXXGate(g.theta), rev(g.qubits))
            elif g.name == "CZ":
                qc.append(CZGate(), rev(g.qubits))
        basis = ["rx", "ry", "rz", "rxx"] if gateset == "ion_trap" else ["rx", "rz", "cz"]
        t = transpile(qc, basis_gates=basis, optimization_level=level)
        ops = t.count_ops()
        total = sum(n for name, n in ops.items() if name not in ("global_phase", "id", "delay"))
        counts = {
            "RX": ops.get("rx", 0), "RY": ops.get("ry", 0), "RZ": ops.get("rz", 0),
            "RXX": ops.get("rxx", 0), "CZ": ops.get("cz", 0),
        }
        u1 = Operator(t).data
        ok = equivalent_up_to_global_phase(circuit_unitary(num_qubits, gates), u1, atol=1e-5)
        return counts, total, ok
    except Exception:
        return None


def _bqskit_compile(gates, num_qubits: int, gateset: str, level: int, timeout_s: float = 240.0):
    """Compile with BQSKit at optimization level L (paper's B-L2..L4).

    Returns (counts, total, ok) or None if unavailable/timed-out/failed.
    """
    import signal

    class _Timeout(Exception):
        pass

    def _handler(signum, frame):
        raise _Timeout()

    try:
        import numpy as np
        from bqskit import Circuit, MachineModel, compile
        from bqskit.compiler.gateset import GateSet
        from bqskit.ir.gates import RXGate, RYGate, RZGate, CZGate, RXXGate

        c = Circuit(num_qubits)
        for g in gates:
            if g.name == "RX":
                c.append_gate(RXGate().with_all_frozen_params([g.theta]), [g.qubits[0]])
            elif g.name == "RY":
                c.append_gate(RYGate().with_all_frozen_params([g.theta]), [g.qubits[0]])
            elif g.name == "RZ":
                c.append_gate(RZGate().with_all_frozen_params([g.theta]), [g.qubits[0]])
            elif g.name == "RXX":
                c.append_gate(RXXGate().with_all_frozen_params([g.theta]), [g.qubits[0], g.qubits[1]])
            elif g.name == "CZ":
                c.append_gate(CZGate(), [g.qubits[0], g.qubits[1]])
        basis = [RXGate(), RYGate(), RZGate(), RXXGate()] if gateset == "ion_trap" \
            else [RXGate(), RZGate(), CZGate()]
        model = MachineModel(num_qubits, gate_set=GateSet(list(basis)))
        old = signal.signal(signal.SIGALRM, _handler)
        signal.alarm(int(timeout_s))
        try:
            out = compile(c, model=model, optimization_level=level)
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
        counts = {}
        for gate, n in dict(out.gate_counts).items():
            name = type(gate).__name__.replace("Gate", "")
            counts[name] = counts.get(name, 0) + int(n)
        u1 = np.asarray(out.get_unitary())
        ok = equivalent_up_to_global_phase(circuit_unitary(num_qubits, gates), u1, atol=1e-5)
        return counts, out.num_operations, ok
    except Exception:
        return None


def _mean_std(values):
    if not values:
        return float("nan"), float("nan")
    return statistics.mean(values), statistics.pstdev(values)


def _rows_for_comparison_csv(gateset: str) -> list[dict]:
    path = COMPARISON_DIR / f"comparison_{gateset}.csv"
    if not path.exists():
        print(f"[warn] {path} missing - run benchmark_comparison.py first")
        return []
    return [r for r in load_rows(path) if r["end"] != "-1" and int(r["end"]) >= 0]


# --------------------------------------------------------------------------- #
# Tables 3-5: single-circuit examples
# --------------------------------------------------------------------------- #


def run_single_circuit(cfg: dict, seed: int) -> dict:
    """Reduce one random circuit with qiskit L1-3, BQSKit L2-4, and our engines."""
    gateset, nq, length, types = cfg["gateset"], cfg["qubits"], cfg["length"], cfg["types"]
    weights = None if gateset == "ion_trap" else {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}
    gates, _ = random_circuit(nq, length, gateset, seed=seed, weights=weights)
    inp = {t: count_gates(gates).get(t, 0) for t in types}

    out: dict[str, dict] = {"original": inp}

    for level in (1, 2, 3):
        res = _qiskit_transpile(gates, nq, gateset, level)
        if res:
            counts, total, ok = res
            out[f"Q-L{level}"] = {t: counts.get(t, 0) for t in types}
        else:
            out[f"Q-L{level}"] = None

    for level in (2, 3, 4):
        res = _bqskit_compile(gates, nq, gateset, level, timeout_s=240.0)
        if res:
            counts, total, ok = res
            out[f"B-L{level}"] = {t: counts.get(t, 0) for t in types}
        else:
            out[f"B-L{level}"] = None

    # our reducers (local 3/4-wire blocks work for any circuit size)
    if gateset == "ion_trap" and nq == 4:
        db = load_or_build_database("ion_trap", {1: 12, 2: 10, 3: 7, 4: 5})
        r, stats = reduce_with_database(gates, nq, db, budget_sec=15.0, seed=seed)
        out["Ours-exact-sweep"] = {t: count_gates(r).get(t, 0) for t in types}
        out["Ours-exact-sweep"]["total"] = len(r)
    else:
        db = load_or_build_database(gateset, {1: 12, 2: 6, 3: 5, 4: 4})
        r, stats = reduce_with_database(gates, nq, db, budget_sec=15.0, seed=seed)
        out["Ours-DB-sweep"] = {t: count_gates(r).get(t, 0) for t in types}
        out["Ours-DB-sweep"]["total"] = len(r)

    return out


def cmd_single_circuits(only: str | None = None) -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    for name, cfg, paper in [
        ("table3_ion4", TABLE3_CFG, PAPER_TABLE3),
        ("table4_nisq6", TABLE4_CFG, PAPER_TABLE4),
        ("table5_nisq15", TABLE5_CFG, PAPER_TABLE5),
    ]:
        if only and name != only:
            continue
        print(f"[single] {name}: {cfg['gateset']} q{cfg['qubits']} len {cfg['length']}")
        out = run_single_circuit(cfg, seed=42)
        rows = []
        for method, counts in out.items():
            if counts is None:
                rows.append({"method": method, **{t: -1 for t in cfg["types"]}, "total": -1})
            else:
                rows.append({"method": method, **{t: counts[t] for t in cfg["types"]},
                             "total": sum(counts[t] for t in cfg["types"])})
        write_csv(RESULT_DIR / f"{name}.csv", rows)
        # also store paper reference
        paper_rows = [{"method": m, **{t: v[i] for i, t in enumerate(cfg["types"])},
                       "total": sum(v)} for m, v in paper.items()]
        write_csv(RESULT_DIR / f"{name}_paper.csv", paper_rows)
        print("   ours:")
        for r in rows:
            print("   ", r["method"], {t: r[t] for t in cfg["types"]}, r["total"])
        print("   paper (reference):")
        for r in paper_rows:
            print("   ", r["method"], {t: r[t] for t in cfg["types"]}, r["total"])


# --------------------------------------------------------------------------- #
# Table 2: timing on 100 length-100 circuits
# --------------------------------------------------------------------------- #


def _reduce_traj_v1(gates, nq, db, budget):
    r, stats = reduce_random_sampling(gates, nq, db, budget_sec=budget, seed=0)
    return len(r), stats.runtime_sec


def _reduce_traj_v2(gates, nq, db, budget):
    r, stats = reduce_with_database(gates, nq, db, budget_sec=budget, seed=0)
    return len(r), stats.runtime_sec


def _reduce_traj_v3(gates, nq, db, budget):
    """V3-RF: random sampling gated by the online classifier (skips lookups)."""
    r, stats, gate = reduce_random_sampling_gated(gates, nq, db, budget_sec=budget, seed=0)
    return len(r), stats.runtime_sec, gate


def _reduce_traj_exact(gates, nq, db, budget):
    r, stats = reduce_circuit_exact(gates, nq, db, budget_s=budget, seed=0, cost_aware=False)
    return len(r), stats.runtime_sec


def cmd_timing(num_circuits: int = 100, v1_circuits: int = 10, budget: float = 20.0) -> None:
    """Table 2 style: reduce length-100 circuits, record time + end length."""
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    db_ion = load_or_build_database("ion_trap", {1: 12, 2: 10, 3: 7, 4: 5})
    db_exact = load_or_build_database("exact_ion_trap", {1: 12, 2: 10, 3: 7, 4: 5}) \
        if False else None  # exact engine loads its own via signature; keep numeric DB here

    methods = {
        "V1-RS (ours)": _reduce_traj_v1,
        "V2-DR (ours)": _reduce_traj_v2,
        "V3-RF (ours)": _reduce_traj_v3,
        "exact (ours)": _reduce_traj_exact,
    }

    # V1 (random sampling) is slow - run on a subset, others on all circuits.
    rows: list[dict] = []
    v3_gate_stats = {"skips": 0, "attempts": 0, "reductions": 0}
    for seed in range(num_circuits):
        gates, _ = random_circuit(4, 100, "ion_trap", seed=seed)
        for label, fn in methods.items():
            n_here = v1_circuits if label.startswith("V1") else num_circuits
            if seed >= n_here:
                continue
            if label.startswith("exact"):
                # exact engine needs its own DB; load via its loader
                from qcr_repro.exact_database import load_or_build_exact
                dbx = load_or_build_exact("ion_trap", {1: 12, 2: 10, 3: 7, 4: 5})
                end, secs = _reduce_traj_exact(gates, 4, dbx, budget)
            else:
                res = fn(gates, 4, db_ion, budget)
                if label.startswith("V3"):
                    end, secs, gate = res
                    v3_gate_stats["skips"] += gate.lookups_skipped
                    v3_gate_stats["attempts"] += gate.lookups_attempted
                    v3_gate_stats["reductions"] += gate.reductions_found
                else:
                    end, secs = res
            rows.append({"method": label, "seed": seed, "start": 100, "end": end,
                         "runtime_sec": round(secs, 3)})
    write_csv(RESULT_DIR / "table2_timing.csv", rows)

    # summary
    summary = []
    for label in methods:
        rs = [r for r in rows if r["method"] == label]
        ends = [r["end"] for r in rs]
        secs = [r["runtime_sec"] for r in rs]
        summary.append({"method": label, "n": len(rs),
                        "mean_len": round(statistics.mean(ends), 2),
                        "mean_time_s": round(statistics.mean(secs), 2),
                        "std_time_s": round(statistics.pstdev(secs), 2)})
    write_csv(RESULT_DIR / "table2_summary.csv", summary)
    for s in summary:
        print(f"   {s['method']:<16} n={s['n']} mean_len={s['mean_len']} "
              f"time={s['mean_time_s']}s (std {s['std_time_s']})")
    if v3_gate_stats["attempts"] + v3_gate_stats["skips"]:
        print(f"   V3-RF gate totals: attempts={v3_gate_stats['attempts']} "
              f"skipped={v3_gate_stats['skips']} "
              f"reductions={v3_gate_stats['reductions']} "
              f"(skip rate {100.0 * v3_gate_stats['skips'] / (v3_gate_stats['attempts'] + v3_gate_stats['skips']):.1f}%)")


# --------------------------------------------------------------------------- #
# Tables 6/7 (100 runs) - use benchmark_comparison.py output
# --------------------------------------------------------------------------- #


def _comparison_summary(gateset: str) -> list[dict]:
    rows = _rows_for_comparison_csv(gateset)
    if not rows:
        return []
    types = GATE_TYPES[gateset]
    out = []
    for method in sorted({r["method"] for r in rows}):
        rs = [r for r in rows if r["method"] == method]
        if not rs:
            continue
        out.append({
            "method": method,
            "n": len(rs),
            **{f"{t}_mean": round(statistics.mean(int(r[t]) for r in rs), 2)
               for t in types},
            **{f"{t}_std": round(statistics.pstdev(int(r[t]) for r in rs), 2)
               for t in types},
            "total_mean": round(statistics.mean(int(r["end"]) for r in rs), 2),
            "total_std": round(statistics.pstdev(int(r["end"]) for r in rs), 2),
            "twq_mean": round(statistics.mean(int(r["twq"]) for r in rs), 2),
            "ok_rate": round(sum(r["ok"] == "True" for r in rs) / len(rs), 3),
        })
    return out


# --------------------------------------------------------------------------- #
# Tables assembly
# --------------------------------------------------------------------------- #


def _md_table(header, rows, label_note=""):
    lines = [f"| {' | '.join(header)} |", f"|{'---|' * len(header)}"]
    for r in rows:
        lines.append(f"| {' | '.join(str(x) for x in r)} |")
    return "\n".join(lines) + "\n"


def cmd_tables() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    md: list[str] = [
        "# Paper results (Rosenhahn et al., NJP 27, 104509, 2025)",
        "",
        "Reference values are transcribed from the paper PDF (Tables 1-7). "
        "Our values come from `results/paper_tables/` and `results/comparison/`. "
        "Figure 10 (hardware) requires authenticated hardware access.",
        "",
    ]

    # Table 2
    md += ["## Table 2: computation-time statistics (100 length-100 circuits)",
           "", _md_table(
               ["method", "n", "mean end length", "mean time (s)", "std time (s)"],
               [[s["method"], s["n"], s["mean_len"], s["mean_time_s"], s["std_time_s"]]
                for s in load_rows(RESULT_DIR / "table2_summary.csv")],
           ), ""]
    md += ["Paper reference (mean/std time): " +
           "; ".join(f"{k} {v['mean']}s/{v['std']}s" for k, v in PAPER_TABLE2.items()),
           ""]

    # Tables 3-5
    for name, cfg, paper in [("table3_ion4", TABLE3_CFG, PAPER_TABLE3),
                             ("table4_nisq6", TABLE4_CFG, PAPER_TABLE4),
                             ("table5_nisq15", TABLE5_CFG, PAPER_TABLE5)]:
        ours = load_rows(RESULT_DIR / f"{name}.csv")
        md += [f"## {name.replace('_', ' ')} ({cfg['gateset']}, q{cfg['qubits']}, "
               f"len {cfg['length']})", "",
               "**Paper (reference):**", "",
               _md_table(["method"] + cfg["types"] + ["total"],
                         [[r["method"]] + [r[t] for t in cfg["types"]] +
                          [sum(r[t] for t in cfg["types"])] for r in
                          load_rows(RESULT_DIR / f"{name}_paper.csv")]),
               "", "**Our instance (seed 42):**", "",
               _md_table(["method"] + cfg["types"] + ["total"],
                         [[r["method"]] + [r[t] for t in cfg["types"]] + [r["total"]]
                          for r in ours]),
               ""]

    # Tables 6/7
    for gateset, paper, label in [("ion_trap", PAPER_TABLE6, "Table 6 (ion trap)"),
                                  ("nisq", PAPER_TABLE7, "Table 7 (NISQ)")]:
        types = GATE_TYPES[gateset]
        ours = _comparison_summary(gateset)
        md += [f"## {label}: 100-run statistics (4 qubits, length 300)", "",
               "**Paper (reference means):**", "",
               _md_table(["method"] + types + ["total"],
                         [[m] + paper[m] + [sum(paper[m])] for m in paper]),
               "", "**Ours (means over identical circuits):**", "",
               _md_table(["method"] + [f"{t}" for t in types] + ["total", "ok"],
                         [[r["method"]] + [r[f"{t}_mean"] for t in types] +
                          [r["total_mean"], r["ok_rate"]] for r in ours]),
               ""]

    (RESULT_DIR / "paper_tables.md").write_text("\n".join(md), encoding="utf-8")
    print(f"Wrote {RESULT_DIR / 'paper_tables.md'}")


# --------------------------------------------------------------------------- #
# Figures 1-10 (paper style)
# --------------------------------------------------------------------------- #

_STYLE = {
    "font.family": "DejaVu Sans", "font.size": 10.5, "axes.grid": True,
    "grid.color": "#d7d7d7", "grid.linestyle": "--", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "figure.dpi": 200, "savefig.dpi": 200,
}
OURS_C = "#2e7d32"
BASE_C = "#c62828"
PAPER_C = "#616161"


def _finish(fig, name):
    fig.tight_layout()
    fig.savefig(FIG_DIR / name)
    plt.close(fig)


def fig1_bell_state():
    """Figure 1: Bell-state mapping - naive vs resource-efficient (schematic)."""
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 4.6))

    def draw(ax, seq, title):
        ax.set_xlim(0, 10)
        ax.set_ylim(-1.2, 1.2)
        ax.axis("off")
        ax.set_title(title, fontsize=10.5)
        ax.hlines([-0.55, 0.55], 0.15, 9.85, color="#444444", linewidths=1.4)
        ax.text(0.05, 0.05, "", fontsize=9)
        ax.text(-0.02, 0.62, "q1", fontsize=9, ha="right", va="center")
        ax.text(-0.02, -0.48, "q0", fontsize=9, ha="right", va="center")
        x = 0.7
        for label, span in seq:
            w = 1.1 if span == 1 else 2.3
            y = 0.55 if span == 1 else 0.0
            rect = plt.Rectangle((x, y - 0.24), w, 0.48, fill=True,
                                 facecolor="#cfd8dc", edgecolor="#333333", linewidth=0.8)
            ax.add_patch(rect)
            ax.text(x + w / 2, y, label, ha="center", va="center", fontsize=7.2)
            x += w + 0.28
        ax.set_xlim(0, x + 0.3)

    draw(axes[0], [("H", 1), ("RX", 1), ("RZ", 1), ("RXX", 2), ("RY", 1),
                   ("RZ", 1), ("RXX", 2), ("RZ", 1)], "Naive: one-gate-at-a-time mapping")
    draw(axes[1], [("RX", 1), ("RXX", 2), ("RY", 1), ("RZ", 1)],
         "Resource-efficient: direct two-qubit synthesis")
    fig.suptitle("Figure 1: Bell-state preparation on an ion-trap gate set (schematic)")
    _finish(fig, "figure1_bell_state.png")


def fig2_database():
    """Figure 2: compute-graph schematic (depth-1 operators around the root)."""
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.axis("off")
    ax.set_title("Figure 2: Compute graph for the operator pool\n"
                 "(depth 1: 14 operator nodes; node counts per depth from Table 1)",
                 fontsize=10.5)

    root = (0.0, 0.0)
    ax.plot(*root, "o", ms=14, color=PAPER_C, zorder=3)
    ax.text(*root, "I", ha="center", va="center", fontsize=8, color="white")
    n = 14
    r = 1.35
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False) + 0.15
    for i, a in enumerate(angles):
        x, y = r * np.cos(a), r * np.sin(a)
        ax.plot([root[0], x], [root[1], y], color="#bbbbbb", linewidth=0.8, zorder=1)
        ax.plot(x, y, "o", ms=9, color=OURS_C, zorder=3)
        ax.text(x, y, f"g{i + 1}", ha="center", va="center", fontsize=6, color="white")

    # expansion to depth 2 (schematic, a few nodes)
    for a in angles[:4]:
        x, y = r * np.cos(a), r * np.sin(a)
        for j in range(2):
            x2 = x * 1.9 + 0.25 * (j - 0.5)
            y2 = y * 1.9 + 0.25 * (j - 0.5)
            ax.plot([x, x2], [y, y2], color="#cccccc", linewidth=0.6)
            ax.plot(x2, y2, ".", ms=3, color="#555555")
    ax.text(3.1, 2.2, "depth 2: 114 nodes", fontsize=8.5, color="#333333")
    ax.text(3.1, 1.85, "depth 3: 584 nodes", fontsize=8.5, color="#333333")
    ax.text(3.1, 1.5, "depth 4: 2024 nodes", fontsize=8.5, color="#333333")
    ax.set_xlim(-4.2, 4.2)
    ax.set_ylim(-3.0, 3.0)
    _finish(fig, "figure2_database.png")


def fig3_variants():
    """Figure 3: V1/V2/V3 optimization variants (schematic)."""
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.axis("off")
    ax.set_title("Figure 3: Optimization variants (schematic)", fontsize=11)

    variants = [
        ("V1: Random search (RS)", ["sample block", "random replacement", "repeat"]),
        ("V2: Database retrieval (DR)", ["sample block", "DB lookup", "replace if shorter"]),
        ("V3: RF-supported lookup", ["sample block", "RF: reducible?", "DB lookup"]),
    ]
    for i, (title, steps) in enumerate(variants):
        x0 = 0.04 + i * 0.33
        ax.add_patch(plt.Rectangle((x0, 0.42), 0.27, 0.36, fill=False,
                                   edgecolor="#2e7d32", linewidth=1.4))
        ax.text(x0 + 0.135, 0.74, title, ha="center", fontsize=9, fontweight="bold")
        for j, s in enumerate(steps):
            ax.text(x0 + 0.135, 0.62 - j * 0.1, f"{j + 1}. {s}", ha="center", fontsize=8)
        if i < 2:
            ax.annotate("", xy=(x0 + 0.31, 0.6), xytext=(x0 + 0.29, 0.6),
                        arrowprops=dict(arrowstyle="->"))
    _finish(fig, "figure3_variants.png")


def _trajectory(gates, nq, reduce_fn, budgets, seed=0):
    """Run fresh reductions at increasing budgets; return (time, length) points."""
    pts = []
    for b in budgets:
        end, secs = reduce_fn(gates, nq, b)
        pts.append((secs, end))
    return pts


def fig4_optimization_steps():
    """Figure 4: optimization steps on a 40-gate example circuit (real data)."""
    gates, _ = random_circuit(4, 40, "ion_trap", seed=7)
    db = load_or_build_database("ion_trap", {1: 12, 2: 10, 3: 7, 4: 5})
    budgets = [0.5, 1, 2, 4, 8, 16]
    pts = _trajectory(gates, 4, lambda g, n, b: _reduce_traj_v2(g, n, db, b), budgets)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    times = [p[0] for p in pts]
    lens = [p[1] for p in pts]
    ax.plot([0] + times, [len(gates)] + lens, marker="o", markersize=5, linewidth=2,
            color=OURS_C)
    ax.set_xlabel("Computation time (s)")
    ax.set_ylabel("Circuit length")
    ax.set_title("Figure 4: Term-replacement steps on a 40-gate example (ours)")
    ax.annotate(f"start: {len(gates)}", (0, len(gates)), xytext=(6, -14),
                textcoords="offset points", fontsize=9)
    ax.annotate(f"end: {lens[-1]}", (times[-1], lens[-1]), xytext=(-34, 8),
                textcoords="offset points", fontsize=9, fontweight="bold", color=OURS_C)
    _finish(fig, "figure4_optimization_steps.png")


def fig5_convergence():
    """Figure 5: length vs time for V1/V2/exact on a 100-gate example (real data)."""
    gates, _ = random_circuit(4, 100, "ion_trap", seed=3)
    db = load_or_build_database("ion_trap", {1: 12, 2: 10, 3: 7, 4: 5})
    from qcr_repro.exact_database import load_or_build_exact
    dbx = load_or_build_exact("ion_trap", {1: 12, 2: 10, 3: 7, 4: 5})
    budgets = [2, 4, 8, 16, 30]

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for label, fn, color in [
        ("V1-RS (random sampling)", lambda g, n, b: _reduce_traj_v1(g, n, db, b), BASE_C),
        ("V2-DR (DB retrieval)", lambda g, n, b: _reduce_traj_v2(g, n, db, b), "#f9a825"),
        ("exact (symplectic)", lambda g, n, b: _reduce_traj_exact(g, n, dbx, b), OURS_C),
    ]:
        pts = _trajectory(gates, 4, fn, budgets)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", markersize=4.5,
                linewidth=1.8, color=color, label=label)
    ax.set_xlabel("Computation time (s)")
    ax.set_ylabel("Circuit length")
    ax.set_title("Figure 5: Reduction of a 100-gate example circuit (ours)")
    ax.legend(loc="upper right")
    _finish(fig, "figure5_convergence.png")


def fig6_method_stats():
    """Figure 6: boxplot of end lengths per method (Table 2 timing run)."""
    path = RESULT_DIR / "table2_timing.csv"
    if not path.exists():
        print("WARNING: figure6 skipped (run --timing first)")
        return
    rows = load_rows(path)
    methods = sorted({r["method"] for r in rows})
    data = [[int(r["end"]) for r in rows if r["method"] == m] for m in methods]

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    bp = ax.boxplot(data, tick_labels=[m.replace(" (ours)", "") for m in methods],
                    patch_artist=True, medianprops=dict(color="#c62828", linewidth=1.6))
    for patch in bp["boxes"]:
        patch.set(facecolor="#dcebe2", edgecolor=OURS_C, linewidth=1.1)
    ax.set_ylabel("End circuit length (start 100)")
    ax.set_title("Figure 6: Method statistics on 100 length-100 circuits (ours)")
    _finish(fig, "figure6_method_stats.png")


def fig7_scaling():
    """Figure 7: scaling principle - local blocks on larger circuits (schematic)."""
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 4.8))

    def wires(ax, n, title, highlight=None):
        ax.set_xlim(0, 10)
        ax.set_ylim(-0.1, 1.1)
        ax.axis("off")
        ax.set_title(title, fontsize=9.5, loc="left")
        y = 0.9
        for q in range(n):
            ax.hlines(y - q * 0.18, 0.3, 9.7, color="#888888", linewidths=1.0)
        if highlight:
            x0, x1, y0, y1 = highlight
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=True,
                                       facecolor="#ffe082", edgecolor="#c77700",
                                       linewidth=1.2, alpha=0.7))

    wires(axes[0], 6, "Full circuit (6 qubits) with a selected 3-qubit block",
          highlight=(2.2, 6.4, 0.9 - 2 * 0.18 - 0.16, 0.9))
    wires(axes[1], 3, "Reduce the local block in 3-qubit space",
          highlight=(0.3, 9.7, 0.9 - 2 * 0.18 - 0.16, 0.9))
    wires(axes[2], 6, "Lift the reduced block back to the full circuit")
    fig.suptitle("Figure 7: Scaling to arbitrary qubit counts via local blocks")
    _finish(fig, "figure7_scaling.png")


def _boxplot_from_comparison(gateset: str, types, title, out_name, paper_rows=None):
    rows = _rows_for_comparison_csv(gateset)
    if not rows:
        print(f"WARNING: {out_name} skipped (comparison CSV missing)")
        return
    order = [m for m in ["paper"] + sorted({r["method"] for r in rows})]
    data = []
    labels = []
    if paper_rows:
        data.append([sum(paper_rows["Ours"])])
        labels.append("Paper\n'Ours'")
    for m in sorted({r["method"] for r in rows}):
        rs = [r for r in rows if r["method"] == m]
        data.append([int(r["end"]) for r in rs])
        labels.append(m.replace("_", " "))
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True,
                    medianprops=dict(color="#c62828", linewidth=1.6))
    for patch, k in zip(bp["boxes"], range(len(data))):
        patch.set(facecolor=PAPER_C if k == 0 else "#dcebe2", edgecolor="#333333", linewidth=1.0)
    ax.set_ylabel("Final gate count (start 300)")
    ax.set_title(title)
    _finish(fig, out_name)


def fig8_ion_boxplot():
    """Figure 8: ion-trap statistical summary (boxplot over the 100-run set)."""
    _boxplot_from_comparison("ion_trap", GATE_TYPES["ion_trap"],
                             "Figure 8: Ion trap (4 qubits, 100 runs) - gate-count distribution",
                             "figure8_ion_boxplot.png")


def fig9_nisq_boxplot():
    """Figure 9: NISQ statistical summary (boxplot over the 100-run set)."""
    _boxplot_from_comparison("nisq", GATE_TYPES["nisq"],
                             "Figure 9: NISQ (4 qubits, 100 runs) - gate-count distribution",
                             "figure9_nisq_boxplot.png")


def fig10_hardware():
    """Figure 10: hardware experiments - placeholder (requires IBM access)."""
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    ax.axis("off")
    ax.set_title("Figure 10: Quantum hardware experiments (IBM Brisbane / Kyiv)", fontsize=11)
    ax.text(0.5, 0.62,
            "Unavailable in this environment:\n"
            "the paper measures equivalent long/short circuits on IBM Eagle r3\n"
            "hardware (Brisbane, Kyiv) and compares outcome distributions with\n"
            "the ideal simulation.  Requires authenticated hardware access.",
            ha="center", va="center", fontsize=10.5, transform=ax.transAxes)
    ax.text(0.5, 0.22,
            "Paper finding: the shorter (reduced) circuit deviates less from the\n"
            "ideal simulation, confirming the noise benefit of reduction.",
            ha="center", va="center", fontsize=9.5, color="#555555",
            transform=ax.transAxes, style="italic")
    _finish(fig, "figure10_hardware.png")


def cmd_figures() -> None:
    plt.rcParams.update(_STYLE)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig1_bell_state()
    fig2_database()
    fig3_variants()
    fig4_optimization_steps()
    fig5_convergence()
    fig6_method_stats()
    fig7_scaling()
    fig8_ion_boxplot()
    fig9_nisq_boxplot()
    fig10_hardware()
    print(f"Saved paper figures to {FIG_DIR}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate the paper's results and figures.")
    parser.add_argument("--single-circuits", action="store_true", help="Tables 3-5 examples")
    parser.add_argument("--single-table", type=str, default=None,
                        choices=["table3_ion4", "table4_nisq6", "table5_nisq15"],
                        help="only run this single-circuit example (BQSKit is slow)")
    parser.add_argument("--timing", action="store_true", help="Table 2 timing run")
    parser.add_argument("--timing-circuits", type=int, default=100)
    parser.add_argument("--timing-v1-circuits", type=int, default=10)
    parser.add_argument("--timing-budget", type=float, default=20.0)
    parser.add_argument("--figures", action="store_true", help="paper-style figures 1-10")
    parser.add_argument("--tables", action="store_true", help="assemble paper_tables.md")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    if args.all:
        args.single_circuits = args.timing = args.figures = args.tables = True

    if args.single_circuits:
        cmd_single_circuits(only=args.single_table)
    if args.timing:
        cmd_timing(num_circuits=args.timing_circuits,
                   v1_circuits=args.timing_v1_circuits, budget=args.timing_budget)
    if args.figures:
        cmd_figures()
    if args.tables:
        cmd_tables()
    if not any([args.single_circuits, args.timing, args.figures, args.tables, args.all]):
        parser.print_help()


if __name__ == "__main__":
    main()
