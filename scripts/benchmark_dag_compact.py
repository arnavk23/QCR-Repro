"""Benchmark: DAG block-compaction (src/dag.py) vs the baseline reducer and the paper's reported numbers.

dag_compact deterministically reorders the circuit before each sweep so every
<=3-wire block becomes contiguous -- a valid reordering of the circuit's
dependency DAG (unitary-preserving, scripts/check_dag_compact.py), exposing
windows to the exhaustive sweep that were previously only found by chance
via the stochastic transport_shuffle. Compares baseline / prepass / dag /
prepass+dag on identical random circuits (Table 6/7 style) under a fixed
time budget, with per-circuit exactness checks. Writes
results/dag_compact/comparison_dag_report.{md,csv,json}.

Usage:
    PYTHONPATH=src python scripts/benchmark_dag_compact.py --gateset nisq --num-circuits 30 --budget 30
    PYTHONPATH=src python scripts/benchmark_dag_compact.py --gateset ion_trap --num-circuits 30 --budget 30
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import time
from pathlib import Path

from qcr_repro.circuits import count_gates, random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.prepass import apply_prepass
from qcr_repro.reducer import reduce_circuit
from qcr_repro.unitary import equivalent_up_to_global_phase

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "dag_compact"

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}

# Paper Table 6 ("Ours", ion trap) / Table 7 (NISQ): per-type means.
PAPER_OURS = {
    "ion_trap": {"RX": 10, "RY": 29, "RZ": 29, "RXX": 43, "total": 111},
    "nisq": {"RX": 45, "RZ": 19, "CZ": 43, "total": 107},
}

METHODS = ["baseline", "prepass", "dag", "prepass_dag"]


def run_one(gates, num_qubits, db, budget, seed, method, gate_set, angles, two):
    t0 = time.perf_counter()
    use_prepass = method in ("prepass", "prepass_dag")
    use_dag = method in ("dag", "prepass_dag")
    if use_prepass:
        pre, prepass_removed = apply_prepass(list(gates), gate_set, angles, two, num_qubits, zx=True)
    else:
        pre, prepass_removed = list(gates), 0
    out, passes, replacements = reduce_circuit(
        pre, num_qubits, db, budget, seed, dag_compact=use_dag
    )
    elapsed = time.perf_counter() - t0
    u0 = circuit_unitary(num_qubits, gates)
    ok = equivalent_up_to_global_phase(u0, circuit_unitary(num_qubits, out), atol=1e-5)
    return {
        "method": method,
        "input_len": len(gates),
        "prepass_len": len(pre),
        "prepass_removed": prepass_removed,
        "final_len": len(out),
        "passes": passes,
        "replacements": replacements,
        "runtime_s": elapsed,
        "exact_ok": ok,
        "counts": count_gates(out),
        "seed": seed,
    }


def summarize(rows):
    final = [r["final_len"] for r in rows]
    totals = [sum(r["counts"].values()) for r in rows]
    twq = [sum(v for k, v in r["counts"].items() if k not in ("RX", "RY", "RZ")) for r in rows]
    return {
        "final_len_mean": statistics.mean(final),
        "final_len_std": statistics.stdev(final) if len(final) > 1 else 0.0,
        "total_mean": statistics.mean(totals),
        "twq_mean": statistics.mean(twq),
        "runtime_mean": statistics.mean(r["runtime_s"] for r in rows),
        "prepass_removed_mean": statistics.mean(r["prepass_removed"] for r in rows),
        "prepass_len_mean": statistics.mean(r["prepass_len"] for r in rows),
        "exact_ok": all(r["exact_ok"] for r in rows),
    }


def _t_test_vs_paper(final_lens: list[float], paper_total: int) -> tuple[float, float]:
    """One-sample t-test mean(final_lens) vs paper_total; returns (t, approx two-sided p via normal approx)."""
    n = len(final_lens)
    if n < 2:
        return float("nan"), float("nan")
    mean = statistics.mean(final_lens)
    sd = statistics.stdev(final_lens)
    if sd == 0:
        return float("inf") if mean != paper_total else 0.0, 0.0
    t = (mean - paper_total) / (sd / (n**0.5))
    # normal approximation to the two-sided p-value (n is large enough here)
    from math import erf, sqrt

    p = 2 * (1 - 0.5 * (1 + erf(abs(t) / sqrt(2))))
    return t, p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--num-circuits", type=int, default=30)
    ap.add_argument("--budget", type=float, default=30.0)
    ap.add_argument("--length", type=int, default=300)
    ap.add_argument("--gateset", default="ion_trap", choices=["ion_trap", "nisq"])
    ap.add_argument("--seed-start", type=int, default=1)
    ap.add_argument("--methods", default=None,
                     help="comma-separated subset of " + ",".join(METHODS) + " (default: all)")
    args = ap.parse_args()
    if args.methods:
        selected = args.methods.split(",")
        for m in selected:
            if m not in METHODS:
                raise SystemExit(f"unknown method {m!r}, choose from {METHODS}")
        METHODS[:] = selected

    gate_set = args.gateset
    depths = ION_DEPTHS if gate_set == "ion_trap" else NISQ_DEPTHS
    weights = None if gate_set == "ion_trap" else {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}
    num_qubits = 4

    print(f"== building/loading {gate_set} DB ==")
    db = load_or_build_database(gate_set, depths, verbose=True)
    angles = db.angles
    two = db.two_qubit_angles

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"{gate_set}_b{int(args.budget)}"
    csv_path = OUT_DIR / f"comparison_dag_{tag}.csv"
    all_rows = []
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["seed", "method", "input_len", "prepass_len", "prepass_removed",
                         "final_len", "passes", "replacements", "runtime_s", "exact_ok",
                         "counts"])
        for seed in range(args.seed_start, args.seed_start + args.num_circuits):
            gates, _ = random_circuit(num_qubits, args.length, gate_set, seed=seed, weights=weights)
            for method in METHODS:
                row = run_one(gates, num_qubits, db, args.budget, seed, method, gate_set, angles, two)
                all_rows.append(row)
                writer.writerow([seed, method, row["input_len"], row["prepass_len"],
                                 row["prepass_removed"], row["final_len"], row["passes"],
                                 row["replacements"], f"{row['runtime_s']:.3f}", int(row["exact_ok"]),
                                 json.dumps(row["counts"])])
                print(f"  seed{seed:>4} {method:<12} {row['input_len']}->{row['prepass_len']}->{row['final_len']} "
                      f"({row['runtime_s']:.1f}s, ok={int(row['exact_ok'])})", flush=True)

    by_method = {m: summarize([r for r in all_rows if r["method"] == m]) for m in METHODS}
    paper = PAPER_OURS[gate_set]
    paper_total = paper["total"]
    paper_twq = paper.get("RXX", paper.get("CZ", 0))

    md = []
    md.append("# DAG block-compaction comparison report\n")
    md.append(f"- Gate set: `{gate_set}` (paper Table {'6' if gate_set == 'ion_trap' else '7'}, "
              f"{num_qubits} qubits, length {args.length})")
    md.append(f"- Circuits per method: {args.num_circuits} (seeds {args.seed_start}-"
              f"{args.seed_start + args.num_circuits - 1}), per-circuit budget: {args.budget}s")
    md.append(f"- Verifier: numeric 1e-5 (input vs output unitary up to global phase)")
    md.append(f"- Paper 'Ours' (mean over 100 runs): total {paper_total}, two-qubit {paper_twq}\n")

    md.append("| method | input | after prepass | final total (mean +/- std) | two-qubit | vs paper | t-stat | p | runtime (s) | exact |")
    md.append("|---|---|---|---|---|---|---|---|---|---|")
    labels = {"baseline": "baseline", "prepass": "prepass", "dag": "**dag_compact**",
              "prepass_dag": "**prepass+dag_compact**"}
    for method in METHODS:
        s = by_method[method]
        rows = [r["final_len"] for r in all_rows if r["method"] == method]
        t, p = _t_test_vs_paper(rows, paper_total)
        vs_paper = (s["final_len_mean"] - paper_total) / paper_total * 100
        verdict = "WIN" if vs_paper < 0 else "LOSE"
        md.append(
            f"| {labels[method]} | {args.length} | {s['prepass_len_mean']:.0f} | "
            f"{s['final_len_mean']:.1f} (+/- {s['final_len_std']:.1f}) | {s['twq_mean']:.1f} | "
            f"{verdict} ({vs_paper:+.1f}%) | {t:.1f} | {p:.2g} | "
            f"{s['runtime_mean']:.1f} | {'OK' if s['exact_ok'] else 'FAIL'} |"
        )

    base_mean = by_method["baseline"]["final_len_mean"]
    dag_mean = by_method["prepass_dag"]["final_len_mean"]
    delta_pct = (dag_mean - base_mean) / base_mean * 100
    md.append(f"\n`prepass+dag_compact` vs `baseline` at equal {args.budget}s budget: "
              f"{dag_mean:.1f} vs {base_mean:.1f} gates ({delta_pct:+.1f}%).\n")

    md.append("## What dag_compact does\n")
    md.append("`src/dag.py` deterministically reorders the circuit so every <=3-wire "
              "block becomes contiguous, using the true per-wire dependency DAG rather "
              "than physical adjacency in the gate list -- generalizing Qiskit-style "
              "2-qubit block collection to k<=3-wire blocks. Any two gates that swap "
              "position under this reordering act on disjoint wires, so it is a valid "
              "topological order and preserves the circuit's unitary exactly "
              "(scripts/check_dag_compact.py). It replaces reliance on the stochastic "
              "transport_shuffle/shuffle_commuting_pairs passes for exposing reducible "
              "windows: the previous pipeline could only find a window if random "
              "shuffling happened to bring its gates physically adjacent; dag_compact "
              "finds every such window in one deterministic pass.")

    report_md = "\n".join(md)
    (OUT_DIR / f"comparison_dag_{tag}_report.md").write_text(report_md)
    (OUT_DIR / f"comparison_dag_{tag}_report.json").write_text(json.dumps(
        {"by_method": by_method, "paper": paper, "args": vars(args)}, indent=2))
    print("\nReport written to", OUT_DIR / f"comparison_dag_{tag}_report.md")
    print(report_md)


if __name__ == "__main__":
    main()
