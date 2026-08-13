"""Benchmark: pre-passes + batched sweep vs the baseline reducer and the paper's reported numbers.

Compares baseline / prepass / prepass_batched pipelines on identical random ion-trap circuits (Table 6 style) under a fixed time budget, with per-circuit exactness checks. Writes results/prepass/comparison_prepass_report.{md,csv,json}.

Usage:
    PYTHONPATH=src python scripts/benchmark_prepass.py [--num-circuits N] [--budget SEC] [--length L] [--gateset ion_trap|nisq]"""

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
OUT_DIR = ROOT / "results" / "prepass"

DEPTHS = {1: 10, 2: 8, 3: 5, 4: 4}
NISQ_DEPTHS = {1: 10, 2: 5, 3: 4, 4: 3}

# Paper Table 6 ("Ours", ion trap) / Table 7 (NISQ): per-type means.
PAPER_OURS = {
    "ion_trap": {"RX": 10, "RY": 29, "RZ": 29, "RXX": 43, "total": 111},
    "nisq": {"RX": 45, "RZ": 19, "CZ": 43, "total": 107},
}
# Paper Table 2: mean computation time (s) of V1/V2/V3 on 100 length-100 circuits.
PAPER_TIMES = {"V1-random": 199.0, "V2-db": 55.0, "V3-rf": 38.0}


def run_one(gates, num_qubits, db, budget, seed, method, gate_set, angles, two):
    """Run one pipeline on one circuit; return a dict of results."""
    t0 = time.perf_counter()
    if method == "baseline":
        out, passes, replacements = reduce_circuit(
            list(gates), num_qubits, db, budget, seed
        )
        prepass_removed = 0
        prepass_len = len(gates)
    else:
        pre, prepass_removed = apply_prepass(list(gates), gate_set, angles, two, num_qubits, zx=True)
        prepass_len = len(pre)
        kwargs = {"algebraic": True, "zx": True}
        if method == "prepass_batched":
            kwargs["use_batched"] = True
        out, passes, replacements = reduce_circuit(
            pre, num_qubits, db, budget, seed, **kwargs
        )
    elapsed = time.perf_counter() - t0
    u0 = circuit_unitary(num_qubits, gates)
    ok = equivalent_up_to_global_phase(u0, circuit_unitary(num_qubits, out), atol=1e-5)
    return {
        "method": method,
        "input_len": len(gates),
        "prepass_len": prepass_len,
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--num-circuits", type=int, default=12)
    ap.add_argument("--budget", type=float, default=10.0)
    ap.add_argument("--length", type=int, default=300)
    ap.add_argument("--gateset", default="ion_trap", choices=["ion_trap", "nisq"])
    args = ap.parse_args()

    gate_set = args.gateset
    depths = DEPTHS if gate_set == "ion_trap" else NISQ_DEPTHS
    weights = None if gate_set == "ion_trap" else {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}
    num_qubits = 4

    print(f"== building/loading {gate_set} DB ==")
    db = load_or_build_database(gate_set, depths, verbose=False)
    gs = db.graphs
    angles = db.angles
    two = db.two_qubit_angles

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "comparison_prepass.csv"
    all_rows = []
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["seed", "method", "input_len", "prepass_len", "prepass_removed",
                         "final_len", "passes", "replacements", "runtime_s", "exact_ok",
                         "counts"])
        for seed in range(1, args.num_circuits + 1):
            gates, _ = random_circuit(num_qubits, args.length, gate_set, seed=seed, weights=weights)
            for method in ("baseline", "prepass", "prepass_batched"):
                row = run_one(gates, num_qubits, db, args.budget, seed, method, gate_set, angles, two)
                all_rows.append(row)
                writer.writerow([seed, method, row["input_len"], row["prepass_len"],
                                 row["prepass_removed"], row["final_len"], row["passes"],
                                 row["replacements"], f"{row['runtime_s']:.3f}", int(row["exact_ok"]),
                                 json.dumps(row["counts"])])
                print(f"  seed{seed:>3} {method:<16} {row['input_len']}->{row['prepass_len']}->{row['final_len']} "
                      f"({row['runtime_s']:.1f}s, ok={int(row['exact_ok'])})", flush=True)

    # Aggregate
    by_method = {m: summarize([r for r in all_rows if r["method"] == m]) for m in
                 ("baseline", "prepass", "prepass_batched")}
    paper = PAPER_OURS[gate_set]
    paper_total = paper["total"]
    paper_twq = paper.get("RXX", paper.get("CZ", 0))

    md = []
    md.append("# Pre-pass + batched sweep comparison report\n")
    md.append(f"- Gate set: `{gate_set}` (paper Table {'6' if gate_set == 'ion_trap' else '7'}, "
              f"{num_qubits} qubits, length {args.length})")
    md.append(f"- Circuits per method: {args.num_circuits}, per-circuit budget: {args.budget}s")
    md.append(f"- Exacter: numeric 1e-5 (input vs output unitary up to global phase)")
    md.append(f"- Paper 'Ours' (mean over 100 runs): total {paper_total}, "
              f"two-qubit {paper_twq}\n")

    md.append("| method | input | after prepass | final total | two-qubit | vs paper Ours | prepass removed | runtime (s) | exact |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    for method, label in (("baseline", "**baseline**"), ("prepass", "**prepass**"),
                          ("prepass_batched", "**prepass+batched**")):
        s = by_method[method]
        vs_paper = (s["final_len_mean"] - paper_total) / paper_total * 100
        verdict = "WIN" if vs_paper < 0 else "LOSE"
        twq_s = s["twq_mean"]
        md.append(
            f"| {label} | {args.length} | {s['prepass_len_mean']:.0f} | "
            f"{s['final_len_mean']:.1f} (+/- {s['final_len_std']:.1f}) | {twq_s:.1f} | "
            f"{verdict} ({vs_paper:+.1f}%) | {s['prepass_removed_mean']:.1f} | "
            f"{s['runtime_mean']:.1f} | {'OK' if s['exact_ok'] else 'FAIL'} |"
        )

    md.append("\n## Effect of the pre-passes\n")
    md.append("- `algebraic_merge`: fuses adjacent same-axis rotations whose sum snaps to a pool angle;"
              "drops exact cancellations (e.g. RZ(a)RZ(-a) = I).")
    md.append("- `zx_cancellations`: drops adjacent same-pair CZ CZ = I; for the diagonal-CZ pools, gathers "
              "RZ gates across CZ so runs become adjacent and fuse.")
    md.append("- `apply_prepass` is applied to fixpoint *before* the database loop; the output stays in the "
              "discrete pool, so the DB loop is unaffected except that the input is shorter.\n")

    md.append("## Wall-clock vs the paper (contextual)\n")
    md.append("- Paper Table 2 (reported, 100 length-100 circuits, 3-4 qubit): "
              + ", ".join(f"{k} {v}s" for k, v in PAPER_TIMES.items()) + ".")
    md.append("- Our pipelines above ran on identical random circuits of length "
              f"{args.length} on this machine; the same budget was given to every method, "
              "so the comparison of interest is the final length at fixed budget.\n")

    md.append("## Batched sweep vs scalar sweep (sweep-only microbenchmark)\n")
    md.append("- The vectorized batched sweep is bit-identical to the scalar sweep "
              "(verified in scripts/check_batched_vs_scalar.py) and measures "
              "~1.5-1.7x faster on the length-300 ion-trap fixpoint in our microbenchmark; "
              "the remaining per-window cost is the SHA-256 digest.")

    report_md = "\n".join(md)
    (OUT_DIR / "comparison_prepass_report.md").write_text(report_md)
    (OUT_DIR / "comparison_prepass_report.json").write_text(json.dumps(
        {"by_method": by_method, "paper": paper, "paper_times": PAPER_TIMES,
         "args": vars(args)}, indent=2))
    print("\nReport written to", OUT_DIR / "comparison_prepass_report.md")
    print(report_md)


if __name__ == "__main__":
    main()
