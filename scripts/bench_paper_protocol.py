"""Full paper-protocol benchmark (Rosenhahn et al., Table 6/7 style).

Reduces ``--num-circuits`` random 4-qubit circuits of length 300 per gate set
under a fixed time budget and reports mean/std of the total and per-type gate
counts, runtimes, and equivalence pass rates.

Usage:
    python scripts/bench_paper_protocol.py --gateset ion_trap --num-circuits 100
    python scripts/bench_paper_protocol.py --gateset nisq --num-circuits 100 \\
        --weights RX:4,RZ:4,CZ:3 --rz-pass
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import os
import statistics
import time
from pathlib import Path

from qcr_repro.circuits import count_gates, random_circuit
from qcr_repro.compute_graph import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.unitary_utils import equivalent_up_to_global_phase

DEPTHS = {
    "ion_trap": {1: 12, 2: 10, 3: 7, 4: 5},
    "nisq": {1: 12, 2: 5, 3: 4, 4: 4},
}


def circuit_depth(num_qubits: int, gates) -> int:
    """Depth in the target gate set (longest dependency chain)."""
    last = [0] * num_qubits
    for gate in gates:
        d = 1 + max(last[q] for q in gate.qubits)
        for q in gate.qubits:
            last[q] = d
    return max(last) if last else 0


def worker(args):
    (gate_set, num_qubits, length, budget_s, seed, weights, rz_pass, db) = args
    gates, _ = random_circuit(num_qubits, length, gate_set, seed=seed, weights=weights)
    u0 = circuit_unitary(num_qubits, gates)
    t0 = time.time()
    r, passes, reds = reduce_circuit(gates, num_qubits, db, budget_s, seed, rz_pass=rz_pass)
    dt = time.time() - t0
    ok = equivalent_up_to_global_phase(u0, circuit_unitary(num_qubits, r), atol=1e-5)
    return {
        "seed": seed,
        "start": len(gates),
        "end": len(r),
        "counts": count_gates(r),
        "depth_in": circuit_depth(num_qubits, gates),
        "depth_out": circuit_depth(num_qubits, r),
        "passes": passes,
        "reds": reds,
        "ok": ok,
        "secs": dt,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gateset", type=str, required=True, choices=["ion_trap", "nisq"])
    parser.add_argument("--num-circuits", type=int, default=100)
    parser.add_argument("--num-qubits", type=int, default=4)
    parser.add_argument("--length", type=int, default=300)
    parser.add_argument("--budget", type=float, default=30.0)
    parser.add_argument("--seed-base", type=int, default=0)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--weights", type=str, default="", help="e.g. RX:4,RZ:4,CZ:3")
    parser.add_argument("--rz-pass", action="store_true")
    parser.add_argument("--outdir", type=str, default="results_paper")
    args = parser.parse_args()

    weights = None
    if args.weights:
        weights = {}
        for item in args.weights.split(","):
            name, value = item.split(":")
            weights[name.strip()] = float(value)

    workers = args.workers or min(args.num_circuits, os.cpu_count() or 4)
    db = load_or_build_database(args.gateset, DEPTHS[args.gateset], verbose=False)

    seeds = [args.seed_base + s for s in range(args.num_circuits)]
    tasks = [
        (args.gateset, args.num_qubits, args.length, args.budget, s, weights, args.rz_pass, db)
        for s in seeds
    ]

    t0 = time.time()
    with mp.get_context("fork").Pool(processes=workers) as pool:
        results = pool.map(worker, tasks)
    wall = time.time() - t0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / f"paper_{args.gateset}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["seed", "start", "end", "runtime_s", "ok", "depth_in", "depth_out",
                         "RX", "RY", "RZ", "RXX", "CZ", "passes", "reds"])
        for r in results:
            c = r["counts"]
            writer.writerow([r["seed"], r["start"], r["end"], round(r["secs"], 3), r["ok"],
                             r["depth_in"], r["depth_out"], c.get("RX", 0), c.get("RY", 0),
                             c.get("RZ", 0), c.get("RXX", 0), c.get("CZ", 0), r["passes"], r["reds"]])

    ends = [r["end"] for r in results]
    secs = [r["secs"] for r in results]
    ok_rate = sum(r["ok"] for r in results) / len(results)
    gate_names = ["RX", "RY", "RZ", "RXX", "CZ"]
    print(f"[{args.gateset}] {len(results)} circuits x {args.length} gates, budget {args.budget}s, "
          f"wall {wall:.1f}s across {workers} workers")
    print(f"  total gates:  mean {statistics.mean(ends):.1f} +/- {statistics.pstdev(ends):.1f}   "
          f"median {statistics.median(ends):.0f}   best {min(ends)}   eq-rate {ok_rate:.3f}")
    for name in gate_names:
        vals = [r["counts"].get(name, 0) for r in results]
        print(f"  {name:>3}: mean {statistics.mean(vals):.1f} +/- {statistics.pstdev(vals):.1f}")
    print(f"  depth: in {statistics.mean(r['depth_in'] for r in results):.1f} -> "
          f"out {statistics.mean(r['depth_out'] for r in results):.1f}")
    print(f"  runtime: mean {statistics.mean(secs):.1f}s  median {statistics.median(secs):.1f}s")
    print(f"  saved CSV: {csv_path}")


if __name__ == "__main__":
    main()
