"""Parallel best-of-N circuit reduction benchmark.

Each worker process reduces an independent random circuit under a fixed time
budget, so N seeds produce N independent reductions whose best/mean/median give
a robust picture of the achievable end-length (paper Table 6/7 style).

Usage:
    python scripts/bench_best_of.py [num_qubits] [length] [budget_s]
        [--gateset ion_trap|nisq] [--workers N] [--seeds 1 2 3 4 5 6 7 8]
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import statistics
import time

import numpy as np

from qcr_repro.circuits import count_gates, random_circuit
from qcr_repro.compute_graph import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.unitary_utils import equivalent_up_to_global_phase

from bench_ion import reduce_circuit

DEPTHS = {
    "ion_trap": {1: 12, 2: 10, 3: 7, 4: 5},
    "nisq": {1: 12, 2: 5, 3: 4, 4: 4},
}


def worker(args):
    num_qubits, length, budget_s, seed, gate_set, db, weights, rz_pass, prefer = args
    gates, _ = random_circuit(num_qubits, length, gate_set, seed=seed, weights=weights)
    u0 = circuit_unitary(num_qubits, gates)
    t0 = time.time()
    r, passes, reduced = reduce_circuit(
        gates, num_qubits, db, budget_s, seed, prefer=prefer, rz_pass=rz_pass
    )
    dt = time.time() - t0
    ok = equivalent_up_to_global_phase(u0, circuit_unitary(num_qubits, r), atol=1e-5)
    return {
        "seed": seed,
        "end": len(r),
        "counts": count_gates(r),
        "passes": passes,
        "reds": reduced,
        "ok": ok,
        "secs": dt,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("num_qubits", type=int, default=4, nargs="?")
    parser.add_argument("length", type=int, default=300, nargs="?")
    parser.add_argument("budget_s", type=float, default=20.0, nargs="?")
    parser.add_argument("--gateset", type=str, default="ion_trap", choices=["ion_trap", "nisq"])
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--seeds", type=int, nargs="*", default=[1, 2, 3, 4, 5, 6, 7, 8])
    parser.add_argument("--weights", type=str, default="", help="Per-gate sampling weights, e.g. RX:4,RZ:4,CZ:3")
    parser.add_argument("--rz-pass", action="store_true", help="Enable the RZ-across-CZ transport pass (NISQ)")
    parser.add_argument("--prefer", type=str, default="", help="Escape preference weights, e.g. CZ:1.0")
    args = parser.parse_args()

    workers = args.workers or min(len(args.seeds), os.cpu_count() or 4)
    depths = DEPTHS[args.gateset]
    db = load_or_build_database(args.gateset, depths, verbose=False)

    weights = None
    if args.weights:
        weights = {}
        for item in args.weights.split(","):
            name, value = item.split(":")
            weights[name.strip()] = float(value)

    prefer = None
    if args.prefer:
        prefer = {}
        for item in args.prefer.split(","):
            name, value = item.split(":")
            prefer[name.strip()] = float(value)

    tasks = [
        (args.num_qubits, args.length, args.budget_s, s, args.gateset, db, weights, args.rz_pass, prefer)
        for s in args.seeds
    ]
    t0 = time.time()
    with mp.get_context("fork").Pool(processes=workers) as pool:
        results = pool.map(worker, tasks)
    wall = time.time() - t0

    ends = [r["end"] for r in results]
    best = min(results, key=lambda r: r["end"])
    print("\n".join(
        f"seed{r['seed']:>2}: {args.length}->{r['end']:>3} "
        f"{r['counts']} ok={r['ok']} ({r['secs']:.1f}s)"
        for r in results
    ))
    print("-" * 72)
    print(f"best   : {args.length}->{best['end']}  {best['counts']}")
    print(f"mean   : {statistics.mean(ends):.1f}")
    print(f"median : {statistics.median(ends):.0f}")
    print(f"wall   : {wall:.1f}s across {workers} workers, {len(results)} seeds [{args.gateset}]")


if __name__ == "__main__":
    main()
