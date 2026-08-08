"""Head-to-head benchmark: numeric vs exact (length) vs exact (cost-aware).

Compares the three reducers on random circuits matching the composition of
Rosenhahn et al. (New J. Phys. 27, 104509, 2025) Tables 3/4/5, reporting mean
gate counts, per-type counts, runtimes, and equivalence rates.

Usage:
    python scripts/benchmark_exact.py --gateset ion_trap --num-qubits 3 \
        --num-circuits 16 --length 300 --budget 10 --outdir results_exact
    python scripts/benchmark_exact.py --gateset nisq_clifford --num-qubits 3 \
        --num-circuits 16 --length 300 --budget 10 --outdir results_exact
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
from qcr_repro.exact_graph import load_or_build_exact
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.exact_reducer import reduce_circuit_exact, verify_exact
from qcr_repro.unitary_utils import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 6, 4: 4}
NISQ_ANGLES = (-3.141592653589793 / 2, 3.141592653589793 / 2)  # Clifford subset


def worker(args):
    (method, gate_set, num_qubits, length, budget_s, seed, weights, rz_pass, db_paths) = args
    # "nisq_clifford" = the NISQ pool restricted to the Clifford angles +-pi/2.
    gen_set = "nisq" if gate_set == "nisq_clifford" else gate_set
    gen_angles = db_paths.get("angles") if gate_set == "nisq_clifford" else None
    gates, _ = random_circuit(num_qubits, length, gen_set, seed=seed, weights=weights, angles=gen_angles)
    t0 = time.time()
    if method == "numeric":
        db = load_or_build_database(gate_set, db_paths["depths"], db_paths.get("angles"))
        r, _, _ = reduce_circuit(gates, num_qubits, db, budget_s, seed, rz_pass=rz_pass)
        ok = equivalent_up_to_global_phase(circuit_unitary(num_qubits, gates), circuit_unitary(num_qubits, r), atol=1e-5)
        verifier = "numeric-1e-5"
    elif method == "exact_len":
        db = load_or_build_exact(gate_set, db_paths["depths"], db_paths.get("angles"))
        r, stats = reduce_circuit_exact(gates, num_qubits, db, budget_s, seed, cost_aware=False, rz_pass=rz_pass)
        ok = verify_exact(gates, r, num_qubits)
        verifier = "exact"
    else:
        db = load_or_build_exact(gate_set, db_paths["depths"], db_paths.get("angles"))
        r, stats = reduce_circuit_exact(gates, num_qubits, db, budget_s, seed, cost_aware=True, rz_pass=rz_pass)
        ok = verify_exact(gates, r, num_qubits)
        verifier = "exact"
    dt = time.time() - t0
    return {
        "seed": seed,
        "method": method,
        "start": len(gates),
        "end": len(r),
        "counts": count_gates(r),
        "twq": sum(1 for g in r if len(g.qubits) == 2),
        "secs": dt,
        "ok": ok,
        "verifier": verifier,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gateset", type=str, default="ion_trap", choices=["ion_trap", "nisq", "nisq_clifford"])
    parser.add_argument("--num-qubits", type=int, default=3)
    parser.add_argument("--num-circuits", type=int, default=16)
    parser.add_argument("--length", type=int, default=300)
    parser.add_argument("--budget", type=float, default=10.0)
    parser.add_argument("--seed-base", type=int, default=0)
    parser.add_argument("--weights", type=str, default="")
    parser.add_argument("--rz-pass", action="store_true")
    parser.add_argument("--methods", type=str, default="numeric,exact_len,exact_cost")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--outdir", type=str, default="results/costaware_quick")
    args = parser.parse_args()

    weights = None
    if args.weights:
        weights = {}
        for item in args.weights.split(","):
            name, value = item.split(":")
            weights[name.strip()] = float(value)

    if args.gateset == "ion_trap":
        db_paths = {"depths": ION_DEPTHS}
    elif args.gateset == "nisq_clifford":
        db_paths = {"depths": NISQ_DEPTHS, "angles": NISQ_ANGLES}
    else:
        db_paths = {"depths": NISQ_DEPTHS}

    methods = [m.strip() for m in args.methods.split(",")]
    seeds = [args.seed_base + s for s in range(args.num_circuits)]
    tasks = []
    for m in methods:
        for s in seeds:
            tasks.append((m, args.gateset, args.num_qubits, args.length, args.budget, s, weights, args.rz_pass, db_paths))

    workers = args.workers or min(len(tasks), os.cpu_count() or 4)
    t0 = time.time()
    with mp.get_context("fork").Pool(processes=workers) as pool:
        results = pool.map(worker, tasks)
    wall = time.time() - t0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / f"exact_{args.gateset}_q{args.num_qubits}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["method", "seed", "start", "end", "runtime_s", "ok", "verifier", "twq",
                         "RX", "RY", "RZ", "RXX", "CZ"])
        for r in results:
            c = r["counts"]
            writer.writerow([r["method"], r["seed"], r["start"], r["end"], round(r["secs"], 3),
                             r["ok"], r["verifier"], r["twq"],
                             c.get("RX", 0), c.get("RY", 0), c.get("RZ", 0), c.get("RXX", 0), c.get("CZ", 0)])

    print(f"[{args.gateset} q{args.num_qubits}] {len(seeds)} circuits x {args.length} gates, budget {args.budget}s, "
          f"wall {wall:.1f}s across {workers} workers")
    for method in methods:
        rs = [r for r in results if r["method"] == method]
        ends = [r["end"] for r in rs]
        secs = [r["secs"] for r in rs]
        twqs = [r["twq"] for r in rs]
        ok_rate = sum(r["ok"] for r in rs) / len(rs)
        print(f"  {method:<10} end mean {statistics.mean(ends):6.1f}  median {statistics.median(ends):4.0f}  "
              f"best {min(ends):3d}  RXX mean {statistics.mean(twqs):5.1f}  time mean {statistics.mean(secs):5.1f}s  "
              f"ok {ok_rate:.3f} ({rs[0]['verifier']})")
    print(f"  saved CSV: {csv_path}")


if __name__ == "__main__":
    main()
