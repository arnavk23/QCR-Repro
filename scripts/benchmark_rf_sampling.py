"""Benchmark: RF-gated lookup in the paper's random-sampling loop (V2 vs V3).

The paper's V3 gates the *database lookup* of a random sub-block sampler
(Sec. 3.2-3.4), not an exhaustive sweep.  This script compares, at equal
per-circuit time budget on NISQ inputs:
  - sampling      (V2): random sub-block sampling + database lookup
  - sampling_gated (V3): the same loop with an RF-gated lookup (RfGatedDatabase)
  - sweep         : the exhaustive-sweep reducer, as a reference

Usage:
    python scripts/benchmark_rf_sampling.py [--num-circuits 6] [--budget 10.0]
        [--gateset nisq] [--outdir results/rf_sampling]
"""
from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import (
    reduce_circuit,
    reduce_random_sampling,
    reduce_random_sampling_gated,
)
from qcr_repro.unitary import equivalent_up_to_global_phase

NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}
NISQ_WEIGHTS = {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}


def verify(gates, reduced, num_qubits: int) -> bool:
    if len(reduced) >= len(gates):
        return True
    u = circuit_unitary(num_qubits, gates)
    v = circuit_unitary(num_qubits, reduced)
    return bool(equivalent_up_to_global_phase(u, v, atol=1e-5))


def run_one(db, method: str, gates, num_qubits: int, budget: float, seed: int):
    t0 = time.time()
    if method == "sampling":
        reduced, stats = reduce_random_sampling(gates, num_qubits, db, budget, seed)
        gate = None
    elif method == "sampling_gated":
        reduced, stats, gate = reduce_random_sampling_gated(gates, num_qubits, db, budget, seed)
    else:
        reduced, _p, _r = reduce_circuit(gates, num_qubits, db, budget, seed)
        stats = None
        gate = None
    return {
        "method": method,
        "seed": seed,
        "end": len(reduced),
        "replacements": getattr(stats, "replacements", -1),
        "iterations": getattr(stats, "iterations", -1),
        "lookups_skipped": gate.lookups_skipped if gate else -1,
        "lookups_attempted": gate.lookups_attempted if gate else -1,
        "ok": verify(gates, reduced, num_qubits),
        "secs": time.time() - t0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--num-circuits", type=int, default=6)
    ap.add_argument("--budget", type=float, default=10.0)
    ap.add_argument("--seed-base", type=int, default=0)
    ap.add_argument("--outdir", default="results/rf_sampling")
    args = ap.parse_args()

    num_qubits, length, gateset = 4, 300, "nisq"
    print("loading NISQ database (cached build)...", flush=True)
    db = load_or_build_database(gateset, NISQ_DEPTHS, backend="ram")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    methods = ["sampling", "sampling_gated", "sweep"]
    rows = []
    for seed in range(args.seed_base, args.seed_base + args.num_circuits):
        gates, _ = random_circuit(num_qubits, length, gateset, seed=seed, weights=NISQ_WEIGHTS)
        for method in methods:
            rows.append(run_one(db, method, gates, num_qubits, args.budget, seed))
            print(f"seed {seed} {method}: end={rows[-1]['end']} ok={rows[-1]['ok']}", flush=True)

    path = outdir / "comparison_rf_sampling.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {path}")
    print("\nmethod            end(mean)  sd     ok   replacements(mean)")
    for m in methods:
        rs = [r for r in rows if r["method"] == m]
        ends = [r["end"] for r in rs]
        reps = [r["replacements"] for r in rs if r["replacements"] >= 0]
        rep_s = f"{statistics.mean(reps):.1f}" if reps else "-"
        ok = sum(r["ok"] for r in rs)
        print(f"{m:16s} {statistics.mean(ends):8.1f}  {statistics.stdev(ends):5.1f}  {ok}/{len(rs)}  {rep_s}")
    gated = [r for r in rows if r["method"] == "sampling_gated"]
    if gated:
        sk = [r["lookups_skipped"] for r in gated]
        at = [r["lookups_attempted"] for r in gated]
        print(f"\nV3 gate: mean lookups attempted={statistics.mean(at):.0f}, skipped={statistics.mean(sk):.0f}")


if __name__ == "__main__":
    main()
