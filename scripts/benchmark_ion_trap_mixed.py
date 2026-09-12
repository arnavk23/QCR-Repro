"""Benchmark on the mixed Clifford/non-Clifford ion-trap pool (RX/RY/RZ at
+/-pi/2, +/-pi/4, +/-pi/8; RXX fixed at the native +pi/2 entangler).

Unlike benchmark_comparison.py, there is no published paper baseline for
this gate set -- it's an extension beyond Rosenhahn et al.'s original two
protocols, so this compares our hybrid (exact-engine-for-Clifford-windows,
numeric-fallback-elsewhere) reducer directly against qiskit and BQSKit.

Usage:
    PYTHONPATH=src python scripts/benchmark_ion_trap_mixed.py --num-circuits 20 --budget 15
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
from qcr_repro.exact_database import load_or_build_exact
from qcr_repro.gates import circuit_unitary
from qcr_repro.hybrid import HybridDatabase
from qcr_repro.reducer import reduce_circuit
from qcr_repro.unitary import equivalent_up_to_global_phase

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_comparison import _bqskit_compile, _qiskit_transpile, _run_bqskit_tasks  # noqa: E402

# Calibrated against this pool's own branching factor (18/37/57/78 tokens at
# k=1..4, vs. plain ion_trap's 6/13/21/30): its BFS explodes much faster per
# depth, so depths here are shallower than ION_DEPTHS despite the smaller
# pool count difference looking modest (scripts/... calibration, see paper
# Section on the mixed ion-trap extension).
MIXED_DEPTHS = {1: 6, 2: 4, 3: 4, 4: 3}
# The Clifford sub-pool is structurally identical to plain ion_trap (same
# gate names, same +/-pi/2-only angles), so it reuses ion_trap's own proven
# depths.
CLIFFORD_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}

GATESET = "ion_trap_mixed"


def _count_twq(gates) -> int:
    return sum(1 for g in gates if len(g.qubits) == 2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--num-circuits", type=int, default=20)
    ap.add_argument("--length", type=int, default=300)
    ap.add_argument("--num-qubits", type=int, default=4)
    ap.add_argument("--budget", type=float, default=15.0)
    ap.add_argument("--seed-base", type=int, default=0)
    ap.add_argument("--bqskit-levels", type=str, default="2,3")
    ap.add_argument("--no-bqskit", action="store_true")
    ap.add_argument("--outdir", type=str, default="results/comparison_ion_trap_mixed")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[{GATESET}] building numeric database {MIXED_DEPTHS} ...", flush=True)
    t0 = time.time()
    numeric_db = load_or_build_database(GATESET, MIXED_DEPTHS)
    print(f"[{GATESET}] numeric database ready ({time.time() - t0:.1f}s)", flush=True)

    print(f"[{GATESET}] building exact Clifford sub-pool database {CLIFFORD_DEPTHS} ...", flush=True)
    t0 = time.time()
    exact_db = load_or_build_exact(f"{GATESET}_clifford", CLIFFORD_DEPTHS)
    print(f"[{GATESET}] exact database ready ({time.time() - t0:.1f}s)", flush=True)

    hybrid_db = HybridDatabase(numeric_db, exact_db)

    seeds = [args.seed_base + s for s in range(args.num_circuits)]
    rows = []

    for method in ("hybrid_len", "hybrid_cost"):
        cost_aware = method == "hybrid_cost"
        for seed in seeds:
            gates, _ = random_circuit(args.num_qubits, args.length, GATESET, seed=seed)
            u0 = circuit_unitary(args.num_qubits, gates)
            t0 = time.time()
            r, _passes, _reds = reduce_circuit(
                list(gates), args.num_qubits, hybrid_db, args.budget, seed,
                cost_aware=cost_aware, use_fast_sweep=True,
            )
            secs = time.time() - t0
            ok = equivalent_up_to_global_phase(u0, circuit_unitary(args.num_qubits, r), atol=1e-5)
            rows.append({
                "method": method, "seed": seed, "start": len(gates), "end": len(r),
                "counts": count_gates(r), "twq": _count_twq(r), "secs": secs, "ok": ok,
                "verifier": "numeric-1e-5",
            })
        print(f"[{GATESET}] {method}: {len(seeds)} circuits done "
              f"(exact_lookups={hybrid_db.exact_lookups}, numeric_lookups={hybrid_db.numeric_lookups})",
              flush=True)

    qiskit_ok = True
    for level in (1, 2, 3):
        for seed in seeds:
            gates, _ = random_circuit(args.num_qubits, args.length, GATESET, seed=seed)
            res = _qiskit_transpile(gates, args.num_qubits, GATESET, level)
            if res is None:
                qiskit_ok = False
                continue
            counts, total, ok = res
            rows.append({
                "method": f"qiskit_l{level}", "seed": seed, "start": len(gates), "end": total,
                "counts": counts, "twq": counts.get("RXX", 0), "secs": 0.0, "ok": ok,
                "verifier": "qiskit-1e-5",
            })
    print(f"[{GATESET}] qiskit: {'ok' if qiskit_ok else 'some circuits unavailable/failed'}", flush=True)

    if not args.no_bqskit:
        bqskit_levels = [int(x) for x in args.bqskit_levels.split(",") if x.strip()]
        bqskit_tasks = [
            (GATESET, f"bqskit_l{level}", args.num_qubits, args.length, args.budget, seed,
             None, False, {}, 8, 1, "ram", False, False, False)
            for level in bqskit_levels for seed in seeds
        ]
        print(f"[{GATESET}] running {len(bqskit_tasks)} bqskit task(s) sequentially...", flush=True)
        # _run_bqskit_tasks calls random_circuit(gateset=...) internally using
        # task[0]; passing "ion_trap" there (not "ion_trap_mixed") would draw
        # the wrong (Clifford-only) angle set, so patch gateset through as-is
        # and let _bqskit_compile's basis selection (fixed to cover
        # ion_trap_mixed) handle the rest.
        rows.extend(_run_bqskit_tasks(bqskit_tasks))

    csv_path = outdir / "comparison_ion_trap_mixed.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "seed", "start", "end", "runtime_s", "ok", "verifier", "twq",
                          "RX", "RY", "RZ", "RXX", "CZ"])
        for r in sorted(rows, key=lambda r: (r["method"], r["seed"])):
            c = r["counts"]
            writer.writerow([r["method"], r["seed"], r["start"], r["end"], round(r["secs"], 3),
                              r["ok"], r["verifier"], r["twq"],
                              c.get("RX", 0), c.get("RY", 0), c.get("RZ", 0), c.get("RXX", 0), c.get("CZ", 0)])

    methods = sorted({r["method"] for r in rows})
    summary = {}
    print(f"\n[{GATESET}] summary (n={args.num_circuits}, budget={args.budget}s):", flush=True)
    for m in methods:
        mrows = [r for r in rows if r["method"] == m and r["end"] >= 0]
        if not mrows:
            continue
        ends = [r["end"] for r in mrows]
        twqs = [r["twq"] for r in mrows]
        ok_rate = sum(r["ok"] for r in mrows) / len(mrows)
        summary[m] = {
            "n": len(mrows), "end_mean": statistics.mean(ends), "end_std": statistics.pstdev(ends),
            "twq_mean": statistics.mean(twqs), "twq_std": statistics.pstdev(twqs), "ok_rate": ok_rate,
        }
        print(f"  {m:<14} end {summary[m]['end_mean']:7.1f} (+- {summary[m]['end_std']:5.1f})  "
              f"twq {summary[m]['twq_mean']:6.1f} (+- {summary[m]['twq_std']:5.1f})  "
              f"ok {ok_rate:.3f}  n={summary[m]['n']}", flush=True)

    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nsaved CSV: {csv_path}")
    print(f"saved summary: {outdir / 'summary.json'}")


if __name__ == "__main__":
    main()
