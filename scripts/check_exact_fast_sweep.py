"""Exact-engine analogue of check_fast_sweep*.py: verifies
_sweep_reduce_len_fast / _sweep_reduce_cost_fast (exact_reducer.py) against
the originals on the ion_trap (Clifford) pool.

Checks:
  1) per-window parity: for many explicit windows, does the incrementally
     built tableau find a candidate at the same length (len-only) and the
     same Pareto-optimal (twq, len) as SymplecticGraph.try_reduce /
     try_reduce_cost via the original from-scratch rebuild?
  2) reduce_circuit_exact A/B under a fixed time budget: does the
     incremental sweep's throughput gain turn into better reduction, like
     it did for the numeric engine?

Usage:
    PYTHONPATH=src python scripts/check_exact_fast_sweep.py
"""

from __future__ import annotations

import sys
import time

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.exact_database import load_or_build_exact
from qcr_repro.exact_reducer import (
    _build_incremental_tableau_candidates,
    reduce_circuit_exact,
    verify_exact,
)
from qcr_repro.symplectic import circuits_equal_exact

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}

_PREFER = {"RX": 0.01, "RY": 0.01, "RZ": 0.01, "CZ": 5.0, "RXX": 5.0}


def check_parity(max_block_len: int = 8) -> int:
    failures = 0
    db = load_or_build_exact("ion_trap", ION_DEPTHS)
    total_windows = 0
    for seed in range(1, 6):
        gates, _ = random_circuit(4, 150, "ion_trap", seed=seed)
        n = len(gates)
        for pos in range(0, n - 1, 3):
            hi = min(max_block_len, n - pos)
            window = gates[pos : pos + hi]
            results, _lwm = _build_incremental_tableau_candidates(window, db, max_block_len)
            for length in range(2, hi + 1):
                total_windows += 1
                block = gates[pos : pos + length]
                wires = sorted({q for g in block for q in g.qubits})
                graph = db.graphs.get(len(wires))
                remapped, _, _ = db._remap(block)
                slow_key = graph.block_key(remapped) if graph is not None else None
                slow_chain = graph.buckets.get(slow_key) if (graph is not None and slow_key is not None) else None
                slow_len = len(slow_chain) if slow_chain is not None else None
                slow_alts = sorted(graph.alts.get(slow_key, [])) if (graph is not None and slow_key is not None) else []

                entry = results.get(length)
                if entry is None:
                    fast_len = None
                    fast_alts = []
                else:
                    fast_chain, fast_alts_raw, _g = entry
                    fast_len = len(fast_chain) if fast_chain is not None else None
                    fast_alts = sorted(fast_alts_raw)

                # Different (but equally valid) local wire embeddings can
                # legitimately realize a window with different specific
                # token chains -- only the achievable (twq, len) Pareto
                # frontier is guaranteed invariant (established for the
                # numeric engine in check_symmetry.py; same argument
                # applies here since the exact graph is built by exploring
                # every pool gate on every local wire/pair uniformly too).
                slow_pareto = sorted({(twq, ln) for twq, ln, _chain in slow_alts})
                fast_pareto = sorted({(twq, ln) for twq, ln, _chain in fast_alts})
                if slow_len != fast_len or slow_pareto != fast_pareto:
                    print(f"  MISMATCH seed{seed} pos{pos} length{length}: "
                          f"slow_len={slow_len} fast_len={fast_len} "
                          f"slow_pareto={slow_pareto} fast_pareto={fast_pareto}")
                    failures += 1
    print(f"  ion_trap: {total_windows} windows probed, {failures} mismatches")
    return failures


def check_ab(num_seeds: int = 10, budget_s: float = 8.0) -> int:
    failures = 0
    db = load_or_build_exact("ion_trap", ION_DEPTHS)
    total_base = 0
    total_fast = 0
    t_base = 0.0
    t_fast = 0.0
    print("-- ion_trap (exact engine, cost_aware) --")
    for seed in range(1, num_seeds + 1):
        gates, _ = random_circuit(4, 300, "ion_trap", seed=seed)

        t0 = time.time()
        base, stats_b = reduce_circuit_exact(
            list(gates), 4, db, budget_s=budget_s, seed=seed, cost_aware=True, prefer=_PREFER
        )
        t_base += time.time() - t0

        t0 = time.time()
        fast, stats_f = reduce_circuit_exact(
            list(gates), 4, db, budget_s=budget_s, seed=seed, cost_aware=True, prefer=_PREFER,
            use_fast_sweep=True,
        )
        t_fast += time.time() - t0

        ok_base = verify_exact(gates, base, 4)
        ok_fast = verify_exact(gates, fast, 4)
        total_base += len(base)
        total_fast += len(fast)
        status = "ok  " if (ok_base and ok_fast) else "FAIL"
        print(f"  {status} seed{seed}: baseline {len(base)} (passes {stats_b.iterations}, ok={ok_base}) vs "
              f"fast_sweep {len(fast)} (passes {stats_f.iterations}, ok={ok_fast})")
        if not (ok_base and ok_fast):
            failures += 1
    print(f"  mean baseline {total_base / num_seeds:.1f}  mean fast_sweep {total_fast / num_seeds:.1f}  "
          f"delta {total_fast / num_seeds - total_base / num_seeds:+.1f}  "
          f"wall-clock base {t_base:.1f}s fast {t_fast:.1f}s")
    return failures


def main() -> None:
    print("== 1) per-window parity: incremental tableau vs from-scratch rebuild (ion_trap) ==")
    failures = check_parity()
    print("== 2) reduce_circuit_exact A/B: baseline vs use_fast_sweep, same seed/budget ==")
    failures += check_ab()
    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
