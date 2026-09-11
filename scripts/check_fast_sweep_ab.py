"""reduce_circuit(use_fast_sweep=True) A/B against the baseline scalar
sweep, same seed/time budget, on both gate sets -- the real test: does the
per-window throughput gain (scripts/check_fast_sweep.py: 1.27x-1.37x wall
clock at full fixpoint) turn into better reduction under a *fixed* time
budget, unlike every previous lookup-side idea this session (CZ-parity,
SymmetricDatabase, SuffixExtendedDatabase), which all made per-lookup work
more expensive and lost to it.

Usage:
    PYTHONPATH=src python scripts/check_fast_sweep_ab.py
"""

from __future__ import annotations

import sys
import time

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}

# rz_pass is NISQ-only (RZ-across-CZ gathering assumes RX is the only
# non-diagonal gate, false for ion_trap's RY/RXX).
CASES = [
    ("ion_trap", None, ION_DEPTHS, False),
    ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}, NISQ_DEPTHS, True),
]


def check_ab(num_seeds: int = 10, budget_s: float = 6.0) -> int:
    failures = 0
    for gate_set, weights, depths, rz_pass in CASES:
        db = load_or_build_database(gate_set, depths)
        total_base = 0
        total_fast = 0
        t_base = 0.0
        t_fast = 0.0
        print(f"-- {gate_set} --")
        for seed in range(1, num_seeds + 1):
            gates, _ = random_circuit(4, 300, gate_set, seed=seed, weights=weights)
            u0 = circuit_unitary(4, gates)

            t0 = time.time()
            base, pb, _ = reduce_circuit(
                list(gates), 4, db, budget_s=budget_s, seed=seed, rz_pass=rz_pass
            )
            t_base += time.time() - t0

            t0 = time.time()
            fast, pf, _ = reduce_circuit(
                list(gates), 4, db, budget_s=budget_s, seed=seed, rz_pass=rz_pass, use_fast_sweep=True
            )
            t_fast += time.time() - t0

            ok_base = equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-5)
            ok_fast = equivalent_up_to_global_phase(u0, circuit_unitary(4, fast), atol=1e-5)
            total_base += len(base)
            total_fast += len(fast)
            status = "ok  " if (ok_base and ok_fast) else "FAIL"
            print(f"  {status} seed{seed}: baseline {len(base)} (passes {pb}, ok={ok_base}) vs "
                  f"fast_sweep {len(fast)} (passes {pf}, ok={ok_fast})")
            if not (ok_base and ok_fast):
                failures += 1
        print(f"  mean baseline {total_base / num_seeds:.1f}  mean fast_sweep {total_fast / num_seeds:.1f}  "
              f"delta {total_fast / num_seeds - total_base / num_seeds:+.1f}  "
              f"wall-clock base {t_base:.1f}s fast {t_fast:.1f}s")
    return failures


def main() -> None:
    print("== reduce_circuit A/B: baseline scalar sweep vs use_fast_sweep, same seed/budget ==")
    failures = check_ab()
    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
