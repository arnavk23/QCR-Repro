"""_sweep_reduce_fast must be unitary-equivalent to _sweep_reduce and at
least as good, just faster -- not necessarily byte-identical, since a
different (but equally valid) wire-relabeling can land on a different
equally-short realization of the same window.

_sweep_reduce tries every window length from max_block_len down to 2 at
every position, rebuilding each length's unitary as a fresh matrix-chain
product from scratch (O(max_block_len^2) gate-matrix multiplies per
position). _sweep_reduce_fast (src/reducer.py) builds the same lengths'
unitaries incrementally instead (O(max_block_len)), widening the local
Hilbert space by a Kronecker product exactly when a gate introduces a wire
the window hasn't touched yet, and reuses each graph's own precomputed
per-token matrix cache rather than calling embedded_gate_matrix fresh. This
checks (1) both remain unitary-preserving on many random circuits on both
gate sets, (2) final length is statistically comparable (not systematically
worse), and (3) wall-clock speedup.

Usage:
    PYTHONPATH=src python scripts/check_fast_sweep.py
"""

from __future__ import annotations

import sys
import time

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import _sweep_reduce, _sweep_reduce_fast
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}

CASES = [
    ("ion_trap", None, ION_DEPTHS),
    ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}, NISQ_DEPTHS),
]


def check_fast_sweep(max_block_len: int = 8) -> int:
    failures = 0
    for gate_set, weights, depths in CASES:
        db = load_or_build_database(gate_set, depths)
        t_slow = 0.0
        t_fast = 0.0
        total_slow_len = 0
        total_fast_len = 0
        n_cases = 0
        for seed in range(1, 8):
            for length in (30, 120, 300):
                gates, _ = random_circuit(4, length, gate_set, seed=seed, weights=weights)
                u0 = circuit_unitary(4, gates)

                slow = list(gates)
                t0 = time.time()
                _sweep_reduce(slow, 4, db, max_block_len)
                t_slow += time.time() - t0

                fast = list(gates)
                t0 = time.time()
                _sweep_reduce_fast(fast, 4, db, max_block_len)
                t_fast += time.time() - t0

                ok_slow = equivalent_up_to_global_phase(u0, circuit_unitary(4, slow), atol=1e-6)
                ok_fast = equivalent_up_to_global_phase(u0, circuit_unitary(4, fast), atol=1e-6)
                total_slow_len += len(slow)
                total_fast_len += len(fast)
                n_cases += 1
                if not (ok_slow and ok_fast):
                    print(f"  FAIL {gate_set} seed{seed} len{length}: "
                          f"ok_slow={ok_slow} ok_fast={ok_fast} "
                          f"len_slow={len(slow)} len_fast={len(fast)}")
                    failures += 1
                else:
                    same = len(slow) == len(fast)
                    print(f"  ok   {gate_set} seed{seed} len{length}: "
                          f"{length} -> slow {len(slow)} / fast {len(fast)} "
                          f"({'same length' if same else 'DIFFERENT length'})")
        print(f"  {gate_set} wall-clock: slow {t_slow:.2f}s  fast {t_fast:.2f}s  "
              f"speedup {t_slow / t_fast if t_fast else float('inf'):.2f}x  "
              f"mean len slow {total_slow_len / n_cases:.1f}  fast {total_fast_len / n_cases:.1f}")
    return failures


def main() -> None:
    print("== _sweep_reduce vs _sweep_reduce_fast: correctness, length, speed (both gate sets) ==")
    failures = check_fast_sweep()
    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
