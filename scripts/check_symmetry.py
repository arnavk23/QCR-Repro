"""SymmetricDatabase is a unitary-preserving, gate-set-agnostic lookup upgrade.

Every ReductionDatabase.try_reduce embeds a block's touched wires into local
graph indices in exactly one way (ascending numeric order). SymmetricDatabase
(src/symmetry.py) tries every one of the k! relabelings against the same
graph and keeps the best hit, decoded back through that relabeling's own
inverse. This checks (1) it never changes the block's unitary, on both
ion_trap and nisq, and (2) an A/B of reduce_circuit with vs without it, same
seed/budget, on both gate sets -- the point being a technique that has to
show up (or honestly not) on both, not one tuned to a single gate set's
algebra.

Usage:
    PYTHONPATH=src python scripts/check_symmetry.py
"""

from __future__ import annotations

import sys
import time

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.symmetry import SymmetricDatabase
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}

CASES = [
    # rz_pass is a NISQ-only pass (RZ-across-CZ gathering): it assumes RX is
    # the only non-diagonal gate, which is false for ion_trap's RY/RXX, so it
    # must stay off there -- matching the paper's own protocol
    # (rz_pass = args.rz_pass or args.gateset == "nisq" in benchmark_comparison.py).
    ("ion_trap", None, ION_DEPTHS, False),
    ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}, NISQ_DEPTHS, True),
]


def check_try_reduce_unitary_preserving() -> int:
    failures = 0
    for gate_set, weights, depths, rz_pass in CASES:
        db = SymmetricDatabase(load_or_build_database(gate_set, depths))
        for seed in range(1, 6):
            gates, _ = random_circuit(4, 150, gate_set, seed=seed, weights=weights)
            # Probe try_reduce directly on every window the sweep would see,
            # so a failure points at the lookup itself rather than the whole
            # search loop.
            n = len(gates)
            checked = 0
            for start in range(0, n, 7):
                for length in range(2, min(8, n - start) + 1):
                    block = gates[start : start + length]
                    u0 = circuit_unitary(4, block)
                    cand = db.try_reduce(block)
                    checked += 1
                    if cand is None:
                        continue
                    u1 = circuit_unitary(4, cand)
                    if not equivalent_up_to_global_phase(u0, u1, atol=1e-8):
                        print(f"  FAIL {gate_set} seed{seed} start{start} len{length}: "
                              f"try_reduce changed the unitary")
                        failures += 1
            print(f"  ok   {gate_set} seed{seed}: {checked} windows probed, "
                  f"identity_hits={db.identity_hits} permuted_hits={db.permuted_hits}")
    return failures


def check_reduce_circuit_ab(num_seeds: int = 8, budget_s: float = 6.0) -> int:
    failures = 0
    for gate_set, weights, depths, rz_pass in CASES:
        base_db = load_or_build_database(gate_set, depths)
        sym_db = SymmetricDatabase(load_or_build_database(gate_set, depths))
        total_base = 0
        total_sym = 0
        t_base = 0.0
        t_sym = 0.0
        print(f"-- {gate_set} --")
        for seed in range(1, num_seeds + 1):
            gates, _ = random_circuit(4, 300, gate_set, seed=seed, weights=weights)
            u0 = circuit_unitary(4, gates)

            t0 = time.time()
            base, pb, _ = reduce_circuit(
                list(gates), 4, base_db, budget_s=budget_s, seed=seed, rz_pass=rz_pass
            )
            t_base += time.time() - t0

            t0 = time.time()
            sym, ps, _ = reduce_circuit(
                list(gates), 4, sym_db, budget_s=budget_s, seed=seed, rz_pass=rz_pass
            )
            t_sym += time.time() - t0

            ok_base = equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-5)
            ok_sym = equivalent_up_to_global_phase(u0, circuit_unitary(4, sym), atol=1e-5)
            status = "ok  " if (ok_base and ok_sym) else "FAIL"
            total_base += len(base)
            total_sym += len(sym)
            print(f"  {status} seed{seed}: baseline {len(base)} (passes {pb}, ok={ok_base}) vs "
                  f"symmetric {len(sym)} (passes {ps}, ok={ok_sym})")
            if not (ok_base and ok_sym):
                failures += 1
        print(f"  mean baseline {total_base / num_seeds:.1f}  mean symmetric {total_sym / num_seeds:.1f}  "
              f"delta {total_sym / num_seeds - total_base / num_seeds:+.1f}  "
              f"identity_hits={sym_db.identity_hits} permuted_hits={sym_db.permuted_hits}  "
              f"wall-clock base {t_base:.1f}s sym {t_sym:.1f}s")
    return failures


def main() -> None:
    print("== 1) try_reduce is unitary-preserving on every probed window (ion_trap + nisq) ==")
    failures = check_try_reduce_unitary_preserving()
    print("== 2) reduce_circuit A/B (baseline vs SymmetricDatabase, same seed/budget) ==")
    failures += check_reduce_circuit_ab()
    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
