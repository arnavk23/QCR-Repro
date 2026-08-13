"""Correctness checks: the batched sweep is bit-identical to the scalar sweep, and the algebraic/ZX pre-passes preserve the unitary with only pool-representable gates.

Usage:
    PYTHONPATH=src python scripts/check_batched_vs_scalar.py"""

from __future__ import annotations

import sys

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.prepass import apply_prepass
from qcr_repro.reducer import _sweep_reduce, reduce_circuit
from qcr_repro.token_pool import TokenPool
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 10, 2: 8, 3: 5, 4: 4}
NISQ_DEPTHS = {1: 10, 2: 5, 3: 4, 4: 3}


def same_gates(a, b) -> bool:
    if len(a) != len(b):
        return False
    for g1, g2 in zip(a, b):
        if g1.name != g2.name or tuple(g1.qubits) != tuple(g2.qubits):
            return False
        if g1.theta is None or g2.theta is None:
            if g1.theta is not None or g2.theta is not None:
                return False
        elif abs(g1.theta - g2.theta) > 1e-12:
            return False
    return True


def is_pool_representable(pool: TokenPool, gates) -> bool:
    """True iff every gate exists verbatim in the discrete token pool."""
    try:
        pool.encode(list(gates))
        return True
    except KeyError:
        return False


def check_batched_vs_scalar(db_ion, db_nisq) -> int:
    failures = 0
    cases = [("ion_trap", None, db_ion), ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}, db_nisq)]
    for gate_set, weights, db in cases:
        for seed in range(1, 4):
            for length in (60, 120):
                gates, _ = random_circuit(4, length, gate_set, seed=seed, weights=weights)
                scalar = list(gates)
                batched = list(gates)
                n_scalar = _sweep_reduce(scalar, 4, db, 8)
                from qcr_repro.batched import batched_sweep
                n_batched = batched_sweep(batched, 4, db, 8)
                identical = same_gates(scalar, batched)
                if n_scalar != n_batched or not identical:
                    print(f"  FAIL {gate_set} seed{seed} len{length}: "
                          f"replacements scalar={n_scalar} batched={n_batched} identical={identical}")
                    failures += 1
                else:
                    print(f"  ok   {gate_set} seed{seed} len{length}: {length}->{len(scalar)} replacements={n_scalar}")
    return failures


def check_prepass() -> int:
    failures = 0
    for gate_set in ("ion_trap", "nisq"):
        for seed in range(1, 4):
            gates, pool = random_circuit(4, 200, gate_set, seed=seed)
            out, removed = apply_prepass(gates, gate_set, pool.angles, pool.two_qubit_angles, 4, zx=True)
            u0 = circuit_unitary(4, gates)
            u1 = circuit_unitary(4, out)
            ok = equivalent_up_to_global_phase(u0, u1, atol=1e-6)
            representable = is_pool_representable(pool, out)
            if not ok or not representable:
                print(f"  FAIL prepass {gate_set} seed{seed}: removed={removed} unitary_ok={ok} representable={representable}")
                failures += 1
            else:
                print(f"  ok   prepass {gate_set} seed{seed}: {len(gates)}->{len(out)} removed={removed}")
    return failures


def check_reduce_circuit_parity(db) -> int:
    failures = 0
    for seed in (1, 2):
        gates, _ = random_circuit(4, 120, "ion_trap", seed=seed)
        u0 = circuit_unitary(4, gates)
        a, pa, _ = reduce_circuit(list(gates), 4, db, budget_s=4.0, seed=seed)
        b, pb, _ = reduce_circuit(list(gates), 4, db, budget_s=4.0, seed=seed, use_batched=True)
        ok_a = equivalent_up_to_global_phase(u0, circuit_unitary(4, a), atol=1e-5)
        ok_b = equivalent_up_to_global_phase(u0, circuit_unitary(4, b), atol=1e-5)
        no_worse = len(b) <= len(a)
        print(f"  {'ok  ' if (ok_a and ok_b and no_worse) else 'FAIL'} reduce_circuit parity seed{seed}: "
              f"scalar {len(a)} (passes {pa}) vs batched {len(b)} (passes {pb}) "
              f"unitary_ok={ok_a and ok_b} no_worse={no_worse}")
        if not (ok_a and ok_b and no_worse):
            failures += 1
    return failures


def main() -> None:
    print("== building/loading DBs (first run may take a few minutes) ==")
    db_ion = load_or_build_database("ion_trap", ION_DEPTHS, verbose=False)
    db_nisq = load_or_build_database("nisq", NISQ_DEPTHS, verbose=False)

    failures = 0
    print("== 1) batched sweep vs scalar sweep (bit-identical) ==")
    failures += check_batched_vs_scalar(db_ion, db_nisq)
    print("== 2) pre-pass unitary preservation + pool representability ==")
    failures += check_prepass()
    print("== 3) reduce_circuit parity (scalar vs batched, same seed/budget) ==")
    failures += check_reduce_circuit_parity(db_ion)

    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
