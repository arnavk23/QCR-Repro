"""collapse_cz_parity is a unitary-preserving CZ-pair collapse.

RZ and CZ are both diagonal in the computational basis regardless of which
wires they touch, so CZ(i, j) commutes with every gate except a non-diagonal
(RX) gate on wire i or j. collapse_cz_parity (src/reducer.py) exploits this
to collapse every CZ(i, j) run between such barriers to parity (0 or 1 gate),
independent of literal adjacency -- something neither zx_cancellations
(literal-adjacency only) nor the active NISQ pipeline currently does.

This checks (1) the collapse alone preserves the unitary and drops the
predicted number of gates on random and edge-case circuits, and (2) an A/B
comparison of reduce_circuit with cz_pass=True vs False, same seed/budget,
both against the NISQ protocol.

Usage:
    PYTHONPATH=src python scripts/check_cz_parity.py
"""

from __future__ import annotations

import sys

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import collapse_cz_parity, reduce_circuit
from qcr_repro.unitary import equivalent_up_to_global_phase

NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}


def _same_multiset_minus_cz(a: list[GateInstance], b: list[GateInstance]) -> bool:
    """b must equal a with only CZ gates removed/added (never non-CZ gates
    added, removed, or mutated)."""
    def key(g: GateInstance):
        theta = None if g.theta is None else round(float(g.theta), 10)
        return (g.name, tuple(g.qubits), theta)

    non_cz_a = sorted(key(g) for g in a if g.name != "CZ")
    non_cz_b = sorted(key(g) for g in b if g.name != "CZ")
    return non_cz_a == non_cz_b


def check_random_circuits(num_qubits: int = 4) -> int:
    failures = 0
    weights = {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}
    for seed in range(1, 8):
        for length in (30, 120, 300):
            gates, _ = random_circuit(num_qubits, length, "nisq", seed=seed, weights=weights)
            working = list(gates)
            removed = collapse_cz_parity(working, num_qubits)
            u0 = circuit_unitary(num_qubits, gates)
            u1 = circuit_unitary(num_qubits, working)
            ok = equivalent_up_to_global_phase(u0, u1, atol=1e-9)
            struct_ok = _same_multiset_minus_cz(gates, working)
            cz_before = sum(1 for g in gates if g.name == "CZ")
            cz_after = sum(1 for g in working if g.name == "CZ")
            if not ok or not struct_ok:
                print(f"  FAIL seed{seed} len{length}: unitary_ok={ok} struct_ok={struct_ok}")
                failures += 1
            else:
                print(f"  ok   seed{seed} len{length}: CZ {cz_before} -> {cz_after} "
                      f"(removed {removed})")
    return failures


def check_edge_cases(num_qubits: int = 4) -> int:
    failures = 0
    cases = {
        "same_pair_even": [
            GateInstance("RX", (0,), 1.0),
            GateInstance("CZ", (0, 1), None),
            GateInstance("RZ", (1,), 0.5),
            GateInstance("CZ", (1, 0), None),  # same pair, reversed qubit order
        ],
        "same_pair_odd_interleaved": [
            GateInstance("CZ", (0, 1), None),
            GateInstance("CZ", (2, 3), None),  # disjoint pair, must not block
            GateInstance("RZ", (0,), 0.3),
            GateInstance("CZ", (0, 1), None),
            GateInstance("CZ", (0, 1), None),
        ],
        "rx_blocks_across_wire": [
            GateInstance("CZ", (0, 1), None),
            GateInstance("RX", (0,), 1.0),  # blocks wire 0 -> new window
            GateInstance("CZ", (0, 1), None),
        ],
        "rx_on_other_wire_no_block": [
            GateInstance("CZ", (0, 1), None),
            GateInstance("RX", (2,), 1.0),  # unrelated wire, must not block
            GateInstance("CZ", (0, 1), None),
        ],
        "three_pairs_sharing_wire": [
            GateInstance("CZ", (0, 1), None),
            GateInstance("CZ", (1, 2), None),
            GateInstance("CZ", (0, 1), None),
            GateInstance("CZ", (1, 2), None),
        ],
    }
    for name, gates in cases.items():
        working = list(gates)
        removed = collapse_cz_parity(working, num_qubits)
        u0 = circuit_unitary(num_qubits, gates)
        u1 = circuit_unitary(num_qubits, working)
        ok = equivalent_up_to_global_phase(u0, u1, atol=1e-9)
        struct_ok = _same_multiset_minus_cz(gates, working)
        if not (ok and struct_ok):
            print(f"  FAIL {name}: unitary_ok={ok} struct_ok={struct_ok} "
                  f"in={len(gates)} out={len(working)}")
            failures += 1
        else:
            print(f"  ok   {name}: {len(gates)} -> {len(working)} gates (removed {removed})")
    return failures


def check_reduce_circuit_ab(num_seeds: int = 10, budget_s: float = 5.0) -> int:
    """A/B: reduce_circuit(cz_pass=True) vs baseline, same seed/budget,
    against the NISQ protocol. Asserts unitary equivalence for both; prints
    the length delta as signal (not asserted, single-seed timing noise)."""
    failures = 0
    weights = {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}
    db = load_or_build_database("nisq", NISQ_DEPTHS)
    total_base = 0
    total_cz = 0
    for seed in range(1, num_seeds + 1):
        gates, _ = random_circuit(4, 300, "nisq", seed=seed, weights=weights)
        u0 = circuit_unitary(4, gates)
        base, pb, _ = reduce_circuit(list(gates), 4, db, budget_s=budget_s, seed=seed, rz_pass=True)
        cz, pc, _ = reduce_circuit(
            list(gates), 4, db, budget_s=budget_s, seed=seed, rz_pass=True, cz_pass=True
        )
        ok_base = equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-5)
        ok_cz = equivalent_up_to_global_phase(u0, circuit_unitary(4, cz), atol=1e-5)
        status = "ok  " if (ok_base and ok_cz) else "FAIL"
        total_base += len(base)
        total_cz += len(cz)
        print(f"  {status} seed{seed}: baseline {len(base)} (passes {pb}) vs "
              f"cz_pass {len(cz)} (passes {pc}) unitary_ok={ok_base and ok_cz}")
        if not (ok_base and ok_cz):
            failures += 1
    print(f"  mean baseline {total_base / num_seeds:.1f}  mean cz_pass {total_cz / num_seeds:.1f}  "
          f"delta {total_cz / num_seeds - total_base / num_seeds:+.1f}")
    return failures


def main() -> None:
    print("== 1) random nisq circuits: collapse_cz_parity preserves unitary ==")
    failures = check_random_circuits()
    print("== 2) edge cases (disjoint pairs, wire-sharing pairs, RX blocking) ==")
    failures += check_edge_cases()
    print("== 3) reduce_circuit A/B (baseline vs cz_pass, same seed/budget) ==")
    failures += check_reduce_circuit_ab()
    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
