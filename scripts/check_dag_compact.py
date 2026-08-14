"""compact_by_blocks is a unitary-preserving reordering.

collect_wire_blocks/compact_by_blocks (src/dag.py) reorder a circuit so every
<=max_wires-wire block becomes contiguous, without touching gate content.
This asserts the reordering is a valid topological order of the circuit's
true per-wire dependency DAG -- i.e. it never moves a gate across another
gate it shares a wire with -- by checking the resulting unitary exactly
matches the original circuit's, on many random ion_trap/nisq circuits and
edge cases (single-wire runs, single-gate circuits, empty circuits, wide
multi-qubit spans).

Usage:
    PYTHONPATH=src python scripts/check_dag_compact.py
"""

from __future__ import annotations

import sys

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.dag import block_size_stats, collect_wire_blocks, compact_by_blocks
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}


def _same_multiset(a: list[GateInstance], b: list[GateInstance]) -> bool:
    """True iff a and b contain exactly the same gates (as a multiset) --
    compaction must be a pure reordering, never adding/dropping/mutating."""
    def key(g: GateInstance):
        theta = None if g.theta is None else round(float(g.theta), 10)
        return (g.name, tuple(g.qubits), theta)

    ka = sorted(map(key, a))
    kb = sorted(map(key, b))
    return ka == kb


def check_random_circuits(num_qubits: int = 4, max_wires: int = 3) -> int:
    failures = 0
    cases = [("ion_trap", None), ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0})]
    for gate_set, weights in cases:
        for seed in range(1, 8):
            for length in (30, 120, 300):
                gates, _ = random_circuit(num_qubits, length, gate_set, seed=seed, weights=weights)
                compacted = compact_by_blocks(gates, max_wires=max_wires)
                u0 = circuit_unitary(num_qubits, gates)
                u1 = circuit_unitary(num_qubits, compacted)
                ok = equivalent_up_to_global_phase(u0, u1, atol=1e-9)
                multiset_ok = _same_multiset(gates, compacted)
                if not ok or not multiset_ok:
                    print(f"  FAIL {gate_set} seed{seed} len{length}: "
                          f"unitary_ok={ok} multiset_ok={multiset_ok}")
                    failures += 1
                else:
                    stats = block_size_stats(gates, max_wires=max_wires)
                    print(f"  ok   {gate_set} seed{seed} len{length}: "
                          f"blocks={stats['num_blocks']} max_block={stats['max_block_len']} "
                          f"mean_block={stats['mean_block_len']:.2f}")
    return failures


def check_edge_cases(num_qubits: int = 4) -> int:
    failures = 0
    cases = {
        "empty": [],
        "single_gate": [GateInstance("RX", (0,), 1.0)],
        "single_wire_run": [GateInstance("RX", (0,), 1.0) for _ in range(10)],
        "two_wire_only": [
            GateInstance("RX", (0,), 1.0),
            GateInstance("CZ", (0, 1), None),
            GateInstance("RZ", (1,), 0.5),
            GateInstance("CZ", (0, 1), None),
        ],
        "four_wire_bridge": [
            GateInstance("RX", (0,), 1.0),
            GateInstance("RX", (1,), 1.0),
            GateInstance("RX", (2,), 1.0),
            GateInstance("RX", (3,), 1.0),
            GateInstance("CZ", (0, 1), None),
            GateInstance("CZ", (2, 3), None),
            GateInstance("CZ", (1, 2), None),  # bridges two 2-wire blocks -> 4 wires, must split
            GateInstance("RZ", (0,), 0.3),
            GateInstance("RZ", (3,), 0.3),
        ],
    }
    for name, gates in cases.items():
        for max_wires in (1, 2, 3):
            compacted = compact_by_blocks(gates, max_wires=max_wires)
            if not gates:
                ok = compacted == []
                multiset_ok = True
            else:
                u0 = circuit_unitary(num_qubits, gates)
                u1 = circuit_unitary(num_qubits, compacted)
                ok = equivalent_up_to_global_phase(u0, u1, atol=1e-9)
                multiset_ok = _same_multiset(gates, compacted)
            blocks = collect_wire_blocks(gates, max_wires=max_wires)
            # A block can't be narrower than its widest single gate, so the
            # effective cap is max(max_wires, that gate's own arity).
            widths_ok = all(
                len({q for i in members for q in gates[i].qubits})
                <= max(max_wires, max(len(gates[i].qubits) for i in members))
                for members in blocks
            )
            if not (ok and multiset_ok and widths_ok):
                print(f"  FAIL {name} max_wires={max_wires}: unitary_ok={ok} "
                      f"multiset_ok={multiset_ok} widths_ok={widths_ok}")
                failures += 1
            else:
                print(f"  ok   {name} max_wires={max_wires}: {len(blocks)} block(s)")
    return failures


def check_reduce_circuit_parity() -> int:
    """reduce_circuit(dag_compact=True) preserves the unitary at equal
    seed/budget.  Length is only reported, not asserted: both paths run
    under a fixed *wall-clock* budget, so which one ends up a gate or two
    shorter on a single short run is timing noise, not signal -- see
    scripts/benchmark_dag_compact.py for the statistically-powered
    (many-seed, longer-budget) comparison that actually settles that
    question."""
    failures = 0
    cases = [("ion_trap", None, load_or_build_database("ion_trap", ION_DEPTHS)),
              ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}, load_or_build_database("nisq", NISQ_DEPTHS))]
    for gate_set, weights, db in cases:
        for seed in (1, 2, 3):
            gates, _ = random_circuit(4, 150, gate_set, seed=seed, weights=weights)
            u0 = circuit_unitary(4, gates)
            base, pb, _ = reduce_circuit(list(gates), 4, db, budget_s=4.0, seed=seed)
            comp, pc, _ = reduce_circuit(list(gates), 4, db, budget_s=4.0, seed=seed, dag_compact=True)
            ok_base = equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-5)
            ok_comp = equivalent_up_to_global_phase(u0, circuit_unitary(4, comp), atol=1e-5)
            status = "ok  " if (ok_base and ok_comp) else "FAIL"
            print(f"  {status} reduce_circuit parity {gate_set} seed{seed}: "
                  f"baseline {len(base)} (passes {pb}) vs dag_compact {len(comp)} (passes {pc}) "
                  f"unitary_ok={ok_base and ok_comp}")
            if not (ok_base and ok_comp):
                failures += 1
    return failures


def main() -> None:
    print("== 1) random ion_trap/nisq circuits: compact_by_blocks preserves unitary ==")
    failures = check_random_circuits()
    print("== 2) edge cases (empty, single-wire, bridging blocks) ==")
    failures += check_edge_cases()
    print("== 3) reduce_circuit parity (baseline vs dag_compact, same seed/budget) ==")
    failures += check_reduce_circuit_parity()
    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
