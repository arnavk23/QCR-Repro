"""Per-window parity check: for many explicit (pos, length) windows, does
_sweep_reduce_fast's internal candidate length always match what
ReductionDatabase.try_reduce finds for the exact same block? If this holds
universally, the final-length differences check_fast_sweep.py reports are
benign greedy-path divergence (both paths are locally optimal at every
single step, they just sometimes discover reductions in a different order,
which is already a known sensitivity of this search -- not something the
fast path introduces). If it does NOT hold, there's a real bug.

Usage:
    PYTHONPATH=src python scripts/check_fast_sweep_parity.py
"""

from __future__ import annotations

import sys

import numpy as np

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.database import load_or_build_database
from qcr_repro.reducer import _fast_gate_matrix

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}

CASES = [
    ("ion_trap", None, ION_DEPTHS),
    ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}, NISQ_DEPTHS),
]


def fast_length_for_window(gates: list[GateInstance], db, max_block_len: int) -> dict[int, int | None]:
    """Replicates _sweep_reduce_fast's per-length candidate-length
    computation for a single fixed window starting at index 0, without
    applying any replacement -- pure probe."""
    hi = min(max_block_len, len(gates))
    local_wire_map: dict[int, int] = {}
    u: np.ndarray | None = None
    lengths: dict[int, int | None] = {}
    for length in range(1, hi + 1):
        gate = gates[length - 1]
        new_wires = [w for w in gate.qubits if w not in local_wire_map]
        for w in new_wires:
            local_wire_map[w] = len(local_wire_map)
        k = len(local_wire_map)
        if u is None:
            u = np.eye(2**k, dtype=complex)
        elif new_wires:
            u = np.kron(u, np.eye(2 ** len(new_wires), dtype=complex))
        graph = db.graphs.get(k)
        local_gate = GateInstance(
            name=gate.name,
            qubits=tuple(sorted(local_wire_map[w] for w in gate.qubits)),
            theta=gate.theta,
        )
        u = _fast_gate_matrix(graph, k, local_gate) @ u
        if length < 2:
            continue
        if graph is None:
            lengths[length] = None
            continue
        chain = graph.lookup(u)
        lengths[length] = len(chain) if chain is not None else None
    return lengths


def check_parity(max_block_len: int = 8) -> int:
    failures = 0
    total_windows = 0
    for gate_set, weights, depths in CASES:
        db = load_or_build_database(gate_set, depths)
        for seed in range(1, 6):
            gates, _ = random_circuit(4, 150, gate_set, seed=seed, weights=weights)
            n = len(gates)
            for pos in range(0, n - 1, 3):
                hi = min(max_block_len, n - pos)
                fast_lengths = fast_length_for_window(gates[pos : pos + hi], db, max_block_len)
                for length in range(2, hi + 1):
                    total_windows += 1
                    block = gates[pos : pos + length]
                    slow_cand = db.try_reduce(block)
                    slow_len = len(slow_cand) if slow_cand is not None else None
                    fast_len = fast_lengths.get(length)
                    if slow_len != fast_len:
                        print(f"  MISMATCH {gate_set} seed{seed} pos{pos} length{length}: "
                              f"slow={slow_len} fast={fast_len}")
                        failures += 1
        print(f"  {gate_set}: {total_windows} windows probed so far, {failures} mismatches")
    return failures


def main() -> None:
    print("== per-window length parity: db.try_reduce vs _sweep_reduce_fast's internal lookup ==")
    failures = check_parity()
    print(f"\n{'ALL CHECKS PASSED (fast always matches slow length, per window)' if failures == 0 else f'{failures} MISMATCH(ES)'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
