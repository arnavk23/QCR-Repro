"""Ion-trap / NISQ circuit reduction benchmark (paper Table 6/7 style).

Usage:
    python scripts/benchmark_basic.py [num_qubits] [length] [budget_s] [seeds...]

Reduction pipeline (implemented in qcr_repro.reducer.reduce_circuit):
  1. Initial exhaustive sweep (uses the wire-count DB).
  2. Single-qubit clustering + sweep.
  3. Transport shuffle + sweep, alternating direction, with periodic
     restart-from-best when stuck, and equivalence-class escape moves that
     resample irreducible windows with structurally different factorizations.
"""

from __future__ import annotations

import sys

from qcr_repro.circuits import count_gates, random_circuit
from qcr_repro.compute_graph import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.unitary_utils import equivalent_up_to_global_phase

DEPTHS = {1: 10, 2: 8, 3: 5, 4: 4}


def main() -> None:
    args = [a for a in sys.argv[1:]]
    num_qubits = int(args[0]) if len(args) > 0 else 4
    length = int(args[1]) if len(args) > 1 else 300
    budget_s = float(args[2]) if len(args) > 2 else 60.0
    seeds = [int(a) for a in args[3:]] or [1, 2, 3]

    db = load_or_build_database("ion_trap", DEPTHS, verbose=False)
    for seed in seeds:
        gates, _ = random_circuit(num_qubits, length, "ion_trap", seed=seed)
        u0 = circuit_unitary(num_qubits, gates)
        r, passes, reduced = reduce_circuit(gates, num_qubits, db, budget_s, seed)
        ok = equivalent_up_to_global_phase(u0, circuit_unitary(num_qubits, r), atol=1e-5)
        print(
            f"seed{seed}: {length}->{len(r)} {count_gates(r)} "
            f"passes={passes} reds={reduced} ok={ok}",
            flush=True,
        )


if __name__ == "__main__":
    main()
