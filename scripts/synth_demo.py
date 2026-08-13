"""Global Clifford synthesis for the QCOptimDemo ion-trap circuits.

Ports each input circuit's binary symplectic tableau (phase-carrying) through
Aaronson-Gottesman row reduction and returns a globally-renormalized circuit,
verifying the result bit-exactly (tableau) and numerically (up to global phase).
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from qcr_repro.clifford import Tableau, circuits_equal
from qcr_repro.config import GateInstance
from qcr_repro.gates import circuit_unitary
from qcr_repro.qasm import parse_qasm_subset, write_qasm_subset
from qcr_repro.synthesis import _SignedTableau, synth_clifford
from qcr_repro.unitary import equivalent_up_to_global_phase

from collections import Counter


def main() -> None:
    parser = argparse.ArgumentParser(description="Global Clifford synthesis for demo circuits.")
    parser.add_argument("--input", required=True, help="Path to input QASM")
    parser.add_argument("--output", required=True, help="Path to output QASM")
    parser.add_argument("--skip-numeric", action="store_true", help="Skip numeric unitary check")
    parser.add_argument("--atol", type=float, default=1e-3,
                        help="Numeric equivalence tolerance (scaled by gate-count FP accumulation)")
    args = parser.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)

    num_qubits, gates = parse_qasm_subset(in_path)
    counts = Counter((g.name, round(g.theta, 4)) for g in gates)
    print(f"Loaded {len(gates)} gates on {num_qubits} qubits from {in_path}")
    print(f"  counts: {dict(sorted(counts.items()))}")

    tab = Tableau(num_qubits)
    tab.apply_circuit(gates)
    signed = _SignedTableau(num_qubits)
    signed.apply_circuit(gates)

    t0 = time.time()
    out = synth_clifford(tab, signed)
    dt = time.time() - t0

    eq_tableau = circuits_equal(gates, out, num_qubits)
    print(f"Synthesized {len(gates)} -> {len(out)} gates in {dt:.2f}s")
    print(f"  tableau bit-exact (up to global phase): {eq_tableau}")
    print(f"  counts: {dict(sorted(Counter((g.name, round(g.theta, 4)) for g in out).items()))}")

    if not eq_tableau:
        raise SystemExit("FAIL: tableau mismatch")

    if not args.skip_numeric:
        u = circuit_unitary(num_qubits, gates)
        v = circuit_unitary(num_qubits, out)
        ok = equivalent_up_to_global_phase(u, v, atol=args.atol)
        print(f"  numeric unitary equiv (up to global phase): {ok}")
        if not ok:
            raise SystemExit("FAIL: numeric mismatch")

    write_qasm_subset(out_path, num_qubits, out)
    print(f"Saved synthesized circuit to {out_path}")


if __name__ == "__main__":
    main()
