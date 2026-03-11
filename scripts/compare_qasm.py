from __future__ import annotations

import argparse

from qcr_repro.gates import circuit_unitary
from qcr_repro.qasm_io import parse_qasm_subset
from qcr_repro.unitary_utils import equivalent_up_to_global_phase


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare two QASM files up to global phase.")
    parser.add_argument("--a", required=True)
    parser.add_argument("--b", required=True)
    parser.add_argument("--atol", type=float, default=1e-5)
    args = parser.parse_args()

    na, ga = parse_qasm_subset(args.a)
    nb, gb = parse_qasm_subset(args.b)
    if na != nb:
        raise ValueError("Qubit count mismatch between input files")

    ua = circuit_unitary(na, ga)
    ub = circuit_unitary(nb, gb)

    ok = equivalent_up_to_global_phase(ua, ub, atol=args.atol)
    print(f"A gates: {len(ga)}")
    print(f"B gates: {len(gb)}")
    print(f"Equivalent (atol={args.atol}): {ok}")


if __name__ == "__main__":
    main()
