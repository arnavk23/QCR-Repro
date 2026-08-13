from __future__ import annotations

import argparse
import time
from pathlib import Path

from qcr_repro.gates import circuit_unitary
from qcr_repro.qasm import parse_qasm_subset, write_qasm_subset
from qcr_repro.reducer import reduce_with_lookup
from qcr_repro.unitary import equivalent_up_to_global_phase


def main() -> None:
    parser = argparse.ArgumentParser(description="Ported Python reducer for QCOptimDemo QASM files.")
    parser.add_argument("--input", required=True, help="Path to input QASM (e.g., longcode10.txt)")
    parser.add_argument("--output", required=True, help="Path to output reduced QASM")
    parser.add_argument("--iters", type=int, default=15000, help="Reduction iterations")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    parser.add_argument("--local-qubits", type=int, default=3, help="Local wire budget for lookup")
    parser.add_argument("--graph-depth", type=int, default=4, help="Compute graph depth")
    parser.add_argument("--max-block", type=int, default=7, help="Max sampled block length")
    args = parser.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)

    num_qubits, gates = parse_qasm_subset(in_path)
    print(f"Loaded {len(gates)} gates on {num_qubits} qubits from {in_path}")

    start_u = circuit_unitary(num_qubits, gates)

    t0 = time.time()
    reduced, stats = reduce_with_lookup(
        gates,
        num_qubits=num_qubits,
        local_qubits=args.local_qubits,
        max_block_len=args.max_block,
        graph_depth=args.graph_depth,
        iterations=args.iters,
        seed=args.seed,
    )
    dt = time.time() - t0

    end_u = circuit_unitary(num_qubits, reduced)
    ok = equivalent_up_to_global_phase(start_u, end_u, atol=1e-5)

    write_qasm_subset(out_path, num_qubits, reduced)

    print(f"Reduced {stats.start_len} -> {stats.end_len} gates")
    print(f"Replacements: {stats.replacements}")
    print(f"Runtime: {dt:.2f}s")
    print(f"Unitary preserved (1e-5, up to global phase): {ok}")
    print(f"Saved reduced circuit to {out_path}")


if __name__ == "__main__":
    main()
