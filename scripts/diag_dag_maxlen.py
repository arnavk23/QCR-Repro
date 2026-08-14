"""Does dag_compact need a larger max_block_len to pay off?  compact_by_blocks 
produces blocks bigger than the sweep's default 8-gate window cap, so the sweep only
able to chew through a fraction of each block.  Compares baseline vs
dag_compact at max_block_len in {8, 12, 16, 24} on a few circuits.

Usage: PYTHONPATH=src python scripts/diag_dag_maxlen.py --gateset ion_trap
"""
from __future__ import annotations

import argparse
import time

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gateset", default="ion_trap", choices=["ion_trap", "nisq"])
    ap.add_argument("--budget", type=float, default=10.0)
    ap.add_argument("--seeds", type=int, default=5)
    args = ap.parse_args()

    depths = ION_DEPTHS if args.gateset == "ion_trap" else NISQ_DEPTHS
    weights = None if args.gateset == "ion_trap" else {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}
    db = load_or_build_database(args.gateset, depths, verbose=False)

    for max_block_len in (8, 12, 16, 24):
        baseline_lens = []
        dag_lens = []
        for seed in range(1, args.seeds + 1):
            gates, _ = random_circuit(4, 300, args.gateset, seed=seed, weights=weights)
            u0 = circuit_unitary(4, gates)
            b, _, _ = reduce_circuit(list(gates), 4, db, args.budget, seed, max_block_len=max_block_len)
            d, _, _ = reduce_circuit(list(gates), 4, db, args.budget, seed,
                                      max_block_len=max_block_len, dag_compact=True)
            assert equivalent_up_to_global_phase(u0, circuit_unitary(4, b), atol=1e-5)
            assert equivalent_up_to_global_phase(u0, circuit_unitary(4, d), atol=1e-5)
            baseline_lens.append(len(b))
            dag_lens.append(len(d))
        bm = sum(baseline_lens) / len(baseline_lens)
        dm = sum(dag_lens) / len(dag_lens)
        print(f"max_block_len={max_block_len:>3}: baseline={bm:.1f} {baseline_lens}  "
              f"dag={dm:.1f} {dag_lens}  delta={(dm - bm) / bm * 100:+.1f}%")


if __name__ == "__main__":
    main()
