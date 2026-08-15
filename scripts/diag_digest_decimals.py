"""Diagnostic (not part of the test suite): does the compute graph's node
key precision (digest_decimals, src/database.py) cause spurious node
fragmentation on NISQ?

The paper (Rosenhahn/Osborne/Hirche, sec 2.2) merges duplicate compute-graph
nodes with a *tolerance* of 1e-5 on the unitary comparison. Our ComputeGraph
instead keys nodes by SHA-256 of the phase-normalized unitary rounded to
digest_decimals=10 -- far tighter than 1e-5. If floating-point noise on
deeper chains pushes two genuinely-equal (up to global phase) unitaries to
round differently at the 10th decimal, they'd get spuriously split into two
graph nodes, silently reducing effective database coverage -- a candidate
explanation for the NISQ "database factorization coverage" gap (report
sec 2.2). This script checks whether coarsening digest_decimals measurably
changes node counts (evidence for fragmentation) at a small, fast-to-build
depth, without touching the correctness-critical lookup path.

Usage: PYTHONPATH=src python scripts/diag_digest_decimals.py
"""
from __future__ import annotations

from qcr_repro.database import ComputeGraph
from qcr_repro.token_pool import TokenPool

NISQ_ANGLES = (-1.5707963267948966, -0.7853981633974483, 0.7853981633974483, 1.5707963267948966)


def main() -> None:
    for wires in (2, 3):
        pool = TokenPool(num_qubits=wires, gate_set="nisq", angles=NISQ_ANGLES)
        for depth in (3, 4):
            counts = {}
            for decimals in (10, 8, 6, 5, 4):
                g = ComputeGraph(pool=pool, max_depth=depth, digest_decimals=decimals)
                counts[decimals] = g.num_nodes
            base = counts[10]
            print(f"wires={wires} depth={depth} pool_size={len(pool.tokens())}: " +
                  "  ".join(f"d{d}={n}({(n - base) / base * 100:+.2f}%)" for d, n in counts.items()))


if __name__ == "__main__":
    main()
