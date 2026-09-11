"""Does raising max_alts (currently hardcoded to 4 in ComputeGraph) surface
genuinely twq-improving equal-length alternatives for NISQ that are
currently being evicted by the cap? Builds a 3-wire NISQ graph fresh at
max_alts=4 and max_alts=32 and compares how often try_reduce_cost-style
Pareto filtering would find something try_reduce's raw shortest chain
doesn't already have optimal twq for.

Usage:
    PYTHONPATH=src python scripts/diag_nisq_max_alts.py
"""

from __future__ import annotations

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.database import ComputeGraph
from qcr_repro.token_pool import TokenPool

WEIGHTS = {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}


def count_twq(gates):
    return sum(1 for g in gates if len(g.qubits) == 2)


def probe(graph: ComputeGraph, num_qubits: int = 4, n_seeds: int = 5, length: int = 200) -> None:
    total = 0
    alts_at_cap = 0
    improving = 0
    max_alt_len_seen = 0
    for seed in range(1, n_seeds + 1):
        gates, _ = random_circuit(num_qubits, length, "nisq", seed=seed, weights=WEIGHTS)
        n = len(gates)
        for pos in range(0, n - 1, 2):
            for wlen in range(2, min(8, n - pos) + 1):
                block = gates[pos : pos + wlen]
                wires = sorted({q for g in block for q in g.qubits})
                if len(wires) != graph.pool.num_qubits:
                    continue
                total += 1
                forward = {w: i for i, w in enumerate(wires)}
                remapped = [
                    GateInstance(name=g.name, qubits=tuple(sorted(forward[q] for q in g.qubits)), theta=g.theta)
                    for g in block
                ]
                u = graph.block_unitary(remapped)
                if u is None:
                    continue
                key = graph._lookup_key(u)
                if key is None:
                    continue
                shortest = graph.buckets[key]
                if len(shortest) >= wlen:
                    continue  # not a real reduction opportunity anyway
                alts = graph._alt_lists().get(key, ())
                max_alt_len_seen = max(max_alt_len_seen, len(alts))
                if len(alts) >= graph.max_alts:
                    alts_at_cap += 1
                shortest_twq = sum(1 for tok in shortest if len(graph.pool.gate_for_token(tok).qubits) == 2)
                for chain in alts:
                    if len(chain) == len(shortest):
                        twq = sum(1 for tok in chain if len(graph.pool.gate_for_token(tok).qubits) == 2)
                        if twq < shortest_twq:
                            improving += 1
                            break
    print(f"  max_alts={graph.max_alts}: {total} matching-wire-count reducible windows probed, "
          f"{alts_at_cap} hit the alts cap, {improving} had a same-length lower-twq alt, "
          f"max alt-list length actually seen: {max_alt_len_seen}")


def main() -> None:
    pool = TokenPool(num_qubits=3, gate_set="nisq")
    print("Building 3-wire NISQ graph at depth 5, max_alts=4 (current default)...")
    g4 = ComputeGraph(pool, max_depth=5, max_alts=4)
    probe(g4)

    print("Building 3-wire NISQ graph at depth 5, max_alts=32...")
    pool2 = TokenPool(num_qubits=3, gate_set="nisq")
    g32 = ComputeGraph(pool2, max_depth=5, max_alts=32)
    probe(g32)


if __name__ == "__main__":
    main()
