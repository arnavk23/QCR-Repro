"""Diagnostic: how often does NISQ's try_reduce_cost actually find a
genuinely different (lower-twq) candidate than try_reduce, versus just
returning the same chain? If it's rarely different, the cost-aware
objective is effectively dormant for NISQ, which would explain why
numeric_cost isn't reliably beating numeric_len on two-qubit count.

Usage:
    PYTHONPATH=src python scripts/diag_nisq_cost_aware.py
"""

from __future__ import annotations

from qcr_repro.circuits import random_circuit
from qcr_repro.config import GateInstance
from qcr_repro.database import load_or_build_database

NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}
WEIGHTS = {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}


def count_twq(gates):
    return sum(1 for g in gates if len(g.qubits) == 2)


def main() -> None:
    db = load_or_build_database("nisq", NISQ_DEPTHS)
    total_windows = 0
    both_none = 0
    only_len_hits = 0
    only_cost_hits = 0
    same_result = 0
    cost_shorter = 0
    cost_fewer_twq_same_len = 0
    cost_worse_or_equal = 0
    alts_nonempty_count = 0
    alts_len_hist = {}

    for seed in range(1, 6):
        gates, _ = random_circuit(4, 200, "nisq", seed=seed, weights=WEIGHTS)
        n = len(gates)
        for pos in range(0, n - 1, 2):
            for length in range(2, min(8, n - pos) + 1):
                total_windows += 1
                block = gates[pos : pos + length]
                raw_len_cand = db.try_reduce(block)
                # try_reduce returns its raw match unconditionally (even
                # same-length or longer); _sweep_reduce only ever *applies*
                # it when strictly shorter. Apply that same filter here so
                # the comparison to try_reduce_cost (which already filters
                # internally) is apples-to-apples.
                len_cand = raw_len_cand if (raw_len_cand is not None and len(raw_len_cand) < length) else None
                cost_cand = db.try_reduce_cost(block)

                # peek at the raw alts list this window's key has, regardless
                # of acceptance criteria, to see how much material is even
                # available to be cost-aware about.
                wires = sorted({q for g in block for q in g.qubits})
                graph = db.graphs.get(len(wires))
                if graph is not None:
                    forward = {w: i for i, w in enumerate(wires)}
                    remapped = [
                        GateInstance(name=g.name, qubits=tuple(sorted(forward[q] for q in g.qubits)), theta=g.theta)
                        for g in block
                    ]
                    u = graph.block_unitary(remapped)
                    key = graph._lookup_key(u) if u is not None else None
                    if key is not None:
                        alts = graph._alt_lists().get(key, ())
                        if len(alts) > 0:
                            alts_nonempty_count += 1
                        alts_len_hist[len(alts)] = alts_len_hist.get(len(alts), 0) + 1

                if len_cand is None and cost_cand is None:
                    both_none += 1
                    continue
                if len_cand is not None and cost_cand is None:
                    only_len_hits += 1
                    continue
                if len_cand is None and cost_cand is not None:
                    only_cost_hits += 1
                    continue

                len_len, len_twq = len(len_cand), count_twq(len_cand)
                cost_len, cost_twq = len(cost_cand), count_twq(cost_cand)
                if len_len == cost_len and len_twq == cost_twq:
                    same_result += 1
                elif cost_len < len_len:
                    cost_shorter += 1
                elif cost_len == len_len and cost_twq < len_twq:
                    cost_fewer_twq_same_len += 1
                else:
                    cost_worse_or_equal += 1

    print(f"total windows probed: {total_windows}")
    print(f"  both None (no reduction found either way): {both_none}")
    print(f"  only try_reduce hit: {only_len_hits}")
    print(f"  only try_reduce_cost hit: {only_cost_hits}")
    print(f"  identical result (len,twq): {same_result}")
    print(f"  cost strictly shorter: {cost_shorter}")
    print(f"  cost same length, fewer twq: {cost_fewer_twq_same_len}")
    print(f"  cost worse/equal-but-not-identical (unexpected): {cost_worse_or_equal}")
    print(f"alts list nonempty for {alts_nonempty_count} windows; length histogram: {sorted(alts_len_hist.items())}")


if __name__ == "__main__":
    main()
