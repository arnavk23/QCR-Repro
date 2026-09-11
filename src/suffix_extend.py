"""Suffix-extension (bidirectional / meet-in-the-middle) database lookup.

A compute graph only records the shortest chain for unitaries reached by
forward BFS from the identity within ``max_depth``. If a block's unitary T
is not itself in that reachable set, undoing one trailing gate from T --
querying M_c^dagger @ T for every pool token c -- may land back inside the
known set: if the graph has chain R for that residual, then R + (c,)
realizes T using one gate more than the graph's own depth, at the cost of
one extra hash lookup per pool token. Iterating gives depth-k extension at
|pool|^k extra lookups. This is standard bidirectional search over the same
forward-only graph; it needs no change to how the graph is built and no
gate-set-specific algebra, so the same code applies to ion_trap, nisq, and
nisq_clifford alike.
"""

from __future__ import annotations

import numpy as np

from .config import GateInstance
from .gates import embedded_gate_matrix


def _phase_key(graph, matrix: np.ndarray) -> bytes:
    flat = matrix.reshape(-1)
    idx = int(np.argmax(np.round(np.abs(flat), 8)))
    phase = np.angle(flat[idx])
    nflat = flat * np.exp(-1j * phase)
    return graph._node_key(nflat, graph.digest_decimals)


def _token_matrices(graph) -> dict[int, np.ndarray]:
    cached = getattr(graph, "_token_matrices", None)
    if cached:
        return cached
    return {
        tok: embedded_gate_matrix(graph.pool.num_qubits, graph.pool.gate_for_token(tok))
        for tok in graph.pool.tokens()
    }


def extended_lookup(graph, target: np.ndarray, max_extra: int = 1) -> tuple[int, ...] | None:
    """Shortest known chain for ``target``, extending up to ``max_extra``
    trailing gates beyond the graph's direct (depth-limited) reach. None if
    nothing is found even after extension."""
    best = graph.buckets.get(_phase_key(graph, target))
    if max_extra < 1:
        return best

    matrices = _token_matrices(graph)
    frontier: list[tuple[tuple[int, ...], np.ndarray]] = [((), target)]
    for _ in range(max_extra):
        next_frontier: list[tuple[tuple[int, ...], np.ndarray]] = []
        for suffix, resid in frontier:
            for tok, mat in matrices.items():
                resid2 = mat.conj().T @ resid
                chain = graph.buckets.get(_phase_key(graph, resid2))
                if chain is not None:
                    candidate = tuple(chain) + (tok,) + suffix
                    if best is None or len(candidate) < len(best):
                        best = candidate
                next_frontier.append(((tok,) + suffix, resid2))
        frontier = next_frontier
    return best


class SuffixExtendedDatabase:
    """Wraps a ReductionDatabase, extending each per-wire-count graph's
    reach with :func:`extended_lookup` (depth ``max_extra`` trailing-gate
    extension) instead of a single direct query.

    Forwards attribute access to the wrapped database so it drops into
    reduce_circuit unchanged. Tracks direct_hits / extended_hits so a
    benchmark can report how often extension actually mattered.
    """

    def __init__(self, db, max_extra: int = 1):
        self._db = db
        self.max_extra = max_extra
        self.direct_hits = 0
        self.extended_hits = 0

    def __getattr__(self, name):
        return getattr(self._db, name)

    def _remap(self, block: list[GateInstance]) -> tuple[list[GateInstance], dict[int, int], object] | None:
        wires = sorted({q for gate in block for q in gate.qubits})
        graph = self._db.graphs.get(len(wires))
        if graph is None:
            return None
        forward = {wire: idx for idx, wire in enumerate(wires)}
        remapped = [
            GateInstance(name=g.name, qubits=tuple(sorted(forward[q] for q in g.qubits)), theta=g.theta)
            for g in block
        ]
        reverse = {idx: wire for wire, idx in forward.items()}
        return remapped, reverse, graph

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        result = self._remap(block)
        if result is None:
            return None
        remapped, reverse, graph = result
        u = graph.block_unitary(remapped)
        if u is None:
            return None
        direct = graph.buckets.get(_phase_key(graph, u))
        chain = extended_lookup(graph, u, self.max_extra)
        if chain is None:
            return None
        self.direct_hits += 1 if (direct is not None and tuple(chain) == tuple(direct)) else 0
        self.extended_hits += 1 if (direct is None or tuple(chain) != tuple(direct)) else 0
        decoded = graph.pool.decode(list(chain))
        return [
            GateInstance(name=g.name, qubits=tuple(sorted(reverse[q] for q in g.qubits)), theta=g.theta)
            for g in decoded
        ]

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        # Length-only extension for now; cost-aware composition of extended
        # chains is left to a follow-up if the length-only result justifies it.
        return self.try_reduce(block)

    def try_reduce_escape(self, block, rng, slack: int = 3, prefer=None):
        return self._db.try_reduce_escape(block, rng, slack, prefer)
