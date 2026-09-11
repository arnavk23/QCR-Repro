"""Wire-permutation-symmetric database lookup.

ReductionDatabase.try_reduce always embeds a block's touched real wires into
local graph indices in one fixed way: ascending numeric order of the real
wire labels. That is one specific choice out of k! possible relabelings of a
k-wire block, and the compute graph itself has no notion of which relabeling
produced a query -- it just answers "is this exact local unitary a node I've
reached, and if so, what's the shortest chain to it."

Two different relabelings of the same block are, in general, two different
local unitaries (single-qubit gates on distinct wires don't commute past a
relabeling the way a swap-symmetric two-qubit gate would), so they can land
on different graph nodes with different shortest-known chains -- one may be
found where another isn't, or found shorter. Querying every relabeling and
decoding the best hit back through that relabeling's own inverse yields a
circuit on the original real wires implementing the exact same unitary as
the input block, since conjugating by a wire permutation and then undoing it
is lossless. This costs up to k! lookups per block instead of one, against
the *same* graph already built -- no bigger database, no gate-set-specific
algebra, applies identically to every pool (ion_trap, nisq, nisq_clifford).
"""

from __future__ import annotations

from itertools import permutations

from .config import GateInstance


def _remap(block: list[GateInstance], forward: dict[int, int]) -> list[GateInstance]:
    return [
        GateInstance(name=g.name, qubits=tuple(sorted(forward[q] for q in g.qubits)), theta=g.theta)
        for g in block
    ]


def _restore(chain: list[GateInstance], reverse: dict[int, int]) -> list[GateInstance]:
    return [
        GateInstance(name=g.name, qubits=tuple(sorted(reverse[q] for q in g.qubits)), theta=g.theta)
        for g in chain
    ]


def _count_twq(chain: list[GateInstance]) -> int:
    return sum(1 for g in chain if len(g.qubits) == 2)


class SymmetricDatabase:
    """Wraps a ReductionDatabase, querying every wire relabeling of a block
    against the same per-wire-count graphs and keeping the best hit.

    Forwards attribute access (graphs, gate_set_name, angles, ...) to the
    wrapped database so it drops into reduce_circuit unchanged. Tracks
    permuted_hits / identity_hits so a benchmark can report how often a
    non-default relabeling was the one that actually found something.
    """

    def __init__(self, db):
        self._db = db
        self.identity_hits = 0
        self.permuted_hits = 0

    def __getattr__(self, name):
        return getattr(self._db, name)

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        wires = sorted({q for gate in block for q in gate.qubits})
        graph = self._db.graphs.get(len(wires))
        if graph is None:
            return None
        best: tuple[list[GateInstance], dict[int, int], int] | None = None
        for i, perm in enumerate(permutations(wires)):
            forward = {wire: idx for idx, wire in enumerate(perm)}
            candidate = graph.try_reduce(_remap(block, forward))
            if candidate is None:
                continue
            if best is None or len(candidate) < len(best[0]):
                reverse = {idx: wire for wire, idx in forward.items()}
                best = (candidate, reverse, i)
        if best is None:
            return None
        self._note(best[2])
        return _restore(best[0], best[1])

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        wires = sorted({q for gate in block for q in gate.qubits})
        graph = self._db.graphs.get(len(wires))
        if graph is None:
            return None
        best_key: tuple[int, int] | None = None
        best: tuple[list[GateInstance], dict[int, int], int] | None = None
        for i, perm in enumerate(permutations(wires)):
            forward = {wire: idx for idx, wire in enumerate(perm)}
            candidate = graph.try_reduce_cost(_remap(block, forward))
            if candidate is None:
                continue
            key = (_count_twq(candidate), len(candidate))
            if best_key is None or key < best_key:
                reverse = {idx: wire for wire, idx in forward.items()}
                best_key, best = key, (candidate, reverse, i)
        if best is None:
            return None
        self._note(best[2])
        return _restore(best[0], best[1])

    def try_reduce_escape(self, block, rng, slack: int = 3, prefer=None):
        # Escape resampling perturbs a block deliberately; permuting the
        # embedding on top would confound its own diversity mechanism, so
        # this forwards to the identity embedding only.
        return self._db.try_reduce_escape(block, rng, slack, prefer)

    def _note(self, perm_index: int) -> None:
        if perm_index == 0:
            self.identity_hits += 1
        else:
            self.permuted_hits += 1
