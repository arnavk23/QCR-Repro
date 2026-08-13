"""Exact compute graphs for Clifford gate pools.

Nodes are signed symplectic tableaux (integer keys): bit-exact, tolerance-free identity and provably correct replacements. Each node stores Pareto-optimal factorizations indexed by (two-qubit count, length) for cost-aware lookups."""

from __future__ import annotations

import pickle
import random
import time
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

from .config import GateInstance
from .symplectic import SignedTableau
from .token_pool import TokenPool

_PHASE = {1 + 0j: 0, -1 + 0j: 1, 1j: 2, -1j: 3}
_PHASE_R = (1 + 0j, -1 + 0j, 1j, -1j)
_LETTER = {"I": 0, "X": 1, "Y": 2, "Z": 3}


def _row_value(factors: tuple[tuple[int, str], ...], coeff: complex, n: int) -> int:
    """Encode one signed-Pauli row as an int (< 4 * 5^n)."""
    base = 5**n
    val = 0
    for q, letter in factors:
        val += _LETTER[letter] * (5**q)
    return _PHASE[coeff] * base + val


def _tableau_key(t: SignedTableau) -> int:
    """Pack the whole tableau into a single int (rows are independent)."""
    n = t.n
    width = 4 * 5**n
    key = 0
    for factors, coeff in t.rows:
        key = key * width + _row_value(factors, coeff, n)
    return key


def _is_clifford_pool(pool: TokenPool) -> bool:
    gs = pool.gate_set
    from math import pi

    allowed = {pi / 2, -pi / 2}
    for gate in pool.gates():
        if gate.name in ("RX", "RY", "RZ"):
            if gate.theta not in allowed:
                return False
        elif gate.name == "RXX":
            if gate.theta not in allowed:
                return False
        elif gate.name != "CZ":
            return False
    return True


def _twq_of_chain(chain: tuple[int, ...], pool: TokenPool) -> int:
    return sum(1 for tok in chain if len(pool.gate_for_token(tok).qubits) == 2)


@dataclass
class SymplecticGraph:
    """Exact compute graph over a Clifford token pool up to ``max_depth``."""

    pool: TokenPool
    max_depth: int
    max_alts: int = 4
    buckets: dict[int, tuple[int, ...]] = field(default_factory=dict)
    twq: dict[int, int] = field(default_factory=dict)
    alts: dict[int, list[tuple[int, int, tuple[int, ...]]]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.buckets:
            self._build()

    def _build(self) -> None:
        if not _is_clifford_pool(self.pool):
            raise ValueError("SymplecticGraph requires a fully Clifford token pool")
        n = self.pool.num_qubits
        tokens = self.pool.tokens()
        root = SignedTableau(n)
        root_key = _tableau_key(root)
        self.buckets = {root_key: ()}
        self.twq = {root_key: 0}
        self.alts = {root_key: [(0, 0, ())]}
        alts = self.alts

        queue: deque[tuple[SignedTableau, tuple[int, ...]]] = deque([(root, ())])
        while queue:
            current, chain = queue.popleft()
            depth = len(chain)
            if depth >= self.max_depth:
                continue
            for token in tokens:
                gate = self.pool.gate_for_token(token)
                nxt = SignedTableau(n)
                nxt.rows = list(current.rows)
                nxt.apply(gate)
                key = _tableau_key(nxt)
                new_chain = chain + (token,)
                new_len = len(new_chain)
                new_twq = _twq_of_chain(new_chain, self.pool)
                existing = self.buckets.get(key)
                if existing is None:
                    self.buckets[key] = new_chain
                    self.twq[key] = new_twq
                    alts[key] = [(new_twq, new_len, new_chain)]
                    queue.append((nxt, new_chain))
                else:
                    # Record a few Pareto-optimal (twq, len) factorizations.
                    entry = (new_twq, new_len, new_chain)
                    bucket = alts.get(key)
                    if bucket is None:
                        bucket = []
                        alts[key] = bucket
                    dominated = any(e[0] <= new_twq and e[1] <= new_len and e != entry for e in bucket)
                    if dominated:
                        continue
                    bucket = [e for e in bucket if not (new_twq <= e[0] and new_len <= e[1] and e != entry)]
                    bucket.append(entry)
                    bucket.sort(key=lambda e: (e[0], e[1]))
                    alts[key] = bucket[: self.max_alts]

    @property
    def num_nodes(self) -> int:
        return len(self.buckets)

    @property
    def num_edges(self) -> int:
        return len(self.pool.tokens()) * self.num_nodes

    # ------------------------------------------------------------------ #
    # lookups
    # ------------------------------------------------------------------ #

    def block_key(self, block: list[GateInstance]) -> int | None:
        """Exact tableau key of ``block``, or None if a gate is not Clifford."""
        try:
            t = SignedTableau(self.pool.num_qubits)
            for gate in block:
                token = self.pool.token_for_gate(gate)
                t.apply(self.pool.gate_for_token(token))
        except KeyError:
            return None
        except ValueError:
            return None
        return _tableau_key(t)

    def _candidates(self, block: list[GateInstance]) -> list[tuple[int, int, tuple[int, ...]]] | None:
        key = self.block_key(block)
        if key is None or key not in self.buckets:
            return None
        cands = [e for e in self.alts.get(key, [])]
        return cands

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        """Shortest factorization of ``block``, or None if not reducible."""
        key = self.block_key(block)
        if key is None:
            return None
        chain = self.buckets.get(key)
        if chain is None or len(chain) >= len(block):
            return None
        return self.pool.decode(list(chain))

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        """Reduce block minimizing (two-qubit count, length) among stored factorizations.

Accepts strictly shorter candidates and equal-length candidates with strictly fewer two-qubit gates."""
        key = self.block_key(block)
        if key is None:
            return None
        cands = self.alts.get(key)
        if not cands:
            return None
        block_twq = _twq_of_chain(tuple(self.pool.token_for_gate(g) for g in block), self.pool)
        block_len = len(block)
        best: tuple[int, int, tuple[int, ...]] | None = None
        for (twq, ln, chain) in cands:
            if ln > block_len:
                continue
            if ln == block_len and twq >= block_twq:
                continue
            if best is None or (twq, ln) < (best[0], best[1]):
                best = (twq, ln, chain)
        if best is None:
            return None
        return self.pool.decode(list(best[2]))

    def try_reduce_escape(
        self,
        block: list[GateInstance],
        rng: random.Random,
        slack: int = 3,
        prefer: dict[str, float] | None = None,
    ) -> list[GateInstance] | None:
        """Resample ``block`` with a structurally different equivalent word."""
        key = self.block_key(block)
        if key is None:
            return None
        shortest = self.buckets.get(key)
        if shortest is None:
            return None
        cands = [e for e in self.alts.get(key, []) if e[2] != shortest and 1 <= e[1] <= len(block) + slack]
        if not cands:
            return None
        # Prefer candidates that change the two-qubit structure the most.
        block_twq = _twq_of_chain(tuple(self.pool.token_for_gate(g) for g in block), self.pool)
        best = max(cands, key=lambda e: abs(e[0] - block_twq))
        ties = [e for e in cands if abs(e[0] - block_twq) == abs(best[0] - block_twq)]
        chosen = rng.choice(ties)
        return self.pool.decode(list(chosen[2]))

    def signature(self) -> tuple:
        return (self.pool.num_qubits, self.pool.gate_set.name, self.pool.angles,
                self.pool.two_qubit_angles, self.max_depth)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, path: Path) -> "SymplecticGraph":
        with path.open("rb") as file:
            return pickle.load(file)


@dataclass
class SymplecticDatabase:
    """Databases of exact graphs indexed by the number of wires a block touches."""

    gate_set_name: str
    depths: dict[int, int] = field(default_factory=dict)
    angles: tuple[float, ...] | None = None
    two_qubit_angles: tuple[float, ...] | None = None
    graphs: dict[int, SymplecticGraph] = field(default_factory=dict)

    def build(self) -> None:
        for wires, depth in sorted(self.depths.items()):
            pool = TokenPool(
                num_qubits=wires,
                gate_set=self.gate_set_name,
                angles=self.angles,
                two_qubit_angles=self.two_qubit_angles,
            )
            self.graphs[wires] = SymplecticGraph(pool, depth)

    def _graph_for(self, block: list[GateInstance]) -> SymplecticGraph | None:
        wires = sorted({q for gate in block for q in gate.qubits})
        return self.graphs.get(len(wires))

    @staticmethod
    def _remap(block: list[GateInstance]) -> tuple[list[GateInstance], dict[int, int], dict[int, int]]:
        wires = sorted({q for gate in block for q in gate.qubits})
        forward = {wire: idx for idx, wire in enumerate(wires)}
        reverse = {idx: wire for wire, idx in forward.items()}
        remapped = [
            GateInstance(name=gate.name, qubits=tuple(sorted(forward[q] for q in gate.qubits)), theta=gate.theta)
            for gate in block
        ]
        return remapped, forward, reverse

    def _restore(self, candidate: list[GateInstance], reverse: dict[int, int]) -> list[GateInstance]:
        return [
            GateInstance(name=gate.name, qubits=tuple(sorted(reverse[q] for q in gate.qubits)), theta=gate.theta)
            for gate in candidate
        ]

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        graph = self._graph_for(block)
        if graph is None:
            return None
        remapped, _, reverse = self._remap(block)
        candidate = graph.try_reduce(remapped)
        if candidate is None:
            return None
        return self._restore(candidate, reverse)

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        graph = self._graph_for(block)
        if graph is None:
            return None
        remapped, _, reverse = self._remap(block)
        candidate = graph.try_reduce_cost(remapped)
        if candidate is None:
            return None
        return self._restore(candidate, reverse)

    def try_reduce_escape(
        self,
        block: list[GateInstance],
        rng: random.Random,
        slack: int = 3,
        prefer: dict[str, float] | None = None,
    ) -> list[GateInstance] | None:
        graph = self._graph_for(block)
        if graph is None:
            return None
        remapped, _, reverse = self._remap(block)
        candidate = graph.try_reduce_escape(remapped, rng, slack, prefer)
        if candidate is None:
            return None
        return self._restore(candidate, reverse)

    def signature(self) -> tuple:
        return (self.gate_set_name, tuple(sorted(self.depths.items())), self.angles, self.two_qubit_angles)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, path: Path) -> "SymplecticDatabase":
        with path.open("rb") as file:
            return pickle.load(file)


_CACHE_DIR = Path(__file__).resolve().parents[2] / ".cache"


def _exact_filename(db: SymplecticDatabase) -> str:
    depths = "_".join(f"w{w}d{d}" for w, d in sorted(db.depths.items()))
    angles = "".join(str(a) for a in db.angles).replace(".", "p").replace("-", "m")
    tq = "".join(str(a) for a in db.two_qubit_angles).replace(".", "p").replace("-", "m")
    return f"exact_{db.gate_set_name}_{depths}_ang{angles}_tq{tq}.pkl"


def load_or_build_exact(
    gate_set_name: str,
    depths: dict[int, int],
    angles: tuple[float, ...] | None = None,
    two_qubit_angles: tuple[float, ...] | None = None,
    cache_dir: Path = _CACHE_DIR,
    verbose: bool = False,
) -> SymplecticDatabase:
    """Build (and cache on disk) a SymplecticDatabase for the given configuration."""
    from .config import gateset_for

    gs = gateset_for(gate_set_name)
    if angles is None:
        angles = gs.angles
    if two_qubit_angles is None:
        two_qubit_angles = gs.two_angles

    db = SymplecticDatabase(gate_set_name=gate_set_name, depths=dict(depths), angles=angles,
                            two_qubit_angles=two_qubit_angles)
    path = cache_dir / _exact_filename(db)
    if path.exists():
        if verbose:
            print(f"[cache] loading {path.name}")
        with path.open("rb") as file:
            return pickle.load(file)
    if verbose:
        print(f"[cache] building {path.name}")
    t0 = time.time()
    db.build()
    elapsed = time.time() - t0
    db.save(path)
    if verbose:
        for wires, graph in sorted(db.graphs.items()):
            print(f"  {wires}-wire exact graph: {graph.num_nodes} nodes (depth {graph.max_depth}, built {elapsed:.1f}s total)")
    return db
