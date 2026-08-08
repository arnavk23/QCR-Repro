from __future__ import annotations

import hashlib
import pickle
import random
import time
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import GateInstance
from .gates import embedded_gate_matrix
from .tokenizer import TokenPool


def _token_type_counts(gates: list[GateInstance]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for gate in gates:
        counts[gate.name] = counts.get(gate.name, 0) + 1
    return counts


def _chain_two_qubit(chain: tuple[int, ...], pool: TokenPool) -> int:
    """Number of two-qubit gates in a stored token chain."""
    return sum(1 for tok in chain if len(pool.gate_for_token(tok).qubits) == 2)


@dataclass
class ComputeGraph:
    """Exhaustive compute graph over a token pool up to ``max_depth``.

    Nodes are phase-normalized unitaries; edges are pool gates.  For every node
    we store the *shortest* pool-token chain implementing it (found by BFS), so a
    lookup returns the optimal factorization of any reachable unitary.

    Keys are 256-bit SHA-256 digests of the rounded, phase-normalized matrix so
    that even large (2^4 x 2^4) unitaries produce tiny in-memory entries.
    """

    pool: TokenPool
    max_depth: int
    digest_decimals: int = 10
    buckets: dict[bytes, tuple[int, ...]] = field(default_factory=dict)
    alts: dict[bytes, list[tuple[int, ...]]] = field(default_factory=dict)
    max_alts: int = 4
    _token_matrices: dict[int, np.ndarray] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        if not self.buckets:
            self._build()

    def _alt_lists(self) -> dict[bytes, list[tuple[int, ...]]]:
        # Tolerate pickles built before the alternative-chains feature existed.
        alts = getattr(self, "alts", None)
        if alts is None:
            alts = {}
            self.alts = alts
        return alts

    @staticmethod
    def _node_key(flat_normalized: np.ndarray, decimals: int) -> bytes:
        rounded = (np.round(flat_normalized, decimals) + 0.0).tobytes()
        return hashlib.sha256(rounded).digest()

    def _build(self) -> None:
        dim = 2**self.pool.num_qubits
        identity = np.eye(dim, dtype=complex)
        cached = {
            token: embedded_gate_matrix(self.pool.num_qubits, self.pool.gate_for_token(token))
            for token in self.pool.tokens()
        }
        self._token_matrices = cached

        root_key = self._node_key(identity.reshape(-1), self.digest_decimals)
        self.buckets = {root_key: ()}

        queue: deque[tuple[np.ndarray, tuple[int, ...]]] = deque([(identity, ())])
        while queue:
            current_u, chain = queue.popleft()
            depth = len(chain)
            if depth >= self.max_depth:
                continue
            for token, gate_u in cached.items():
                next_u = gate_u @ current_u
                flat = next_u.reshape(-1)
                idx = int(np.argmax(np.round(np.abs(flat), 8)))
                phase = np.angle(flat[idx])
                nflat = flat * np.exp(-1j * phase)
                key = self._node_key(nflat, self.digest_decimals)
                existing = self.buckets.get(key)
                if existing is not None:
                    # The compute graph contains cycles: the same unitary admits
                    # many different gate words.  Keep a few *structurally
                    # different* factorizations per node (without expanding
                    # them further) to enable equivalence-class resampling
                    # during reduction (see try_reduce_escape).
                    alts = self._alt_lists()
                    bucket = alts.get(key)
                    if bucket is None:
                        bucket = []
                        alts[key] = bucket
                    candidate = chain + (token,)
                    if len(bucket) < self.max_alts and candidate != existing and candidate not in bucket:
                        bucket.append(candidate)
                    continue
                self.buckets[key] = chain + (token,)
                queue.append((next_u, chain + (token,)))

    @property
    def num_nodes(self) -> int:
        return len(self.buckets)

    @property
    def num_edges(self) -> int:
        return self.pool.num_qubits * len(self.pool.tokens()) * self.num_nodes

    def lookup(self, unitary: np.ndarray) -> tuple[int, ...] | None:
        """Return the shortest pool-token chain for ``unitary`` or ``None``.

        Phase normalization is performed inline on the flattened matrix to
        avoid allocation-heavy intermediate arrays on the hot lookup loop.
        """
        flat = unitary.reshape(-1)
        idx = int(np.argmax(np.round(np.abs(flat), 8)))
        phase = np.angle(flat[idx])
        nflat = flat * np.exp(-1j * phase)
        key = self._node_key(nflat, self.digest_decimals)
        return self.buckets.get(key)

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        """Shortest factorization of ``block`` in this graph's wire space, or None."""
        u = self.block_unitary(block)
        if u is None:
            return None
        chain = self.lookup(u)
        if chain is None:
            return None
        return self.pool.decode(list(chain))

    def _lookup_key(self, unitary: np.ndarray) -> bytes | None:
        flat = unitary.reshape(-1)
        idx = int(np.argmax(np.round(np.abs(flat), 8)))
        phase = np.angle(flat[idx])
        key = self._node_key(flat * np.exp(-1j * phase), self.digest_decimals)
        return key if key in self.buckets else None

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        """Reduce ``block`` minimizing (two-qubit count, length).

        Considers every stored factorization of the block's unitary (the BFS
        shortest word plus the alternative chains found through graph cycles)
        and returns the best one that is strictly shorter or equal-length with
        strictly fewer two-qubit gates -- a numeric analogue of the cost-aware
        objective the exact engine implements for Clifford pools.
        """
        u = self.block_unitary(block)
        if u is None:
            return None
        key = self._lookup_key(u)
        if key is None:
            return None
        block_twq = sum(1 for g in block if len(g.qubits) == 2)
        block_len = len(block)
        candidates = [self.buckets[key]] + list(self._alt_lists().get(key, ()))
        seen: set[tuple[int, ...]] = set()
        best: tuple[int, int, tuple[int, ...]] | None = None
        for chain in candidates:
            if chain in seen:
                continue
            seen.add(chain)
            if len(chain) > block_len:
                continue
            if len(chain) == block_len and _chain_two_qubit(chain, self.pool) >= block_twq:
                continue
            twq = _chain_two_qubit(chain, self.pool)
            if best is None or (twq, len(chain)) < (best[0], best[1]):
                best = (twq, len(chain), chain)
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
        """Resample ``block`` with a structurally different equivalent word.

        The compute graph contains cycles, so a unitary usually admits several
        distinct gate words of comparable length.  Replacing an irreducible
        block with a *different* (possibly longer) equivalent word perturbs the
        circuit so that sweeps can find reductions that were previously out of
        reach.  Returns None if no sufficiently different word exists.
        """
        u = self.block_unitary(block)
        if u is None:
            return None
        key = self._lookup_key(u)
        if key is None:
            return None
        shortest = self.buckets[key]
        alts = self._alt_lists().get(key, ())
        cap = len(block) + slack
        candidates = [c for c in alts if 1 <= len(c) <= cap and c != shortest]
        if not candidates:
            return None

        block_counts = _token_type_counts(block)
        prefer = prefer or {}

        def score(candidate: tuple[int, ...]) -> float:
            cand_counts = _token_type_counts(self.pool.decode(list(candidate)))
            diff = sum(abs(cand_counts.get(k, 0) - block_counts.get(k, 0)) for k in set(cand_counts) | set(block_counts))
            pref = sum(prefer.get(g, 0.0) * cand_counts.get(g, 0) for g in prefer)
            return diff + pref

        best_score = max(score(c) for c in candidates)
        ties = [c for c in candidates if score(c) == best_score]
        chosen = rng.choice(ties)
        return self.pool.decode(list(chosen))

    def block_unitary(self, block: list[GateInstance]) -> np.ndarray | None:
        """Unitary of ``block`` via precomputed pool token matrices, or None.

        Returns None if any gate is not representable by this graph's pool.
        """
        dim = 2**self.pool.num_qubits
        u = np.eye(dim, dtype=complex)
        try:
            for gate in block:
                token = self.pool.token_for_gate(gate)
                u = self._token_matrices[token] @ u
        except KeyError:
            return None
        return u

    def signature(self) -> tuple:
        return (
            self.pool.num_qubits,
            self.pool.gate_set.name,
            self.pool.angles,
            self.pool.two_qubit_angles,
            self.max_depth,
            self.digest_decimals,
        )

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, path: Path) -> "ComputeGraph":
        with path.open("rb") as file:
            return pickle.load(file)


@dataclass
class ReductionDatabase:
    """Databases indexed by the number of wires a sampled block touches."""

    gate_set_name: str
    depths: dict[int, int] = field(default_factory=dict)
    angles: tuple[float, ...] | None = None
    two_qubit_angles: tuple[float, ...] | None = None
    digest_decimals: int = 10
    graphs: dict[int, ComputeGraph] = field(default_factory=dict)

    def build(self) -> None:
        for wires, depth in sorted(self.depths.items()):
            pool = TokenPool(
                num_qubits=wires,
                gate_set=self.gate_set_name,
                angles=self.angles,
                two_qubit_angles=self.two_qubit_angles,
            )
            graph = ComputeGraph(pool, depth, self.digest_decimals)
            self.graphs[wires] = graph

    def signature(self) -> tuple:
        return (
            self.gate_set_name,
            tuple(sorted(self.depths.items())),
            self.angles,
            self.two_qubit_angles,
            self.digest_decimals,
        )

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, path: Path) -> "ReductionDatabase":
        with path.open("rb") as file:
            return pickle.load(file)

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        """Cost-aware reduction (minimizes two-qubit count, then length)."""
        wires = sorted({q for gate in block for q in gate.qubits})
        wire_count = len(wires)
        graph = self.graphs.get(wire_count)
        if graph is None:
            return None

        forward = {wire: idx for idx, wire in enumerate(wires)}
        reverse = {idx: wire for wire, idx in forward.items()}
        remapped: list[GateInstance] = []
        for gate in block:
            qubits = tuple(sorted(forward[q] for q in gate.qubits))
            remapped.append(GateInstance(name=gate.name, qubits=qubits, theta=gate.theta))

        candidate = graph.try_reduce_cost(remapped)
        if candidate is None:
            return None

        return [
            GateInstance(
                name=gate.name,
                qubits=tuple(sorted(reverse[q] for q in gate.qubits)),
                theta=gate.theta,
            )
            for gate in candidate
        ]

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        """Reduce ``block`` (gates on arbitrary wires) if a shorter equivalent exists.

        Blocks touching at most 3 wires are mapped into a low-dimensional space,
        looked up, and mapped back.  Returns None if no reduction is found or the
        block touches more wires than any available graph.
        """
        wires = sorted({q for gate in block for q in gate.qubits})
        wire_count = len(wires)
        graph = self.graphs.get(wire_count)
        if graph is None:
            return None

        forward = {wire: idx for idx, wire in enumerate(wires)}
        reverse = {idx: wire for wire, idx in forward.items()}
        remapped: list[GateInstance] = []
        for gate in block:
            qubits = tuple(sorted(forward[q] for q in gate.qubits))
            remapped.append(GateInstance(name=gate.name, qubits=qubits, theta=gate.theta))

        candidate = graph.try_reduce(remapped)
        if candidate is None:
            return None

        restored = [
            GateInstance(
                name=gate.name,
                qubits=tuple(sorted(reverse[q] for q in gate.qubits)),
                theta=gate.theta,
            )
            for gate in candidate
        ]
        return restored

    def try_reduce_escape(
        self,
        block: list[GateInstance],
        rng: random.Random,
        slack: int = 3,
        prefer: dict[str, float] | None = None,
    ) -> list[GateInstance] | None:
        """Database-level escape resampling (see ComputeGraph.try_reduce_escape)."""
        wires = sorted({q for gate in block for q in gate.qubits})
        wire_count = len(wires)
        graph = self.graphs.get(wire_count)
        if graph is None:
            return None

        forward = {wire: idx for idx, wire in enumerate(wires)}
        reverse = {idx: wire for wire, idx in forward.items()}
        remapped = [
            GateInstance(name=gate.name, qubits=tuple(sorted(forward[q] for q in gate.qubits)), theta=gate.theta)
            for gate in block
        ]
        candidate = graph.try_reduce_escape(remapped, rng, slack, prefer)
        if candidate is None:
            return None
        return [
            GateInstance(
                name=gate.name,
                qubits=tuple(sorted(reverse[q] for q in gate.qubits)),
                theta=gate.theta,
            )
            for gate in candidate
        ]


_CACHE_DIR = Path(__file__).resolve().parents[2] / ".cache"


def _db_filename(db: ReductionDatabase) -> str:
    depths = "_".join(f"w{w}d{d}" for w, d in sorted(db.depths.items()))
    angles = "".join(str(a) for a in db.angles).replace(".", "p").replace("-", "m")
    tq = "".join(str(a) for a in db.two_qubit_angles).replace(".", "p").replace("-", "m")
    return f"db_{db.gate_set_name}_{depths}_ang{angles}_tq{tq}.pkl"


def load_or_build_database(
    gate_set_name: str,
    depths: dict[int, int],
    angles: tuple[float, ...] | None = None,
    two_qubit_angles: tuple[float, ...] | None = None,
    cache_dir: Path = _CACHE_DIR,
    verbose: bool = False,
) -> ReductionDatabase:
    """Build (and cache on disk) a ReductionDatabase for the given configuration."""
    from .config import gateset_for

    gs = gateset_for(gate_set_name)
    if angles is None:
        angles = gs.angles
    if two_qubit_angles is None:
        two_qubit_angles = gs.two_angles

    db = ReductionDatabase(
        gate_set_name=gate_set_name,
        depths=dict(depths),
        angles=angles,
        two_qubit_angles=two_qubit_angles,
    )
    path = cache_dir / _db_filename(db)

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
            print(f"  {wires}-wire graph: {graph.num_nodes} nodes (depth {graph.max_depth}, built {elapsed:.1f}s total)")
    return db
