from __future__ import annotations

import hashlib
import json
import pickle
import random
import shutil
import sqlite3
import time
from collections import OrderedDict, deque
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import GateInstance
from .gates import embedded_gate_matrix
from .token_pool import TokenPool


def _token_type_counts(gates: list[GateInstance]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for gate in gates:
        counts[gate.name] = counts.get(gate.name, 0) + 1
    return counts


def _chain_two_qubit(chain: tuple[int, ...], pool: TokenPool) -> int:
    """Number of two-qubit gates in a stored token chain."""
    return sum(1 for tok in chain if len(pool.gate_for_token(tok).qubits) == 2)


_MISSING = object()


class SqliteDict:
    """Dict-like view over a SQLite table (BLOB key -> pickled BLOB value).

    Batched writes in transactions and a bounded LRU read cache; preserves
    value identity so callers may mutate and re-assign.  Backs buckets/alts
    when backend="sqlite", trading build RAM for disk."""

    def __init__(
        self,
        conn: sqlite3.Connection,
        table: str,
        cache_size: int = 65536,
        write_buffer: int = 8192,
    ) -> None:
        self._conn = conn
        self._table = table
        self._cache: OrderedDict[bytes, object] = OrderedDict()
        self._buffer: dict[bytes, object] = {}
        self._cache_size = cache_size
        self._write_buffer = write_buffer

    def flush(self) -> None:
        """Write all buffered entries to the table in one transaction.

The buffer is cleared only after commit, so a failed write never drops entries."""
        if not self._buffer:
            return
        items = [(k, pickle.dumps(v)) for k, v in self._buffer.items()]
        with self._conn:
            self._conn.executemany(
                f"INSERT OR REPLACE INTO {self._table} (k, v) VALUES (?, ?)", items
            )
        self._buffer.clear()

    def get(self, key: bytes, default=None):
        if key in self._buffer:
            return self._buffer[key]
        value = self._cache.get(key, _MISSING)
        if value is _MISSING:
            row = self._conn.execute(
                f"SELECT v FROM {self._table} WHERE k = ?", (key,)
            ).fetchone()
            if row is None:
                return default
            value = pickle.loads(row[0])
            self._cache[key] = value
            if len(self._cache) > self._cache_size:
                self._cache.popitem(last=False)
        return value

    def __getitem__(self, key: bytes):
        value = self.get(key, _MISSING)
        if value is _MISSING:
            raise KeyError(key)
        return value

    def __setitem__(self, key: bytes, value) -> None:
        self._buffer[key] = value
        self._cache[key] = value
        if len(self._cache) > self._cache_size:
            self._cache.popitem(last=False)
        if len(self._buffer) >= self._write_buffer:
            self.flush()

    def __contains__(self, key: bytes) -> bool:
        if key in self._buffer or key in self._cache:
            return True
        return (
            self._conn.execute(
                f"SELECT 1 FROM {self._table} WHERE k = ?", (key,)
            ).fetchone()
            is not None
        )

    def __len__(self) -> int:
        row = self._conn.execute(f"SELECT COUNT(*) FROM {self._table}").fetchone()
        return int(row[0])

    def __iter__(self):
        for key in self.keys():
            yield key

    def keys(self):
        """Iterate all stored keys (committed table rows, then buffered ones)."""
        seen: set[bytes] = set()
        for (key,) in self._conn.execute(f"SELECT k FROM {self._table}"):
            seen.add(key)
            yield key
        for key in self._buffer:
            if key not in seen:
                yield key

    def items(self):
        for key in self.keys():
            yield key, self.get(key)


@dataclass
class ComputeGraph:
    """Exhaustive compute graph over a token pool up to max_depth.

    Nodes are phase-normalized unitaries keyed by SHA-256 digests; edges are
    pool gates; each node stores the shortest pool-token chain (BFS).  With
    backend="sqlite" the tables live in a SQLite file instead of RAM."""

    pool: TokenPool
    max_depth: int
    digest_decimals: int = 10
    buckets: dict[bytes, tuple[int, ...]] = field(default_factory=dict)
    alts: dict[bytes, list[tuple[int, ...]]] = field(default_factory=dict)
    max_alts: int = 4
    _token_matrices: dict[int, np.ndarray] = field(default_factory=dict, repr=False)
    _conn: sqlite3.Connection | None = field(default=None, repr=False, init=False)
    backend: str = "ram"
    store_path: Path | None = None

    def __post_init__(self) -> None:
        if self.backend == "sqlite":
            self._attach_sqlite()
        if not self.buckets:
            self._build()

    def _attach_sqlite(self) -> None:
        """Open (or create) the SQLite tables backing buckets/alts."""
        if self.store_path is None:
            raise ValueError("store_path is required with backend='sqlite'")
        self.store_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(self.store_path))
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute(
            "CREATE TABLE IF NOT EXISTS buckets (k BLOB PRIMARY KEY, v BLOB NOT NULL)"
        )
        conn.execute(
            "CREATE TABLE IF NOT EXISTS alts (k BLOB PRIMARY KEY, v BLOB NOT NULL)"
        )
        conn.commit()
        self._conn = conn
        self.buckets = SqliteDict(conn, "buckets")
        self.alts = SqliteDict(conn, "alts")
        self._init_token_matrices()

    def _init_token_matrices(self) -> None:
        self._token_matrices = {
            token: embedded_gate_matrix(self.pool.num_qubits, self.pool.gate_for_token(token))
            for token in self.pool.tokens()
        }

    def flush(self) -> None:
        """Flush pending writes to disk (no-op in RAM mode)."""
        if self.backend != "sqlite":
            return
        if isinstance(self.buckets, SqliteDict):
            self.buckets.flush()
        if isinstance(self.alts, SqliteDict):
            self.alts.flush()
        if self._conn is not None:
            self._conn.commit()

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
        if self.backend == "ram":
            self.buckets = {root_key: ()}
        else:
            self.buckets[root_key] = ()

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
                        # write the mutation back so the SQLite store persists it
                        alts[key] = bucket
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
        """Return the shortest pool-token chain for unitary, or None.

Phase normalization runs inline on the flattened matrix for speed."""
        flat = unitary.reshape(-1)
        idx = int(np.argmax(np.round(np.abs(flat), 8)))
        phase = np.angle(flat[idx])
        nflat = flat * np.exp(-1j * phase)
        key = self._node_key(nflat, self.digest_decimals)
        return self.buckets.get(key)

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        """Reduce block (gates on arbitrary wires) if a shorter equivalent exists.

Maps <=3-wire blocks into a low-dimensional space, looks up, maps back; None if no reduction is found."""
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
        """Reduce block minimizing (two-qubit count, length).

Considers every stored factorization and keeps the best strictly-shorter or equal-length-with-fewer-two-qubit candidate."""
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
        """Resample block with a structurally different equivalent word.

Perturbs irreducible blocks so later sweeps can find new reductions; None if no sufficiently different word exists."""
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
        if self.backend == "sqlite":
            # The graph data lives in its SQLite file; it is persisted through
            # ReductionDatabase.save, not pickling (SqliteDict is not picklable).
            raise TypeError(
                "sqlite-backed graphs persist via ReductionDatabase.save(); "
                "call flush() and save the store directory instead"
            )
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, path: Path) -> "ComputeGraph":
        if path.is_dir():
            raise TypeError(
                "sqlite-backed graphs load via ReductionDatabase.load(store_dir)"
            )
        with path.open("rb") as file:
            return pickle.load(file)


@dataclass
class DiskComputeGraph:
    """Compute graph backed by SQLite instead of in-RAM dicts.

    The BFS is level-synchronous: only the current frontier's unitaries are
    held in memory; node tables (key -> shortest token chain) stream to a
    SQLite file, with an in-RAM key set for O(1) membership checks.  Deep
    (3-wire, depth 6, ~3.2M nodes) graphs build on a laptop where the
    all-RAM graph would not fit.  The public interface mirrors
    :class:`ComputeGraph`, so the rest of the pipeline works unchanged.
    """

    pool: TokenPool
    max_depth: int
    db_path: Path
    digest_decimals: int = 10
    max_alts: int = 4
    _keys: set[bytes] = field(default_factory=set, repr=False)
    _conn: sqlite3.Connection | None = field(default=None, repr=False)
    _token_matrices: dict[int, np.ndarray] = field(default_factory=dict, repr=False)
    _num_nodes: int = 0

    # ------------------------------------------------------------------ #
    # construction / persistence
    # ------------------------------------------------------------------ #

    def __post_init__(self) -> None:
        if not self._keys:
            self._build()
        self._connect()

    def _connect(self) -> sqlite3.Connection:
        if self._conn is None:
            conn = sqlite3.connect(self.db_path)
            conn.execute("PRAGMA journal_mode=OFF")
            conn.execute("PRAGMA synchronous=OFF")
            conn.execute("PRAGMA cache_size=-200000")
            self._conn = conn
        return self._conn

    @staticmethod
    def _create_schema(conn: sqlite3.Connection) -> None:
        conn.execute("CREATE TABLE IF NOT EXISTS nodes (key BLOB PRIMARY KEY, chain BLOB, twq INT)")
        conn.execute("CREATE TABLE IF NOT EXISTS alts (key BLOB, chain BLOB, twq INT)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_alts_key ON alts (key)")

    def _build(self) -> None:
        db_path = Path(self.db_path)
        db_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        self._create_schema(conn)

        dim = 2 ** self.pool.num_qubits
        identity = np.eye(dim, dtype=complex)
        cached = {
            token: embedded_gate_matrix(self.pool.num_qubits, self.pool.gate_for_token(token))
            for token in self.pool.tokens()
        }
        self._token_matrices = cached

        root_key = self._node_key(identity.reshape(-1), self.digest_decimals)
        self._keys = {root_key}
        self._num_nodes = 1
        conn.execute("INSERT OR IGNORE INTO nodes VALUES (?, ?, ?)", (root_key, b"", 0))

        # During the build the shortest chains and alternative chains are held
        # in RAM (they are needed for O(1) dedup / alt bookkeeping); they are
        # streamed to SQLite level by level and dropped before return so the
        # frontier matrices are the only large resident structure.
        shortest: dict[bytes, tuple[int, ...]] = {root_key: ()}
        alts: dict[bytes, list[tuple[int, ...]]] = {root_key: []}
        twq_by_token = {
            token: 1 if len(self.pool.gate_for_token(token).qubits) == 2 else 0
            for token in self.pool.tokens()
        }

        frontier: list[tuple[np.ndarray, tuple[int, ...]]] = [(identity, ())]
        level = 0
        t0 = time.time()
        while frontier and level < self.max_depth:
            level += 1
            next_frontier: list[tuple[np.ndarray, tuple[int, ...]]] = []
            rows: list[tuple[bytes, bytes, int]] = []
            alt_rows: list[tuple[bytes, bytes, int]] = []
            for u, chain in frontier:
                for token, gate_u in cached.items():
                    next_u = gate_u @ u
                    flat = next_u.reshape(-1)
                    idx = int(np.argmax(np.round(np.abs(flat), 8)))
                    phase = np.angle(flat[idx])
                    nflat = flat * np.exp(-1j * phase)
                    key = self._node_key(nflat, self.digest_decimals)
                    new_chain = chain + (token,)
                    existing = shortest.get(key)
                    if existing is not None:
                        # The graph contains cycles: the same unitary admits
                        # many words.  Keep a few structurally different
                        # factorizations for escape / cost-aware resampling.
                        bucket = alts[key]
                        if (
                            len(bucket) < self.max_alts
                            and new_chain != existing
                            and new_chain not in bucket
                        ):
                            bucket.append(new_chain)
                            alt_rows.append((key, pickle.dumps(new_chain), 0))
                        continue
                    shortest[key] = new_chain
                    alts[key] = []
                    self._keys.add(key)
                    self._num_nodes += 1
                    rows.append((key, pickle.dumps(new_chain), twq_by_token[token]))
                    next_frontier.append((next_u, new_chain))
            conn.executemany("INSERT OR IGNORE INTO nodes VALUES (?, ?, ?)", rows)
            if alt_rows:
                conn.executemany("INSERT INTO alts VALUES (?, ?, ?)", alt_rows)
            conn.commit()
            frontier = next_frontier
            print(f"    [disk graph] level {level}: +{len(rows)} nodes "
                  f"(frontier {len(frontier)}, total {self._num_nodes}, "
                  f"{time.time() - t0:.1f}s)", flush=True)
        conn.commit()
        conn.close()
        # Free the build-time tables; lookups read from SQLite from here on.
        del shortest, alts

    @staticmethod
    def _node_key(flat_normalized: np.ndarray, decimals: int) -> bytes:
        rounded = (np.round(flat_normalized, decimals) + 0.0).tobytes()
        return hashlib.sha256(rounded).digest()

    def save(self, path: Path) -> None:
        # Persist the manifest; the SQLite file lives at self.db_path.
        path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = None  # drop the live connection before pickling
        with path.open("wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, path: Path) -> "DiskComputeGraph":
        with path.open("rb") as file:
            obj = pickle.load(file)
        obj._connect()
        obj._rebuild_keys()
        return obj

    def _rebuild_keys(self) -> None:
        conn = self._connect()
        self._keys = set()
        self._num_nodes = 0
        for (key,) in conn.execute("SELECT key FROM nodes"):
            self._keys.add(key)
            self._num_nodes += 1

    def __getstate__(self):
        state = dict(self.__dict__)
        state["_conn"] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._connect()
        self._rebuild_keys()

    # ------------------------------------------------------------------ #
    # lookups
    # ------------------------------------------------------------------ #

    @property
    def num_nodes(self) -> int:
        return self._num_nodes

    @property
    def num_edges(self) -> int:
        return self.pool.num_qubits * len(self.pool.tokens()) * self._num_nodes

    @property
    def buckets(self) -> "_BucketView":
        return _BucketView(self)

    def lookup(self, unitary: np.ndarray) -> tuple[int, ...] | None:
        flat = unitary.reshape(-1)
        idx = int(np.argmax(np.round(np.abs(flat), 8)))
        phase = np.angle(flat[idx])
        nflat = flat * np.exp(-1j * phase)
        key = self._node_key(nflat, self.digest_decimals)
        return self._chain_for(key)

    def _chain_for(self, key: bytes) -> tuple[int, ...] | None:
        if key not in self._keys:
            return None
        conn = self._connect()
        row = conn.execute("SELECT chain FROM nodes WHERE key = ?", (key,)).fetchone()
        if row is None:
            return None
        return pickle.loads(row[0]) if row[0] else ()

    def _alts_for(self, key: bytes) -> list[tuple[int, ...]]:
        if key not in self._keys:
            return []
        conn = self._connect()
        out = []
        for (chain_blob,) in conn.execute("SELECT chain FROM alts WHERE key = ?", (key,)):
            if chain_blob:
                out.append(pickle.loads(chain_blob))
        return out

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
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
        return key if key in self._keys else None

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        u = self.block_unitary(block)
        if u is None:
            return None
        key = self._lookup_key(u)
        if key is None:
            return None
        block_twq = sum(1 for g in block if len(g.qubits) == 2)
        block_len = len(block)
        shortest = self._chain_for(key)
        candidates = [shortest] + self._alts_for(key)
        seen: set[tuple[int, ...]] = set()
        best: tuple[int, int, tuple[int, ...]] | None = None
        for chain in candidates:
            if chain is None or chain in seen:
                continue
            seen.add(chain)
            if len(chain) > block_len:
                continue
            twq = _chain_two_qubit(chain, self.pool)
            if len(chain) == block_len and twq >= block_twq:
                continue
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
        u = self.block_unitary(block)
        if u is None:
            return None
        key = self._lookup_key(u)
        if key is None:
            return None
        shortest = self._chain_for(key)
        if shortest is None:
            return None
        alts = [c for c in self._alts_for(key) if 1 <= len(c) <= len(block) + slack and c != shortest]
        if not alts:
            return None

        block_counts = _token_type_counts(block)
        prefer = prefer or {}

        def score(candidate: tuple[int, ...]) -> float:
            cand_counts = _token_type_counts(self.pool.decode(list(candidate)))
            diff = sum(abs(cand_counts.get(k, 0) - block_counts.get(k, 0)) for k in set(cand_counts) | set(block_counts))
            pref = sum(prefer.get(g, 0.0) * cand_counts.get(g, 0) for g in prefer)
            return diff + pref

        best_score = max(score(c) for c in alts)
        ties = [c for c in alts if score(c) == best_score]
        chosen = rng.choice(ties)
        return self.pool.decode(list(chosen))

    def block_unitary(self, block: list[GateInstance]) -> np.ndarray | None:
        dim = 2 ** self.pool.num_qubits
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


class _BucketView:
    """Read-only ``dict``-like view over a :class:`DiskComputeGraph` node table.

    Exists so the batched sweep (``batched.lookup_batch`` does
    ``graph.buckets.get(key)``) works unchanged against disk-backed graphs.
    """

    def __init__(self, graph: DiskComputeGraph):
        self._graph = graph

    def get(self, key: bytes, default=None):
        chain = self._graph._chain_for(key)
        return chain if chain is not None else default

    def keys(self):
        for key in self._graph._keys:
            yield key

    def values(self):
        for key in self._graph._keys:
            yield self._graph._chain_for(key)

    def items(self):
        for key in self._graph._keys:
            yield key, self._graph._chain_for(key)

    def __iter__(self):
        return self.keys()

    def __getitem__(self, key: bytes):
        chain = self._graph._chain_for(key)
        if chain is None:
            raise KeyError(key)
        return chain

    def __contains__(self, key: bytes) -> bool:
        return key in self._graph._keys

    def __len__(self) -> int:
        return self._graph._num_nodes


@dataclass
class ReductionDatabase:
    """Databases indexed by the number of wires a sampled block touches."""

    gate_set_name: str
    depths: dict[int, int] = field(default_factory=dict)
    angles: tuple[float, ...] | None = None
    two_qubit_angles: tuple[float, ...] | None = None
    digest_decimals: int = 10
    graphs: dict[int, ComputeGraph] = field(default_factory=dict)
    backend: str = "ram"
    store_dir: Path | None = None

    def build(self, storage: str = "ram", disk_dir: Path | None = None) -> None:
        for wires, depth in sorted(self.depths.items()):
            pool = TokenPool(
                num_qubits=wires,
                gate_set=self.gate_set_name,
                angles=self.angles,
                two_qubit_angles=self.two_qubit_angles,
            )
            if storage == "disk":
                db_path = disk_dir / f"graph_{self.gate_set_name}_{wires}w_d{depth}.sqlite"
                self.graphs[wires] = DiskComputeGraph(pool, depth, db_path, self.digest_decimals)
            else:
                self.graphs[wires] = ComputeGraph(pool, depth, self.digest_decimals)

    def signature(self) -> tuple:
        return (
            self.gate_set_name,
            tuple(sorted(self.depths.items())),
            self.angles,
            self.two_qubit_angles,
            self.digest_decimals,
        )

    def save(self, path: Path) -> None:
        if self.graphs and isinstance(next(iter(self.graphs.values())), DiskComputeGraph):
            # DiskComputeGraph data lives in per-wire SQLite files; persist a
            # JSON meta record plus a small manifest pickle per graph (the
            # pickle only carries the key set + path; chains stay in SQLite).
            path.mkdir(parents=True, exist_ok=True)
            meta = {
                "gate_set_name": self.gate_set_name,
                "depths": dict(self.depths),
                "angles": self.angles,
                "two_qubit_angles": self.two_qubit_angles,
                "digest_decimals": self.digest_decimals,
                "storage": "disk",
            }
            with (path / "meta.json").open("w", encoding="utf-8") as file:
                json.dump(meta, file)
            for wires, graph in self.graphs.items():
                graph.save(path / f"graph_{self.gate_set_name}_{wires}w_d{self.depths[wires]}.pkl")
            return
        if self.backend == "sqlite":
            # The graph data already lives in SQLite files under ``path``;
            # persist a small JSON meta record and flush any pending writes.
            for graph in self.graphs.values():
                graph.flush()
            path.mkdir(parents=True, exist_ok=True)
            meta = {
                "gate_set_name": self.gate_set_name,
                "depths": dict(self.depths),
                "angles": self.angles,
                "two_qubit_angles": self.two_qubit_angles,
                "digest_decimals": self.digest_decimals,
            }
            with (path / "meta.json").open("w", encoding="utf-8") as file:
                json.dump(meta, file)
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, path: Path) -> "ReductionDatabase":
        if path.is_dir():
            with (path / "meta.json").open("r", encoding="utf-8") as file:
                meta = json.load(file)
            db = cls(
                gate_set_name=meta["gate_set_name"],
                depths={int(k): int(v) for k, v in meta["depths"].items()},
                angles=tuple(meta["angles"]) if meta.get("angles") else None,
                two_qubit_angles=tuple(meta["two_qubit_angles"])
                if meta.get("two_qubit_angles")
                else None,
                digest_decimals=meta.get("digest_decimals", 10),
                backend="sqlite",
                store_dir=path,
            )
            if meta.get("storage") == "disk":
                for wires, depth in sorted(meta["depths"].items()):
                    wires = int(wires)
                    graph = DiskComputeGraph.load(
                        path / f"graph_{meta['gate_set_name']}_{wires}w_d{depth}.pkl"
                    )
                    db.graphs[wires] = graph
                return db
            # legacy backend: build() attaches to the existing SQLite files
            # (no rebuild: the buckets tables are non-empty) and initializes
            # token matrices.
            db.build()
            return db
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
        """Reduce ``block`` (touching <=3 wires) to a shorter equivalent, or None.

        Blocks touching more wires than any available graph are not reducible.
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


def _repo_root() -> Path:
    """Locate the repository root (contains pyproject.toml) from any layout."""
    p = Path(__file__).resolve().parent
    for _ in range(4):
        if (p / "pyproject.toml").exists():
            return p
        p = p.parent
    return Path(__file__).resolve().parents[2]


_CACHE_DIR = _repo_root() / ".cache"


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
    backend: str = "ram",
) -> ReductionDatabase:
    """Build (and disk-cache) a ReductionDatabase for the given configuration.

backend="sqlite" stores bucket tables in SQLite files so deep graphs build within a laptop's memory."""
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

    if backend == "sqlite":
        store_dir = cache_dir / (_db_filename(db) + "_disk")
        if (store_dir / "meta.json").exists():
            if verbose:
                print(f"[cache] loading {store_dir.name}")
            return ReductionDatabase.load(store_dir)
        # No meta record -> any leftover directory is a partial build from an
        # interrupted run; wipe it so the rebuild starts from a clean slate.
        if store_dir.exists():
            shutil.rmtree(store_dir, ignore_errors=True)
        if verbose:
            print(f"[cache] building {store_dir.name}")
        db.backend = "sqlite"
        db.store_dir = store_dir
        t0 = time.time()
        db.build(storage="disk", disk_dir=store_dir)
        elapsed = time.time() - t0
        db.save(store_dir)
        if verbose:
            for wires, graph in sorted(db.graphs.items()):
                print(f"  {wires}-wire disk graph: {graph.num_nodes} nodes (depth {graph.max_depth}, built {elapsed:.1f}s total)")
        return db

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
