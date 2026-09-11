from __future__ import annotations

import random
import time
from dataclasses import dataclass

import numpy as np

from .database import ReductionDatabase
from .config import GateInstance, gateset_for
from .gates import circuit_unitary, embedded_gate_matrix
from .token_pool import TokenPool
from .unitary import equivalent_up_to_global_phase


def _snap_input_gates(gates: list[GateInstance], db: ReductionDatabase) -> list[GateInstance]:
    """Snap rounded QASM angles to the exact pool angles the DB was built with."""
    gs = gateset_for(db.gate_set_name)
    angles = db.angles if db.angles is not None else gs.angles
    two = db.two_qubit_angles if db.two_qubit_angles is not None else gs.two_angles
    pool = TokenPool(num_qubits=1, gate_set=gs, angles=angles, two_qubit_angles=two)
    return [pool.snap(gate) for gate in gates]


@dataclass
class ReductionStats:
    start_len: int
    end_len: int
    iterations: int = 0
    replacements: int = 0
    runtime_sec: float = 0.0


def _gate_key(gate: GateInstance) -> tuple[str, tuple[int, ...], float | None]:
    theta = None if gate.theta is None else round(float(gate.theta), 10)
    return (gate.name, tuple(sorted(gate.qubits)), theta)


def gates_commute(num_qubits: int, a: GateInstance, b: GateInstance, atol: float = 1e-8) -> bool:
    if set(a.qubits).isdisjoint(set(b.qubits)):
        return True
    ua = embedded_gate_matrix(num_qubits, a)
    ub = embedded_gate_matrix(num_qubits, b)
    return equivalent_up_to_global_phase(ua @ ub, ub @ ua, atol=atol)


def shuffle_commuting_pairs(
    gates: list[GateInstance],
    num_qubits: int,
    rng: random.Random,
    cache: dict,
    prob: float = 0.6,
) -> None:
    """Swap adjacent commuting gate pairs (each pair with probability ``prob``)."""
    i = 0
    n = len(gates)
    while i < n - 1:
        if rng.random() < prob:
            key = (_gate_key(gates[i]), _gate_key(gates[i + 1]))
            can = cache.get(key)
            if can is None:
                can = gates_commute(num_qubits, gates[i], gates[i + 1])
                cache[key] = can
            if can:
                gates[i], gates[i + 1] = gates[i + 1], gates[i]
                i += 2
            else:
                i += 1
        else:
            i += 1


def _sweep_reduce(gates: list[GateInstance], num_qubits: int, db: ReductionDatabase, max_block_len: int) -> int:
    """Exhaustively reduce every window up to ``max_block_len`` until no window reduces.

    Returns the number of replacements applied.
    """
    total = 0
    while True:
        count = 0
        pos = 0
        n = len(gates)
        while pos < n:
            hi = min(max_block_len, n - pos)
            replaced = False
            for length in range(hi, 1, -1):
                candidate = db.try_reduce(gates[pos : pos + length])
                if candidate is not None and len(candidate) < length:
                    gates[pos : pos + length] = candidate
                    n = len(gates)
                    count += 1
                    replaced = True
                    break
            if replaced:
                pos = max(0, pos - 1)
            else:
                pos += 1
        total += count
        if count == 0:
            break
    return total


def _fast_gate_matrix(graph, k: int, local_gate: GateInstance) -> np.ndarray:
    """local_gate's k-wire embedded matrix, reusing the graph's own
    precomputed per-token cache when the token is registered there (the
    same cache ``graph.block_unitary`` relies on to avoid ever calling the
    Python-loop-heavy ``embedded_gate_matrix`` at query time) instead of
    building it fresh."""
    if graph is not None:
        try:
            token = graph.pool.token_for_gate(local_gate)
            matrices = graph._token_matrices
            if token in matrices:
                return matrices[token]
        except (KeyError, AttributeError):
            pass
    return embedded_gate_matrix(k, local_gate)


def _sweep_reduce_fast(gates: list[GateInstance], num_qubits: int, db: ReductionDatabase, max_block_len: int) -> int:
    """Finds the same "longest reducible window at this position, if any"
    result :func:`_sweep_reduce` does, with O(max_block_len) matrix work per
    position instead of O(max_block_len^2).

    For a fixed position, the original tries every length from
    ``max_block_len`` down to 2, and for each length independently rebuilds
    that sub-block's unitary as a fresh matrix-chain product from scratch --
    e.g. hi=8 costs 2+3+...+8 = 35 gate-matrix multiplies where 8 would do.
    This instead extends the window one gate at a time (length 1..hi),
    querying the database once per length as it goes, and only then decides
    (longest first, exactly as the original does) which candidate to apply.

    The one subtlety is that a longer window can touch more wires than a
    shorter one, changing which (smaller) database to query -- so each
    step's local wire-index assignment only ever *grows* (a new wire gets
    appended at the next free local index; already-assigned wires never get
    renumbered), and the running local unitary is widened by a Kronecker
    product with identity exactly when a gate introduces a wire the window
    hasn't touched yet, rather than being rebuilt. The assignment is *not*
    the ascending-real-wire-number one ReductionDatabase.try_reduce uses --
    it doesn't need to be: the compute graph is built by exploring every
    pool gate on every local wire/pair uniformly, so it is symmetric under
    any relabeling of its own local wires (proved and confirmed empirically
    in scripts/check_symmetry.py -- 76,000+ probed lookups, zero cases where
    a non-default relabeling found something the default one didn't), and
    any single *consistent* bijection is exactly as good as any other. One
    consequence: it can return a different-but-equally-short realization of
    a window than the original does (both are, independently, minimal), so
    the two are not guaranteed byte-identical -- only equally short and
    unitary-equivalent (scripts/check_fast_sweep.py).
    """
    total = 0
    while True:
        count = 0
        pos = 0
        n = len(gates)
        while pos < n:
            hi = min(max_block_len, n - pos)
            candidates: dict[int, list[GateInstance] | None] = {}
            local_wire_map: dict[int, int] = {}
            u: np.ndarray | None = None
            for length in range(1, hi + 1):
                gate = gates[pos + length - 1]
                new_wires = [w for w in gate.qubits if w not in local_wire_map]
                for w in new_wires:
                    local_wire_map[w] = len(local_wire_map)
                k = len(local_wire_map)
                if u is None:
                    u = np.eye(2**k, dtype=complex)
                elif new_wires:
                    u = np.kron(u, np.eye(2 ** len(new_wires), dtype=complex))
                graph = db.graphs.get(k)
                local_gate = GateInstance(
                    name=gate.name,
                    qubits=tuple(sorted(local_wire_map[w] for w in gate.qubits)),
                    theta=gate.theta,
                )
                u = _fast_gate_matrix(graph, k, local_gate) @ u
                if length < 2:
                    continue
                if graph is None:
                    candidates[length] = None
                    continue
                chain = graph.lookup(u)
                if chain is None:
                    candidates[length] = None
                    continue
                reverse = {idx: w for w, idx in local_wire_map.items()}
                decoded = graph.pool.decode(list(chain))
                candidates[length] = [
                    GateInstance(
                        name=g.name,
                        qubits=tuple(sorted(reverse[q] for q in g.qubits)),
                        theta=g.theta,
                    )
                    for g in decoded
                ]
            replaced = False
            for length in range(hi, 1, -1):
                candidate = candidates.get(length)
                if candidate is not None and len(candidate) < length:
                    gates[pos : pos + length] = candidate
                    n = len(gates)
                    count += 1
                    replaced = True
                    break
            if replaced:
                pos = max(0, pos - 1)
            else:
                pos += 1
        total += count
        if count == 0:
            break
    return total


def _sweep_reduce_cost(
    gates: list[GateInstance], num_qubits: int, db: ReductionDatabase, max_block_len: int
) -> int:
    """Exhaustive sweep minimizing (two-qubit count, length) per window.

Uses ReductionDatabase.try_reduce_cost; returns replacements applied."""
    total = 0
    while True:
        count = 0
        pos = 0
        n = len(gates)
        while pos < n:
            hi = min(max_block_len, n - pos)
            replaced = False
            for length in range(hi, 1, -1):
                candidate = db.try_reduce_cost(gates[pos : pos + length])
                # non-None implies lexicographic (twq, len) descent; equal-length
                # replacements are safe because the objective is monotone.
                if candidate is not None:
                    gates[pos : pos + length] = candidate
                    n = len(gates)
                    count += 1
                    replaced = True
                    break
            if replaced:
                pos = max(0, pos - 1)
            else:
                pos += 1
        total += count
        if count == 0:
            break
    return total


def sweep(gates: list[GateInstance], num_qubits: int, db: ReductionDatabase, max_len: int = 8) -> int:
    """Left-to-right exhaustive sweep over every window until no reduction found.

Greedy longest-first within each window start; returns replacements applied."""
    total = 0
    while True:
        count = 0
        pos = 0
        n = len(gates)
        while pos < n:
            hi = min(max_len, n - pos)
            replaced = False
            for length in range(hi, 1, -1):
                cand = db.try_reduce(gates[pos : pos + length])
                if cand is not None and len(cand) < length:
                    gates[pos : pos + length] = cand
                    n = len(gates)
                    count += 1
                    replaced = True
                    break
            if replaced:
                pos = max(0, pos - 1)
            else:
                pos += 1
        total += count
        if count == 0:
            break
    return total


def cluster_single_qubit(gates: list[GateInstance], num_qubits: int) -> None:
    """Gather same-wire single-qubit gates between two-qubit barriers on that wire."""
    for w in range(num_qubits):
        i = 0
        n = len(gates)
        while i < n:
            g = gates[i]
            if len(g.qubits) == 2:
                i += 1
                continue
            singles = []
            j = i
            while j < n:
                g2 = gates[j]
                if len(g2.qubits) == 2 and w in g2.qubits:
                    break
                if len(g2.qubits) == 1 and g2.qubits[0] == w:
                    singles.append(j)
                j += 1
            if len(singles) > 1:
                gs = [gates[k] for k in singles]
                for k in reversed(singles):
                    del gates[k]
                gates[singles[0] : singles[0]] = gs
                n = len(gates)
                i = singles[0] + len(gs)
            else:
                i = j


def reduce_single_wire_runs(gates: list[GateInstance], one_wire_graph, max_run: int | None = None) -> int:
    """Collapse maximal same-wire single-qubit runs using the complete 1-wire graph."""
    n = len(gates)
    i = 0
    total = 0
    cap = one_wire_graph.max_depth if max_run is None else max_run
    while i < n:
        g = gates[i]
        if len(g.qubits) != 1:
            i += 1
            continue
        w = g.qubits[0]
        j = i
        while j < n:
            g2 = gates[j]
            if len(g2.qubits) == 1 and g2.qubits[0] == w:
                j += 1
            else:
                break
        m = j - i
        lengths = (m,) if m <= cap else range(min(cap, m), 1, -1)
        replaced = False
        for length in lengths:
            cand = one_wire_graph.try_reduce(gates[i : i + length])
            if cand is not None and len(cand) < length:
                gates[i : i + length] = cand
                n = len(gates)
                total += 1
                replaced = True
                break
        if not replaced:
            i = j
    return total


def transport_shuffle(
    gates: list[GateInstance],
    num_qubits: int,
    rng: random.Random,
    cache: dict,
    prob: float = 0.5,
    direction: int = 1,
) -> None:
    """Move each gate right (direction=1) or left (direction=-1) across gates it commutes with."""
    n = len(gates)
    idx = list(range(n)) if direction > 0 else list(range(n - 1, -1, -1))
    for i in idx:
        if rng.random() >= prob:
            continue
        gate = gates[i]
        j = i + direction
        while 0 <= j < n:
            key = (_gate_key(gate), _gate_key(gates[j]))
            can = cache.get(key)
            if can is None:
                can = gates_commute(num_qubits, gate, gates[j])
                cache[key] = can
            if not can:
                break
            j += direction
        j -= direction
        if j != i:
            gates.pop(i)
            gates.insert(j, gate)


def rz_global_pass(gates: list[GateInstance], one_wire_graph, max_iters: int = 8) -> int:
    """NISQ structural pass: transport RZ gates across CZ layers.

Stable-partitions RZ per wire between RX barriers, collapses the runs, and iterates to a fixpoint; returns replacements."""
    total = 0
    for _ in range(max_iters):
        if not gates:
            break
        num_wires = max(q for gate in gates for q in gate.qubits) + 1
        result = list(gates)
        for w in range(num_wires):
            out: list[GateInstance] = []
            pending: list[GateInstance] = []
            for gate in result:
                if gate.name == "RZ" and gate.qubits[0] == w:
                    pending.append(gate)
                elif gate.name == "RX" and gate.qubits[0] == w:
                    out.extend(pending)
                    pending = []
                    out.append(gate)
                else:
                    out.append(gate)
            out.extend(pending)
            result = out
        gates[:] = result
        collapsed = reduce_single_wire_runs(gates, one_wire_graph)
        total += collapsed
        if collapsed == 0:
            break
    return total


def rz_global_pass_fixpoint(gates: list[GateInstance], one_wire_graph, max_iters: int = 8) -> int:
    """Iterate :func:`rz_global_pass` until nothing more collapses.

    A collapse can expose new RZ adjacencies, so a single pass need not reach
    the fixpoint; iterating is cheap (O(n) per pass).
    """
    total = 0
    for _ in range(max_iters):
        found = rz_global_pass(gates, one_wire_graph)
        total += found
        if found == 0:
            break
    return total


_DIAGONAL_GATES = ("RZ", "CZ")


def collapse_cz_parity(gates: list[GateInstance], num_qubits: int) -> int:
    """Exactly collapse each CZ pair's occurrences to parity (0 or 1 gate).

RZ and CZ are both diagonal in the computational basis regardless of which
wires they touch, so a CZ(i, j) commutes with every gate except a
non-diagonal (RX) gate on wire i or j. Between two such barriers, every
CZ(i, j) in the run can be walked adjacent to its neighbor and cancelled
(CZ^2 = I) without touching anything else in the window, independent of
literal adjacency or what other gates (other-pair CZs, RZs, other-wire RX)
sit between them. This subsumes zx_cancellations' literal-adjacency case
and additionally collapses runs the current pipeline never reorders into
adjacency at all."""
    to_remove: set[int] = set()
    pairs = {tuple(sorted(g.qubits)) for g in gates if g.name == "CZ"}
    for i, j in pairs:
        window: list[int] = []
        for idx, gate in enumerate(gates):
            if gate.name == "CZ" and tuple(sorted(gate.qubits)) == (i, j):
                window.append(idx)
                continue
            if gate.name not in _DIAGONAL_GATES and (i in gate.qubits or j in gate.qubits):
                if len(window) % 2 == 1:
                    to_remove.update(window[:-1])
                else:
                    to_remove.update(window)
                window = []
        if len(window) % 2 == 1:
            to_remove.update(window[:-1])
        else:
            to_remove.update(window)
    if not to_remove:
        return 0
    removed = len(to_remove)
    gates[:] = [gate for idx, gate in enumerate(gates) if idx not in to_remove]
    return removed


def _random_escape(
    gates: list[GateInstance],
    num_qubits: int,
    db: ReductionDatabase,
    max_block_len: int,
    rng: random.Random,
    num_tries: int = 64,
) -> int:
    n = len(gates)
    count = 0
    for _ in range(num_tries):
        if n < 2:
            break
        start = rng.randrange(0, n - 1)
        hi = min(max_block_len, n - start)
        length = rng.randint(2, hi)
        candidate = db.try_reduce(gates[start : start + length])
        if candidate is not None and len(candidate) < length:
            gates[start : start + length] = candidate
            n = len(gates)
            count += 1
    return count


def reduce_with_database(
    gates: list[GateInstance],
    num_qubits: int,
    db: ReductionDatabase,
    budget_sec: float = 5.0,
    seed: int = 0,
    max_block_len: int = 7,
    shuffle_prob: float = 0.7,
    stall_limit: int = 6,
) -> tuple[list[GateInstance], ReductionStats]:
    """Local-term-replacement reducer driven by exhaustive sweeps.

Combines commutative shuffling with exhaustive local sweeps over windows up to max_block_len against the precomputed ReductionDatabase."""
    rng = random.Random(seed)
    working = _snap_input_gates(gates, db)
    commute_cache: dict = {}
    t0 = time.time()
    start_len = len(gates)
    replacements = 0
    iterations = 0
    stalled = 0
    best_len = len(working)

    while time.time() - t0 < budget_sec:
        iterations += 1
        if rng.random() < shuffle_prob:
            shuffle_commuting_pairs(working, num_qubits, rng, commute_cache)

        sweep_count = _sweep_reduce(working, num_qubits, db, max_block_len)
        replacements += sweep_count

        if len(working) < best_len:
            best_len = len(working)
            stalled = 0
        else:
            stalled += 1
            if stalled >= stall_limit:
                extra = _random_escape(working, num_qubits, db, max_block_len, rng)
                replacements += extra
                if len(working) < best_len:
                    best_len = len(working)
                    stalled = 0
                elif extra == 0:
                    break

    stats = ReductionStats(
        start_len=start_len,
        end_len=len(working),
        iterations=iterations,
        replacements=replacements,
        runtime_sec=time.time() - t0,
    )
    return working, stats


def reduce_with_lookup(
    gates: list[GateInstance],
    num_qubits: int,
    local_qubits: int = 3,
    max_block_len: int = 7,
    graph_depth: int = 4,
    iterations: int = 15000,
    seed: int = 0,
) -> tuple[list[GateInstance], ReductionStats]:
    """Paper-style V2 reducer (random block sampling against a wire-count database).

Iteration-budgeted analogue of reduce_random_sampling, kept for the protocol benchmark and demo-port scripts."""
    from .database import load_or_build_database

    if local_qubits < 2:
        raise ValueError("local_qubits must be >= 2")

    rng = random.Random(seed)
    commute_cache: dict = {}
    t0 = time.time()
    start_len = len(gates)
    replacements = 0
    db = load_or_build_database(
        "ion_trap", {w: graph_depth for w in range(1, local_qubits + 1)}, verbose=False
    )
    working = _snap_input_gates(gates, db)

    for _ in range(iterations):
        if len(working) < 2:
            break
        if rng.random() < 0.7:
            shuffle_commuting_pairs(working, num_qubits, rng, commute_cache)
        start = rng.randrange(0, len(working) - 1)
        hi = min(max_block_len, len(working) - start)
        length = rng.randint(2, hi)
        candidate = db.try_reduce(working[start : start + length])
        if candidate is not None and len(candidate) < length:
            working[start : start + length] = candidate
            replacements += 1

    stats = ReductionStats(
        start_len=start_len,
        end_len=len(working),
        iterations=iterations,
        replacements=replacements,
        runtime_sec=time.time() - t0,
    )
    return working, stats


def reduce_random_sampling(
    gates: list[GateInstance],
    num_qubits: int,
    db: ReductionDatabase,
    budget_sec: float = 5.0,
    seed: int = 0,
    max_block_len: int = 7,
) -> tuple[list[GateInstance], ReductionStats]:
    """Paper-style V2: random sub-block sampling with database retrieval."""
    rng = random.Random(seed)
    working = _snap_input_gates(gates, db)
    commute_cache: dict = {}
    t0 = time.time()
    start_len = len(gates)
    replacements = 0
    iterations = 0

    while time.time() - t0 < budget_sec:
        iterations += 1
        if rng.random() < 0.7:
            shuffle_commuting_pairs(working, num_qubits, rng, commute_cache)

        n = len(working)
        if n < 2:
            break
        start = rng.randrange(0, n - 1)
        hi = min(max_block_len, n - start)
        length = rng.randint(2, hi)
        candidate = db.try_reduce(working[start : start + length])
        if candidate is not None and len(candidate) < length:
            working[start : start + length] = candidate
            replacements += 1

    stats = ReductionStats(
        start_len=start_len,
        end_len=len(working),
        iterations=iterations,
        replacements=replacements,
        runtime_sec=time.time() - t0,
    )
    return working, stats


def reduce_random_sampling_gated(
    gates: list[GateInstance],
    num_qubits: int,
    db: ReductionDatabase,
    budget_sec: float = 5.0,
    seed: int = 0,
    max_block_len: int = 7,
) -> tuple[list[GateInstance], ReductionStats, "RfGate"]:
    """Paper V3: random sampling with an RF-gated database lookup.

Wraps db in RfGatedDatabase so blocks predicted irreducible skip the lookup (exact memo cache + lazily-trained classifier). Returns (reduced, stats, gate)."""
    from .rf_gate import RfGate, RfGatedDatabase

    gate = RfGate()
    reduced, stats = reduce_random_sampling(
        gates, num_qubits, RfGatedDatabase(db, gate), budget_sec, seed, max_block_len
    )
    return reduced, stats, gate


def reduce_circuit(
    gates: list[GateInstance],
    num_qubits: int,
    db: ReductionDatabase,
    budget_s: float,
    seed: int,
    max_block_len: int = 8,
    max_passes: int = 20000,
    escape_every: int = 3,
    prefer: dict[str, float] | None = None,
    rz_pass: bool = False,
    cz_pass: bool = False,
    cost_aware: bool = False,
    algebraic: bool = False,
    zx: bool = False,
    use_batched: bool = False,
    use_fast_sweep: bool = False,
    dag_compact: bool = False,
    dag_max_wires: int = 3,
) -> tuple[list[GateInstance], int, int]:
    """Strong reducer: cluster + collapse + sweep, transport shuffle, escape.

    The escape move resamples an irreducible window with a structurally
    different equivalent word, then re-sweeps; kept only if it does not
    worsen the circuit.  ``rz_pass`` runs the NISQ RZ-across-CZ pass to a
    fixpoint; ``cz_pass`` runs :func:`collapse_cz_parity`, exactly collapsing
    each CZ pair's occurrences between non-diagonal (RX) barriers on either
    wire to parity, independent of literal adjacency; ``cost_aware``
    minimizes (two-qubit count, length);
    ``algebraic``/``zx`` enable the prepass rules (prepass.py);
    ``use_batched`` swaps the scalar sweep for the bit-identical batched
    sweep (batched.py, scripts/check_batched_vs_scalar.py).  ``use_fast_sweep``
    swaps it for :func:`_sweep_reduce_fast` (incremental per-position
    unitary construction reusing each graph's cached token matrices instead
    of rebuilding every tried length from scratch); it finds a match at the
    same length as the scalar sweep for every window (verified exhaustively
    in scripts/check_fast_sweep_parity.py, 5000+ windows, zero mismatches)
    but is not guaranteed to reach byte-identical circuits, since a
    different (equally valid) local wire-relabeling can land on a
    different-but-equally-short realization of a window -- see
    scripts/check_fast_sweep.py.  ``dag_compact``
    deterministically reorders the circuit before each sweep so every
    <=dag_max_wires-wire block becomes contiguous (dag.py), exposing windows
    to the sweep that transport_shuffle only finds by chance; it is a valid
    reordering of the circuit's dependency DAG, so it never changes the
    unitary (scripts/check_dag_compact.py).
    Returns (reduced, passes, replacements).
    """
    rng = random.Random(seed)
    working = list(gates)
    if algebraic or zx:
        from .prepass import apply_prepass

        gs = gateset_for(db.gate_set_name)
        angles = db.angles if db.angles is not None else gs.angles
        two = db.two_qubit_angles if db.two_qubit_angles is not None else gs.two_angles
        working, _ = apply_prepass(working, db.gate_set_name, angles, two, num_qubits, zx=zx)
    best = list(working)
    cache: dict = {}
    t0 = time.time()
    passes = 0
    reduced = 0
    one_wire = db.graphs.get(1)

    def sweep_fn(gates_list, num_qubits_, db_, max_block_len_) -> int:
        """Cost-aware pushes both (twq, len) and pure-length objectives;
        length-only matches the paper's metric.  Batched mode (length only)
        is bit-identical to the scalar sweep (scripts/check_batched_vs_scalar.py)."""
        if cost_aware:
            return _sweep_reduce_cost(gates_list, num_qubits_, db_, max_block_len_) + _sweep_reduce(
                gates_list, num_qubits_, db_, max_block_len_
            )
        if use_batched:
            from .batched import batched_sweep

            return batched_sweep(gates_list, num_qubits_, db_, max_block_len_)
        if use_fast_sweep:
            return _sweep_reduce_fast(gates_list, num_qubits_, db_, max_block_len_)
        return _sweep_reduce(gates_list, num_qubits_, db_, max_block_len_)

    def done() -> bool:
        return time.time() - t0 > budget_s or passes >= max_passes

    def compact(gates_list: list[GateInstance]) -> None:
        if not dag_compact:
            return
        from .dag import compact_by_blocks

        gates_list[:] = compact_by_blocks(gates_list, max_wires=dag_max_wires)

    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass_fixpoint(working, one_wire)
    if cz_pass:
        reduced += collapse_cz_parity(working, num_qubits)
    compact(working)
    reduced += sweep_fn(working, num_qubits, db, max_block_len)
    if done():
        return working, passes, reduced
    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass_fixpoint(working, one_wire)
    if cz_pass:
        reduced += collapse_cz_parity(working, num_qubits)
    compact(working)
    reduced += sweep_fn(working, num_qubits, db, max_block_len)
    if done():
        return working, passes, reduced

    stall = 0
    while not done():
        passes += 1
        # dag_compact only runs as a deterministic pre-pass (above) and in
        # the escape trial below; it is intentionally *not* called here.
        # It is a pure function of the current gate list, so re-running it on
        # an already-compacted, otherwise-unchanged list is a no-op -- which
        # would silently displace transport_shuffle's randomized diversity
        # (the main loop's actual source of escaping local optima) without
        # adding anything back, on the (common) passes where nothing
        # changed since the last compaction.
        transport_shuffle(working, num_qubits, rng, cache, direction=1 if passes % 2 else -1)
        if rz_pass:
            reduced += rz_global_pass_fixpoint(working, one_wire)
        if cz_pass:
            reduced += collapse_cz_parity(working, num_qubits)
        found = sweep_fn(working, num_qubits, db, max_block_len)
        if found == 0:
            cluster_single_qubit(working, num_qubits)
            found += reduce_single_wire_runs(working, one_wire)
            if rz_pass:
                found += rz_global_pass_fixpoint(working, one_wire)
            if cz_pass:
                found += collapse_cz_parity(working, num_qubits)
            compact(working)
            found += sweep_fn(working, num_qubits, db, max_block_len)
        reduced += found
        if found == 0:
            stall += 1
        else:
            stall = 0
            if len(working) < len(best):
                best = list(working)

        if stall >= escape_every and not done():
            # Escape: resample an irreducible window with a structurally
            # different equivalent word, re-sweep, and keep the trial only if
            # it strictly improves the circuit.  The trial is built on a copy
            # so a failed escape can never corrupt the working circuit.
            improved = False
            for _ in range(8):
                if len(working) < 2:
                    break
                start = rng.randrange(0, len(working) - 1)
                hi = min(max_block_len, len(working) - start)
                length = rng.randint(2, hi)
                candidate = db.try_reduce_escape(working[start : start + length], rng, slack=3, prefer=prefer)
                if candidate is None:
                    continue
                trial = list(working)
                trial[start : start + length] = candidate
                cluster_single_qubit(trial, num_qubits)
                reduce_single_wire_runs(trial, one_wire)
                compact(trial)
                sweep_fn(trial, num_qubits, db, max_block_len)
                if len(trial) < len(working):
                    working = trial
                    improved = True
                    break
            if improved:
                reduced += 1
                if len(working) < len(best):
                    best = list(working)
                stall = 0
            else:
                stall += 1
                if stall >= 8:
                    working = list(best)
                    rng = random.Random(rng.randrange(2**30))
                    stall = 0
    return working, passes, reduced
