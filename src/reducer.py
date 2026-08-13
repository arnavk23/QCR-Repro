from __future__ import annotations

import random
import time
from dataclasses import dataclass

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
    """Iterate :func:`rz_global_pass` to a fixpoint.

    A single pass gathers RZ runs per wire and collapses them, but a collapse
    can expose new RZ adjacencies (the sweep that follows can also merge gates
    that leave fresh RZ runs).  Iterating to a fixpoint is cheap (O(n) per
    pass) and guarantees the pass no longer finds anything to remove.
    """
    total = 0
    for _ in range(max_iters):
        found = rz_global_pass(gates, one_wire_graph)
        total += found
        if found == 0:
            break
    return total


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
    working = list(gates)
    commute_cache: dict = {}
    t0 = time.time()
    start_len = len(gates)
    replacements = 0
    db = load_or_build_database(
        "ion_trap", {w: graph_depth for w in range(1, local_qubits + 1)}, verbose=False
    )

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
    cost_aware: bool = False,
    algebraic: bool = False,
    zx: bool = False,
    use_batched: bool = False,
) -> tuple[list[GateInstance], int, int]:
    """Strong reducer: cluster + collapse + sweep, transport shuffle, escape.

    The escape move resamples an irreducible window with a structurally
    different (possibly longer) equivalent word from the compute graph's cycle
    structure, then re-sweeps; the change is kept only if it does not worsen
    the circuit.  ``rz_pass`` enables the NISQ-specific RZ-across-CZ transport
    pass (iterated to a fixpoint).  ``cost_aware=True`` minimizes
    (two-qubit count, length) per window replacement instead of length alone
    (matches the exact engine's objective).

    ``algebraic`` / ``zx`` enable the cheap pre-passes (qcr_repro.prepass) that
    shrink the input *before* the database loop using exact rotation-fusion and
    ZX-cancellation rules (fused angles are always pool-representable).
    ``use_batched`` swaps the per-window scalar sweep for the vectorized batched
    sweep (qcr_repro.batched), which is bit-identical in results (see
    scripts/check_batched_vs_scalar.py).
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
        """Cost-aware mode pushes both the (twq, len) and the pure-length
        objective on every pass; length-only mode matches the paper's metric.
        Batched mode (length objective only) is bit-identical to the scalar
        sweep (see scripts/check_batched_vs_scalar.py)."""
        if cost_aware:
            return _sweep_reduce_cost(gates_list, num_qubits_, db_, max_block_len_) + _sweep_reduce(
                gates_list, num_qubits_, db_, max_block_len_
            )
        if use_batched:
            from .batched import batched_sweep

            return batched_sweep(gates_list, num_qubits_, db_, max_block_len_)
        return _sweep_reduce(gates_list, num_qubits_, db_, max_block_len_)

    def done() -> bool:
        return time.time() - t0 > budget_s or passes >= max_passes

    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass_fixpoint(working, one_wire)
    reduced += sweep_fn(working, num_qubits, db, max_block_len)
    if done():
        return working, passes, reduced
    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass_fixpoint(working, one_wire)
    reduced += sweep_fn(working, num_qubits, db, max_block_len)
    if done():
        return working, passes, reduced

    stall = 0
    while not done():
        passes += 1
        transport_shuffle(working, num_qubits, rng, cache, direction=1 if passes % 2 else -1)
        if rz_pass:
            reduced += rz_global_pass_fixpoint(working, one_wire)
        found = sweep_fn(working, num_qubits, db, max_block_len)
        if found == 0:
            cluster_single_qubit(working, num_qubits)
            found += reduce_single_wire_runs(working, one_wire)
            if rz_pass:
                found += rz_global_pass_fixpoint(working, one_wire)
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
