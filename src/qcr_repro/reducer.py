from __future__ import annotations

import random
import time
from dataclasses import dataclass

from .compute_graph import ReductionDatabase
from .config import GateInstance, gateset_for
from .gates import circuit_unitary, embedded_gate_matrix
from .tokenizer import TokenPool
from .unitary_utils import equivalent_up_to_global_phase


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


def sweep(gates: list[GateInstance], num_qubits: int, db: ReductionDatabase, max_len: int = 8) -> int:
    """Left-to-right exhaustive sweep over every window until no reduction found.

    Greedy longest-first within each window start; returns the number of
    replacements applied.
    """
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


def rz_global_pass(gates: list[GateInstance], one_wire_graph) -> int:
    """NISQ-specific structural pass: transport RZ gates across CZ layers.

    RZ commutes with CZ (both diagonal) and with RZ on any wire, so every RZ
    gate may be moved freely across CZ gates and across single-qubit gates on
    other wires; the only obstruction on a wire is an RX on that same wire.
    This pass performs a stable partition per wire: within each maximal region
    delimited by RX_w barriers, all RZ_w gates are gathered together, forming
    per-wire RZ runs which then collapse exactly.  Returns replacements.
    """
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
    return reduce_single_wire_runs(gates, one_wire_graph)


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

    Combines commutative shuffling (global sampling) with exhaustive local
    sweeps over all windows up to ``max_block_len`` (local optimization) using the
    precomputed ReductionDatabase.
    """
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
) -> tuple[list[GateInstance], int, int]:
    """Strong reducer: cluster + collapse + sweep, transport shuffle, and escape.

    The escape move resamples an irreducible window with a structurally
    different (possibly longer) equivalent word from the compute graph's cycle
    structure, then re-sweeps; the change is kept only if it does not worsen
    the circuit.  ``rz_pass`` enables the NISQ-specific RZ-across-CZ transport
    pass.  Returns (reduced, passes, replacements).
    """
    rng = random.Random(seed)
    working = list(gates)
    best = list(gates)
    cache: dict = {}
    t0 = time.time()
    passes = 0
    reduced = 0
    one_wire = db.graphs.get(1)

    def done() -> bool:
        return time.time() - t0 > budget_s or passes >= max_passes

    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass(working, one_wire)
    reduced += sweep(working, num_qubits, db, max_block_len)
    if done():
        return working, passes, reduced
    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass(working, one_wire)
    reduced += sweep(working, num_qubits, db, max_block_len)
    if done():
        return working, passes, reduced

    stall = 0
    while not done():
        passes += 1
        transport_shuffle(working, num_qubits, rng, cache, direction=1 if passes % 2 else -1)
        if rz_pass:
            reduced += rz_global_pass(working, one_wire)
        found = sweep(working, num_qubits, db, max_block_len)
        if found == 0:
            cluster_single_qubit(working, num_qubits)
            found += reduce_single_wire_runs(working, one_wire)
            if rz_pass:
                found += rz_global_pass(working, one_wire)
            found += sweep(working, num_qubits, db, max_block_len)
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
                sweep(trial, num_qubits, db, max_block_len)
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
