"""Exact, cost-aware circuit reduction for Clifford gate pools.

Drives SymplecticDatabase with exact tableau identity, a (two-qubit, length) cost objective, and exact verification -- no numeric tolerance anywhere."""

from __future__ import annotations

import random
import time
from dataclasses import dataclass

from .config import GateInstance
from .exact_database import SymplecticDatabase
from .reducer import (
    _gate_key,
    cluster_single_qubit,
    gates_commute,
    reduce_single_wire_runs,
    rz_global_pass,
)
from .symplectic import SignedTableau, circuits_equal_exact


@dataclass
class ExactStats:
    start_len: int
    end_len: int
    start_two_qubit: int
    end_two_qubit: int
    iterations: int = 0
    replacements: int = 0
    runtime_sec: float = 0.0


def _count_two_qubit(gates: list[GateInstance]) -> int:
    return sum(1 for g in gates if len(g.qubits) == 2)


def _improves_cost(candidate: list[GateInstance], block: list[GateInstance]) -> bool:
    """True if ``candidate`` is lexicographically better in (twq, len)."""
    if len(candidate) < len(block):
        return True
    return len(candidate) == len(block) and _count_two_qubit(candidate) < _count_two_qubit(block)


def _cost_weight(gate_name: str) -> float:
    """Return a weight for a gate type that models hardware cost.

    Lower is cheaper: single-qubit gates < two-qubit gates.
    """
    weights = {
        "RX": 0.01,
        "RY": 0.01,
        "RZ": 0.01,
        "CZ": 5.0,
        "RXX": 5.0,
        "RXY": 5.0,
        "RYY": 5.0,
    }
    return weights.get(gate_name, 1.0)


def _cost_weighted_reduction(
    block: list[GateInstance],
    db: SymplecticDatabase,
    prefer: dict[str, float] | None = None,
) -> tuple[list[GateInstance], float | None]:
    """Reduce a block minimizing (two-qubit count, length, weighted cost).

Returns (reduced_block, cost) or (None, None) if no improvement."""
    # Get the SymplecticGraph for this block
    graph = db._graph_for(block)
    if graph is None:
        return None, None

    # Remap to graph-local wire indices
    remapped, _, reverse = db._remap(block)
    key = graph.block_key(remapped)
    if key is None or key not in graph.buckets:
        return None, None

    # Get all Pareto-optimal candidates
    cands = graph.alts.get(key, [])
    if not cands:
        return None, None

    block_len = len(block)
    block_twq = _count_two_qubit(block)
    if prefer is None:
        prefer = {}

    best = None
    best_key = None
    for (twq, ln, chain) in cands:
        # Accept strictly shorter candidates, and equal-length candidates that
        # strictly reduce the two-qubit count (so cost-aware mode can rewrite
        # e.g. CZ-heavy words into CZ-sparse ones of the same length, which the
        # paper's length-only objective never attempts).
        if ln > block_len:
            continue
        if ln == block_len and twq >= block_twq:
            continue
        # Decode and compute cost using prefer weights
        decoded = graph.pool.decode(list(chain))
        cand_cost = sum(prefer.get(g.name, 1.0) for g in decoded)
        # Pareto-optimal by (two-qubit count, length) then cost
        comparison = (twq, ln, cand_cost)
        if best is None or comparison < best_key:
            best = decoded
            best_key = comparison

    if best is None:
        return None, None

    # Restore to original wire indices
    restored = db._restore(best, reverse)
    return restored, best_key[2]


def _sweep_reduce_cost(
    gates: list[GateInstance],
    db: SymplecticDatabase,
    max_block_len: int,
    prefer: dict[str, float] | None = None,
) -> int:
    """Exhaustive sweep minimizing (two-qubit count, length, cost).

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
                candidate, cost = _cost_weighted_reduction(
                    gates[pos : pos + length], db, prefer
                )
                # candidate is non-None only when it strictly improves
                # (two-qubit count, length); equal-length replacements are safe
                # because the lexicographic objective is monotone and bounded.
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


def _sweep_reduce_len(
    gates: list[GateInstance],
    db: SymplecticDatabase,
    max_block_len: int,
    prefer: dict[str, float] | None = None,
) -> int:
    """Exhaustive sweep minimizing length only (ignores prefer)."""
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


def _random_escape_cost(
    gates: list[GateInstance],
    db: SymplecticDatabase,
    max_block_len: int,
    rng: random.Random,
    num_tries: int = 64,
    prefer: dict[str, float] | None = None,
) -> int:
    """Escape: resample irreducible windows with structurally different words.

Tries num_tries random windows, accepting any candidate that reduces (two-qubit, length) under the hardware-cost model."""
    n = len(gates)
    count = 0

    # Precompute token indices for quick two-qubit count estimation.
    for _ in range(num_tries):
        if n < 2:
            break
        start = rng.randrange(0, n - 1)
        hi = min(max_block_len, n - start)
        length = rng.randint(2, hi)

        block = gates[start : start + length]
        candidate = db.try_reduce_cost(block)

        # Accept cost-aware reductions that strictly improve (twq, len); only
        # apply escape resampling when it is strictly shorter (the escape move
        # is meant to perturb, not grow, the circuit in this hot path).
        if candidate is not None and _improves_cost(candidate, block):
            gates[start : start + length] = candidate
            n = len(gates)
            count += 1
            continue
        candidate = db.try_reduce_escape(block, rng, prefer=prefer)
        if candidate is not None and len(candidate) < length:
            gates[start : start + length] = candidate
            n = len(gates)
            count += 1
    return count


def reduce_circuit_exact(
    gates: list[GateInstance],
    num_qubits: int,
    db: SymplecticDatabase,
    budget_s: float = 5.0,
    seed: int = 0,
    max_block_len: int = 8,
    cost_aware: bool = True,
    rz_pass: bool = False,
    max_passes: int = 20000,
    prefer: dict[str, float] | None = None,
) -> tuple[list[GateInstance], ExactStats]:
    """Strong exact reducer: cluster + collapse + cost-aware sweep + transport.

cost_aware=True minimizes (two-qubit count, length); prefer maps gate names to cost weights (lower = more preferred). Equivalence is bit-exact up to global phase."""
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

    # Cost-aware preferences: RXX/CZ are expensive, single-qubit gates are cheap
    if cost_aware and prefer is None:
        prefer = {
            "RX": 0.01,
            "RY": 0.01,
            "RZ": 0.01,
            "CZ": 5.0,
            "RXX": 5.0,
            "RXY": 5.0,
            "RYY": 5.0,
        }

    def do_sweeps(gates_list: list[GateInstance], db_, max_block_len_) -> int:
        """Cost-aware mode runs the (twq, len) sweep and then the pure-length
        sweep, so both objectives are pushed on every pass; length-only mode
        keeps the paper's single objective."""
        if cost_aware:
            return _sweep_reduce_cost(gates_list, db_, max_block_len_, prefer) + _sweep_reduce_len(
                gates_list, db_, max_block_len_
            )
        return _sweep_reduce_len(gates_list, db_, max_block_len_)

    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass(working, one_wire)
    reduced += do_sweeps(working, db, max_block_len)
    if done():
        return working, _stats(gates, working, passes, reduced, t0)
    cluster_single_qubit(working, num_qubits)
    reduced += reduce_single_wire_runs(working, one_wire)
    if rz_pass:
        reduced += rz_global_pass(working, one_wire)
    reduced += do_sweeps(working, db, max_block_len)
    if done():
        return working, _stats(gates, working, passes, reduced, t0)

    stall = 0
    while not done():
        passes += 1
        # transport shuffle using exact commutation via memoized numeric checks
        _transport_shuffle(working, num_qubits, rng, cache, direction=1 if passes % 2 else -1)
        if rz_pass:
            reduced += rz_global_pass(working, one_wire)
        found = do_sweeps(working, db, max_block_len)
        if found == 0:
            cluster_single_qubit(working, num_qubits)
            found += reduce_single_wire_runs(working, one_wire)
            if rz_pass:
                found += rz_global_pass(working, one_wire)
            found += do_sweeps(working, db, max_block_len)
        reduced += found
        if found == 0:
            stall += 1
        else:
            stall = 0
            if len(working) < len(best):
                best = list(working)

        if stall >= 3 and not done():
            # Escape with a structurally different equivalent word; keep only
            # strictly improving trials (built on a copy).
            improved = False
            for _ in range(8):
                if len(working) < 2:
                    break
                start = rng.randrange(0, len(working) - 1)
                hi = min(max_block_len, len(working) - start)
                length = rng.randint(2, hi)
                candidate = db.try_reduce_escape(working[start : start + length], rng, prefer=prefer)
                if candidate is None:
                    continue
                trial = list(working)
                trial[start : start + length] = candidate
                cluster_single_qubit(trial, num_qubits)
                reduce_single_wire_runs(trial, one_wire)
                do_sweeps(trial, db, max_block_len)
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
    return working, _stats(gates, working, passes, reduced, t0)


def _stats(
    start: list[GateInstance],
    end: list[GateInstance],
    passes: int,
    replacements: int,
    t0: float,
) -> ExactStats:
    return ExactStats(
        start_len=len(start),
        end_len=len(end),
        start_two_qubit=_count_two_qubit(start),
        end_two_qubit=_count_two_qubit(end),
        iterations=passes,
        replacements=replacements,
        runtime_sec=time.time() - t0,
    )


def _transport_shuffle(
    gates: list[GateInstance],
    num_qubits: int,
    rng: random.Random,
    cache: dict,
    prob: float = 0.5,
    direction: int = 1,
) -> None:
    """Move each gate across gates it commutes with (exact commutation checks)."""
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


def verify_exact(
    original: list[GateInstance], reduced: list[GateInstance], num_qubits: int
) -> bool:
    """Bit-exact verification up to global phase via signed tableaux."""
    return circuits_equal_exact(original, reduced, num_qubits)
