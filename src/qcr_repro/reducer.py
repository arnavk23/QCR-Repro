from __future__ import annotations

import random
from dataclasses import dataclass
from math import pi

from .compute_graph import ComputeGraphBuilder
from .config import GateInstance
from .gates import circuit_unitary, embedded_gate_matrix
from .tokenizer import TokenPool
from .unitary_utils import equivalent_up_to_global_phase, unitary_key


@dataclass
class ReductionStats:
    iterations: int
    replacements: int
    start_len: int
    end_len: int


def _commute_pair(num_qubits: int, a: GateInstance, b: GateInstance, atol: float = 1e-8) -> bool:
    if set(a.qubits).isdisjoint(set(b.qubits)):
        return True
    ua = embedded_gate_matrix(num_qubits, a)
    ub = embedded_gate_matrix(num_qubits, b)
    return equivalent_up_to_global_phase(ua @ ub, ub @ ua, atol=atol)


def _gate_cache_key(gate: GateInstance) -> tuple[str, tuple[int, ...], float | None]:
    theta = None if gate.theta is None else round(float(gate.theta), 8)
    return (gate.name, gate.qubits, theta)


def _try_shuffle_commuting_pairs(
    gates: list[GateInstance],
    num_qubits: int,
    rng: random.Random,
    commute_cache: dict[tuple[tuple[str, tuple[int, ...], float | None], tuple[str, tuple[int, ...], float | None]], bool],
) -> None:
    i = 0
    while i < len(gates) - 1:
        if rng.random() < 0.6:
            k1 = _gate_cache_key(gates[i])
            k2 = _gate_cache_key(gates[i + 1])
            pair_key = (k1, k2)
            if pair_key not in commute_cache:
                commute_cache[pair_key] = _commute_pair(num_qubits, gates[i], gates[i + 1])
            can_swap = commute_cache[pair_key]
        else:
            can_swap = False

        if can_swap:
            gates[i], gates[i + 1] = gates[i + 1], gates[i]
            i += 2
        else:
            i += 1


def _remap_block(block: list[GateInstance]) -> tuple[list[GateInstance], dict[int, int], dict[int, int]]:
    wires = sorted({q for gate in block for q in gate.qubits})
    forward = {wire: idx for idx, wire in enumerate(wires)}
    reverse = {idx: wire for wire, idx in forward.items()}

    remapped: list[GateInstance] = []
    for gate in block:
        new_qubits = tuple(forward[q] for q in gate.qubits)
        if len(new_qubits) == 2 and new_qubits[0] > new_qubits[1]:
            new_qubits = (new_qubits[1], new_qubits[0])
        remapped.append(GateInstance(name=gate.name, qubits=new_qubits, theta=gate.theta))

    return remapped, forward, reverse


def _restore_block(block: list[GateInstance], reverse: dict[int, int]) -> list[GateInstance]:
    restored: list[GateInstance] = []
    for gate in block:
        qubits = tuple(reverse[q] for q in gate.qubits)
        if len(qubits) == 2 and qubits[0] > qubits[1]:
            qubits = (qubits[1], qubits[0])
        restored.append(GateInstance(name=gate.name, qubits=qubits, theta=gate.theta))
    return restored


def reduce_with_lookup(
    gates: list[GateInstance],
    num_qubits: int,
    local_qubits: int = 3,
    max_block_len: int = 7,
    graph_depth: int = 4,
    iterations: int = 15000,
    seed: int = 0,
) -> tuple[list[GateInstance], ReductionStats]:
    if local_qubits < 2:
        raise ValueError("local_qubits must be >= 2")

    rng = random.Random(seed)
    working = list(gates)
    commute_cache: dict[
        tuple[tuple[str, tuple[int, ...], float | None], tuple[str, tuple[int, ...], float | None]],
        bool,
    ] = {}

    pool = TokenPool(
        num_qubits=local_qubits,
        gate_set="ion_trap",
        angles=(-pi / 2, pi / 2),
        rxx_angles=(pi / 2,),
    )
    graph = ComputeGraphBuilder(pool, max_depth=graph_depth).build()

    replacements = 0
    for _ in range(iterations):
        if len(working) < 2:
            break

        if rng.random() < 0.7:
            _try_shuffle_commuting_pairs(working, num_qubits, rng, commute_cache)

        start = rng.randrange(0, len(working) - 1)
        length = rng.randint(2, min(max_block_len, len(working) - start))
        block = working[start : start + length]

        remapped, _forward, reverse = _remap_block(block)
        wire_count = len(reverse)
        if wire_count > local_qubits:
            continue

        block_unitary = circuit_unitary(local_qubits, remapped)
        key = unitary_key(block_unitary)
        node = graph.get(key)
        if node is None:
            continue

        candidate_tokens = list(node.token_chain)
        if len(candidate_tokens) >= len(block):
            continue

        candidate_remapped = pool.decode(candidate_tokens)
        candidate = _restore_block(candidate_remapped, reverse)

        if not equivalent_up_to_global_phase(
            circuit_unitary(num_qubits, block),
            circuit_unitary(num_qubits, candidate),
            atol=1e-5,
        ):
            continue

        working = working[:start] + candidate + working[start + length :]
        replacements += 1

    stats = ReductionStats(
        iterations=iterations,
        replacements=replacements,
        start_len=len(gates),
        end_len=len(working),
    )
    return working, stats
