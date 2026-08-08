from __future__ import annotations

import random
from dataclasses import dataclass

from .config import GateInstance, GateSetName, gateset_for
from .tokenizer import TokenPool


def build_pool(
    num_qubits: int,
    gate_set: GateSetName,
    angles: tuple[float, ...] | None = None,
    two_qubit_angles: tuple[float, ...] | None = None,
) -> TokenPool:
    gs = gateset_for(gate_set)
    if angles is None:
        angles = gs.angles
    if two_qubit_angles is None:
        two_qubit_angles = gs.two_angles
    return TokenPool(
        num_qubits=num_qubits,
        gate_set=gate_set,
        angles=angles,
        two_qubit_angles=two_qubit_angles,
    )


def random_circuit(
    num_qubits: int,
    length: int,
    gate_set: GateSetName,
    angles: tuple[float, ...] | None = None,
    two_qubit_angles: tuple[float, ...] | None = None,
    seed: int = 0,
    weights: dict[str, float] | None = None,
) -> tuple[list[GateInstance], TokenPool]:
    """Random circuit drawn i.i.d. from the gate pool.

    ``weights`` optionally overrides the sampling weight per gate name (default
    1.0 for every pool gate).  The NISQ inputs in Rosenhahn et al. (Table 7)
    exhibit RX:RZ:CZ ~ 108:109:82, which is exactly what a pool with CZ
    weighted twice as strongly as the single-qubit gates produces; the ion-trap
    inputs (Table 6) match uniform weights.
    """
    pool = build_pool(num_qubits, gate_set, angles, two_qubit_angles)
    rng = random.Random(seed)
    if weights is None:
        tokens = [rng.choice(pool.tokens()) for _ in range(length)]
    else:
        pool_gates = pool.gates()
        wlist = [weights.get(g.name, 1.0) for g in pool_gates]
        indices = rng.choices(range(len(pool_gates)), weights=wlist, k=length)
        tokens = [idx + 1 for idx in indices]
    return pool.decode(tokens), pool


def count_gates(gates: list[GateInstance]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for gate in gates:
        counts[gate.name] = counts.get(gate.name, 0) + 1
    return counts


def gate_totals(gates: list[GateInstance]) -> int:
    return len(gates)
