from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Optional

from .config import ANGLE_EPS, DEFAULT_ANGLES, GateInstance, GateSet, GateSetName, gateset_for


@dataclass
class TokenPool:
    num_qubits: int
    gate_set: GateSetName | GateSet
    angles: Optional[tuple[float, ...]] = None
    two_qubit_angles: Optional[tuple[float, ...]] = None

    def __post_init__(self) -> None:
        if isinstance(self.gate_set, str):
            self.gate_set = gateset_for(self.gate_set)
        if self.angles is None:
            self.angles = self.gate_set.angles
        if self.two_qubit_angles is None:
            self.two_qubit_angles = self.gate_set.two_angles
        self._token_to_gate: dict[int, GateInstance] = {}
        self._gate_to_token: dict[GateInstance, int] = {}
        self._build()

    def _add(self, gate: GateInstance) -> None:
        token = len(self._token_to_gate) + 1
        self._token_to_gate[token] = gate
        self._gate_to_token[gate] = token

    def _build(self) -> None:
        gs = self.gate_set
        angles = self.angles
        for name in gs.single_qubit:
            for qubit in range(self.num_qubits):
                for theta in angles:
                    self._add(GateInstance(name=name, qubits=(qubit,), theta=theta))

        two_angles = self.two_qubit_angles if self.two_qubit_angles is not None else gs.two_angles
        for name in gs.two_qubit:
            if name == "CZ":
                for q0, q1 in combinations(range(self.num_qubits), 2):
                    self._add(GateInstance(name=name, qubits=(q0, q1), theta=None))
            else:
                for q0, q1 in combinations(range(self.num_qubits), 2):
                    for theta in two_angles:
                        self._add(GateInstance(name=name, qubits=(q0, q1), theta=theta))

    def gate_for_token(self, token: int) -> GateInstance:
        return self._token_to_gate[token]

    def token_for_gate(self, gate: GateInstance) -> int:
        return self._gate_to_token[gate]

    def tokens(self) -> list[int]:
        return list(self._token_to_gate.keys())

    def gates(self) -> list[GateInstance]:
        return [self._token_to_gate[token] for token in self.tokens()]

    def decode(self, token_chain: list[int]) -> list[GateInstance]:
        return [self.gate_for_token(token) for token in token_chain]

    def encode(self, gates: list[GateInstance]) -> list[int]:
        return [self.token_for_gate(gate) for gate in gates]

    def snap(self, gate: GateInstance) -> GateInstance:
        """Return a pool gate equivalent to ``gate`` if one exists within ANGLE_EPS.

        Used to restore exactness for QASM inputs whose angles were rounded to a
        few decimals but were generated from a discrete pool.
        """
        gs = self.gate_set
        if gate.name in gs.single_qubit and gate.theta is not None:
            best = min(gs.angles, key=lambda a: abs(a - gate.theta))
            if abs(best - gate.theta) <= ANGLE_EPS:
                return GateInstance(name=gate.name, qubits=tuple(sorted(gate.qubits)), theta=best)
        if gate.name in gs.two_qubit:
            qubits = tuple(sorted(gate.qubits))
            if gate.name == "CZ":
                return GateInstance(name=gate.name, qubits=qubits, theta=None)
            if gate.theta is not None:
                best = min(self.two_qubit_angles if self.two_qubit_angles is not None else gs.two_angles,
                           key=lambda a: abs(a - gate.theta))
                if abs(best - gate.theta) <= ANGLE_EPS:
                    return GateInstance(name=gate.name, qubits=qubits, theta=best)
        return gate
