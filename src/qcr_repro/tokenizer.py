from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Optional

from .config import DEFAULT_ANGLES, GateInstance, GateSetName


@dataclass
class TokenPool:
    num_qubits: int
    gate_set: GateSetName
    angles: tuple[float, ...] = DEFAULT_ANGLES
    rxx_angles: Optional[tuple[float, ...]] = None

    def __post_init__(self) -> None:
        self._token_to_gate: dict[int, GateInstance] = {}
        self._gate_to_token: dict[GateInstance, int] = {}
        self._build()

    def _add(self, gate: GateInstance) -> None:
        token = len(self._token_to_gate) + 1
        self._token_to_gate[token] = gate
        self._gate_to_token[gate] = token

    def _build(self) -> None:
        single_names = ["RX", "RZ"]
        two_qubit_names = ["CZ"]
        if self.gate_set == "ion_trap":
            single_names = ["RX", "RY", "RZ"]
            two_qubit_names = ["RXX"]

        for name in single_names:
            for qubit in range(self.num_qubits):
                for theta in self.angles:
                    self._add(GateInstance(name=name, qubits=(qubit,), theta=theta))

        for name in two_qubit_names:
            for q0, q1 in combinations(range(self.num_qubits), 2):
                if name == "CZ":
                    self._add(GateInstance(name=name, qubits=(q0, q1), theta=None))
                else:
                    thetas = self.rxx_angles if self.rxx_angles is not None else self.angles
                    for theta in thetas:
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
