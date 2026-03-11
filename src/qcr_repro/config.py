from __future__ import annotations

from dataclasses import dataclass
from math import pi
from typing import Literal

DEFAULT_ANGLES = (-pi / 2, -pi / 4, pi / 4, pi / 2)
GateSetName = Literal["nisq", "ion_trap"]


@dataclass(frozen=True)
class GateInstance:
    name: str
    qubits: tuple[int, ...]
    theta: float | None = None

    @property
    def arity(self) -> int:
        return len(self.qubits)

    def label(self) -> str:
        theta_part = "" if self.theta is None else f", {self.theta:.6f}"
        qubit_part = ", ".join(str(q) for q in self.qubits)
        return f"{self.name}({qubit_part}{theta_part})"
