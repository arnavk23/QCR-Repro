from __future__ import annotations

from dataclasses import dataclass
from math import pi
from typing import Literal

DEFAULT_ANGLES = (-pi / 2, -pi / 4, pi / 4, pi / 2)
PAPER_ION_ANGLES = (-pi / 2, pi / 2)
GateSetName = Literal["nisq", "ion_trap"]

ANGLE_EPS = 1e-3


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


@dataclass(frozen=True)
class GateSet:
    """An available-gate pool definition (angle-discretized)."""

    name: GateSetName
    single_qubit: tuple[str, ...]
    two_qubit: tuple[str, ...]
    angles: tuple[float, ...]
    two_qubit_angles: tuple[float, ...] | None = None

    @property
    def two_angles(self) -> tuple[float, ...]:
        """Angles used for parameterized two-qubit gates (RXX)."""
        return self.two_qubit_angles if self.two_qubit_angles is not None else self.angles


ION_TRAP = GateSet(
    name="ion_trap",
    single_qubit=("RX", "RY", "RZ"),
    two_qubit=("RXX",),
    angles=PAPER_ION_ANGLES,
    two_qubit_angles=(pi / 2,),
)

NISQ = GateSet(
    name="nisq",
    single_qubit=("RX", "RZ"),
    two_qubit=("CZ",),
    angles=DEFAULT_ANGLES,
)

GATE_SETS: dict[GateSetName, GateSet] = {
    "ion_trap": ION_TRAP,
    "nisq": NISQ,
}


def gateset_for(name: GateSetName) -> GateSet:
    return GATE_SETS[name]
