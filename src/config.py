from __future__ import annotations

from dataclasses import dataclass
from math import pi
from typing import Literal

DEFAULT_ANGLES = (-pi / 2, -pi / 4, pi / 4, pi / 2)
PAPER_ION_ANGLES = (-pi / 2, pi / 2)
ION_TRAP_FINE_ANGLES = (-pi / 2, -pi / 4, -pi / 8, pi / 8, pi / 4, pi / 2)
GateSetName = Literal["nisq", "ion_trap", "nisq_clifford", "ion_trap_mixed", "ion_trap_mixed_clifford"]

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

NISQ_CLIFFORD_GATE_SET = GateSet(
    name="nisq_clifford",
    single_qubit=("RX", "RZ"),
    two_qubit=("CZ",),
    angles=(pi / 2, -pi / 2),  # Clifford subset: only ±π/2 rotations
)

# Ion-trap with a finer single-qubit angle grid (±π/2, ±π/4, ±π/8): mixes
# Clifford (±π/2) and non-Clifford (±π/4, ±π/8) rotations in the same pool,
# unlike ION_TRAP above which is entirely Clifford. RXX stays fixed at π/2,
# matching real ion-trap hardware's single tunable native entangler; only
# the single-qubit rotations get the finer grid.
ION_TRAP_MIXED = GateSet(
    name="ion_trap_mixed",
    single_qubit=("RX", "RY", "RZ"),
    two_qubit=("RXX",),
    angles=ION_TRAP_FINE_ANGLES,
    two_qubit_angles=(pi / 2,),
)

# The Clifford sub-pool of ION_TRAP_MIXED (±π/2 only), for the exact engine
# via the hybrid database (see hybrid.py): windows entirely within this
# sub-pool get bit-exact tableau reduction; anything touching ±π/4 or ±π/8
# falls back to the numeric engine.
ION_TRAP_MIXED_CLIFFORD = GateSet(
    name="ion_trap_mixed_clifford",
    single_qubit=("RX", "RY", "RZ"),
    two_qubit=("RXX",),
    angles=PAPER_ION_ANGLES,
    two_qubit_angles=(pi / 2,),
)

GATE_SETS: dict[GateSetName, GateSet] = {
    "ion_trap": ION_TRAP,
    "nisq": NISQ,
    "nisq_clifford": NISQ_CLIFFORD_GATE_SET,
    "ion_trap_mixed": ION_TRAP_MIXED,
    "ion_trap_mixed_clifford": ION_TRAP_MIXED_CLIFFORD,
}


def gateset_for(name: GateSetName) -> GateSet:
    return GATE_SETS[name]
