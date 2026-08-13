"""Cheap algebraic and ZX-style pre-passes run before the database loop.

Rotation-spider fusion, CZ-pair cancellation, and RZ-gathering across diagonal CZ -- all free rules applied by exact angle arithmetic, shrinking the input before expensive search."""

from __future__ import annotations

import math

from .config import ANGLE_EPS, GateInstance

# Parameterized gates whose angles compose algebraically.
_MERGEABLE = ("RX", "RY", "RZ", "RXX")


def _canonical(theta: float) -> float:
    """Map an angle into (-pi, pi]."""
    two_pi = 2.0 * math.pi
    t = math.fmod(theta, two_pi)
    if t > math.pi:
        t -= two_pi
    elif t <= -math.pi:
        t += two_pi
    return t


def _fuse(
    g1: GateInstance,
    g2: GateInstance,
    angles: tuple[float, ...],
    two_angles: tuple[float, ...],
    atol: float,
) -> tuple[str, GateInstance | None]:
    """Try to fuse incoming g2 with top-of-stack g1.

Returns (status, replacement) with status in {"drop", "replace", "keep"}."""
    if g1.name != g2.name or tuple(sorted(g1.qubits)) != tuple(sorted(g2.qubits)):
        return "keep", g2
    total = _canonical(g1.theta + g2.theta)
    if abs(total) <= atol:
        return "drop", None
    pool = angles if len(g1.qubits) == 1 else two_angles
    for p in pool:
        if abs(_canonical(total - p)) <= atol:
            return "replace", GateInstance(name=g1.name, qubits=g1.qubits, theta=p)
    return "keep", g2


def algebraic_merge(
    gates: list[GateInstance],
    angles: tuple[float, ...],
    two_angles: tuple[float, ...],
    atol: float = 1e-6,
) -> tuple[list[GateInstance], int]:
    """One left-to-right pass of adjacent same-axis rotation fusion.

Returns (out, removed); iterate apply_prepass to reach the fixpoint."""
    out: list[GateInstance] = []
    removed = 0
    for gate in gates:
        if (
            out
            and gate.name in _MERGEABLE
            and out[-1].name == gate.name
            and out[-1].theta is not None
            and gate.theta is not None
        ):
            status, replacement = _fuse(out[-1], gate, angles, two_angles, atol)
            if status == "drop":
                out.pop()
                removed += 2
                continue
            if status == "replace":
                out[-1] = replacement
                removed += 1
                continue
        out.append(gate)
    return out, removed


def zx_cancellations(gates: list[GateInstance], atol: float = 1e-9) -> tuple[list[GateInstance], int]:
    """Cheap ZX rule: two adjacent CZ gates on the same pair are the identity."""
    out: list[GateInstance] = []
    removed = 0
    for gate in gates:
        if (
            gate.name == "CZ"
            and out
            and out[-1].name == "CZ"
            and tuple(sorted(out[-1].qubits)) == tuple(sorted(gate.qubits))
        ):
            out.pop()
            removed += 2
            continue
        out.append(gate)
    return out, removed


def gather_rz_across_cz(gates: list[GateInstance], num_qubits: int) -> list[GateInstance]:
    """Move RZ gates rightward across gates they commute with (CZ, other-wire gates).

Valid only for diagonal-CZ pools (nisq / nisq_clifford)."""
    result = list(gates)
    for w in range(num_qubits):
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
    return result


def apply_prepass(
    gates: list[GateInstance],
    gate_set_name: str,
    angles: tuple[float, ...],
    two_angles: tuple[float, ...],
    num_qubits: int,
    zx: bool = True,
    atol: float = 1e-6,
    max_iters: int = 8,
) -> tuple[list[GateInstance], int]:
    """Apply the pre-passes to fixpoint; zx=True adds CZ cancellation and RZ gathering.

Returns (out, removed); output contains only pool-representable gates."""
    working = list(gates)
    total = 0
    for _ in range(max_iters):
        removed = 0
        if zx:
            if gate_set_name != "ion_trap":
                working = gather_rz_across_cz(working, num_qubits)
            working, r = zx_cancellations(working)
            removed += r
        working, r = algebraic_merge(working, angles, two_angles, atol)
        removed += r
        total += removed
        if removed == 0:
            break
    return working, total
