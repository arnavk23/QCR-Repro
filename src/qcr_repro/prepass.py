"""Cheap algebraic and ZX-style pre-passes applied *before* the expensive
database-driven reduction loop.

These passes implement reductions that require no search at all:

* ``algebraic_merge`` -- rotation spider fusion.  Adjacent same-axis rotations
  on the same qubit(s) compose exactly (RZ(a)RZ(b) = RZ(a+b), and the same for
  RX, RY and RXX).  The fused angle is kept only when it is the identity
  (dropped) or snaps to a pool angle, so the output stays representable in the
  database's discrete angle pool.
* ``zx_cancellations`` -- the two cheapest ZX-calculus rules on this gate set:
  adjacent same-pair CZ CZ = identity, and (NISQ only) RZ-gathering across the
  diagonal CZ gate -- RZ commutes with CZ -- so that same-axis RZ runs become
  adjacent and can then be fused/cancelled by ``algebraic_merge``.

Rosenhahn et al. state that ZX-calculus alone underperformed their method; the
point here is the opposite -- use only the *free* rules as a cheap pre-pass to
shrink the input, and spend the expensive DB/random search only on the
residual.  Every rule is applied by exact angle arithmetic, so the pre-pass is
deterministic and preserves the circuit unitary exactly (verified in
``scripts/check_batched_matches_scalar.py``).
"""

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
    """Try to fuse incoming ``g2`` with the top-of-stack ``g1``.

    Returns ``(status, replacement)`` where status is one of:

    * ``"drop"``    -- the two gates cancel to the identity (both removed)
    * ``"replace"`` -- ``g1`` is replaced by the fused gate (a pool angle)
    * ``"keep"``    -- no fusion possible; ``g2`` is appended unchanged
    """
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

    Returns ``(out, removed)`` where ``removed`` is the number of gates deleted
    (a fused pair reports 1 removed, a cancelling pair 2).  Fixpoint over the
    whole circuit is achieved by iterating :func:`apply_prepass`; a single pass
    is enough for most adjacent pairs because dropping a gate exposes the
    gates around it to the next iteration.
    """
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
    """Move every RZ gate on wire ``w`` rightward across gates it commutes with.

    RZ commutes with CZ (both diagonal) and with any gate on a different wire;
    the only obstruction on a wire is an RX on the same wire.  This is the
    paper's RZ-across-CZ structural pass applied as a pre-pass.  It is valid
    only when the two-qubit pool gate is diagonal (the NISQ CZ pools), so call
    it only for ``nisq`` / ``nisq_clifford``.
    """
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
    """Apply the algebraic / ZX-style pre-passes to fixpoint.

    ``zx=True`` enables the ZX extras (CZ-pair cancellation and, for the
    diagonal-CZ pools, RZ gathering) on top of the rotation fusion; rotation
    fusion is always applied.  Returns ``(out, removed)``.  The result is
    guaranteed to contain only pool-representable gates (fused angles snap to
    pool angles or vanish), so it can feed the database loop directly.
    """
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
