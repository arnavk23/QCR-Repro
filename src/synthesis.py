"""Global Clifford synthesis for the ion-trap gate pool.

Given any circuit over the Clifford pool {RX, RY, RZ, RXX} with +-pi/2 angles,
build a phase-tracking binary-symplectic tableau and reduce it back with the
Aaronson & Gottesman (2004) row-reduction algorithm (ported from Qiskit's
qiskit/synthesis/clifford/clifford_decompose_ag.py plus the phase-update rules
of qiskit/quantum_info/operators/symplectic/clifford_circuits.py, both
Apache-2.0).  The output is exact up to global phase and matches the input
tableau bit-exactly, phase bits included.

Two-qubit interactions are emitted as the native pool gate RXX (1 CNOT is 1
RXX plus single-qubit dressings), and all per-wire single-qubit runs are then
re-canonicalized to at most 3 pool rotations, which is optimal for one wire.

Reference: S. Aaronson, D. Gottesman, Phys. Rev. A 70, 052328 (2004).
"""

from __future__ import annotations

import numpy as np

from .clifford import Tableau, circuits_equal
from .config import GateInstance
from .gates import I2, gate_matrix, rz

H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)

# ---------------------------------------------------------------------------
# Phase-tracking tableau (Qiskit-style _append_* rules on integer bitmasks)
# ---------------------------------------------------------------------------


class _SignedTableau:
    """Rows 0..n-1 hold images of X_i, rows n..2n-1 of Z_i; phases[i] in {0,1}."""

    __slots__ = ("n", "xrows", "zrows", "phases")

    def __init__(self, n: int) -> None:
        self.n = n
        self.xrows = [1 << i for i in range(n)] + [0] * n
        self.zrows = [0] * n + [1 << i for i in range(n)]
        self.phases = [0] * (2 * n)

    def h(self, a: int) -> None:
        bit = 1 << a
        for i in range(2 * self.n):
            x = (self.xrows[i] >> a) & 1
            z = (self.zrows[i] >> a) & 1
            if x & z:
                self.phases[i] ^= 1
            if x != z:
                self.xrows[i] ^= bit
                self.zrows[i] ^= bit

    def s(self, a: int) -> None:
        bit = 1 << a
        for i in range(2 * self.n):
            if (self.xrows[i] >> a) & 1:
                if (self.zrows[i] >> a) & 1:
                    self.phases[i] ^= 1
                self.zrows[i] ^= bit

    def sdg(self, a: int) -> None:
        bit = 1 << a
        for i in range(2 * self.n):
            if (self.xrows[i] >> a) & 1:
                if not ((self.zrows[i] >> a) & 1):
                    self.phases[i] ^= 1
                self.zrows[i] ^= bit

    def cx(self, c: int, t: int) -> None:
        tb, cb = 1 << t, 1 << c
        xc = [(self.xrows[i] >> c) & 1 for i in range(2 * self.n)]
        zc = [(self.zrows[i] >> c) & 1 for i in range(2 * self.n)]
        xt = [(self.xrows[i] >> t) & 1 for i in range(2 * self.n)]
        zt = [(self.zrows[i] >> t) & 1 for i in range(2 * self.n)]
        for i in range(2 * self.n):
            if (xt[i] ^ zc[i] ^ 1) and zt[i] and xc[i]:
                self.phases[i] ^= 1
            if xc[i]:
                self.xrows[i] ^= tb
            if zt[i]:
                self.zrows[i] ^= cb

    def swap(self, a: int, b: int) -> None:
        xa, xb = 1 << a, 1 << b
        for i in range(2 * self.n):
            x = self.xrows[i]
            z = self.zrows[i]
            xt = ((x >> a) & 1) ^ ((x >> b) & 1)
            zt = ((z >> a) & 1) ^ ((z >> b) & 1)
            self.xrows[i] = x ^ (xt << a) ^ (xt << b)
            self.zrows[i] = z ^ (zt << a) ^ (zt << b)

    def x(self, a: int) -> None:
        for i in range(2 * self.n):
            if (self.zrows[i] >> a) & 1:
                self.phases[i] ^= 1

    def z(self, a: int) -> None:
        for i in range(2 * self.n):
            if (self.xrows[i] >> a) & 1:
                self.phases[i] ^= 1

    # -- pool gates as composites of the primitives (mirror clifford.Tableau)

    def rz(self, a: int, theta: float) -> None:
        if theta > 0:
            self.s(a)
        else:
            self.sdg(a)

    def rx(self, a: int, theta: float) -> None:
        self.h(a)
        self.rz(a, theta)
        self.h(a)

    def ry(self, a: int, theta: float) -> None:
        # RY(theta) = S RX(theta) S3 (matrix product), so the circuit applies
        # [S3, RX, S] in order and the tableau is built in that same order.
        self.sdg(a)
        self.rx(a, theta)
        self.s(a)

    def rxx(self, a: int, b: int) -> None:
        self.h(a)
        self.h(b)
        self.cx(a, b)
        self.s(b)
        self.cx(a, b)
        self.h(a)
        self.h(b)

    def cz(self, a: int, b: int) -> None:
        self.h(b)
        self.cx(a, b)
        self.h(b)

    def apply(self, gate: GateInstance) -> None:
        name = gate.name
        q = gate.qubits
        if name == "RX":
            self.rx(q[0], gate.theta)
        elif name == "RY":
            self.ry(q[0], gate.theta)
        elif name == "RZ":
            self.rz(q[0], gate.theta)
        elif name == "RXX":
            self.rxx(q[0], q[1])
        elif name == "CZ":
            self.cz(q[0], q[1])
        else:
            raise ValueError(f"Clifford engine does not support {name}")

    def apply_circuit(self, gates: list[GateInstance]) -> None:
        for gate in gates:
            self.apply(gate)

    def key(self) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
        return (tuple(self.xrows), tuple(self.zrows), tuple(self.phases))


# ---------------------------------------------------------------------------
# 1-qubit Clifford -> pool-gate canonical table (up to global phase)
# ---------------------------------------------------------------------------

_POOL_1Q: tuple[tuple[str, float], ...] = (
    ("RX", np.pi / 2), ("RX", -np.pi / 2),
    ("RY", np.pi / 2), ("RY", -np.pi / 2),
    ("RZ", np.pi / 2), ("RZ", -np.pi / 2),
)


def _single_key(m: np.ndarray) -> bytes:
    """Canonical byte key of a 2x2 unitary up to global phase."""
    idx = np.unravel_index(np.argmax(np.abs(m)), m.shape)
    m = m * (m[idx].conjugate() / abs(m[idx]))
    m = np.round(m, 8)
    m.real[np.abs(m.real) < 1e-9] = 0.0
    m.imag[np.abs(m.imag) < 1e-9] = 0.0
    return m.tobytes()


def _build_single_qubit_table() -> dict[bytes, list[GateInstance]]:
    """Minimal pool decompositions of every single-qubit Clifford."""
    from itertools import product

    table: dict[bytes, list[GateInstance]] = {_single_key(I2): []}
    for length in (1, 2, 3, 4):
        for seq in product(_POOL_1Q, repeat=length):
            m = I2
            for name, theta in seq:
                m = gate_matrix(GateInstance(name, (0,), theta)) @ m
            key = _single_key(m)
            if key not in table:
                table[key] = [GateInstance(name, (0,), theta) for name, theta in seq]

    # Verify coverage: the 24 single-qubit Cliffords generated by {S, H}.
    s = rz(np.pi / 2)
    frontier = [np.eye(2, dtype=complex)]
    elements = [np.eye(2, dtype=complex)]
    for _ in range(6):
        nxt = []
        for m in frontier:
            for gm in (s, H):
                p = gm @ m
                pk = _single_key(p)
                if all(pk != _single_key(e) for e in elements):
                    elements.append(p)
                    nxt.append(p)
        frontier = nxt
    missing = [e for e in elements if _single_key(e) not in table]
    if missing:
        raise RuntimeError(f"single-qubit Clifford table is incomplete ({len(missing)} missing)")
    return table


_SINGLE_TABLE = _build_single_qubit_table()


def _single_decompose(m: np.ndarray) -> list[GateInstance] | None:
    return _SINGLE_TABLE.get(_single_key(m))


# ---------------------------------------------------------------------------
# Aaronson-Gottesman row reduction
# ---------------------------------------------------------------------------


def _row_reduce(tab: _SignedTableau) -> list[tuple]:
    """Reduce tableau to identity, recording the applied primitives.

Items are ('H', a), ('S', a), ('S3', a), ('CX', c, t), ('SWAP', a, b),
('X', a), ('Z', a); the inverse of this sequence is the synthesized circuit."""
    n = tab.n
    rec: list[tuple] = []

    def h(a: int) -> None:
        tab.h(a)
        rec.append(("H", a))

    def s(a: int) -> None:
        tab.s(a)
        rec.append(("S", a))

    def s3(a: int) -> None:
        tab.sdg(a)
        rec.append(("S3", a))

    def cx(c: int, t: int) -> None:
        tab.cx(c, t)
        rec.append(("CX", c, t))

    def sw(a: int, b: int) -> None:
        tab.swap(a, b)
        rec.append(("SWAP", a, b))

    for i in range(n):
        # -- _set_qubit_x_true: put a 1 in position i of the X row
        if not ((tab.xrows[i] >> i) & 1):
            found = False
            for j in range(i + 1, n):
                if (tab.xrows[i] >> j) & 1:
                    sw(i, j)
                    found = True
                    break
            if not found:
                for j in range(i, n):
                    if (tab.zrows[i] >> j) & 1:
                        h(j)
                        if j != i:
                            sw(i, j)
                        break

        # -- _set_row_x_zero: clear X entries j > i, then Z entries
        for j in range(i + 1, n):
            if (tab.xrows[i] >> j) & 1:
                cx(i, j)
        if any((tab.zrows[i] >> k) & 1 for k in range(i, n)):
            if not ((tab.zrows[i] >> i) & 1):
                s(i)
            for j in range(i + 1, n):
                if (tab.zrows[i] >> j) & 1:
                    cx(j, i)
            s(i)

        # -- _set_row_z_zero: clear stabilizer Z entries j > i, then X entries
        if any((tab.zrows[n + i] >> j) & 1 for j in range(i + 1, n)):
            for j in range(i + 1, n):
                if (tab.zrows[n + i] >> j) & 1:
                    cx(j, i)
        if any((tab.xrows[n + i] >> j) & 1 for j in range(i, n)):
            h(i)
            for j in range(i + 1, n):
                if (tab.xrows[n + i] >> j) & 1:
                    cx(i, j)
            if (tab.zrows[n + i] >> i) & 1:
                s(i)
            h(i)

    for i in range(n):
        if tab.phases[i]:
            tab.z(i)
            rec.append(("Z", i))
        if tab.phases[n + i]:
            tab.x(i)
            rec.append(("X", i))

    return rec


# ---------------------------------------------------------------------------
# Primitive translation to the ion pool
# ---------------------------------------------------------------------------


def _cnot_pool(c: int, t: int) -> list[GateInstance]:
    """CNOT_{c->t} over {RZ, H(virtual), RXX}: 1 RXX + single-qubit dressings."""
    return [
        GateInstance("H", (t,)),
        GateInstance("RZ", (t,), -np.pi / 2),
        GateInstance("RZ", (c,), -np.pi / 2),
        GateInstance("H", (c,)),
        GateInstance("H", (t,)),
        GateInstance("RXX", (c, t), np.pi / 2),
        GateInstance("H", (t,)),
        GateInstance("H", (c,)),
        GateInstance("H", (t,)),
    ]


def _translate(primitives: list[tuple]) -> list[GateInstance]:
    """Expand primitives into pool gates (with virtual 'H'/'X'/'Z' singles)."""
    out: list[GateInstance] = []
    for item in primitives:
        op = item[0]
        if op in ("H", "X", "Z"):
            out.append(GateInstance(op, (item[1],)))
        elif op == "S":
            out.append(GateInstance("RZ", (item[1],), np.pi / 2))
        elif op == "S3":
            out.append(GateInstance("RZ", (item[1],), -np.pi / 2))
        elif op == "CX":
            out.extend(_cnot_pool(item[1], item[2]))
        elif op == "SWAP":
            a, b = item[1], item[2]
            out.extend(_cnot_pool(a, b) + _cnot_pool(b, a) + _cnot_pool(a, b))
        else:
            raise ValueError(f"unknown primitive {item}")
    return out


def _compress_single_qubits(gates: list[GateInstance], num_qubits: int) -> list[GateInstance]:
    """Re-canonicalize maximal per-wire single-qubit runs to <=3 pool gates."""
    runs: list[list[GateInstance]] = [[] for _ in range(num_qubits)]
    out: list[GateInstance] = []
    mats = {"H": H, "X": _X, "Z": _Z}

    def flush(w: int) -> None:
        if not runs[w]:
            return
        m = I2
        for gt in runs[w]:
            if gt.name in mats:
                m = mats[gt.name] @ m
            else:
                m = gate_matrix(gt) @ m
        seq = _single_decompose(m)
        if seq is None:
            out.extend(runs[w])
        else:
            out.extend(GateInstance(g.name, (w,), g.theta) for g in seq)
        runs[w] = []

    for gate in gates:
        if gate.name == "RXX":
            flush(gate.qubits[0])
            flush(gate.qubits[1])
            out.append(gate)
        else:
            runs[gate.qubits[0]].append(gate)
    for w in range(num_qubits):
        flush(w)
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def synth_clifford(tableau: Tableau, signed: _SignedTableau | None = None) -> list[GateInstance]:
    """Synthesize an n-qubit Clifford into ion-pool gates (exact up to global phase).

``signed`` optionally supplies the phase-carrying tableau (from the input
circuit); if omitted, phases are assumed zero."""
    n = tableau.n
    tab = _SignedTableau(n)
    tab.xrows = list(tableau.xrows)
    tab.zrows = list(tableau.zrows)
    if signed is not None:
        tab.phases = list(signed.phases)

    rec = _row_reduce(tab)
    primitives = []
    for item in reversed(rec):
        op = item[0]
        if op == "H":
            primitives.append(("H", item[1]))
        elif op == "S":
            primitives.append(("S3", item[1]))
        elif op == "S3":
            primitives.append(("S", item[1]))
        elif op == "CX":
            primitives.append(("CX", item[1], item[2]))
        elif op == "SWAP":
            primitives.append(("SWAP", item[1], item[2]))
        elif op == "X":
            primitives.append(("X", item[1]))
        elif op == "Z":
            primitives.append(("Z", item[1]))

    return _compress_single_qubits(_translate(primitives), n)


def synth_circuit(gates: list[GateInstance], num_qubits: int, verify: bool = True) -> list[GateInstance]:
    """Synthesize a full Clifford-pool circuit; optionally assert equivalence."""
    tab = Tableau(num_qubits)
    tab.apply_circuit(gates)
    signed = _SignedTableau(num_qubits)
    signed.apply_circuit(gates)
    out = synth_clifford(tab, signed)
    if verify:
        if not circuits_equal(gates, out, num_qubits):
            raise AssertionError("synthesis output is not equivalent to the input circuit")
    return out


def self_test(num_qubits: int = 5, trials: int = 60, seed: int = 0) -> bool:
    """Cross-validate synthesis against random Clifford-pool circuits."""
    import random

    from .gates import circuit_unitary
    from .unitary import equivalent_up_to_global_phase

    rng = random.Random(seed)
    names = ["RX", "RY", "RZ"]
    thetas = [np.pi / 2, -np.pi / 2]
    for _ in range(trials):
        gates: list[GateInstance] = []
        for _ in range(rng.randint(1, 18)):
            if num_qubits > 1 and rng.random() < 0.4:
                i, j = rng.sample(range(num_qubits), 2)
                gates.append(GateInstance(name="RXX", qubits=(i, j), theta=np.pi / 2))
            else:
                q = rng.randrange(num_qubits)
                gates.append(GateInstance(name=rng.choice(names), qubits=(q,), theta=rng.choice(thetas)))
        out = synth_circuit(gates, num_qubits, verify=True)
        u = circuit_unitary(num_qubits, gates)
        v = circuit_unitary(num_qubits, out)
        if not equivalent_up_to_global_phase(u, v):
            print(f"SELF-TEST FAILURE: numeric unitary mismatch on trial {_}")
            print("input:", [g.label() for g in gates])
            print("synth:", [g.label() for g in out])
            return False
    return True


if __name__ == "__main__":
    ok = self_test()
    print("synthesis self_test:", "PASS" if ok else "FAIL")
