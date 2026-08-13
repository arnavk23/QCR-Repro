"""Exact Clifford (tableau) engine for the ion-trap gate pool.

The pool is all-Clifford, so binary symplectic tableaux capture circuits exactly up to global phase: bit-exact equivalence and O(n^2) integer simulation."""

from __future__ import annotations

import random

import numpy as np

from .config import GateInstance
from .gates import circuit_unitary, gate_matrix


class Tableau:
    """Binary symplectic tableau of an n-qubit Clifford circuit.

Rows 0..n-1 hold images of X_i, rows n..2n-1 of Z_i, as (xmask, zmask) int pairs; phases are dropped (global-phase equivalence)."""

    __slots__ = ("n", "xrows", "zrows")

    def __init__(self, n: int) -> None:
        self.n = n
        # identity tableau
        self.xrows = [1 << i for i in range(n)] + [0] * n
        self.zrows = [0] * n + [1 << i for i in range(n)]

    # -- elementary Clifford operations (row transforms on the symplectic matrix)

    def h(self, a: int) -> None:
        """Hadamard on qubit a: swap X and Z columns."""
        bit = 1 << a
        for i in range(2 * self.n):
            x = self.xrows[i]
            z = self.zrows[i]
            if ((x >> a) & 1) != ((z >> a) & 1):
                self.xrows[i] = x ^ bit
                self.zrows[i] = z ^ bit

    def s(self, a: int) -> None:
        """Phase gate S = RZ(pi/2) on qubit a: X_a -> Y_a, Z_a -> Z_a."""
        bit = 1 << a
        for i in range(2 * self.n):
            if (self.xrows[i] >> a) & 1:
                self.zrows[i] ^= bit

    def s3(self, a: int) -> None:
        """S^3 = RZ(-pi/2): X_a -> -Y_a, Z_a -> Z_a."""
        self.s(a)
        self.s(a)
        self.s(a)

    def cnot(self, c: int, t: int) -> None:
        """CNOT with control c, target t (conjugation action)."""
        xb = 1 << t
        zb = 1 << c
        for i in range(2 * self.n):
            if (self.xrows[i] >> c) & 1:
                self.xrows[i] ^= xb
            if (self.zrows[i] >> t) & 1:
                self.zrows[i] ^= zb

    # -- pool gates (compositions of the primitives, up to global phase)

    def rz(self, a: int, theta: float) -> None:
        if theta > 0:
            self.s(a)
        else:
            self.s3(a)

    def rx(self, a: int, theta: float) -> None:
        # RX(theta) = H RZ(theta) H
        self.h(a)
        self.rz(a, theta)
        self.h(a)

    def ry(self, a: int, theta: float) -> None:
        # RY(theta) = S RX(theta) S^3  (conjugation by S maps X -> Y)
        self.s(a)
        self.rx(a, theta)
        self.s3(a)

    def rxx(self, a: int, b: int) -> None:
        # RXX(theta) = (H_a H_b) CNOT_{a->b} RZ_b(theta) CNOT_{a->b} (H_a H_b)
        # (up to global phase), derived from the standard identities
        #   RXX(theta) = (H otimes H) RZZ(theta) (H otimes H)
        #   RZZ(theta) = CNOT (I otimes RZ_target(theta)) CNOT
        self.h(a)
        self.h(b)
        self.cnot(a, b)
        self.rz(b, np.pi / 2)
        self.cnot(a, b)
        self.h(a)
        self.h(b)

    # -- whole-circuit application

    def apply(self, gate: GateInstance) -> None:
        """Apply one pool gate (RX/RY/RZ at +/-pi/2 or RXX at pi/2)."""
        name = gate.name
        q = gate.qubits
        theta = gate.theta
        if name == "RX":
            self.rx(q[0], theta)
        elif name == "RY":
            self.ry(q[0], theta)
        elif name == "RZ":
            self.rz(q[0], theta)
        elif name == "RXX":
            self.rxx(q[0], q[1])
        elif name == "CZ":
            # CZ = H_b CNOT_{a->b} H_b
            self.h(q[1])
            self.cnot(q[0], q[1])
            self.h(q[1])
        else:
            raise ValueError(f"Clifford engine does not support {name}")

    def apply_circuit(self, gates: list[GateInstance]) -> None:
        for gate in gates:
            self.apply(gate)

    # -- comparisons / hashing

    def key(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        return (tuple(self.xrows), tuple(self.zrows))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Tableau):
            return NotImplemented
        return self.xrows == other.xrows and self.zrows == other.zrows

    def __hash__(self) -> int:
        return hash((tuple(self.xrows), tuple(self.zrows)))

    # -- numeric cross-validation helpers

    def pauli_images(self) -> dict[str, str]:
        """Return {Pauli_label: image_label} for X_i, Z_i, Y_i (sign-free)."""
        n = self.n
        out: dict[str, str] = {}

        def label(xm: int, zm: int) -> str:
            parts = []
            for i in range(n):
                x = (xm >> i) & 1
                z = (zm >> i) & 1
                if x and z:
                    parts.append(f"Y{i}")
                elif x:
                    parts.append(f"X{i}")
                elif z:
                    parts.append(f"Z{i}")
            return "*".join(parts) if parts else "I"

        for i in range(n):
            out[f"X{i}"] = label(self.xrows[i], self.zrows[i])
            out[f"Z{i}"] = label(self.xrows[n + i], self.zrows[n + i])
            # Y_i = i X_i Z_i  ->  image is the product of the two images
            xm = self.xrows[i] ^ self.xrows[n + i]
            zm = self.zrows[i] ^ self.zrows[n + i]
            out[f"Y{i}"] = label(xm, zm)
        return out


def circuits_equal(c1: list[GateInstance], c2: list[GateInstance], num_qubits: int) -> bool:
    """Bit-exact check that two Clifford circuits are equal up to global phase."""
    t1 = Tableau(num_qubits)
    t1.apply_circuit(c1)
    t2 = Tableau(num_qubits)
    t2.apply_circuit(c2)
    return t1 == t2


def pauli_matrix(label: str, n: int) -> np.ndarray:
    """Numeric n-qubit Pauli matrix from a sign-free label like 'Y0*X2'."""
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.array([[1, 0], [0, -1]], dtype=complex)
    mats = {"X": x, "Y": y, "Z": z}
    total = np.eye(1, dtype=complex)
    factors: dict[int, np.ndarray] = {}
    if label == "I":
        return np.eye(2**n, dtype=complex)
    for term in label.split("*"):
        factors[int(term[1])] = mats[term[0]]
    # Qubit 0 is the most significant bit (matches circuit_unitary's convention).
    for i in range(n - 1, -1, -1):
        total = np.kron(factors.get(i, np.eye(2, dtype=complex)), total)
    return total


def _is_scalar_multiple(u: np.ndarray, v: np.ndarray, atol: float = 1e-8) -> bool:
    a = u.reshape(-1)
    b = v.reshape(-1)
    scale = np.vdot(a, b) / np.vdot(a, a)
    return np.linalg.norm(b - scale * a) <= atol


def self_test(num_qubits: int = 3, trials: int = 40, seed: int = 0) -> bool:
    """Cross-validate the tableau engine against numeric unitaries."""
    rng = random.Random(seed)
    names = ["RX", "RY", "RZ"]
    thetas = [np.pi / 2, -np.pi / 2]
    for _ in range(trials):
        gates: list[GateInstance] = []
        for _ in range(rng.randint(1, 14)):
            if rng.random() < 0.35:
                i, j = rng.sample(range(num_qubits), 2)
                gates.append(GateInstance(name="RXX", qubits=(i, j), theta=np.pi / 2))
            else:
                q = rng.randrange(num_qubits)
                gates.append(GateInstance(name=rng.choice(names), qubits=(q,), theta=rng.choice(thetas)))
        tab = Tableau(num_qubits)
        tab.apply_circuit(gates)
        u = circuit_unitary(num_qubits, gates)
        images = tab.pauli_images()
        for pauli_label, image_label in images.items():
            p = pauli_matrix(pauli_label, num_qubits)
            conj = u @ p @ u.conj().T
            target = pauli_matrix(image_label, num_qubits)
            if not _is_scalar_multiple(conj, target):
                print(f"SELF-TEST FAILURE: {pauli_label} -> {image_label} but numeric says otherwise")
                print("gates:", [g.label() for g in gates])
                return False
    return True


if __name__ == "__main__":
    ok = self_test()
    print("clifford self_test:", "PASS" if ok else "FAIL")
