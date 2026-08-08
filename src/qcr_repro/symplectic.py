"""Exact phase-preserving symplectic (signed-tableau) engine for Clifford pools.

For any Clifford gate set (every gate is a Clifford unitary, e.g. the ion-trap
pool {RX, RY, RZ at +-pi/2, RXX(pi/2)} or the Clifford-NISQ pool {RX(+-pi/2),
RZ(+-pi/2), CZ}), the unitary implemented by a circuit is, up to global phase,
*fully* characterized by the conjugation action on the Pauli group:

    U X_i U^dagger  and  U Z_i U^dagger   (i = 0 .. n-1).

Two Cliffords are equal up to a global phase *iff* all 2n images coincide
(if W = V^dagger U commutes with every X_i, Z_i then W = c I).

Representation of an image: a signed Pauli, stored as a list of per-qubit
factors (qubit, letter) with letter in {X, Y, Z} plus a complex coefficient in
{1, -1, i, -i}.  Conjugation by a pool gate is applied factor-by-factor using
precomputed maps (validated numerically at import time against direct matrix
conjugation), and the product is canonicalized afterwards.  This yields:

  * bit-exact, tolerance-free circuit equality checks (no 1e-5 fuzz),
  * no floating point in the hot loop (the factor maps are exact),
  * hashable keys that are invariant under global phase by construction.
"""

from __future__ import annotations

import random
from itertools import product

import numpy as np

from .config import GateInstance
from .gates import circuit_unitary, gate_matrix
from .unitary_utils import equivalent_up_to_global_phase

# --------------------------------------------------------------------------- #
# signed Pauli arithmetic
# --------------------------------------------------------------------------- #

# letter -> (k) with Pauli (i, j) multiplication:  i*j = s * k, s in {1,i,-i,-1}
# Encode letters: I=0, X=1, Y=2, Z=3  (standard index i -> sigma_i)
# sigma_a sigma_b = delta_ab I + i eps_abc sigma_c
_LETTERS = {"I": 0, "X": 1, "Y": 2, "Z": 3}
_LETTERS_R = {0: "I", 1: "X", 2: "Y", 3: "Z"}


def _mul_table() -> dict[tuple[int, int], tuple[int, complex]]:
    """Pauli letter multiplication table: (a, b) -> (c, phase) with a*b = phase*c."""
    eps = {
        (1, 2): (3, 1j), (2, 1): (3, -1j),
        (2, 3): (1, 1j), (3, 2): (1, -1j),
        (3, 1): (2, 1j), (1, 3): (2, -1j),
    }
    table: dict[tuple[int, int], tuple[int, complex]] = {}
    for a in range(4):
        for b in range(4):
            if a == 0:
                table[(a, b)] = (b, 1 + 0j)
            elif b == 0:
                table[(a, b)] = (a, 1 + 0j)
            elif a == b:
                table[(a, b)] = (0, 1 + 0j)
            else:
                table[(a, b)] = eps[(a, b)]
    return table


_MUL = _mul_table()


def _canonicalize(factors: list[tuple[int, int]], coeff: complex) -> tuple[tuple[tuple[int, int], ...], complex]:
    """Combine per-qubit Pauli factors into canonical (letter-per-qubit) form.

    ``factors`` is a list of (qubit, letter) tuples (order = product order).
    Returns ((qubit, letter) pairs for non-identity qubits, coeff).
    """
    n = max((q for q, _ in factors), default=-1) + 1
    # per-qubit current letter (0 = I) and accumulated coefficient
    cur = [0] * n
    for q, letter in factors:
        c, phase = _MUL[(cur[q], letter)]
        cur[q] = c
        coeff *= phase
    return tuple((q, _LETTERS_R[cur[q]]) for q in range(n) if cur[q] != 0), coeff


def _factors_of_label(label: str) -> list[tuple[int, int]]:
    if label == "I":
        return []
    return [(int(term[1:]), _LETTERS[term[0]]) for term in label.split("*")]


# --------------------------------------------------------------------------- #
# per-gate conjugation maps (precomputed and numerically validated)
# --------------------------------------------------------------------------- #

def _single_qubit_map(name: str, theta: float) -> dict[tuple[int, int], tuple[list[tuple[int, int]], complex]]:
    """Conjugation map of a single-qubit pool gate on qubit a.

    Returns {(0, letter): (replacement factors on qubit a, coeff)}.
    """
    from .gates import X, Y, Z, I2

    mats = {"X": X, "Y": Y, "Z": Z}
    g = gate_matrix(GateInstance(name=name, qubits=(0,), theta=theta))
    out = {}
    table = _build_signed_pauli_table(1)
    for letter, p in mats.items():
        conj = g @ p @ g.conj().T
        label, c = _round_to_signed_pauli(conj, table)
        out[(0, _LETTERS[letter])] = (_factors_of_label(label), c)
    return out


def _two_qubit_map(name: str, theta: float | None) -> dict[tuple[int, int], tuple[list[tuple[int, int]], complex]]:
    """Conjugation map of a two-qubit pool gate on qubits (0, 1)."""
    from .gates import X, Y, Z, I2

    mats = {"X": X, "Y": Y, "Z": Z}
    g = gate_matrix(GateInstance(name=name, qubits=(0, 1), theta=theta))
    out = {}
    table = _build_signed_pauli_table(2)
    for q in (0, 1):
        for letter, p in mats.items():
            pauli = np.eye(4, dtype=complex)
            if q == 0:
                pauli = np.kron(p, I2)
            else:
                pauli = np.kron(I2, p)
            conj = g @ pauli @ g.conj().T
            label, c = _round_to_signed_pauli(conj, table)
            out[(q, _LETTERS[letter])] = (_factors_of_label(label), c)
    return out


# global cache of validated maps
_MAP_CACHE: dict[tuple, dict] = {}


def _conjugation_map(name: str, theta: float | None) -> dict:
    key = (name, theta)
    cached = _MAP_CACHE.get(key)
    if cached is not None:
        return cached
    if name == "RXX":
        cached = _two_qubit_map("RXX", theta)
    elif name == "CZ":
        cached = _two_qubit_map("CZ", None)
    else:
        cached = _single_qubit_map(name, theta)
    _MAP_CACHE[key] = cached
    return cached


# --------------------------------------------------------------------------- #
# SignedTableau
# --------------------------------------------------------------------------- #

class SignedTableau:
    """Phase-preserving tableau of an n-qubit Clifford unitary (up to global phase).

    Rows 0..n-1 hold the images of X_i, rows n..2n-1 the images of Z_i.
    A row is ``(factors, coeff)`` where ``factors`` is a canonical tuple of
    (qubit, letter) pairs and ``coeff`` in {1, -1, i, -i}.
    """

    __slots__ = ("n", "rows")

    def __init__(self, n: int) -> None:
        self.n = n
        rows: list[tuple[tuple[tuple[int, str], ...], complex]] = []
        for i in range(n):
            rows.append((((i, "X"),), 1 + 0j))
        for i in range(n):
            rows.append((((i, "Z"),), 1 + 0j))
        self.rows = rows

    def _apply_map(self, a: int, gate_map: dict) -> None:
        """Apply a conjugation map acting on qubit ``a`` to every row."""
        rows = self.rows
        out = [None] * len(rows)
        for r, (factors, coeff) in enumerate(rows):
            new_factors: list[tuple[int, int]] = []
            new_coeff = coeff
            for q, letter in factors:
                if q == a:
                    repl, c = gate_map[(0, _LETTERS[letter])]
                    new_factors.extend((a, l) for (_, l) in repl)
                    new_coeff *= c
                else:
                    new_factors.append((q, _LETTERS[letter]))
            canon, coeff2 = _canonicalize(new_factors, new_coeff)
            out[r] = (canon, coeff2)
        self.rows = out

    def apply(self, gate: GateInstance) -> None:
        name = gate.name
        q = gate.qubits
        theta = gate.theta
        if name in ("RX", "RY", "RZ"):
            gate_map = _conjugation_map(name, theta)
            self._apply_map(q[0], gate_map)
        elif name in ("RXX", "CZ"):
            gate_map = _conjugation_map(name, theta)
            # apply on qubit 0 (local coordinate), with qubit 1 the other wire
            self._apply_two_qubit(q[0], q[1], gate_map)
        else:
            raise ValueError(f"SignedTableau cannot apply non-Clifford gate {name}")

    def _apply_two_qubit(self, a: int, b: int, gate_map: dict) -> None:
        """Apply a two-qubit conjugation map (acting on local qubits 0,1)."""
        rows = self.rows
        out = [None] * len(rows)
        for r, (factors, coeff) in enumerate(rows):
            new_factors: list[tuple[int, int]] = []
            new_coeff = coeff
            for qubit, letter in factors:
                if qubit == a:
                    repl, c = gate_map[(0, _LETTERS[letter])]
                    new_factors.extend((a if lq == 0 else b, l) for (lq, l) in repl)
                    new_coeff *= c
                elif qubit == b:
                    repl, c = gate_map[(1, _LETTERS[letter])]
                    new_factors.extend((a if lq == 0 else b, l) for (lq, l) in repl)
                    new_coeff *= c
                else:
                    new_factors.append((qubit, _LETTERS[letter]))
            canon, coeff2 = _canonicalize(new_factors, new_coeff)
            out[r] = (canon, coeff2)
        self.rows = out

    def apply_circuit(self, gates: list[GateInstance]) -> None:
        for gate in gates:
            self.apply(gate)

    # ------------------------------------------------------------------ #
    # keys / comparison
    # ------------------------------------------------------------------ #

    def key(self) -> tuple[tuple[tuple[tuple[int, str], ...], complex], ...]:
        """Hashable exact key (global-phase invariant)."""
        return tuple((factors, coeff) for (factors, coeff) in self.rows)

    def copy(self) -> "SignedTableau":
        t = SignedTableau(self.n)
        t.rows = [tuple(r) for r in self.rows]
        return t

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SignedTableau):
            return NotImplemented
        return self.rows == other.rows

    def __hash__(self) -> int:
        return hash(self.key())


def _factors_of_label_simple(factors: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Convert (qubit, letter-index) factors to (qubit, letter) int pairs."""
    return [(q, _LETTERS_R[l]) for q, l in factors]


# --------------------------------------------------------------------------- #
# numeric helpers (used for validation and precomputation)
# --------------------------------------------------------------------------- #

def pauli_matrix(label: str, n: int) -> np.ndarray:
    """Numeric n-qubit Pauli matrix from a label like 'Y0*X2' (qubit 0 = MSB)."""
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    z = np.array([[1, 0], [0, -1]], dtype=complex)
    mats = {"X": x, "Y": y, "Z": z}
    if label == "I":
        return np.eye(2**n, dtype=complex)
    factors: dict[int, np.ndarray] = {}
    for term in label.split("*"):
        factors[int(term[1])] = mats[term[0]]
    total = np.eye(1, dtype=complex)
    for i in range(n - 1, -1, -1):
        total = np.kron(factors.get(i, np.eye(2, dtype=complex)), total)
    return total


def _build_signed_pauli_table(n: int) -> list[tuple[np.ndarray, tuple[str, complex]]]:
    """All 4^n x 4 signed Paulis on n qubits as (matrix, (label, coeff))."""
    letters = ["I", "X", "Y", "Z"]
    table = []
    for combo in product(letters, repeat=n):
        label = "*".join(f"{l}{i}" for i, l in enumerate(combo) if l != "I") or "I"
        p = pauli_matrix(label, n)
        for c in (1, -1, 1j, -1j):
            table.append((c * p, (label, c)))
    return table


def _round_to_signed_pauli(
    matrix: np.ndarray, table: list[tuple[np.ndarray, tuple[str, complex]]]
) -> tuple[str, complex]:
    """Round a numeric Pauli image matrix to (label, coeff) using a precomputed table."""
    a = matrix.reshape(-1)
    a2 = np.vdot(a, a)
    for p, (label, c) in table:
        b = p.reshape(-1)
        scale = np.vdot(a, b) / a2
        if abs(scale - 1.0) < 1e-6 and np.linalg.norm(b - scale * a) <= 1e-6:
            return label, c
    raise RuntimeError(f"Could not round matrix to signed Pauli:\n{matrix}")


def circuits_equal_exact(c1: list[GateInstance], c2: list[GateInstance], num_qubits: int) -> bool:
    """Bit-exact check that two Clifford circuits are equal up to global phase."""
    t1 = SignedTableau(num_qubits)
    t1.apply_circuit(c1)
    t2 = SignedTableau(num_qubits)
    t2.apply_circuit(c2)
    return t1 == t2


def self_test(num_qubits: int = 3, trials: int = 60, seed: int = 0) -> bool:
    """Cross-validate the exact engine against numeric 2^n x 2^n conjugation."""
    rng = random.Random(seed)
    names = ["RX", "RY", "RZ"]
    thetas = [np.pi / 2, -np.pi / 2]
    table = _build_signed_pauli_table(num_qubits)
    for trial in range(trials):
        gates: list[GateInstance] = []
        for _ in range(rng.randint(1, 18)):
            if rng.random() < 0.35:
                i, j = rng.sample(range(num_qubits), 2)
                gates.append(GateInstance(name="RXX", qubits=(i, j), theta=np.pi / 2))
            elif rng.random() < 0.1:
                i, j = rng.sample(range(num_qubits), 2)
                gates.append(GateInstance(name="CZ", qubits=(i, j), theta=None))
            else:
                q = rng.randrange(num_qubits)
                gates.append(GateInstance(name=rng.choice(names), qubits=(q,), theta=rng.choice(thetas)))
        # exact images
        t = SignedTableau(num_qubits)
        t.apply_circuit(gates)
        exact = t.rows
        # numeric images
        u = circuit_unitary(num_qubits, gates)
        for i in range(2 * num_qubits):
            pauli = pauli_matrix(f"X{i}" if i < num_qubits else f"Z{i - num_qubits}", num_qubits)
            conj = u @ pauli @ u.conj().T
            label, c = _round_to_signed_pauli(conj, table)
            e_factors, e_c = exact[i]
            numeric_factors = tuple(sorted((q, _LETTERS_R[l]) for q, l in _factors_of_label(label)))
            if numeric_factors != e_factors or abs(c - e_c) > 1e-9:
                print(f"SELF-TEST FAILURE trial {trial} image {i}: exact {e_factors}{e_c} vs numeric {numeric_factors}{c}")
                print("gates:", [g.label() for g in gates])
                return False
    return True


if __name__ == "__main__":
    ok = self_test()
    print("symplectic self_test:", "PASS" if ok else "FAIL")
