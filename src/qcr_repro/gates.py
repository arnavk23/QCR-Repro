from __future__ import annotations

from math import cos, sin

import numpy as np

from .config import GateInstance

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def rx(theta: float) -> np.ndarray:
    return cos(theta / 2) * I2 - 1j * sin(theta / 2) * X


def ry(theta: float) -> np.ndarray:
    return cos(theta / 2) * I2 - 1j * sin(theta / 2) * Y


def rz(theta: float) -> np.ndarray:
    return np.array(
        [
            [np.exp(-1j * theta / 2), 0],
            [0, np.exp(1j * theta / 2)],
        ],
        dtype=complex,
    )


def cz() -> np.ndarray:
    return np.diag([1, 1, 1, -1]).astype(complex)


def rxx(theta: float) -> np.ndarray:
    xx = np.kron(X, X)
    return cos(theta / 2) * np.eye(4, dtype=complex) - 1j * sin(theta / 2) * xx


SINGLE_QUBIT = {
    "RX": rx,
    "RY": ry,
    "RZ": rz,
}

TWO_QUBIT = {
    "CZ": lambda _theta=None: cz(),
    "RXX": rxx,
}


def gate_matrix(gate: GateInstance) -> np.ndarray:
    if gate.name in SINGLE_QUBIT:
        if gate.theta is None:
            raise ValueError(f"{gate.name} requires theta")
        return SINGLE_QUBIT[gate.name](gate.theta)
    if gate.name in TWO_QUBIT:
        return TWO_QUBIT[gate.name](gate.theta)
    raise ValueError(f"Unsupported gate: {gate.name}")


def apply_local_gate(num_qubits: int, local_matrix: np.ndarray, qubits: tuple[int, ...]) -> np.ndarray:
    if len(set(qubits)) != len(qubits):
        raise ValueError("Gate qubits must be distinct")
    if any(q < 0 or q >= num_qubits for q in qubits):
        raise ValueError("Gate qubits out of range")

    dim = 2**num_qubits
    gate_dim = 2 ** len(qubits)
    if local_matrix.shape != (gate_dim, gate_dim):
        raise ValueError("Local matrix shape does not match gate arity")

    ordered = tuple(qubits)
    others = tuple(q for q in range(num_qubits) if q not in ordered)
    basis_order = ordered + others

    full = np.zeros((dim, dim), dtype=complex)

    for col in range(dim):
        bits = [(col >> (num_qubits - 1 - i)) & 1 for i in range(num_qubits)]

        local_col = 0
        for q in ordered:
            local_col = (local_col << 1) | bits[q]

        other_bits = [bits[q] for q in others]

        for local_row in range(gate_dim):
            amp = local_matrix[local_row, local_col]
            if abs(amp) == 0:
                continue

            out_bits_map: dict[int, int] = {}
            local_row_bits = [(local_row >> (len(ordered) - 1 - i)) & 1 for i in range(len(ordered))]
            for q, bit in zip(ordered, local_row_bits):
                out_bits_map[q] = bit
            for q, bit in zip(others, other_bits):
                out_bits_map[q] = bit

            row = 0
            for i in range(num_qubits):
                row = (row << 1) | out_bits_map[i]
            full[row, col] += amp

    return full


def embedded_gate_matrix(num_qubits: int, gate: GateInstance) -> np.ndarray:
    return apply_local_gate(num_qubits, gate_matrix(gate), gate.qubits)


def circuit_unitary(num_qubits: int, gates: list[GateInstance]) -> np.ndarray:
    result = np.eye(2**num_qubits, dtype=complex)
    for gate in gates:
        result = embedded_gate_matrix(num_qubits, gate) @ result
    return result
