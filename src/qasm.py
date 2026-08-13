from __future__ import annotations

import re
from pathlib import Path

from .config import GateInstance

_SINGLE_RE = re.compile(r"^(rx|ry|rz)\(([-+0-9.eE]+)\)\s+q\[(\d+)\];$")
_RXX_RE = re.compile(r"^rxx\(([-+0-9.eE]+)\)\s+q\[(\d+)\],\s*q\[(\d+)\];$")
_CZ_RE = re.compile(r"^cz\s+q\[(\d+)\],\s*q\[(\d+)\];$")
_QREG_RE = re.compile(r"^qreg\s+q\[(\d+)\];$")


def parse_qasm_subset(path: str | Path) -> tuple[int, list[GateInstance]]:
    num_qubits = None
    gates: list[GateInstance] = []

    with open(path, "r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if not line or line.startswith("//"):
                continue

            qreg_match = _QREG_RE.match(line)
            if qreg_match:
                num_qubits = int(qreg_match.group(1))
                continue

            single_match = _SINGLE_RE.match(line)
            if single_match:
                name = single_match.group(1).upper()
                theta = float(single_match.group(2))
                qubit = int(single_match.group(3))
                gates.append(GateInstance(name=name, qubits=(qubit,), theta=theta))
                continue

            cz_match = _CZ_RE.match(line)
            if cz_match:
                q0 = int(cz_match.group(1))
                q1 = int(cz_match.group(2))
                if q0 > q1:
                    q0, q1 = q1, q0
                gates.append(GateInstance(name="CZ", qubits=(q0, q1), theta=None))
                continue

            rxx_match = _RXX_RE.match(line)
            if rxx_match:
                theta = float(rxx_match.group(1))
                q0 = int(rxx_match.group(2))
                q1 = int(rxx_match.group(3))
                if q0 == q1:
                    raise ValueError(f"Invalid two-qubit gate in line: {line}")
                if q0 > q1:
                    q0, q1 = q1, q0
                gates.append(GateInstance(name="RXX", qubits=(q0, q1), theta=theta))
                continue

    if num_qubits is None:
        raise ValueError("Could not find qreg declaration in QASM file")

    return num_qubits, gates


def snap_to_pool(gates: list[GateInstance], pool) -> list[GateInstance]:
    """Restore exact pool angles for QASM inputs generated from a discrete pool."""
    return [pool.snap(gate) for gate in gates]


def write_qasm_subset(path: str | Path, num_qubits: int, gates: list[GateInstance]) -> None:
    lines = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f" qreg q[{num_qubits}]; ",
        " creg c[1]; ",
    ]

    for gate in gates:
        if gate.name in {"RX", "RY", "RZ"}:
            lines.append(f"{gate.name.lower()}({gate.theta:.4f}) q[{gate.qubits[0]}];")
        elif gate.name == "RXX":
            lines.append(f"rxx({gate.theta:.4f}) q[{gate.qubits[0]}], q[{gate.qubits[1]}];")
        elif gate.name == "CZ":
            lines.append(f"cz q[{gate.qubits[0]}], q[{gate.qubits[1]}];")
        else:
            raise ValueError(f"Unsupported gate for writer: {gate.name}")

    with open(path, "w", encoding="utf-8") as file:
        file.write("\n".join(lines) + "\n")
