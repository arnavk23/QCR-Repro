"""Vectorized window-unitary computation for the exhaustive sweep.

Same greedy protocol and results as the scalar sweep (see scripts/check_batched_vs_scalar.py); only the wall-clock time differs."""

from __future__ import annotations

import hashlib
from typing import Optional

import numpy as np

from .config import GateInstance
from .gates import embedded_gate_matrix

_MISSING = object()


def _wire_mask(gates: list[GateInstance]) -> np.ndarray:
    """Per-gate bitmask of the wires the gate touches (int64)."""
    n = len(gates)
    masks = np.zeros(n, dtype=np.int64)
    for i, gate in enumerate(gates):
        m = 0
        for q in gate.qubits:
            m |= 1 << q
        masks[i] = m
    return masks


def _bits_to_list(mask: int) -> list[int]:
    return [q for q in range(mask.bit_length()) if (mask >> q) & 1]


def lookup_batch(graph, unitaries: np.ndarray) -> list[Optional[tuple[int, ...]]]:
    """ComputeGraph.lookup for a batch of (B, d, d) unitaries.

Keys are bit-identical to scalar lookups; returns token chains (or None)."""
    batch = unitaries.shape[0]
    dim = unitaries.shape[1]
    flat = unitaries.reshape(batch, dim * dim)
    mag = np.abs(flat)
    idx = np.argmax(np.round(mag, 8), axis=1)
    phase = np.angle(flat[np.arange(batch), idx])
    norm = flat * np.exp(-1j * phase)[:, None]
    rounded = np.round(norm, graph.digest_decimals) + 0.0
    keys = [hashlib.sha256(rounded[i].tobytes()).digest() for i in range(batch)]
    return [graph.buckets.get(key) for key in keys]


def _lookup_group(gates, positions: np.ndarray, length: int, graph, mat_table: np.ndarray) -> list:
    """Look up windows [pos, pos+length) for every pos in positions via batched matmuls."""
    d = 2 ** graph.pool.num_qubits
    gate_idx = positions[:, None] + np.arange(length)
    tensor = mat_table[gate_idx]  # (B, L, d, d), vectorized gather
    u = tensor[:, 0, :, :]
    for k in range(1, length):
        u = tensor[:, k, :, :] @ u
    return lookup_batch(graph, u)


def _compute_windows(
    gates: list[GateInstance],
    db,
    max_block_len: int,
) -> dict[tuple[int, int], tuple]:
    """Vectorized computation of every window unitary of the current circuit.

Returns {(pos, length): (token_chain, wires)}."""
    n = len(gates)
    masks = _wire_mask(gates)
    results: dict[tuple[int, int], tuple] = {}
    # Per-wire-set embedded-matrix tables are reused across lengths.
    tables: dict[int, np.ndarray] = {}
    for length in range(max_block_len, 1, -1):
        batch = n - length + 1
        if batch <= 0:
            continue
        # Sliding OR of the per-gate masks over each length-L window.
        acc = np.zeros(batch, dtype=np.int64)
        for k in range(length):
            acc |= masks[k : k + batch]
        # Group window positions by their wire count.
        for w in np.unique(acc):
            pos_idx = np.nonzero(acc == w)[0]
            graph = db.graphs.get(int(np.bitwise_count(w)))
            if graph is None:
                for pos in pos_idx:
                    results[(int(pos), length)] = (None, _bits_to_list(int(w)))
                continue
            mat_table = tables.get(int(w))
            if mat_table is None:
                wires = _bits_to_list(int(w))
                forward = {wire: i for i, wire in enumerate(wires)}
                dim = 2 ** graph.pool.num_qubits
                table = np.zeros((n, dim, dim), dtype=complex)
                token_mats = graph._token_matrices
                wire_set = set(wires)
                for i, gate in enumerate(gates):
                    if not set(gate.qubits).issubset(wire_set):
                        continue
                    qubits = tuple(sorted(forward[q] for q in gate.qubits))
                    remapped = GateInstance(name=gate.name, qubits=qubits, theta=gate.theta)
                    try:
                        token = graph.pool.token_for_gate(remapped)
                    except KeyError:
                        # Not representable in this pool: the window can never
                        # be in the DB (scalar path returns None for these).
                        continue
                    table[i] = token_mats[token]
                tables[int(w)] = table
                mat_table = table
            wires = _bits_to_list(int(w))
            chains = _lookup_group(gates, pos_idx, length, graph, mat_table)
            for pos, chain in zip(pos_idx, chains):
                results[(int(pos), length)] = (chain, wires)
    return results


def _candidate_from_entry(
    entry: tuple,
    db,
) -> list[GateInstance] | None:
    """Reconstruct the reduced gate list for a stored (chain, wires) entry."""
    chain, wires = entry
    if chain is None:
        return None
    graph = db.graphs.get(len(wires))
    if graph is None:
        return None
    decoded = graph.pool.decode(list(chain))
    reverse = {i: wire for i, wire in enumerate(wires)}
    return [
        GateInstance(
            name=gate.name,
            qubits=tuple(sorted(reverse[q] for q in gate.qubits)),
            theta=gate.theta,
        )
        for gate in decoded
    ]


def batched_sweep(
    gates: list[GateInstance],
    num_qubits: int,
    db,
    max_block_len: int,
) -> int:
    """Exhaustive sweep with vectorized window unitaries.

Same protocol and results as reducer._sweep_reduce; returns replacements applied."""
    total = 0
    while True:
        count = 0
        n = len(gates)
        windows = _compute_windows(gates, db, max_block_len)
        pos = 0
        while pos < n:
            hi = min(max_block_len, n - pos)
            replaced = False
            for length in range(hi, 1, -1):
                entry = windows.get((pos, length), _MISSING)
                if entry is _MISSING:
                    # Window overlapped a previous replacement in this pass:
                    # re-derive it with the scalar path.
                    candidate = db.try_reduce(gates[pos : pos + length])
                else:
                    candidate = _candidate_from_entry(entry, db)
                if candidate is not None and len(candidate) < length:
                    shrink = length - len(candidate)
                    # Shift stored windows to the right of the replaced span
                    # (their gate content is unchanged) and drop those that
                    # overlap it.
                    for p in range(pos + length, n):
                        for l in range(2, max_block_len + 1):
                            key = (p, l)
                            if key in windows:
                                windows[(p - shrink, l)] = windows.pop(key)
                    lo = max(0, pos - max_block_len + 1)
                    for p in range(lo, pos + length):
                        for l in range(2, max_block_len + 1):
                            windows.pop((p, l), None)
                    gates[pos : pos + length] = candidate
                    n = len(gates)
                    count += 1
                    replaced = True
                    break
            if replaced:
                pos = max(0, pos - 1)
            else:
                pos += 1
        total += count
        if count == 0:
            break
    return total
