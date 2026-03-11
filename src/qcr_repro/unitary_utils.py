from __future__ import annotations

import numpy as np


def remove_global_phase(unitary: np.ndarray, atol: float = 1e-12) -> np.ndarray:
    flat = unitary.reshape(-1)
    idx = None
    for i, value in enumerate(flat):
        if abs(value) > atol:
            idx = i
            break
    if idx is None:
        return unitary.copy()
    phase = np.angle(flat[idx])
    return unitary * np.exp(-1j * phase)


def equivalent_up_to_global_phase(u: np.ndarray, v: np.ndarray, atol: float = 1e-5) -> bool:
    if u.shape != v.shape:
        return False
    return np.allclose(remove_global_phase(u), remove_global_phase(v), atol=atol, rtol=0)


def unitary_key(unitary: np.ndarray, decimals: int = 5) -> tuple[tuple[float, float], ...]:
    normalized = remove_global_phase(unitary)
    rounded = np.round(normalized.real, decimals=decimals) + 1j * np.round(normalized.imag, decimals=decimals)
    return tuple((float(value.real), float(value.imag)) for value in rounded.reshape(-1))
