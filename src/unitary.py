from __future__ import annotations

import numpy as np


def _reference_index(unitary: np.ndarray, arground: int = 8) -> int:
    """Index of the phase-normalization reference entry.

Argmax of the *rounded* magnitudes, so ties resolve deterministically and the key is stable for equivalent unitaries."""
    mag = np.abs(unitary.reshape(-1))
    return int(np.argmax(np.round(mag, arground)))


def remove_global_phase(unitary: np.ndarray) -> np.ndarray:
    """Return a copy of `unitary` with a robust, deterministic global phase."""
    idx = _reference_index(unitary)
    phase = np.angle(unitary.reshape(-1)[idx])
    return unitary * np.exp(-1j * phase)


def normalized_stack(unitary: np.ndarray) -> np.ndarray:
    """Phase-normalized (real, imag) interleaved vector of ``unitary``."""
    normalized = remove_global_phase(unitary)
    return np.stack((normalized.real, normalized.imag), axis=-1).reshape(-1)


def coarse_key(unitary: np.ndarray, decimals: int = 4) -> bytes:
    return key_from_stack(normalized_stack(unitary), decimals)


def fine_key(unitary: np.ndarray, decimals: int = 10) -> bytes:
    return key_from_stack(normalized_stack(unitary), decimals)


def key_from_stack(stack: np.ndarray, decimals: int) -> bytes:
    rounded = np.round(stack, decimals=decimals) + 0.0
    return rounded.tobytes()


def equivalent_up_to_global_phase(u: np.ndarray, v: np.ndarray, atol: float = 1e-5) -> bool:
    """True if u and v are equal up to a global phase within atol.

Optimal-global-scale Frobenius residual, as in the reference MATLAB compare.m."""
    if u.shape != v.shape:
        return False
    a = u.reshape(-1)
    b = v.reshape(-1)
    scale = np.vdot(a, b) / np.vdot(a, a)
    residual = np.linalg.norm(b - scale * a)
    return residual <= atol
