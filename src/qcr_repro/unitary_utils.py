from __future__ import annotations

import numpy as np


def _reference_index(unitary: np.ndarray, arground: int = 8) -> int:
    """Index of the phase-normalization reference entry.

    The reference is the argmax of the *rounded* magnitudes.  For unitaries
    whose entries have repeated maximal magnitudes (e.g. all entries of
    RX(pi/2) equal 1/sqrt(2)), a raw argmax is unstable under floating-point
    perturbations and would produce different keys for the same unitary.  After
    rounding the magnitudes the tie is resolved deterministically (first index),
    which keeps the key stable for every equivalent computation.
    """
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
    """Check ``u`` and ``v`` are equal up to a global phase within ``atol``.

    Uses the optimal-global-scale residual of the reference MATLAB ``compare.m``:
    the Frobenius norm of ``v - u * scale`` after fitting the scalar ``scale``.
    """
    if u.shape != v.shape:
        return False
    a = u.reshape(-1)
    b = v.reshape(-1)
    scale = np.vdot(a, b) / np.vdot(a, a)
    residual = np.linalg.norm(b - scale * a)
    return residual <= atol
