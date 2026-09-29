"""Finite-depth linear dispersion utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def intrinsic_omega(
    k: ArrayLike, *, depth: float, gravity: float = 9.81
) -> NDArray[np.float64]:
    """Return ``sqrt(g k tanh(k h))``."""
    values = np.asarray(k, dtype=float)
    if depth <= 0 or gravity <= 0 or np.any(values < 0):
        raise ValueError("depth and gravity must be positive and k must be nonnegative")
    return np.sqrt(gravity * values * np.tanh(depth * values))


def dimensionless_wavenumber_from_frequency(
    omega: ArrayLike, *, depth: float, gravity: float = 9.81
) -> NDArray[np.float64]:
    """Solve ``q tanh(q) = omega**2 h/g`` for nonnegative ``q=k h``."""
    values = np.asarray(omega, dtype=float)
    if depth <= 0 or gravity <= 0 or np.any(values < 0):
        raise ValueError("depth and gravity must be positive and omega nonnegative")

    target = values**2 * depth / gravity
    q = np.where(target < 1.0, np.sqrt(target), target)
    q = np.maximum(q, np.finfo(float).tiny)
    for _ in range(40):
        tanh_q = np.tanh(q)
        residual = q * tanh_q - target
        derivative = tanh_q + q * (1.0 - tanh_q**2)
        update = residual / np.maximum(derivative, np.finfo(float).tiny)
        candidate = np.maximum(q - update, 0.5 * q)
        if np.max(np.abs(candidate - q) / np.maximum(1.0, q)) < 5.0e-15:
            q = candidate
            break
        q = candidate
    if not np.allclose(q * np.tanh(q), target, rtol=2.0e-13, atol=2.0e-15):
        raise ArithmeticError("finite-depth dispersion solve did not converge")
    return q

