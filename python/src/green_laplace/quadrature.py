"""Quadrature rules used by the Green--Laplace operators."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray


def gauss_laguerre_rule(rank: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return the MATLAB-compatible standard Gauss--Laguerre rule."""
    if not isinstance(rank, int) or rank < 1:
        raise ValueError("rank must be a positive integer")

    diagonal = 2.0 * np.arange(1, rank + 1) - 1.0
    off_diagonal = np.arange(1, rank, dtype=float)
    jacobi = np.diag(diagonal)
    if rank > 1:
        jacobi += np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
    nodes, vectors = np.linalg.eigh(jacobi)
    weights = vectors[0, :] ** 2

    for degree in range(2 * rank):
        exact = math.factorial(degree)
        relative_error = abs(np.sum(weights * nodes**degree) - exact) / max(1, exact)
        if relative_error > 3.0e-11:
            raise ArithmeticError(
                f"Gauss--Laguerre moment check failed at rank={rank}, degree={degree}"
            )
    return nodes, weights

