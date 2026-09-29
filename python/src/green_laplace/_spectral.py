"""Shared FFT conventions and numerical support operations."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def fft_wavenumber_grid(
    lx: float, ly: float, nx: int, ny: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if lx <= 0 or ly <= 0 or nx < 2 or ny < 2 or nx % 2 or ny % 2:
        raise ValueError("positive domain lengths and positive even grid sizes are required")
    kx_axis = 2.0 * np.pi * np.fft.fftfreq(nx, d=lx / nx)
    ky_axis = 2.0 * np.pi * np.fft.fftfreq(ny, d=ly / ny)
    return np.meshgrid(kx_axis, ky_axis)


def analytic_spectrum_from_real_field(
    eta1: NDArray[np.float64],
    *,
    depth: float,
    kx: NDArray[np.float64],
    tolerance: float = 1.0e-11,
) -> tuple[NDArray[np.complex128], float]:
    """Extract the strict-forward analytic spectrum from a real snapshot."""
    if eta1.ndim != 2 or eta1.shape != kx.shape or not np.all(np.isfinite(eta1)):
        raise ValueError("eta1 must be a finite two-dimensional grid")
    if depth <= 0:
        raise ValueError("depth must be positive")

    raw = np.fft.fft2(eta1 / depth)
    scale = max(float(np.max(np.abs(raw))), np.finfo(float).tiny)
    zero_line = np.abs(kx) <= 32.0 * np.finfo(float).eps * max(1.0, np.max(np.abs(kx)))
    if np.max(np.abs(raw[zero_line]), initial=0.0) > tolerance * scale:
        raise ValueError("eta1 contains a non-negligible kx=0 component")

    spectrum = np.zeros_like(raw, dtype=np.complex128)
    spectrum[kx > 0] = 2.0 * raw[kx > 0]
    reconstructed = np.real(depth * np.fft.ifft2(spectrum))
    relative_error = float(
        np.linalg.norm(reconstructed - eta1)
        / max(np.linalg.norm(eta1), np.finfo(float).tiny)
    )
    if relative_error > tolerance:
        raise ValueError(
            "eta1 is not consistent with the declared strict-forward x propagation"
        )
    return spectrum, relative_error


def suppress_fft_roundoff(
    spectrum: NDArray[np.complex128], support: NDArray[np.bool_]
) -> NDArray[np.complex128]:
    result = spectrum.copy()
    active = np.abs(result[support])
    if active.size == 0:
        return result
    scale = float(np.max(active))
    if scale == 0.0:
        return result
    floor = 16.0 * np.finfo(float).eps * max(1.0, np.log2(result.size))
    result[support & (np.abs(result) < floor * scale)] = 0.0
    return result


def assert_quadratic_support_is_alias_safe(support: NDArray[np.bool_]) -> None:
    ny, nx = support.shape
    mode_x = np.fft.fftfreq(nx) * nx
    mode_y = np.fft.fftfreq(ny) * ny
    mx, my = np.meshgrid(mode_x, mode_y)
    if not np.any(support):
        raise ValueError("eta1 has empty analytic support")
    if np.max(2.0 * np.abs(mx[support])) >= nx / 2.0 or np.max(
        2.0 * np.abs(my[support])
    ) >= ny / 2.0:
        raise ValueError("quadratic output support reaches a Nyquist boundary")

