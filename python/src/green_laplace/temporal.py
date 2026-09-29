"""Single-direction Green--Laplace reconstruction from ``eta1(t)``."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .dispersion import dimensionless_wavenumber_from_frequency
from .quadrature import gauss_laguerre_rule
from ._eta20 import unidirectional_eta20


@dataclass(frozen=True)
class TimeSeriesResult:
    """Time histories and execution information."""

    time: NDArray[np.float64]
    eta: NDArray[np.float64]
    components: dict[str, NDArray[np.float64] | NDArray[np.complex128]]
    audit: dict[str, Any]


def reconstruct_unidirectional_timeseries(
    eta1: ArrayLike,
    time: ArrayLike,
    *,
    depth: float,
    gravity: float = 9.81,
    omega_max: float | None = None,
    energy_fraction: float = 0.999,
    quadrature_rank: int = 8,
    order: int = 2,
    domain_lengths: Sequence[float] = (25.0, 50.0, 100.0, 200.0),
    relative_tolerance: float = 5.0e-4,
    spatial_points: Sequence[int] | None = None,
    include_eta20: bool = True,
    eta20_rank: int = 16,
    eta20_grid: tuple[int, int] = (128, 16),
) -> TimeSeriesResult:
    """Return ``eta11(t)`` and positive sum-frequency bound-wave components.

    The input is a uniformly sampled real first-order elevation record.  The
    GL calculation uses filtered fields and FFTs; it does not enumerate wave
    pairs.
    """
    values = np.asarray(eta1, dtype=float).reshape(-1)
    times = np.asarray(time, dtype=float).reshape(-1)
    if values.size != times.size or values.size < 9:
        raise ValueError("eta1 and time must have the same length of at least nine")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(times)):
        raise ValueError("eta1 and time must be finite")
    if depth <= 0 or gravity <= 0:
        raise ValueError("depth and gravity must be positive")
    if not 0.0 < energy_fraction <= 1.0:
        raise ValueError("energy_fraction must be in (0, 1]")
    if quadrature_rank not in (4, 6, 8, 12, 16):
        raise ValueError("quadrature_rank must be one of 4, 6, 8, 12, 16")
    if order not in (2, 3):
        raise ValueError("order must be 2 or 3")
    if order == 3 and quadrature_rank not in (4, 6, 8):
        raise ValueError("third-order time reconstruction supports rank 4, 6, or 8")

    dt = float(np.mean(np.diff(times)))
    expected = times[0] + np.arange(times.size) * dt
    deviation = float(np.max(np.abs(times - expected)))
    if dt <= 0 or deviation > max(1.0e-10, 1024.0 * np.finfo(float).eps * max(1.0, times[-1] - times[0])):
        raise ValueError("a uniformly increasing time record is required")

    count = values.size
    period = count * dt
    positive_max = (count - 1) // 2
    all_bins = np.arange(1, positive_max + 1)
    all_omega = 2.0 * np.pi * all_bins / period
    coefficients = np.fft.fft(values) / count
    positive_energy = np.abs(coefficients[all_bins]) ** 2
    total_positive_energy = float(np.sum(positive_energy))
    if total_positive_energy == 0.0:
        raise ValueError("eta1 has no positive-frequency content")

    cumulative = np.cumsum(positive_energy) / total_positive_energy
    energy_index = min(int(np.searchsorted(cumulative, energy_fraction)), all_bins.size - 1)
    energy_bin = int(all_bins[energy_index])
    representable_parent_bin = positive_max // order
    selected_max_bin = min(energy_bin, representable_parent_bin)
    if omega_max is not None:
        if omega_max <= 0:
            raise ValueError("omega_max must be positive")
        selected_max_bin = min(
            selected_max_bin, int(np.floor(omega_max * period / (2.0 * np.pi)))
        )
    if selected_max_bin < 1:
        raise ValueError("the selected eta1 frequency band is empty")

    bins = np.arange(1, selected_max_bin + 1)
    omega = 2.0 * np.pi * bins / period
    amplitudes = 2.0 * np.conj(coefficients[bins])
    base = np.zeros(count, dtype=np.complex128)
    base[bins] = amplitudes
    eta1_used = np.real(np.fft.fft(base))
    retained_energy = float(np.sum(positive_energy[:selected_max_bin]) / total_positive_energy)

    q = dimensionless_wavenumber_from_frequency(omega, depth=depth, gravity=gravity)
    nu = omega * np.sqrt(depth / gravity)
    peak_index = int(np.argmax(np.abs(amplitudes)))
    peak_q = float(q[peak_index])
    scale = float(
        np.sqrt(4.0 * peak_q * np.tanh(peak_q) - 2.0 * peak_q * np.tanh(2.0 * peak_q))
    )
    nodes, weights = gauss_laguerre_rule(quadrature_rank)

    work_count = min(count, _next_power_of_two(order * selected_max_bin + 1))
    if work_count <= order * selected_max_bin:
        raise ValueError("the time grid cannot represent the selected nonlinear sums")
    output_bins = np.arange(2, min(2 * selected_max_bin, positive_max) + 1)
    s = 2.0 * np.pi * output_bins / period * np.sqrt(depth / gravity)

    lower = np.maximum(1, output_bins - selected_max_bin)
    upper = output_bins - lower
    q_max = q[lower - 1] + q[upper - 1]
    q_min = q[np.floor(output_bins / 2).astype(int) - 1] + q[
        np.ceil(output_bins / 2).astype(int) - 1
    ]

    levels: list[dict[str, float | int | None]] = []
    previous: NDArray[np.complex128] | None = None
    analytic: NDArray[np.complex128] | None = None
    analytic3: NDArray[np.complex128] | None = None
    lengths = tuple(float(value) for value in domain_lengths)
    if not lengths or any(value <= 0 for value in lengths):
        raise ValueError("domain_lengths must contain positive values")
    if spatial_points is not None and len(spatial_points) != len(lengths):
        raise ValueError("spatial_points must match domain_lengths")

    for level_index, length in enumerate(lengths):
        if spatial_points is None:
            minimum = max(16, int(np.floor(2 * np.max(q) * length / np.pi + 2)) + 1)
            nx = _next_power_of_two(minimum)
        else:
            nx = int(spatial_points[level_index])
        if nx < 8 or nx % 2:
            raise ValueError("every auxiliary spatial grid size must be positive and even")

        mode = np.fft.fftfreq(nx) * nx
        x = mode * length / nx
        output_k = mode * 2.0 * np.pi / length
        input_spectrum = np.zeros((nx, work_count), dtype=np.complex128)
        input_spectrum[:, bins] = np.exp(1j * np.outer(x, q)) * amplitudes
        filters = np.zeros((6, work_count), dtype=float)
        filters[:, bins] = np.vstack(
            (np.ones_like(q), nu, nu**2, q / nu, q**2 / nu, q)
        )
        v = np.fft.fft(input_spectrum * filters[0], axis=1)
        v_nu = np.fft.fft(input_spectrum * filters[1], axis=1)
        v_nu2 = np.fft.fft(input_spectrum * filters[2], axis=1)
        hx = np.fft.fft(input_spectrum * filters[3], axis=1)
        radial = np.fft.fft(input_spectrum * filters[4], axis=1)
        jx = np.fft.fft(input_spectrum * filters[5], axis=1)
        if order == 3:
            extra = np.zeros((3, work_count), dtype=float)
            extra[:, bins] = np.vstack((nu * q, q**2 * nu, q**2))
            nu_jx = np.fft.fft(input_spectrum * extra[0], axis=1)
            q2_nu = np.fft.fft(input_spectrum * extra[1], axis=1)
            q2_field = np.fft.fft(input_spectrum * extra[2], axis=1)
        source_d = 2.0 * v * v_nu2 + v_nu**2 - hx**2
        source_k = 2.0 * v * radial + 2.0 * hx * jx
        source_d_spectrum = np.fft.fft(
            np.fft.ifft(source_d, axis=1)[:, output_bins], axis=0
        )
        source_k_spectrum = np.fft.fft(
            np.fft.ifft(source_k, axis=1)[:, output_bins], axis=0
        )

        response_q = np.clip(output_k[:, None], q_min[None, :], q_max[None, :])
        response_a = np.sqrt(response_q * np.tanh(response_q))
        response_d = np.zeros_like(response_a)
        response_k = np.zeros_like(response_a)
        for node, weight in zip(nodes, weights, strict=True):
            tau = node / scale
            plus = weight / scale * np.exp(node - (s[None, :] - response_a) * tau)
            minus = weight / scale * np.exp(node - (s[None, :] + response_a) * tau)
            response_d -= response_a * (plus - minus) / 2.0
            response_k += (plus + minus) / 2.0

        output_coefficients = np.zeros(count, dtype=np.complex128)
        output_coefficients[output_bins] = np.sum(
            response_d * source_d_spectrum + response_k * source_k_spectrum,
            axis=0,
        ) / (4.0 * depth * nx)
        current = np.fft.fft(output_coefficients)
        current3 = None
        if order == 3:
            lambda2 = 2.0 * np.sqrt(peak_q * np.tanh(peak_q)) - np.sqrt(
                2.0 * peak_q * np.tanh(2.0 * peak_q)
            )
            lambda3 = 3.0 * np.sqrt(peak_q * np.tanh(peak_q)) - np.sqrt(
                3.0 * peak_q * np.tanh(3.0 * peak_q)
            )
            response_s2, response_c2 = _node_response(
                s, response_a, lambda2, nodes, weights
            )
            eta2_spectrum = (
                -response_a**2 * response_s2 * source_d_spectrum
                + response_c2 * source_k_spectrum
            ) / 4.0
            phi2_spectrum = 1j * (
                response_c2 * source_d_spectrum
                - response_s2 * source_k_spectrum
            ) / 4.0
            lower = _restore(eta2_spectrum, output_bins, work_count)
            forcing_k3 = 1j * lower * radial
            forcing_d3 = -lower * v_nu2
            lower = _restore(1j * response_q * eta2_spectrum, output_bins, work_count)
            forcing_k3 += hx * lower
            lower = _restore(1j * response_q * phi2_spectrum, output_bins, work_count)
            forcing_k3 += 1j * lower * jx
            forcing_d3 += hx * lower
            lower = _restore(response_q**2 * phi2_spectrum, output_bins, work_count)
            forcing_k3 -= v * lower
            lower = _restore(
                -1j * s[None, :] * response_a**2 * phi2_spectrum,
                output_bins,
                work_count,
            )
            forcing_d3 += v * lower
            lower = _restore(response_a**2 * phi2_spectrum, output_bins, work_count)
            forcing_d3 -= 1j * v_nu * lower
            forcing_k3 = forcing_k3 / 2.0 + (
                1j * v * nu_jx * jx + 0.5j * v**2 * q2_nu
            ) / 4.0
            forcing_d3 = forcing_d3 / 2.0 + (
                -0.5 * v**2 * q2_field + v * (hx * nu_jx - v_nu * radial)
            ) / 4.0
            output_bins3 = np.arange(3, 3 * selected_max_bin + 1)
            source_d3 = np.fft.fft(
                np.fft.ifft(forcing_d3, axis=1)[:, output_bins3], axis=0
            )
            source_k3 = np.fft.fft(
                np.fft.ifft(forcing_k3, axis=1)[:, output_bins3], axis=0
            )
            q3_low, q3_high = _cubic_support(output_bins3, bins, q)
            response_q3 = np.clip(output_k[:, None], q3_low[None, :], q3_high[None, :])
            response_a3 = np.sqrt(response_q3 * np.tanh(response_q3))
            s3 = 2.0 * np.pi * output_bins3 / period * np.sqrt(depth / gravity)
            response_s3, response_c3 = _node_response(
                s3, response_a3, lambda3, nodes, weights
            )
            output_coefficients3 = np.zeros(count, dtype=np.complex128)
            output_coefficients3[output_bins3] = np.sum(
                response_a3**2 * response_s3 * source_d3
                - 1j * response_c3 * source_k3,
                axis=0,
            ) / (nx * depth**2)
            current3 = np.fft.fft(output_coefficients3)
        change: float | None = None
        packed = current if current3 is None else np.column_stack((current, current3))
        if previous is not None:
            change = float(
                np.max(
                    np.linalg.norm(packed - previous, axis=0)
                    / np.maximum(np.linalg.norm(packed, axis=0), np.finfo(float).tiny)
                )
            )
        levels.append({"length_over_depth": length, "points": nx, "relative_change": change})
        analytic = current
        analytic3 = current3
        previous = packed
        if change is not None and change <= relative_tolerance:
            break

    assert analytic is not None
    eta22 = np.real(analytic)
    components: dict[str, NDArray[np.float64] | NDArray[np.complex128]] = {
        "eta11": eta1_used,
        "eta22": eta22,
        "eta22_analytic": analytic,
    }
    eta_total = eta1_used + eta22
    if include_eta20:
        eta20, eta20_audit = unidirectional_eta20(
            amplitudes,
            omega,
            q,
            times,
            depth=depth,
            gravity=gravity,
            rank=eta20_rank,
            nx=eta20_grid[0],
            ny=eta20_grid[1],
        )
        components["eta20"] = eta20
        eta_total = eta_total + eta20
    if order == 3:
        assert analytic3 is not None
        components["eta33"] = np.real(analytic3)
        components["eta33_analytic"] = analytic3
        eta_total = eta_total + np.real(analytic3)
    audit = {
        "implementation": "green-laplace-python-unidirectional-time-eta22",
        "input": "real eta1(t)",
        "quadrature_rank": quadrature_rank,
        "order": order,
        "peak_depth_wavenumber": peak_q,
        "parent_bins": bins.tolist(),
        "parent_omega": omega.tolist(),
        "parent_wavenumber": (q / depth).tolist(),
        "retained_positive_frequency_energy": retained_energy,
        "input_projection_relative": float(
            np.linalg.norm(eta1_used - values)
            / max(np.linalg.norm(values), np.finfo(float).tiny)
        ),
        "excluded_dc": float(np.real(coefficients[0])),
        "time_grid_deviation": deviation,
        "interaction_enumeration": False,
        "stokes_correction": False,
        "levels": levels,
    }
    if include_eta20:
        audit["eta20"] = eta20_audit
    return TimeSeriesResult(
        time=times.copy(),
        eta=eta_total,
        components=components,
        audit=audit,
    )


def _next_power_of_two(value: int) -> int:
    return 1 << (int(value) - 1).bit_length()


def _node_response(s, response_a, scale, nodes, weights):
    if scale <= 0:
        raise ValueError("invalid Green--Laplace response scale")
    response_s = np.zeros_like(response_a)
    response_c = np.zeros_like(response_a)
    for node, weight in zip(nodes, weights, strict=True):
        tau = node / scale
        plus = weight / scale * np.exp(node - (s[None, :] - response_a) * tau)
        minus = weight / scale * np.exp(node - (s[None, :] + response_a) * tau)
        response_s += (plus - minus) / (2.0 * response_a)
        response_c += (plus + minus) / 2.0
    return response_s, response_c


def _restore(coefficients, bins, count):
    spectrum = np.zeros((coefficients.shape[0], count), dtype=np.complex128)
    spectrum[:, bins] = np.fft.ifft(coefficients, axis=0)
    return np.fft.fft(spectrum, axis=1)


def _cubic_support(output_bins, parent_bins, q):
    minimum = int(parent_bins[0])
    maximum = int(parent_bins[-1])
    base = np.floor(output_bins / 3).astype(int)
    extra = output_bins % 3
    low = (3 - extra) * q[base - minimum] + extra * q[
        np.minimum(base + 1, maximum) - minimum
    ]
    remaining = output_bins.copy()
    high = np.zeros_like(output_bins, dtype=float)
    for leaf in range(3):
        index = np.minimum(maximum, remaining - (2 - leaf) * minimum)
        high += q[index - minimum]
        remaining -= index
    return low, high
