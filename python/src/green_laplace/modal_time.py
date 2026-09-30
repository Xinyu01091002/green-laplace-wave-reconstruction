"""Single-direction time reconstruction on a one-dimensional modal grid."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._eta20 import unidirectional_eta20
from .dispersion import dimensionless_wavenumber_from_frequency
from .quadrature import gauss_laguerre_rule


@dataclass(frozen=True)
class ModalTimeSeriesResult:
    time: NDArray[np.float64]
    eta: NDArray[np.float64]
    components: dict[str, NDArray[np.float64] | NDArray[np.complex128]]
    audit: dict[str, Any]


def reconstruct_unidirectional_modal_timeseries(
    eta1: ArrayLike,
    time: ArrayLike,
    *,
    depth: float,
    gravity: float = 9.81,
    omega_max: float | None = None,
    energy_fraction: float = 0.999,
    quadrature_rank: int = 8,
    order: int = 3,
    modal_points: int | Sequence[int] = (2048, 4096, 8192, 16384, 32768, 65536),
    modal_limit: float | None = None,
    relative_tolerance: float = 1.0e-3,
    include_eta20: bool = True,
    eta20_rank: int = 16,
    eta20_grid: tuple[int, int] = (128, 16),
) -> ModalTimeSeriesResult:
    """Reconstruct a single-direction record with a direct modal resolvent."""
    started = perf_counter()
    values = np.asarray(eta1, dtype=float).reshape(-1)
    times = np.asarray(time, dtype=float).reshape(-1)
    if values.size != times.size or values.size < 9:
        raise ValueError("eta1 and time must have the same length of at least nine")
    if depth <= 0 or gravity <= 0 or order not in (2, 3):
        raise ValueError("depth and gravity must be positive and order must be 2 or 3")
    if isinstance(modal_points, int):
        modal_ladder = (modal_points,)
    else:
        modal_ladder = tuple(int(value) for value in modal_points)
    if (
        not modal_ladder
        or any(value < 16 or value % 2 for value in modal_ladder)
        or any(right <= left for left, right in zip(modal_ladder, modal_ladder[1:]))
    ):
        raise ValueError("modal_points must be an increasing sequence of even integers")
    if relative_tolerance <= 0:
        raise ValueError("relative_tolerance must be positive")
    if order == 3 and quadrature_rank not in (4, 6, 8, 10):
        raise ValueError("third-order modal reconstruction supports rank 4, 6, 8, or 10")

    count = values.size
    dt = float(np.mean(np.diff(times)))
    expected = times[0] + np.arange(count) * dt
    deviation = float(np.max(np.abs(times - expected)))
    if dt <= 0 or deviation > max(1e-10, 1024 * np.finfo(float).eps * max(1.0, times[-1] - times[0])):
        raise ValueError("a uniformly increasing time record is required")
    period = count * dt
    positive_max = (count - 1) // 2
    all_bins = np.arange(1, positive_max + 1)
    coefficients = np.fft.fft(values) / count
    energy = np.abs(coefficients[all_bins]) ** 2
    total_energy = float(np.sum(energy))
    if total_energy == 0:
        raise ValueError("eta1 has no positive-frequency content")
    cumulative = np.cumsum(energy) / total_energy
    energy_index = min(int(np.searchsorted(cumulative, energy_fraction)), all_bins.size - 1)
    maximum_bin = int(all_bins[energy_index])
    maximum_bin = min(maximum_bin, positive_max // order)
    if omega_max is not None:
        maximum_bin = min(maximum_bin, int(np.floor(omega_max * period / (2 * np.pi))))
    if maximum_bin < 1:
        raise ValueError("the selected eta1 frequency band is empty")

    bins = np.arange(1, maximum_bin + 1)
    omega = 2 * np.pi * bins / period
    amplitudes = 2 * np.conj(coefficients[bins])
    base = np.zeros(count, dtype=np.complex128)
    base[bins] = amplitudes
    eta11 = np.real(np.fft.fft(base))
    q = dimensionless_wavenumber_from_frequency(omega, depth=depth, gravity=gravity)
    peak_q = float(q[int(np.argmax(np.abs(amplitudes)))])

    previous = None
    levels = []
    converged = False
    eta22 = None
    eta33 = None
    modal_audit = None
    for points in modal_ladder:
        level_started = perf_counter()
        eta22, eta33, modal_audit = _modal_sum_components(
            amplitudes,
            bins,
            q,
            count=count,
            dt=dt,
            depth=depth,
            gravity=gravity,
            peak_q=peak_q,
            rank=quadrature_rank,
            order=order,
            modal_points=points,
            modal_limit=modal_limit,
        )
        current = (eta22,) if eta33 is None else (eta22, eta33)
        changes = None
        if previous is not None:
            changes = tuple(
                float(np.linalg.norm(value - old) / np.linalg.norm(value))
                for value, old in zip(current, previous, strict=True)
            )
            converged = max(changes) <= relative_tolerance
        levels.append(
            {
                "modal_points": points,
                "seconds": perf_counter() - level_started,
                "eta22_change": None if changes is None else changes[0],
                "eta33_change": None if changes is None or eta33 is None else changes[1],
                "wavevector_projection_rms_relative": modal_audit[
                    "wavevector_projection_rms_relative"
                ],
            }
        )
        previous = current
        if converged:
            break
    assert eta22 is not None and modal_audit is not None
    components: dict[str, NDArray[np.float64] | NDArray[np.complex128]] = {
        "eta11": eta11,
        "eta22": np.real(eta22),
        "eta22_analytic": eta22,
    }
    eta_total = eta11 + np.real(eta22)
    if eta33 is not None:
        components["eta33"] = np.real(eta33)
        components["eta33_analytic"] = eta33
        eta_total += np.real(eta33)
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
        eta_total += eta20
        modal_audit["eta20"] = eta20_audit
    modal_audit.update(
        {
            "implementation": "one-dimensional-modal-green-laplace-time",
            "retained_positive_frequency_energy": float(np.sum(energy[:maximum_bin]) / total_energy),
            "input_projection_relative": float(
                np.linalg.norm(eta11 - values) / max(np.linalg.norm(values), np.finfo(float).tiny)
            ),
            "solve_seconds": perf_counter() - started,
            "levels": levels,
            "converged": converged,
            "relative_tolerance": relative_tolerance,
            "auxiliary_space_domain": False,
            "interaction_enumeration": False,
            "stokes_correction": False,
        }
    )
    return ModalTimeSeriesResult(times.copy(), eta_total, components, modal_audit)


def _modal_sum_components(
    amplitudes,
    bins,
    q,
    *,
    count,
    dt,
    depth,
    gravity,
    peak_q,
    rank,
    order,
    modal_points,
    modal_limit,
):
    base_bin = int(np.min(bins))
    work_bins = bins - base_bin
    maximum_work_bin = int(np.max(work_bins))
    work_count = order * maximum_work_bin + 1
    if work_count > count:
        raise ValueError("the nonlinear output band does not fit in the time record")
    q_axis, dq = _modal_axis(modal_points, q, order, modal_limit)
    joint_shifted, projection = _deposit_modal(amplitudes, work_bins, q, q_axis, dq, work_count)
    joint = np.fft.ifftshift(joint_shifted, axes=0)
    q_fft = np.fft.ifftshift(q_axis)
    output_q = np.abs(q_fft)
    signed_q = q_fft[:, None]
    q2 = signed_q**2
    frequency_step = 2 * np.pi / (count * dt) * np.sqrt(depth / gravity)
    input_nu_axis = (base_bin + np.arange(work_count)) * frequency_step
    nu = input_nu_axis[None, :]
    safe_nu = np.maximum(nu, np.finfo(float).tiny)

    v = _modal_field(joint, 1.0, modal_points)
    vn = _modal_field(joint, nu, modal_points)
    vn2 = _modal_field(joint, nu**2, modal_points)
    hx = _modal_field(joint, signed_q / safe_nu, modal_points)
    jx = _modal_field(joint, signed_q, modal_points)
    radial = _modal_field(joint, q2 / safe_nu, modal_points)
    source_d = 2 * v * vn2 + vn**2 - hx**2
    source_k = 2 * v * radial + 2 * hx * jx
    source_d_hat = _modal_spectrum(source_d, modal_points)
    source_k_hat = _modal_spectrum(source_k, modal_points)

    nodes, weights = gauss_laguerre_rule(rank)
    peak_nu = np.sqrt(peak_q * np.tanh(peak_q))
    standalone_scale = np.sqrt(4 * peak_q * np.tanh(peak_q) - 2 * peak_q * np.tanh(2 * peak_q))
    output_work_bins2 = np.arange(0, 2 * maximum_work_bin + 1)
    output_bins2 = 2 * base_bin + output_work_bins2
    eta22_coefficients = np.zeros(count, dtype=np.complex128)
    invalid2 = 0.0
    total2 = 0.0
    response_a = np.sqrt(output_q * np.tanh(output_q))
    for work_bin, output_bin in zip(output_work_bins2, output_bins2, strict=True):
        s = output_bin * frequency_step
        sd = source_d_hat[:, work_bin]
        sk = source_k_hat[:, work_bin]
        source_energy = np.abs(sd) ** 2 + np.abs(sk) ** 2
        valid = s > response_a
        total2 += float(np.sum(source_energy))
        invalid2 += float(np.sum(source_energy[~valid]))
        S, C = _node_response(s, response_a, valid, standalone_scale, nodes, weights)
        response = np.zeros(modal_points, dtype=np.complex128)
        aq = output_q * np.tanh(output_q)
        response[valid] = (-aq[valid] * S[valid] * sd[valid] + C[valid] * sk[valid]) / (4 * depth)
        eta22_coefficients[output_bin] = np.sum(response)
    eta22 = np.fft.fft(eta22_coefficients)

    eta33 = None
    invalid3 = None
    if order == 3:
        aq = output_q * np.tanh(output_q)
        lambda2 = 2 * peak_nu - np.sqrt(2 * peak_q * np.tanh(2 * peak_q))
        lambda3 = 3 * peak_nu - np.sqrt(3 * peak_q * np.tanh(3 * peak_q))
        inner_nu_axis = (2 * base_bin + np.arange(work_count)) * frequency_step
        e2 = np.zeros((modal_points, work_count), dtype=np.complex128)
        p2 = np.zeros_like(e2)
        for work_bin in output_work_bins2:
            s = inner_nu_axis[work_bin]
            sd = source_d_hat[:, work_bin]
            sk = source_k_hat[:, work_bin]
            valid = s > response_a
            S, C = _node_response(s, response_a, valid, lambda2, nodes, weights)
            e2[valid, work_bin] = (-aq[valid] * S[valid] * sd[valid] + C[valid] * sk[valid]) / 4
            p2[valid, work_bin] = 1j * (C[valid] * sd[valid] - S[valid] * sk[valid]) / 4

        first = {
            "eta": v,
            "etax": _modal_field(joint, 1j * signed_q, modal_points),
            "phix": hx,
            "phiz": _modal_field(joint, -1j * nu, modal_points),
            "phixz": _modal_field(joint, nu * signed_q, modal_points),
            "phizz": _modal_field(joint, -1j * q2 / safe_nu, modal_points),
            "phizzz": _modal_field(joint, -1j * q2 * nu, modal_points),
            "phitz": _modal_field(joint, -nu**2, modal_points),
            "phitzz": _modal_field(joint, -q2, modal_points),
        }
        aq2 = aq[:, None]
        pair = {
            "eta": _spectral_field(e2, modal_points),
            "etax": _spectral_field(e2 * (1j * signed_q), modal_points),
            "phix": _spectral_field(p2 * (1j * signed_q), modal_points),
            "phiz": _spectral_field(p2 * aq2, modal_points),
            "phizz": _spectral_field(p2 * q2, modal_points),
            "phitz": _spectral_field(p2 * (-1j * inner_nu_axis[None, :]) * aq2, modal_points),
        }
        forcing_k_pair = (
            first["phix"] * pair["etax"]
            + pair["phix"] * first["etax"]
            - first["eta"] * pair["phizz"]
            - pair["eta"] * first["phizz"]
        )
        forcing_k_direct = first["eta"] * first["phixz"] * first["etax"] - 0.5 * first["eta"] ** 2 * first["phizzz"]
        forcing_d_pair = (
            first["eta"] * pair["phitz"]
            + pair["eta"] * first["phitz"]
            + first["phix"] * pair["phix"]
            + first["phiz"] * pair["phiz"]
        )
        forcing_d_direct = 0.5 * first["eta"] ** 2 * first["phitzz"] + first["eta"] * (
            first["phix"] * first["phixz"] + first["phiz"] * first["phizz"]
        )
        forcing_k_hat = _modal_spectrum(forcing_k_direct / 4 + forcing_k_pair / 2, modal_points)
        forcing_d_hat = _modal_spectrum(forcing_d_direct / 4 + forcing_d_pair / 2, modal_points)
        output_work_bins3 = np.arange(0, 3 * maximum_work_bin + 1)
        output_bins3 = 3 * base_bin + output_work_bins3
        eta33_coefficients = np.zeros(count, dtype=np.complex128)
        invalid_energy3 = 0.0
        total_energy3 = 0.0
        for work_bin, output_bin in zip(output_work_bins3, output_bins3, strict=True):
            s = output_bin * frequency_step
            fk = forcing_k_hat[:, work_bin]
            fd = forcing_d_hat[:, work_bin]
            source_energy = np.abs(fk) ** 2 + np.abs(fd) ** 2
            valid = s > response_a
            total_energy3 += float(np.sum(source_energy))
            invalid_energy3 += float(np.sum(source_energy[~valid]))
            S, C = _node_response(s, response_a, valid, lambda3, nodes, weights)
            response = np.zeros(modal_points, dtype=np.complex128)
            response[valid] = (aq[valid] * S[valid] * fd[valid] - 1j * C[valid] * fk[valid]) / depth**2
            eta33_coefficients[output_bin] = np.sum(response)
        eta33 = np.fft.fft(eta33_coefficients)
        invalid3 = invalid_energy3 / max(total_energy3, np.finfo(float).tiny)

    audit = {
        "modal_points": modal_points,
        "modal_limit": float(abs(q_axis[0])),
        "quadrature_rank": rank,
        "base_bin": base_bin,
        "work_frequency_count": work_count,
        "wavevector_projection_rms_relative": projection,
        "eta22_invalid_source_energy_fraction": invalid2 / max(total2, np.finfo(float).tiny),
        "eta33_invalid_source_energy_fraction": invalid3,
    }
    return eta22, eta33, audit


def _modal_axis(count, parent, order, requested_limit):
    if requested_limit is None:
        limit = (order + 0.05) * float(np.max(np.abs(parent))) / (1 - 2 / count)
    else:
        limit = float(requested_limit)
        if limit <= order * float(np.max(np.abs(parent))):
            raise ValueError("modal_limit must exceed the requested nonlinear output support")
    spacing = 2 * limit / count
    return np.arange(-count / 2, count / 2) * spacing, spacing


def _deposit_modal(amplitudes, bins, q, axis, spacing, work_count):
    position = (q - axis[0]) / spacing
    lower = np.clip(np.floor(position).astype(int), 0, axis.size - 2)
    upper = lower + 1
    upper_weight = np.clip(position - lower, 0.0, 1.0)
    lower_weight = 1.0 - upper_weight
    joint = np.zeros((axis.size, work_count), dtype=np.complex128)
    np.add.at(joint, (lower, bins), amplitudes * lower_weight)
    np.add.at(joint, (upper, bins), amplitudes * upper_weight)
    mass = np.abs(amplitudes) ** 2
    errors = np.concatenate((np.abs(axis[lower] - q) / q, np.abs(axis[upper] - q) / q))
    weights = np.concatenate((mass * lower_weight, mass * upper_weight))
    projection = float(np.sqrt(np.sum(weights * errors**2) / max(np.sum(weights), np.finfo(float).tiny)))
    return joint, projection


def _modal_field(joint, multiplier, points):
    return np.fft.fft(np.fft.ifft(joint * multiplier, axis=0) * points, axis=1)


def _spectral_field(spectrum, points):
    return np.fft.fft(np.fft.ifft(spectrum, axis=0) * points, axis=1)


def _modal_spectrum(field, points):
    return np.fft.fft(np.fft.ifft(field, axis=1), axis=0) / points


def _node_response(s, response_a, valid, scale, nodes, weights):
    S = np.zeros_like(response_a)
    C = np.zeros_like(response_a)
    nonzero = valid & (response_a > 0)
    zero = valid & (response_a == 0)
    av = response_a[nonzero]
    for node, weight in zip(nodes, weights, strict=True):
        tau = node / scale
        plus = weight / scale * np.exp(node - (s - av) * tau)
        minus = weight / scale * np.exp(node - (s + av) * tau)
        S[nonzero] += (plus - minus) / (2 * av)
        C[nonzero] += (plus + minus) / 2
        if np.any(zero):
            base = weight / scale * np.exp(node - s * tau)
            S[zero] += base * tau
            C[zero] += base
    return S, C
