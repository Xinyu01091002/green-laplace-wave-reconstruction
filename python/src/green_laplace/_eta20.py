"""Nonzero-wavenumber single-direction difference-frequency elevation."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from .quadrature import gauss_laguerre_rule


def unidirectional_eta20(
    amplitudes: NDArray[np.complex128],
    omega: NDArray[np.float64],
    q: NDArray[np.float64],
    time: NDArray[np.float64],
    *,
    depth: float,
    gravity: float,
    rank: int = 16,
    nx: int = 128,
    ny: int = 16,
) -> tuple[NDArray[np.float64], dict[str, Any]]:
    if rank not in (6, 12, 16):
        raise ValueError("eta20 rank must be 6, 12, or 16")
    if nx < 16 or ny < 16 or nx % 2 or ny % 2:
        raise ValueError("eta20 modal-grid sizes must be even and at least 16")
    count = time.size
    dt = float(np.mean(np.diff(time)))
    bins = np.rint(omega * count * dt / (2.0 * np.pi)).astype(int)
    if np.any(bins < 1) or not np.allclose(
        omega, 2.0 * np.pi * bins / (count * dt), rtol=1e-10, atol=1e-12
    ):
        raise ValueError("eta20 input frequencies must lie on temporal DFT bins")
    base_bin = int(np.min(bins))
    work_bins = bins - base_bin
    span = int(np.max(work_bins))
    work_count = 2 * span + 1
    qx = q
    qy = np.zeros_like(q)
    nu = omega * np.sqrt(depth / gravity)

    qx_axis, dqx = _modal_axis(nx, qx, 2)
    qy_limit = (2.05 * float(np.max(q))) / (1.0 - 2.0 / ny)
    qy_axis, dqy = _modal_axis(ny, qy, 2, requested_limit=qy_limit)
    joint_shifted = _deposit_joint(
        amplitudes, work_bins, qx, qy, qx_axis, qy_axis, dqx, dqy, work_count
    )
    joint = np.fft.ifftshift(np.fft.ifftshift(joint_shifted, axes=0), axes=1)
    qx_fft = np.fft.ifftshift(qx_axis)
    qy_fft = np.fft.ifftshift(qy_axis)
    output_qx, output_qy = np.meshgrid(qx_fft, qy_fft)
    output_q = np.hypot(output_qx, output_qy)
    response_a = np.sqrt(output_q * np.tanh(output_q))
    qx3 = output_qx[:, :, None]
    qy3 = output_qy[:, :, None]
    q23 = qx3**2 + qy3**2
    frequency_step = 2.0 * np.pi / (count * dt) * np.sqrt(depth / gravity)
    input_nu = (base_bin + np.arange(work_count)) * frequency_step
    nu3 = input_nu[None, None, :]
    safe_nu = np.maximum(nu3, np.finfo(float).tiny)

    one = _joint_field(joint, 1.0, nx, ny)
    dno = _joint_field(joint, nu3**2, nx, ny)
    nu_field = _joint_field(joint, nu3, nx, ny)
    bx = _joint_field(joint, qx3 / safe_nu, nx, ny)
    by = _joint_field(joint, qy3 / safe_nu, nx, ny)
    qx_field = _joint_field(joint, qx3, nx, ny)
    qy_field = _joint_field(joint, qy3, nx, ny)
    radial = _joint_field(joint, q23 / safe_nu, nx, ny)
    forcing_d = (
        -dno * np.conj(one)
        - one * np.conj(dno)
        + nu_field * np.conj(nu_field)
        + bx * np.conj(bx)
        + by * np.conj(by)
    )
    forcing_k = (
        1j * radial * np.conj(one)
        - 1j * bx * np.conj(qx_field)
        - 1j * by * np.conj(qy_field)
        + 1j * qx_field * np.conj(bx)
        + 1j * qy_field * np.conj(by)
        - 1j * one * np.conj(radial)
    )
    forcing_d_spectrum = _joint_spectrum(forcing_d, nx, ny)
    forcing_k_spectrum = _joint_spectrum(forcing_k, nx, ny)

    energy = np.abs(amplitudes) ** 2
    energy_sum = float(np.sum(energy))
    mean_q = float(np.sum(q * energy) / energy_sum)
    delta_q = float(np.sqrt(2.0 * np.sum((q - mean_q) ** 2 * energy) / energy_sum))
    if delta_q <= 64.0 * np.finfo(float).eps * max(1.0, float(np.max(q))):
        return np.zeros(count), {
            "implementation": "modal-grid-green-laplace-eta20",
            "status": "degenerate-zero-difference-output",
            "quadrature_rank": rank,
            "interaction_enumeration": False,
        }

    nodes, weights = gauss_laguerre_rule(rank)
    signed_bins = np.concatenate((np.arange(0, span + 1), np.arange(-span, 0)))
    coefficients = np.zeros(count, dtype=np.complex128)
    for index, signed_bin in enumerate(signed_bins):
        sigma = signed_bin * frequency_step
        fd = forcing_d_spectrum[:, :, index]
        fk = forcing_k_spectrum[:, :, index]
        valid = (output_q > 0) & (response_a > abs(sigma))
        md = np.zeros_like(output_q)
        mk = np.zeros(output_q.shape, dtype=np.complex128)
        av = response_a[valid]
        for node, weight in zip(nodes, weights, strict=True):
            tau = node / delta_q
            plus = weight / delta_q * np.exp(node - (av - sigma) * tau)
            minus = weight / delta_q * np.exp(node - (av + sigma) * tau)
            md[valid] -= av * (plus + minus) / 4.0
            mk[valid] += 1j * (plus - minus) / 4.0
        response = np.zeros(output_q.shape, dtype=np.complex128)
        response[valid] = (md[valid] * fd[valid] + mk[valid] * fk[valid]) / (
            2.0 * depth
        )
        coefficients[signed_bin % count] = np.sum(response)
    coefficients = _hermitian_project(coefficients)
    return np.real(np.fft.fft(coefficients)), {
        "implementation": "modal-grid-green-laplace-eta20",
        "quadrature_rank": rank,
        "modal_grid": [ny, nx],
        "delta_q_rms_pair": delta_q,
        "strict_zero_spatial_mode": "excluded",
        "interaction_enumeration": False,
    }


def _joint_field(joint, multiplier, nx, ny):
    return np.fft.fft(
        np.fft.ifft(np.fft.ifft(joint * multiplier, axis=0), axis=1) * (nx * ny),
        axis=2,
    )


def _joint_spectrum(field, nx, ny):
    return np.fft.fft2(np.fft.ifft(field, axis=2), axes=(0, 1)) / (nx * ny)


def _modal_axis(count, parent, order, requested_limit=0.0):
    if requested_limit == 0.0:
        parent_max = float(np.max(np.abs(parent)))
        if parent_max == 0.0:
            raise ValueError("a nonzero modal-axis limit is required")
        requested_limit = (order + 0.05) * parent_max / (1.0 - 2.0 / count)
    dq = 2.0 * requested_limit / count
    return np.arange(-count / 2, count / 2) * dq, dq


def _deposit_joint(amplitudes, bins, qx, qy, qx_axis, qy_axis, dqx, dqy, count):
    ix0, ix1, wx0, wx1 = _bracket(qx, qx_axis, dqx)
    iy0, iy1, wy0, wy1 = _bracket(qy, qy_axis, dqy)
    joint = np.zeros((qy_axis.size, qx_axis.size, count), dtype=np.complex128)
    for ix, wx in ((ix0, wx0), (ix1, wx1)):
        for iy, wy in ((iy0, wy0), (iy1, wy1)):
            np.add.at(joint, (iy, ix, bins), amplitudes * wx * wy)
    return joint


def _bracket(value, axis, spacing):
    position = (value - axis[0]) / spacing
    lower = np.floor(position).astype(int)
    lower = np.clip(lower, 0, axis.size - 2)
    upper = lower + 1
    upper_weight = np.clip(position - lower, 0.0, 1.0)
    return lower, upper, 1.0 - upper_weight, upper_weight


def _hermitian_project(coefficients):
    result = coefficients.copy()
    count = result.size
    for bin_index in range(1, (count - 1) // 2 + 1):
        average = (result[bin_index] + np.conj(result[-bin_index])) / 2.0
        result[bin_index] = average
        result[-bin_index] = np.conj(average)
    result[0] = np.real(result[0])
    if count % 2 == 0:
        result[count // 2] = np.real(result[count // 2])
    return result

