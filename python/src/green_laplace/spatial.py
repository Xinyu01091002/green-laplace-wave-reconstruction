"""Directional spatial Green--Laplace reconstruction from ``eta1(x, y)``."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._spectral import (
    analytic_spectrum_from_real_field,
    assert_quadratic_support_is_alias_safe,
    fft_wavenumber_grid,
    suppress_fft_roundoff,
)
from .quadrature import gauss_laguerre_rule
from ._order3 import green_laplace_order3


@dataclass(frozen=True)
class DirectionalFieldResult:
    """Fields and execution information returned by a spatial reconstruction."""

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    eta: NDArray[np.float64]
    psi: NDArray[np.float64]
    components: dict[str, NDArray[np.float64] | NDArray[np.complex128]]
    audit: dict[str, Any]


def reconstruct_directional_field(
    eta1: ArrayLike,
    *,
    lx: float,
    ly: float,
    depth: float,
    gravity: float = 9.81,
    time: float = 0.0,
    order: int = 2,
    eta22_rank: int = 6,
    allow_below_model_parent_domain: bool = False,
) -> DirectionalFieldResult:
    """Reconstruct a strict-forward directional field from a real ``eta1`` grid.

    ``eta1`` is the physical first-order elevation at the time origin.  The
    current implementation supports the linear field and the positive
    sum-frequency second-order elevation.  Surface ``psi11`` is returned with
    the same definition used by the MATLAB implementation.
    """
    field = np.asarray(eta1, dtype=float)
    if field.ndim != 2:
        raise ValueError("eta1 must be a two-dimensional array")
    if not np.isfinite(time) or depth <= 0 or gravity <= 0:
        raise ValueError("time must be finite; depth and gravity must be positive")
    if order not in (1, 2, 3):
        raise ValueError("order must be 1, 2, or 3")
    if not isinstance(eta22_rank, int) or eta22_rank < 1 or eta22_rank > 16:
        raise ValueError("eta22_rank must be an integer from 1 to 16")

    ny, nx = field.shape
    kx, ky = fft_wavenumber_grid(lx, ly, nx, ny)
    initial_spectrum, input_error = analytic_spectrum_from_real_field(
        field, depth=depth, kx=kx
    )
    qx = depth * kx
    qy = depth * ky
    q = np.hypot(qx, qy)
    nu = np.sqrt(q * np.tanh(q))
    phase = np.exp(-1j * np.sqrt(gravity / depth) * nu * time)
    spectrum = initial_spectrum * phase

    eta11_plus = depth * np.fft.ifft2(spectrum)
    psi11_spectrum = np.zeros_like(spectrum)
    nonzero = nu > 0
    psi11_spectrum[nonzero] = -1j * spectrum[nonzero] / nu[nonzero]
    potential_scale = depth * np.sqrt(gravity * depth)
    psi11_plus = potential_scale * np.fft.ifft2(psi11_spectrum)

    components: dict[str, NDArray[np.float64] | NDArray[np.complex128]] = {
        "eta11_plus": eta11_plus,
        "psi11_plus": psi11_plus,
        "eta11": np.real(eta11_plus),
        "psi11": np.real(psi11_plus),
    }
    eta_plus = eta11_plus.copy()
    psi_plus = psi11_plus.copy()
    component_audits: dict[str, Any] = {}

    if order >= 2:
        eta22_dimensionless, eta22_audit = _green_laplace_eta22(
            spectrum, qx, qy, eta22_rank
        )
        eta22_plus = depth * eta22_dimensionless
        psi22_dimensionless, psi22_audit = _green_laplace_psi22(
            spectrum,
            qx,
            qy,
            allow_below_model_parent_domain=allow_below_model_parent_domain,
        )
        psi22_plus = potential_scale * psi22_dimensionless
        components["eta22_plus"] = eta22_plus
        components["eta22"] = np.real(eta22_plus)
        components["psi22_plus"] = psi22_plus
        components["psi22"] = np.real(psi22_plus)
        eta_plus += eta22_plus
        psi_plus += psi22_plus
        component_audits["eta22"] = eta22_audit
        component_audits["psi22"] = psi22_audit

    if order >= 3:
        eta33_dimensionless, psi33_dimensionless, order3_audit = green_laplace_order3(
            spectrum, qx, qy
        )
        eta33_plus = depth * eta33_dimensionless
        psi33_plus = potential_scale * psi33_dimensionless
        components["eta33_plus"] = eta33_plus
        components["psi33_plus"] = psi33_plus
        components["eta33"] = np.real(eta33_plus)
        components["psi33"] = np.real(psi33_plus)
        eta_plus += eta33_plus
        psi_plus += psi33_plus
        component_audits["order3"] = order3_audit

    x_axis = np.arange(nx, dtype=float) * lx / nx
    y_axis = np.arange(ny, dtype=float) * ly / ny
    x, y = np.meshgrid(x_axis, y_axis)
    components["eta_plus"] = eta_plus
    components["psi_plus"] = psi_plus
    components["eta"] = np.real(eta_plus)
    components["psi"] = np.real(psi_plus)

    audit = {
        "implementation": "green-laplace-python-directional-spatial",
        "order": order,
        "grid": [ny, nx],
        "domain": [ly, lx],
        "time": float(time),
        "input": "real eta1(x,y) strict-forward in x",
        "input_reconstruction_relative": input_error,
        "phase_convention": "exp(i*(kx*x+ky*y)-i*omega*t)",
        "stokes_correction": False,
        "component_audits": component_audits,
    }
    return DirectionalFieldResult(
        x=x,
        y=y,
        eta=np.real(eta_plus),
        psi=np.real(psi_plus),
        components=components,
        audit=audit,
    )


def _green_laplace_eta22(
    eta11_spectrum: NDArray[np.complex128],
    qx: NDArray[np.float64],
    qy: NDArray[np.float64],
    rank: int,
) -> tuple[NDArray[np.complex128], dict[str, Any]]:
    q = np.hypot(qx, qy)
    threshold = 1.0e-13 * max(1.0, float(np.max(np.abs(eta11_spectrum))))
    support = np.abs(eta11_spectrum) > threshold
    if np.any(qx[support] <= 0) or np.any(q[support] <= 0):
        raise ValueError("strict-forward nonzero analytic eta1 support is required")
    assert_quadratic_support_is_alias_safe(support)

    support_count = np.real(np.fft.ifft2(np.fft.fft2(support.astype(float)) ** 2))
    output_support = support_count > 0.5
    nu = np.sqrt(q * np.tanh(q))
    response_a = q * np.tanh(q)
    sqrt_a = np.sqrt(response_a)
    safe_nu = np.maximum(nu, np.finfo(float).tiny)

    support_amplitudes = np.abs(eta11_spectrum[support]) ** 2
    peak_flat = np.flatnonzero(support)[int(np.argmax(support_amplitudes))]
    peak_q = float(q.flat[peak_flat])
    scale = np.sqrt(
        (2.0 * np.sqrt(peak_q * np.tanh(peak_q))) ** 2
        - 2.0 * peak_q * np.tanh(2.0 * peak_q)
    )
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("invalid determinant-root Green--Laplace scale")

    nodes, weights = gauss_laguerre_rule(rank)
    minimum_parent_nu = float(np.min(nu[support]))
    maximum_pair_a = float(np.max(sqrt_a[output_support]))
    balance_shift = 0.5 * (minimum_parent_nu + maximum_pair_a / 2.0)
    candidate = np.zeros_like(eta11_spectrum)

    for node, weight in zip(nodes, weights, strict=True):
        tau = node / scale
        damping = np.exp(-tau * (nu - balance_shift))
        damped = eta11_spectrum * damping
        v = np.fft.ifft2(damped)
        v_nu = np.fft.ifft2(damped * nu)
        v_nu2 = np.fft.ifft2(damped * nu**2)
        hx = np.fft.ifft2(damped * qx / safe_nu)
        hy = np.fft.ifft2(damped * qy / safe_nu)
        radial_over_nu = np.fft.ifft2(damped * q**2 / safe_nu)
        jx = np.fft.ifft2(damped * qx)
        jy = np.fft.ifft2(damped * qy)
        source_d = 2.0 * v * v_nu2 + v_nu**2 - hx**2 - hy**2
        source_k = 2.0 * v * radial_over_nu + 2.0 * (hx * jx + hy * jy)
        coefficient = (
            weight * np.exp(node) / scale / 4.0 * np.exp(-2.0 * balance_shift * tau)
        )
        sinh_over_a = np.full_like(sqrt_a, tau)
        nonzero = sqrt_a > 0
        sinh_over_a[nonzero] = np.sinh(sqrt_a[nonzero] * tau) / sqrt_a[nonzero]
        source_d_spectrum = np.fft.fft2(source_d)
        source_k_spectrum = np.fft.fft2(source_k)
        source_d_spectrum[~output_support] = 0.0
        source_k_spectrum[~output_support] = 0.0
        source_d_spectrum = suppress_fft_roundoff(source_d_spectrum, output_support)
        source_k_spectrum = suppress_fft_roundoff(source_k_spectrum, output_support)
        candidate += coefficient * (
            -response_a * sinh_over_a * source_d_spectrum
            + np.cosh(sqrt_a * tau) * source_k_spectrum
        )

    return np.fft.ifft2(candidate), {
        "implementation": "prescribed-rank-green-laplace-eta22",
        "quadrature_rank": rank,
        "quadrature_exact_moment_degree": 2 * rank - 1,
        "peak_depth_wavenumber": peak_q,
        "fft_ifft_count": 10 * rank + 5,
        "pointwise_product_count": 7 * rank + 1,
        "interaction_enumeration": False,
    }


def _green_laplace_psi22(
    eta11_spectrum: NDArray[np.complex128],
    qx: NDArray[np.float64],
    qy: NDArray[np.float64],
    *,
    allow_below_model_parent_domain: bool,
) -> tuple[NDArray[np.complex128], dict[str, Any]]:
    q = np.hypot(qx, qy)
    nu = np.sqrt(q * np.tanh(q))
    safe_nu = np.maximum(nu, np.finfo(float).tiny)
    threshold = 1.0e-13 * max(1.0, float(np.max(np.abs(eta11_spectrum))))
    support = np.abs(eta11_spectrum) > threshold
    if not np.any(support) or np.any(qx[support] <= 0) or np.any(q[support] <= 0):
        raise ValueError("strict-forward nonzero analytic eta1 support is required")
    below_domain = bool(np.any(q[support] < 0.3))
    if below_domain and not allow_below_model_parent_domain:
        raise ValueError("psi22 requires parent kh >= 0.3")
    below_mass = float(
        np.sum(np.abs(eta11_spectrum[q < 0.3]) ** 2)
        / max(np.sum(np.abs(eta11_spectrum) ** 2), np.finfo(float).tiny)
    )
    assert_quadratic_support_is_alias_safe(support)
    output_support = np.real(
        np.fft.ifft2(np.fft.fft2(support.astype(float)) ** 2)
    ) > 0.5

    masks, shifts, pair_supports = _radial_balance_channels(support, nu, 2)
    sqrt_a = nu
    safe_a = safe_nu
    support_amplitudes = np.abs(eta11_spectrum[support]) ** 2
    peak_flat = np.flatnonzero(support)[int(np.argmax(support_amplitudes))]
    peak_q = float(q.flat[peak_flat])
    peak_nu = np.sqrt(peak_q * np.tanh(peak_q))
    peak_output_nu = np.sqrt(2.0 * peak_q * np.tanh(2.0 * peak_q))
    branch_scales = (2.0 * peak_nu - peak_output_nu, 2.0 * peak_nu + peak_output_nu)
    nodes = np.array([2.0 - np.sqrt(2.0), 2.0 + np.sqrt(2.0)])
    weights = np.array(
        [
            nodes[0] / (4.0 * (np.sqrt(2.0) - 1.0) ** 2),
            nodes[1] / (4.0 * (np.sqrt(2.0) + 1.0) ** 2),
        ]
    )
    resolvent = np.zeros_like(eta11_spectrum)
    for branch_name, branch_scale, source_k_sign in zip(
        ("slow", "fast"), branch_scales, (-1.0, 1.0), strict=True
    ):
        if branch_scale <= 0 or not np.isfinite(branch_scale):
            raise ValueError("invalid dual-branch Green--Laplace scale")
        for node, weight in zip(nodes, weights, strict=True):
            tau = node / branch_scale
            coefficient = weight * np.exp(node) / branch_scale / 8.0
            response = _radial_balanced_branch_response(
                eta11_spectrum,
                qx,
                qy,
                q,
                nu,
                safe_nu,
                sqrt_a,
                safe_a,
                tau,
                branch_name,
                source_k_sign,
                masks,
                shifts,
                pair_supports,
            )
            resolvent += 1j * coefficient * response

    eta11 = np.fft.ifft2(eta11_spectrum)
    phi1z = -1j * np.fft.ifft2(eta11_spectrum * nu)
    taylor_contact = 0.5 * eta11 * phi1z
    psi_spectrum = resolvent + np.fft.fft2(taylor_contact)
    psi_spectrum[~output_support] = 0.0
    psi22 = np.fft.ifft2(psi_spectrum)
    return psi22, {
        "implementation": "dual-branch-green-laplace-psi22",
        "quadrature_rank_total": 4,
        "quadrature_rank_per_branch": 2,
        "lambda_minus": float(branch_scales[0]),
        "lambda_plus": float(branch_scales[1]),
        "radial_balance_channel_count": len(masks),
        "below_model_parent_domain": below_domain,
        "below_model_spectral_mass": below_mass,
        "interaction_enumeration": False,
    }


def _radial_balance_channels(
    support: NDArray[np.bool_], nu: NDArray[np.float64], count: int
) -> tuple[list[NDArray[np.bool_]], list[float], list[list[NDArray[np.bool_]]]]:
    active = nu[support]
    edges = np.linspace(float(np.min(active)), float(np.max(active)), count + 1)
    masks: list[NDArray[np.bool_]] = []
    shifts: list[float] = []
    for channel in range(count):
        if channel < count - 1:
            mask = support & (nu >= edges[channel]) & (nu < edges[channel + 1])
        else:
            mask = support & (nu >= edges[channel]) & (nu <= edges[channel + 1])
        if np.any(mask):
            values = nu[mask]
            masks.append(mask)
            shifts.append(0.5 * (float(np.min(values)) + float(np.max(values))))
    transforms = [np.fft.fft2(mask.astype(float)) for mask in masks]
    pair_supports = [
        [
            np.real(np.fft.ifft2(transforms[first] * transforms[second])) > 0.5
            for second in range(len(masks))
        ]
        for first in range(len(masks))
    ]
    return masks, shifts, pair_supports


def _radial_balanced_branch_response(
    spectrum: NDArray[np.complex128],
    qx: NDArray[np.float64],
    qy: NDArray[np.float64],
    q: NDArray[np.float64],
    nu: NDArray[np.float64],
    safe_nu: NDArray[np.float64],
    sqrt_a: NDArray[np.float64],
    safe_a: NDArray[np.float64],
    tau: float,
    branch_name: str,
    source_k_sign: float,
    masks: list[NDArray[np.bool_]],
    shifts: list[float],
    pair_supports: list[list[NDArray[np.bool_]]],
) -> NDArray[np.complex128]:
    fields = []
    for mask, shift in zip(masks, shifts, strict=True):
        channel_spectrum = spectrum * mask * np.exp(-tau * (nu - shift))
        fields.append(_psi22_source_fields(channel_spectrum, qx, qy, q, nu, safe_nu))
    response = np.zeros_like(spectrum)
    for first, left in enumerate(fields):
        for second, right in enumerate(fields):
            source_d = (
                2.0 * left["v"] * right["v_nu2"]
                + left["v_nu"] * right["v_nu"]
                - left["hx"] * right["hx"]
                - left["hy"] * right["hy"]
            )
            source_k = 2.0 * left["v"] * right["radial"] + 2.0 * (
                left["hx"] * right["jx"] + left["hy"] * right["jy"]
            )
            source_d_spectrum = np.fft.fft2(source_d)
            source_k_spectrum = np.fft.fft2(source_k)
            pair_support = pair_supports[first][second]
            shift_sum = shifts[first] + shifts[second]
            if branch_name == "slow":
                output_factor = np.exp((sqrt_a - shift_sum) * tau)
            else:
                output_factor = np.exp(-(sqrt_a + shift_sum) * tau)
            pair_response = output_factor * (
                source_d_spectrum + source_k_sign * source_k_spectrum / safe_a
            )
            pair_response[~pair_support] = 0.0
            response += pair_response
    return response


def _psi22_source_fields(
    spectrum: NDArray[np.complex128],
    qx: NDArray[np.float64],
    qy: NDArray[np.float64],
    q: NDArray[np.float64],
    nu: NDArray[np.float64],
    safe_nu: NDArray[np.float64],
) -> dict[str, NDArray[np.complex128]]:
    return {
        "v": np.fft.ifft2(spectrum),
        "v_nu": np.fft.ifft2(spectrum * nu),
        "v_nu2": np.fft.ifft2(spectrum * nu**2),
        "hx": np.fft.ifft2(spectrum * qx / safe_nu),
        "hy": np.fft.ifft2(spectrum * qy / safe_nu),
        "radial": np.fft.ifft2(spectrum * q**2 / safe_nu),
        "jx": np.fft.ifft2(spectrum * qx),
        "jy": np.fft.ifft2(spectrum * qy),
    }
