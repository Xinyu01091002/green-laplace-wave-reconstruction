"""No-Stokes third-order directional spatial Green--Laplace graph."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from ._spectral import suppress_fft_roundoff
from .quadrature import gauss_laguerre_rule


Complex = NDArray[np.complex128]
Real = NDArray[np.float64]


def green_laplace_order3(
    spectrum: Complex,
    qx: Real,
    qy: Real,
    *,
    inner_rank: int = 4,
    outer_rank: int = 4,
) -> tuple[Complex, Complex, dict[str, Any]]:
    if inner_rank not in (4, 6, 8, 10, 12) or outer_rank not in (3, 4, 6, 8, 10):
        raise ValueError("unsupported inner or outer third-order GL rank")
    q = np.hypot(qx, qy)
    nu = np.sqrt(q * np.tanh(q))
    response_a = q * np.tanh(q)
    sqrt_a = np.sqrt(response_a)
    safe_nu = np.maximum(nu, np.finfo(float).tiny)
    threshold = 1.0e-13 * max(1.0, float(np.max(np.abs(spectrum))))
    support = np.abs(spectrum) > threshold
    if not np.any(support) or np.any(qx[support] <= 0) or np.any(q[support] <= 0.5):
        raise ValueError("third-order reconstruction requires strict-forward parent kh > 0.5")
    _assert_cubic_support_is_alias_safe(support)
    support_transform = np.fft.fft2(support.astype(float))
    pair_support = np.real(np.fft.ifft2(support_transform**2)) > 0.5
    triple_support = np.real(np.fft.ifft2(support_transform**3)) > 0.5

    minimum_parent_nu = float(np.min(nu[support]))
    inner_shift = 0.5 * (minimum_parent_nu + float(np.max(sqrt_a[pair_support])) / 2.0)
    outer_shift = 0.5 * (minimum_parent_nu + float(np.max(sqrt_a[triple_support])) / 3.0)
    support_amplitudes = np.abs(spectrum[support]) ** 2
    peak_flat = np.flatnonzero(support)[int(np.argmax(support_amplitudes))]
    peak_q = float(q.flat[peak_flat])
    peak_nu = np.sqrt(peak_q * np.tanh(peak_q))
    lambda2 = 2.0 * peak_nu - np.sqrt(2.0 * peak_q * np.tanh(2.0 * peak_q))
    lambda3 = 3.0 * peak_nu - np.sqrt(3.0 * peak_q * np.tanh(3.0 * peak_q))
    if lambda2 <= 0 or lambda3 <= 0:
        raise ValueError("invalid third-order Green--Laplace scales")
    inner_nodes, inner_weights = gauss_laguerre_rule(inner_rank)
    outer_nodes, outer_weights = gauss_laguerre_rule(outer_rank)
    multipliers = _multipliers(qx, qy, q, nu, response_a, safe_nu)

    eta_spectrum = np.zeros_like(spectrum)
    psi_spectrum = np.zeros_like(spectrum)
    for node, weight in zip(outer_nodes, outer_weights, strict=True):
        tau = node / lambda3
        damped = spectrum * np.exp(-tau * (nu - outer_shift))
        first = _first_order_state(damped, multipliers)
        pair = _inner_pair_fields(
            damped,
            q,
            response_a,
            sqrt_a,
            multipliers,
            pair_support,
            inner_nodes,
            inner_weights,
            lambda2,
            inner_shift,
        )

        forcing_k_pair = (
            first["phix"] * pair["etax"]
            + first["phiy"] * pair["etay"]
            + pair["phix"] * first["etax"]
            + pair["phiy"] * first["etay"]
            - first["eta"] * pair["phizz"]
            - pair["eta"] * first["phizz"]
        )
        forcing_k_direct = first["eta"] * (
            first["phixz"] * first["etax"] + first["phiyz"] * first["etay"]
        ) - 0.5 * first["eta"] ** 2 * first["phizzz"]
        forcing_d_pair = (
            first["eta"] * pair["phitz"]
            + pair["eta"] * first["phitz"]
            + first["phix"] * pair["phix"]
            + first["phiy"] * pair["phiy"]
            + first["phiz"] * pair["phiz"]
        )
        forcing_d_direct = 0.5 * first["eta"] ** 2 * first["phitzz"] + first[
            "eta"
        ] * (
            first["phix"] * first["phixz"]
            + first["phiy"] * first["phiyz"]
            + first["phiz"] * first["phizz"]
        )
        forcing_k = forcing_k_direct / 4.0 + forcing_k_pair / 2.0
        forcing_d = forcing_d_direct / 4.0 + forcing_d_pair / 2.0
        forcing_k_spectrum = np.fft.fft2(forcing_k)
        forcing_d_spectrum = np.fft.fft2(forcing_d)
        forcing_k_spectrum[~triple_support] = 0.0
        forcing_d_spectrum[~triple_support] = 0.0
        forcing_k_spectrum = suppress_fft_roundoff(forcing_k_spectrum, triple_support)
        forcing_d_spectrum = suppress_fft_roundoff(forcing_d_spectrum, triple_support)

        contact, contact_t = _surface_taylor_contact(first, pair)
        contact_spectrum = np.fft.fft2(contact)
        contact_t_spectrum = np.fft.fft2(contact_t)
        contact_spectrum[~triple_support] = 0.0
        contact_t_spectrum[~triple_support] = 0.0
        contact_spectrum = suppress_fft_roundoff(contact_spectrum, triple_support)
        contact_t_spectrum = suppress_fft_roundoff(contact_t_spectrum, triple_support)
        forcing_k_psi = suppress_fft_roundoff(
            forcing_k_spectrum + response_a * contact_spectrum, triple_support
        )
        forcing_d_psi = suppress_fft_roundoff(
            forcing_d_spectrum - contact_t_spectrum, triple_support
        )

        sinh_over_a = np.full_like(sqrt_a, tau)
        nonzero = sqrt_a > 0
        sinh_over_a[nonzero] = np.sinh(sqrt_a[nonzero] * tau) / sqrt_a[nonzero]
        cosh_factor = np.cosh(sqrt_a * tau)
        coefficient = weight * np.exp(node) / lambda3 * np.exp(-3.0 * outer_shift * tau)
        eta_spectrum += coefficient * (
            response_a * sinh_over_a * forcing_d_spectrum
            - 1j * cosh_factor * forcing_k_spectrum
        )
        psi_spectrum += coefficient * (
            -1j * cosh_factor * forcing_d_psi - sinh_over_a * forcing_k_psi
        )

    eta_spectrum[~triple_support] = 0.0
    psi_spectrum[~triple_support] = 0.0
    return np.fft.ifft2(eta_spectrum), np.fft.ifft2(psi_spectrum), {
        "implementation": "nested-no-stokes-green-laplace-order3",
        "inner_quadrature_rank": inner_rank,
        "outer_quadrature_rank": outer_rank,
        "inner_scale": float(lambda2),
        "outer_scale": float(lambda3),
        "inner_balance_shift": inner_shift,
        "outer_balance_shift": outer_shift,
        "interaction_enumeration": False,
    }


def _multipliers(qx: Real, qy: Real, q: Real, nu: Real, a: Real, safe_nu: Real):
    q2 = q**2
    return {
        "qx": qx,
        "qy": qy,
        "q2": q2,
        "nu": nu,
        "nu2": nu**2,
        "nu3": nu**3,
        "dx": 1j * qx,
        "dy": 1j * qy,
        "qx_over_nu": qx / safe_nu,
        "qy_over_nu": qy / safe_nu,
        "q2_over_nu": q2 / safe_nu,
        "minus_i_nu": -1j * nu,
        "nu_qx": nu * qx,
        "nu_qy": nu * qy,
        "minus_i_q2_over_nu": -1j * q2 / safe_nu,
        "minus_i_q2_nu": -1j * q2 * nu,
        "minus_a": -a,
        "minus_q2": -q2,
    }


def _first_order_state(spectrum: Complex, m) -> dict[str, Complex]:
    state = {
        "eta": np.fft.ifft2(spectrum),
        "etax": np.fft.ifft2(m["dx"] * spectrum),
        "etay": np.fft.ifft2(m["dy"] * spectrum),
        "phix": np.fft.ifft2(m["qx_over_nu"] * spectrum),
        "phiy": np.fft.ifft2(m["qy_over_nu"] * spectrum),
        "phiz": np.fft.ifft2(m["minus_i_nu"] * spectrum),
        "phixz": np.fft.ifft2(m["nu_qx"] * spectrum),
        "phiyz": np.fft.ifft2(m["nu_qy"] * spectrum),
        "phizz": np.fft.ifft2(m["minus_i_q2_over_nu"] * spectrum),
        "phizzz": np.fft.ifft2(m["minus_i_q2_nu"] * spectrum),
        "phitz": np.fft.ifft2(m["minus_a"] * spectrum),
        "phitzz": np.fft.ifft2(m["minus_q2"] * spectrum),
    }
    state["eta_t"] = state["phiz"]
    state["phiz_t"] = state["phitz"]
    state["phizz_t"] = state["phitzz"]
    return state


def _inner_pair_fields(
    spectrum: Complex,
    q: Real,
    response_a: Real,
    sqrt_a: Real,
    m,
    output_support: NDArray[np.bool_],
    nodes: Real,
    weights: Real,
    scale: float,
    shift: float,
) -> dict[str, Complex]:
    eta_spectrum = np.zeros_like(spectrum)
    eta_t_spectrum = np.zeros_like(spectrum)
    phi_spectrum = np.zeros_like(spectrum)
    phi_t_spectrum = np.zeros_like(spectrum)
    nonzero = q > 0
    for node, weight in zip(nodes, weights, strict=True):
        tau = node / scale
        damped = spectrum * np.exp(-tau * (m["nu"] - shift))
        v = np.fft.ifft2(damped)
        v_nu = np.fft.ifft2(damped * m["nu"])
        v_nu2 = np.fft.ifft2(damped * m["nu2"])
        v_nu3 = np.fft.ifft2(damped * m["nu3"])
        hx = np.fft.ifft2(damped * m["qx_over_nu"])
        hy = np.fft.ifft2(damped * m["qy_over_nu"])
        jx = np.fft.ifft2(damped * m["qx"])
        jy = np.fft.ifft2(damped * m["qy"])
        radial = np.fft.ifft2(damped * m["q2_over_nu"])
        q2_field = np.fft.ifft2(damped * m["q2"])
        nu_jx = np.fft.ifft2(damped * m["nu_qx"])
        nu_jy = np.fft.ifft2(damped * m["nu_qy"])
        source_d = 2.0 * v * v_nu2 + v_nu**2 - hx**2 - hy**2
        source_k = 2.0 * v * radial + 2.0 * (hx * jx + hy * jy)
        source_d_t = 2.0 * v * v_nu3 + 4.0 * v_nu * v_nu2 - 2.0 * hx * jx - 2.0 * hy * jy
        source_k_t = 2.0 * (v_nu * radial + v * q2_field) + 2.0 * (
            jx**2 + hx * nu_jx + jy**2 + hy * nu_jy
        )
        sd, sk, sd_t, sk_t = (
            np.fft.fft2(source_d),
            np.fft.fft2(source_k),
            np.fft.fft2(source_d_t),
            np.fft.fft2(source_k_t),
        )
        for item in (sd, sk, sd_t, sk_t):
            item[~output_support] = 0.0
        sd = suppress_fft_roundoff(sd, output_support)
        sk = suppress_fft_roundoff(sk, output_support)
        sd_t = suppress_fft_roundoff(sd_t, output_support)
        sk_t = suppress_fft_roundoff(sk_t, output_support)
        sinh_over_a = np.full_like(sqrt_a, tau)
        sinh_over_a[nonzero] = np.sinh(sqrt_a[nonzero] * tau) / sqrt_a[nonzero]
        cosh_factor = np.cosh(sqrt_a * tau)
        coefficient = weight * np.exp(node) / scale * np.exp(-2.0 * shift * tau) / 4.0
        eta_spectrum += coefficient * (-response_a * sinh_over_a * sd + cosh_factor * sk)
        eta_t_spectrum -= 1j * coefficient * (
            -response_a * sinh_over_a * sd_t + cosh_factor * sk_t
        )
        phi_spectrum += 1j * coefficient * (cosh_factor * sd - sinh_over_a * sk)
        phi_t_spectrum += coefficient * (cosh_factor * sd_t - sinh_over_a * sk_t)

    for item in (eta_spectrum, eta_t_spectrum, phi_spectrum, phi_t_spectrum):
        item[~output_support] = 0.0
    eta_t_spectrum = suppress_fft_roundoff(eta_t_spectrum, output_support)
    pair = {
        "eta": np.fft.ifft2(eta_spectrum),
        "etax": np.fft.ifft2(m["dx"] * eta_spectrum),
        "etay": np.fft.ifft2(m["dy"] * eta_spectrum),
        "phix": np.fft.ifft2(m["dx"] * phi_spectrum),
        "phiy": np.fft.ifft2(m["dy"] * phi_spectrum),
        "phiz": np.fft.ifft2(response_a * phi_spectrum),
        "phizz": np.fft.ifft2(m["q2"] * phi_spectrum),
        "phitz": np.fft.ifft2(response_a * phi_t_spectrum),
        "eta_t": np.fft.ifft2(eta_t_spectrum),
    }
    pair["phiz_t"] = pair["phitz"]
    return pair


def _surface_taylor_contact(first, pair) -> tuple[Complex, Complex]:
    contact = (
        0.5 * first["eta"] * pair["phiz"]
        + 0.5 * pair["eta"] * first["phiz"]
        + 0.125 * first["eta"] ** 2 * first["phizz"]
    )
    contact_t = (
        0.5 * (first["eta_t"] * pair["phiz"] + first["eta"] * pair["phiz_t"])
        + 0.5 * (pair["eta_t"] * first["phiz"] + pair["eta"] * first["phiz_t"])
        + 0.25 * first["eta"] * first["eta_t"] * first["phizz"]
        + 0.125 * first["eta"] ** 2 * first["phizz_t"]
    )
    return contact, contact_t


def _assert_cubic_support_is_alias_safe(support: NDArray[np.bool_]) -> None:
    ny, nx = support.shape
    mx, my = np.meshgrid(np.fft.fftfreq(nx) * nx, np.fft.fftfreq(ny) * ny)
    if np.max(3.0 * np.abs(mx[support])) >= nx / 2.0 or np.max(
        3.0 * np.abs(my[support])
    ) >= ny / 2.0:
        raise ValueError("cubic output support reaches a Nyquist boundary")

