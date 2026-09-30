import numpy as np

from green_laplace import (
    reconstruct_directional_field,
    reconstruct_unidirectional_modal_timeseries,
)


def test_on_grid_modal_time_matches_spatial_gl():
    q1 = _commensurate_q()
    q2 = 2 * q1
    gravity = 9.81
    depth = 1.0
    omega1 = np.sqrt(gravity / depth) * _nu(q1)
    omega2 = np.sqrt(gravity / depth) * _nu(q2)
    sample_count = 256
    period = 4 * np.pi / omega1
    time = np.arange(sample_count) * period / sample_count
    amplitudes = np.array([0.01 * np.exp(0.2j), 0.006 * np.exp(-0.4j)])
    eta1_time = np.real(
        amplitudes[0] * np.exp(-1j * omega1 * time)
        + amplitudes[1] * np.exp(-1j * omega2 * time)
    )

    nx, ny = 128, 16
    lx, ly = 2 * np.pi / q1, 2 * np.pi
    x = np.arange(nx) * lx / nx
    y = np.arange(ny) * ly / ny
    X, _ = np.meshgrid(x, y)
    eta1_space = np.real(
        amplitudes[0] * np.exp(1j * q1 * X)
        + amplitudes[1] * np.exp(1j * q2 * X)
    )
    spatial = reconstruct_directional_field(
        eta1_space,
        lx=lx,
        ly=ly,
        depth=depth,
        gravity=gravity,
        order=3,
        eta22_rank=4,
    )
    modal_points = 256
    modal = reconstruct_unidirectional_modal_timeseries(
        eta1_time,
        time,
        depth=depth,
        gravity=gravity,
        omega_max=omega2 + 1e-12,
        energy_fraction=1.0,
        quadrature_rank=4,
        order=3,
        modal_points=(128, modal_points),
        modal_limit=modal_points * (q1 / 4) / 2,
        include_eta20=False,
    )
    for component in ("eta11", "eta22", "eta33"):
        assert np.isclose(
            modal.components[component][0],
            spatial.components[component][0, 0],
            rtol=3e-14,
            atol=1e-16,
        )
    assert modal.audit["wavevector_projection_rms_relative"] < 1e-14
    assert modal.audit["converged"]
    assert modal.audit["modal_points"] == modal_points
    assert not modal.audit["auxiliary_space_domain"]


def test_off_grid_modal_ladder_converges():
    sample_count = 512
    dt = 0.2
    time = np.arange(sample_count) * dt
    bins = np.array([20, 22, 24])
    amplitudes = np.array([0.01, 0.006, 0.004])
    phases = np.array([0.2, -0.7, 1.1])
    eta1 = sum(
        amplitude
        * np.cos(2 * np.pi * bin_index * np.arange(sample_count) / sample_count + phase)
        for amplitude, bin_index, phase in zip(amplitudes, bins, phases, strict=True)
    )
    result = reconstruct_unidirectional_modal_timeseries(
        eta1,
        time,
        depth=1.0,
        omega_max=2 * np.pi * bins[-1] / (sample_count * dt) + 1e-12,
        energy_fraction=1.0,
        quadrature_rank=4,
        order=3,
        modal_points=(256, 512, 1024, 2048),
        relative_tolerance=1e-2,
        include_eta20=False,
    )
    assert result.audit["converged"]
    assert result.audit["modal_points"] == 1024
    assert result.audit["levels"][-1]["eta22_change"] < 1e-2
    assert result.audit["levels"][-1]["eta33_change"] < 1e-2
    projection = [
        level["wavevector_projection_rms_relative"]
        for level in result.audit["levels"]
    ]
    assert all(right < left for left, right in zip(projection, projection[1:]))


def _nu(q):
    return np.sqrt(q * np.tanh(q))


def _commensurate_q():
    lower, upper = 1e-6, 20.0
    for _ in range(100):
        middle = 0.5 * (lower + upper)
        if _nu(2 * middle) / _nu(middle) > 1.5:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)
