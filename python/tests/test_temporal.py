import json
from pathlib import Path

import numpy as np

from green_laplace import reconstruct_unidirectional_timeseries


FIXTURE = Path(__file__).parent / "fixtures" / "unidirectional_time_eta22.json"
ORDER3_FIXTURE = Path(__file__).parent / "fixtures" / "unidirectional_time_order3.json"
ETA20_FIXTURE = Path(__file__).parent / "fixtures" / "unidirectional_time_eta20.json"
FIXED_PHASE_FIXTURE = (
    Path(__file__).parent / "fixtures" / "unidirectional_time_fixed_phase.json"
)


def test_unidirectional_time_matches_matlab_fft_gl():
    data = json.loads(FIXTURE.read_text(encoding="utf-8"))
    eta1 = np.asarray(data["eta1"], dtype=float)
    time = np.asarray(data["time"], dtype=float)
    period = len(time) * (time[1] - time[0])
    omega_max = 2 * np.pi * 8 / period + 1e-12
    result = reconstruct_unidirectional_timeseries(
        eta1,
        time,
        depth=data["depth"],
        gravity=data["gravity"],
        omega_max=omega_max,
        energy_fraction=1.0,
        quadrature_rank=data["rank"],
        domain_lengths=data["domain_lengths"],
        spatial_points=data["spatial_points"],
        relative_tolerance=1e-12,
        include_eta20=False,
    )
    expected_eta22 = np.asarray(data["eta22"], dtype=float)
    assert _relative(result.components["eta11"], eta1) < 2e-13
    assert _relative(result.components["eta22"], expected_eta22) < 3e-11


def test_unidirectional_order3_matches_matlab_fft_gl():
    data = json.loads(ORDER3_FIXTURE.read_text(encoding="utf-8"))
    eta1 = np.asarray(data["eta1"], dtype=float)
    time = np.asarray(data["time"], dtype=float)
    period = len(time) * (time[1] - time[0])
    result = reconstruct_unidirectional_timeseries(
        eta1,
        time,
        depth=data["depth"],
        gravity=data["gravity"],
        omega_max=2 * np.pi * 8 / period + 1e-12,
        energy_fraction=1.0,
        quadrature_rank=data["rank"],
        order=3,
        domain_lengths=data["domain_lengths"],
        spatial_points=data["spatial_points"],
        relative_tolerance=1e-12,
        include_eta20=False,
    )
    assert _relative(result.components["eta22"], data["eta22"]) < 3e-11
    assert _relative(result.components["eta33"], data["eta33"]) < 8e-11


def test_unidirectional_eta20_matches_matlab_fft_gl():
    data = json.loads(ETA20_FIXTURE.read_text(encoding="utf-8"))
    eta1 = np.asarray(data["eta1"], dtype=float)
    time = np.asarray(data["time"], dtype=float)
    period = len(time) * (time[1] - time[0])
    result = reconstruct_unidirectional_timeseries(
        eta1,
        time,
        depth=data["depth"],
        gravity=data["gravity"],
        omega_max=2 * np.pi * 8 / period + 1e-12,
        energy_fraction=1.0,
        order=2,
        domain_lengths=(25.0,),
        spatial_points=(64,),
        include_eta20=True,
        eta20_rank=data["rank"],
        eta20_grid=tuple(data["grid"]),
    )
    assert _relative(result.components["eta20"], data["eta20"]) < 8e-11


def test_unidirectional_fixed_phase_record_matches_matlab_fft_gl():
    data = json.loads(FIXED_PHASE_FIXTURE.read_text(encoding="utf-8"))
    result = reconstruct_unidirectional_timeseries(
        np.asarray(data["eta1"], dtype=float),
        np.asarray(data["time"], dtype=float),
        depth=data["depth"],
        gravity=data["gravity"],
        omega_max=data["omega_max"],
        energy_fraction=1.0,
        quadrature_rank=data["rank"],
        order=3,
        domain_lengths=data["domain_lengths"],
        spatial_points=data["spatial_points"],
        relative_tolerance=1e-12,
        include_eta20=True,
        eta20_rank=data["eta20_rank"],
        eta20_grid=tuple(data["eta20_grid"]),
    )
    for component in ("eta20", "eta22", "eta33"):
        assert _relative(result.components[component], data[component]) < 1e-10


def _relative(candidate, reference):
    candidate = np.asarray(candidate)
    reference = np.asarray(reference)
    return np.linalg.norm(candidate - reference) / np.linalg.norm(reference)
