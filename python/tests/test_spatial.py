import json
from pathlib import Path

import numpy as np

from green_laplace import reconstruct_directional_field


FIXTURE = Path(__file__).parent / "fixtures" / "directional_spatial_eta22.json"
ORDER3_FIXTURE = Path(__file__).parent / "fixtures" / "directional_spatial_order3.json"
NONZERO_TIME_FIXTURE = (
    Path(__file__).parent / "fixtures" / "directional_spatial_nonzero_time.json"
)


def test_directional_spatial_matches_matlab_fft_gl():
    data = json.loads(FIXTURE.read_text(encoding="utf-8"))
    eta1 = np.asarray(data["eta1"], dtype=float)
    result = reconstruct_directional_field(
        eta1,
        lx=data["lx"],
        ly=data["ly"],
        depth=data["depth"],
        gravity=data["gravity"],
        order=2,
        eta22_rank=data["rank"],
    )
    expected_eta22 = np.asarray(data["eta22"], dtype=float)
    expected_psi11 = np.asarray(data["psi11"], dtype=float)
    expected_psi22 = np.asarray(data["psi22"], dtype=float)
    assert _relative(result.components["eta11"], eta1) < 2e-13
    assert _relative(result.components["eta22"], expected_eta22) < 2e-12
    assert _relative(result.components["psi11"], expected_psi11) < 2e-13
    assert _relative(result.components["psi22"], expected_psi22) < 3e-12


def test_directional_order3_matches_matlab_fft_gl():
    data = json.loads(ORDER3_FIXTURE.read_text(encoding="utf-8"))
    result = reconstruct_directional_field(
        np.asarray(data["eta1"], dtype=float),
        lx=data["lx"],
        ly=data["ly"],
        depth=data["depth"],
        gravity=data["gravity"],
        order=3,
        eta22_rank=6,
    )
    assert _relative(result.components["eta33"], data["eta33"]) < 5e-11
    assert _relative(result.components["psi33"], data["psi33"]) < 5e-11


def test_directional_nonzero_time_matches_matlab_fft_gl():
    data = json.loads(NONZERO_TIME_FIXTURE.read_text(encoding="utf-8"))
    result = reconstruct_directional_field(
        np.asarray(data["eta1"], dtype=float),
        lx=data["lx"],
        ly=data["ly"],
        depth=data["depth"],
        gravity=data["gravity"],
        time=data["time"],
        order=3,
        eta22_rank=data["rank"],
    )
    for component in ("eta11", "psi11", "eta22", "psi22", "eta33", "psi33"):
        assert _relative(result.components[component], data[component]) < 8e-11


def _relative(candidate, reference):
    candidate = np.asarray(candidate)
    reference = np.asarray(reference)
    return np.linalg.norm(candidate - reference) / np.linalg.norm(reference)
