import json
from pathlib import Path

import numpy as np

from green_laplace import reconstruct_unidirectional_timeseries


FIXTURES = Path(__file__).parent / "fixtures"


def test_default_time_interface_uses_modal_resolvent():
    data = json.loads(
        (FIXTURES / "unidirectional_modal_order3.json").read_text(encoding="utf-8")
    )
    eta1 = np.asarray(data["eta1"], dtype=float)
    time = np.asarray(data["time"], dtype=float)
    result = reconstruct_unidirectional_timeseries(
        eta1,
        time,
        depth=data["depth"],
        gravity=data["gravity"],
        omega_max=data["omega_max"],
        energy_fraction=1.0,
        quadrature_rank=data["rank"],
        order=3,
        modal_points=data["modal_points"],
        modal_limit=data["modal_limit"],
        include_eta20=False,
    )
    assert _relative(result.components["eta22"], data["eta22"]) < 3e-13
    assert _relative(result.components["eta33"], data["eta33"]) < 3e-13
    assert not result.audit["auxiliary_space_domain"]


def test_unidirectional_eta20_matches_matlab_modal_gl():
    data = json.loads(
        (FIXTURES / "unidirectional_time_eta20.json").read_text(encoding="utf-8")
    )
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
        modal_points=512,
        include_eta20=True,
        eta20_rank=data["rank"],
        eta20_grid=tuple(data["grid"]),
    )
    assert _relative(result.components["eta20"], data["eta20"]) < 8e-11


def _relative(candidate, reference):
    candidate = np.asarray(candidate)
    reference = np.asarray(reference)
    return np.linalg.norm(candidate - reference) / np.linalg.norm(reference)
