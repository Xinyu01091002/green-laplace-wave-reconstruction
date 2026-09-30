"""Public single-direction time-series interface."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from numpy.typing import ArrayLike

from .modal_time import (
    ModalTimeSeriesResult,
    reconstruct_unidirectional_modal_timeseries,
)

TimeSeriesResult = ModalTimeSeriesResult


def reconstruct_unidirectional_timeseries(
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
    **unknown: Any,
) -> TimeSeriesResult:
    """Reconstruct a single-direction record with the modal GL resolvent.

    The input is the separated real first-order elevation ``eta1(t)``. The
    implementation uses a one-dimensional dimensionless-wavenumber modal grid
    and does not create an auxiliary spatial domain.
    """
    if unknown:
        names = ", ".join(sorted(unknown))
        raise TypeError(f"unknown time-reconstruction options: {names}")
    return reconstruct_unidirectional_modal_timeseries(
        eta1,
        time,
        depth=depth,
        gravity=gravity,
        omega_max=omega_max,
        energy_fraction=energy_fraction,
        quadrature_rank=quadrature_rank,
        order=order,
        modal_points=modal_points,
        modal_limit=modal_limit,
        relative_tolerance=relative_tolerance,
        include_eta20=include_eta20,
        eta20_rank=eta20_rank,
        eta20_grid=eta20_grid,
    )
