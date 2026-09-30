"""Compare modal and auxiliary-domain eta33 on frozen random-wave HOS cases."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from green_laplace import reconstruct_unidirectional_modal_timeseries


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--modal-points", type=int, default=8192)
    parser.add_argument("--modal-ladder", type=str, default="")
    parser.add_argument("--relative-tolerance", type=float, default=1e-3)
    args = parser.parse_args()
    cases = json.loads(args.input.read_text(encoding="utf-8"))["cases"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    names = [
        "seed",
        "nominal_Akp",
        "actual_kp_Hs_over_2",
        "modal_points",
        "modal_seconds",
        "auxiliary_seconds",
        "speedup",
        "modal_HOS_relative",
        "auxiliary_HOS_relative",
        "modal_auxiliary_relative",
        "modal_HOS_cosine",
        "modal_HOS_norm_ratio",
        "modal_q_projection_rms",
        "modal_last_change",
        "modal_converged",
        "auxiliary_last_domain_change",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        for case in cases:
            eta1 = np.asarray(case["first"], dtype=float)
            dt = float(case["dt"])
            time_axis = np.arange(eta1.size) * dt
            wp = np.sqrt(
                case["g"] * case["kp"] * np.tanh(case["kp"] * case["h"])
            )
            started = time.perf_counter()
            modal_points = (
                tuple(int(value) for value in args.modal_ladder.split(","))
                if args.modal_ladder
                else args.modal_points
            )
            result = reconstruct_unidirectional_modal_timeseries(
                eta1,
                time_axis,
                depth=case["h"],
                gravity=case["g"],
                omega_max=4 * wp,
                energy_fraction=1.0,
                quadrature_rank=8,
                order=3,
                modal_points=modal_points,
                relative_tolerance=args.relative_tolerance,
                include_eta20=False,
            )
            elapsed = time.perf_counter() - started
            score = (time_axis >= 15 * case["Tp"] - 1e-8) & (
                time_axis <= 35 * case["Tp"] + 1e-8
            )
            modal = np.asarray(result.components["eta33"])[score]
            auxiliary = np.asarray(case["auxiliary_GL_eta33"], dtype=float)
            hos = np.asarray(case["HOS_third"], dtype=float)
            auxiliary_seconds = float(case["auxiliary_seconds"])
            final_changes = [
                value
                for value in (
                    result.audit["levels"][-1]["eta22_change"],
                    result.audit["levels"][-1]["eta33_change"],
                )
                if value is not None
            ]
            row = {
                "seed": case["seed"],
                "nominal_Akp": case["Akp"],
                "actual_kp_Hs_over_2": case["actual_steepness"],
                "modal_points": result.audit["modal_points"],
                "modal_seconds": elapsed,
                "auxiliary_seconds": auxiliary_seconds,
                "speedup": auxiliary_seconds / elapsed,
                "modal_HOS_relative": relative(modal, hos),
                "auxiliary_HOS_relative": relative(auxiliary, hos),
                "modal_auxiliary_relative": relative(modal, auxiliary),
                "modal_HOS_cosine": float(
                    np.dot(modal, hos) / (np.linalg.norm(modal) * np.linalg.norm(hos))
                ),
                "modal_HOS_norm_ratio": float(
                    np.linalg.norm(modal) / np.linalg.norm(hos)
                ),
                "modal_q_projection_rms": result.audit[
                    "wavevector_projection_rms_relative"
                ],
                "modal_last_change": max(final_changes) if final_changes else "",
                "modal_converged": result.audit["converged"],
                "auxiliary_last_domain_change": case["auxiliary_levels"][-1][
                    "relative_change"
                ],
            }
            writer.writerow(row)
            stream.flush()
            print(
                f"seed={case['seed']} Akp={case['Akp']:.2f} "
                f"modal={elapsed:.3f}s HOS={100*row['modal_HOS_relative']:.3f}%"
            )


def relative(candidate: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(candidate - reference) / np.linalg.norm(reference))


if __name__ == "__main__":
    main()
