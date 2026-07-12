"""Evaluate a frozen-charge object wrench and vertical release path."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from beach import AdhesionProfile, Beach, finite_shell_convergence


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--target-mesh-id", type=int, required=True)
    parser.add_argument(
        "--periodic-model",
        choices=("configured", "infinite-physical"),
        default="infinite-physical",
    )
    parser.add_argument("--z-max-m", type=float, required=True)
    parser.add_argument("--z-points", type=int, default=65)
    parser.add_argument("--mass-kg", type=float, required=True)
    parser.add_argument(
        "--gravity-m-s2",
        type=float,
        default=9.80665,
        help="gravity for release mechanics (default: Earth 9.80665 m/s^2)",
    )
    parser.add_argument("--adhesion-force-n", type=float, default=0.0)
    parser.add_argument("--adhesion-range-m", type=float, default=0.0)
    parser.add_argument("--eta-translation", type=float, default=1.0)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--library", type=Path)
    parser.add_argument(
        "--shell-max-layers",
        type=int,
        help=(
            "run corrected E_bottom=0 finite-shell convergence for a periodic2 "
            "snapshot"
        ),
    )
    return parser


def _adhesion(args: argparse.Namespace) -> AdhesionProfile:
    force = float(args.adhesion_force_n)
    extent = float(args.adhesion_range_m)
    if not math.isfinite(force) or not math.isfinite(extent):
        raise ValueError("adhesion force and range must be finite")
    if force == 0.0 and extent == 0.0:
        return AdhesionProfile.none()
    if force < 0.0 or extent <= 0.0:
        raise ValueError(
            "finite-range adhesion requires force >= 0 and range > 0"
        )
    return AdhesionProfile.finite_range_constant(
        force_N=force,
        range_m=extent,
    )


def _validate_args(args: argparse.Namespace) -> AdhesionProfile:
    for name in ("z_max_m", "mass_kg"):
        value = float(getattr(args, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if args.z_points < 2:
        raise ValueError("z_points must be >= 2")
    gravity = float(args.gravity_m_s2)
    if not math.isfinite(gravity) or gravity < 0.0:
        raise ValueError("gravity_m_s2 must be finite and non-negative")
    eta = float(args.eta_translation)
    if not math.isfinite(eta) or not 0.0 <= eta <= 1.0:
        raise ValueError("eta_translation must be finite and between 0 and 1")
    if args.shell_max_layers is not None:
        if args.shell_max_layers < 0:
            raise ValueError("shell_max_layers must be non-negative")
    return _adhesion(args)


def _finite_shell_summary(shells) -> dict[str, object]:
    def optional_list(values) -> list[object] | None:
        return None if values is None else values.tolist()

    def path_rows(paths) -> list[dict[str, object]]:
        rows = []
        for layer, path in zip(shells.image_layers, paths):
            rows.append(
                {
                    "image_layers": int(layer),
                    "path_status": path.status,
                    "endpoint_force_N": path.force_N[-1].tolist(),
                    "endpoint_torque_Nm": path.torque_Nm[-1].tolist(),
                    "endpoint_electrostatic_work_J": float(
                        path.electrostatic_work_J[-1]
                    ),
                    "work_relative_mismatch": path.work_relative_mismatch,
                }
            )
        return rows

    return {
        "status": shells.status,
        "selected_closure": "e_bottom_zero",
        "selected_image_layers": shells.selected_image_layers,
        "image_layers": shells.image_layers.tolist(),
        "increment_converged": shells.increment_converged.tolist(),
        "force_increment_error_N": shells.force_increment_error_N.tolist(),
        "work_increment_error_J": shells.work_increment_error_J.tolist(),
        "force_tail_proxy_N": shells.force_tail_proxy_N.tolist(),
        "work_tail_proxy_J": shells.work_tail_proxy_J.tolist(),
        "reference_model": shells.reference_model,
        "reference_force_error_N": optional_list(
            shells.reference_force_error_N
        ),
        "reference_work_error_J": optional_list(shells.reference_work_error_J),
        "reference_converged": optional_list(shells.reference_converged),
        "raw_symmetric": path_rows(shells.symmetric_paths),
        "corrected_e_bottom_zero": path_rows(shells.corrected_paths),
    }


def main() -> None:
    args = _parser().parse_args()
    adhesion = _validate_args(args)

    run = Beach(args.output_dir, config_path=args.config)
    displacement = np.linspace(0.0, args.z_max_m, args.z_points)
    shell_summary: dict[str, object] | None = None
    with run.object_interaction_snapshot(
        periodic_model=args.periodic_model.replace("-", "_"),
        cache_dir=args.cache_dir,
        library_path=args.library,
    ) as snapshot:
        probe = snapshot.object_probe(args.target_mesh_id)
        wrench = probe.wrench()
        path = probe.vertical_path(displacement)
        if args.shell_max_layers is not None:
            shells = finite_shell_convergence(
                snapshot,
                probe,
                path.displacement_m,
                max_layers=args.shell_max_layers,
            )
            shell_summary = _finite_shell_summary(shells)

    release = path.evaluate_release(
        mass_kg=args.mass_kg,
        gravity_m_s2=args.gravity_m_s2,
        adhesion=adhesion,
        eta_translation=args.eta_translation,
    )
    summary = {
        "force_N": wrench.force_N.tolist(),
        "torque_Nm": wrench.torque_Nm.tolist(),
        "path_status": path.status,
        "work_relative_mismatch": path.work_relative_mismatch,
        "barrier_free_from_rest": release.barrier_free_from_rest,
        "first_inaccessible_displacement_m": (
            release.first_inaccessible_displacement_m
        ),
        "endpoint_speed_m_s": release.endpoint_speed_m_s,
        "numerically_qualified": release.numerically_qualified,
        "finite_shell": shell_summary,
    }
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
