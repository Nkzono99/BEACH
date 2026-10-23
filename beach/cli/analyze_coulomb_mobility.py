"""CLI for Coulomb mobility analysis."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from beach import Beach

from . import _options
from ._output import emit


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``mobility`` options."""

    _options.add_output_dir(parser)
    _options.add_config(parser)
    _options.add_step(parser)
    _options.add_library(parser)
    parser.add_argument(
        "--gravity",
        nargs=3,
        type=float,
        metavar=("GX", "GY", "GZ"),
        default=(0.0, 0.0, -9.81),
        help="gravity vector [m/s^2]",
    )
    parser.add_argument(
        "--support-normal",
        nargs=3,
        type=float,
        metavar=("NX", "NY", "NZ"),
        default=None,
        help="support normal (default: opposite of gravity, else +z)",
    )
    parser.add_argument(
        "--support-kinds",
        default="plane",
        help="comma-separated object kinds treated as fixed supports (default: plane)",
    )
    parser.add_argument(
        "--target-kinds",
        default=None,
        help="comma-separated object kinds to analyze (default: all non-support objects)",
    )
    parser.add_argument(
        "--density-kg-m3",
        type=float,
        default=None,
        help="bulk density used for weight estimation of volumetric templates",
    )
    parser.add_argument(
        "--mu-static",
        type=float,
        default=None,
        help="static-friction coefficient for slide-ratio estimation",
    )
    parser.add_argument(
        "--mu-roll",
        type=float,
        default=None,
        help="rolling-resistance coefficient for roll-ratio estimation",
    )
    parser.add_argument(
        "--adhesion-force-N",
        type=float,
        default=0.0,
        help="additional normal adhesion resisting lift [N]",
    )
    _options.add_save(parser, default="<output_dir>/mobility_summary.csv")
    _options.add_json(parser)


def _parse_kind_list(value: str | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    return tuple(part.strip().lower() for part in value.split(",") if part.strip())


def _write_csv(path: Path, analysis) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "mesh_id",
                "label",
                "kind",
                "center_x_m",
                "center_y_m",
                "center_z_m",
                "force_x_N",
                "force_y_N",
                "force_z_N",
                "force_normal_N",
                "force_tangent_N",
                "torque_x_Nm",
                "torque_y_Nm",
                "torque_z_Nm",
                "torque_normal_Nm",
                "torque_tangent_Nm",
                "mass_kg",
                "characteristic_radius_m",
                "weight_support_N",
                "resisting_normal_N",
                "effective_normal_load_N",
                "lift_ratio",
                "slide_ratio",
                "roll_ratio",
                "notes",
            ],
        )
        writer.writeheader()
        for record in analysis.records:
            writer.writerow(
                {
                    "mesh_id": record.mesh_id,
                    "label": record.label,
                    "kind": record.kind,
                    "center_x_m": record.center_m[0],
                    "center_y_m": record.center_m[1],
                    "center_z_m": record.center_m[2],
                    "force_x_N": record.force_N[0],
                    "force_y_N": record.force_N[1],
                    "force_z_N": record.force_N[2],
                    "force_normal_N": record.force_normal_N,
                    "force_tangent_N": record.force_tangent_N,
                    "torque_x_Nm": record.torque_Nm[0],
                    "torque_y_Nm": record.torque_Nm[1],
                    "torque_z_Nm": record.torque_Nm[2],
                    "torque_normal_Nm": record.torque_normal_Nm,
                    "torque_tangent_Nm": record.torque_tangent_Nm,
                    "mass_kg": record.mass_kg,
                    "characteristic_radius_m": record.characteristic_radius_m,
                    "weight_support_N": record.weight_support_N,
                    "resisting_normal_N": record.resisting_normal_N,
                    "effective_normal_load_N": record.effective_normal_load_N,
                    "lift_ratio": record.lift_ratio,
                    "slide_ratio": record.slide_ratio,
                    "roll_ratio": record.roll_ratio,
                    "notes": "; ".join(record.notes),
                }
            )


def run(args: argparse.Namespace, beach: Beach) -> None:
    """Execute the Coulomb mobility command."""

    save_path = (
        args.save
        if args.save is not None
        else args.output_dir / "mobility_summary.csv"
    )
    analysis = beach.analyze_coulomb_mobility(
        step=args.step,
        library_path=args.library,
        gravity=args.gravity,
        support_normal=args.support_normal,
        support_kinds=_parse_kind_list(args.support_kinds),
        target_kinds=_parse_kind_list(args.target_kinds),
        density_kg_m3=args.density_kg_m3,
        mu_static=args.mu_static,
        mu_roll=args.mu_roll,
        adhesion_force_N=args.adhesion_force_N,
    )

    _write_csv(save_path, analysis)
    emit(
        args,
        {
            "saved": save_path,
            "support_normal": analysis.support_normal_m,
            "objects": [
                {
                    "mesh_id": record.mesh_id,
                    "label": record.label,
                    "kind": record.kind,
                    "force_norm_N": float(np.linalg.norm(record.force_N)),
                    "torque_norm_Nm": float(np.linalg.norm(record.torque_Nm)),
                    "lift_ratio": record.lift_ratio,
                    "slide_ratio": record.slide_ratio,
                    "roll_ratio": record.roll_ratio,
                    "notes": list(record.notes),
                }
                for record in analysis.records
            ],
        },
    )
