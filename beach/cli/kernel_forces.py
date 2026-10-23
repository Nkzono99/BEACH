"""CLI for object-wise forces from the Fortran field kernel."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from beach import Beach

from . import _options
from ._output import emit


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``kernel-forces`` options."""

    _options.add_output_dir(parser)
    _options.add_config(parser)
    _options.add_step(parser)
    parser.add_argument(
        "--target-mesh-ids",
        default=None,
        help="comma-separated mesh ids to evaluate (default: all objects)",
    )
    parser.add_argument(
        "--theta",
        type=_options.positive_float,
        default=None,
        help="FMM theta (default: read sim.tree_theta or 0.5)",
    )
    parser.add_argument(
        "--leaf-max",
        type=_options.positive_int,
        default=None,
        help="FMM leaf_max (default: read sim.tree_leaf_max or 16)",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=None,
        help="FMM Cartesian expansion order (default: receipt, or 4 for legacy output)",
    )
    _options.add_library(parser)
    _options.add_save(parser, default="<output_dir>/object_forces_kernel.csv")
    _options.add_json(parser)


def run(args: argparse.Namespace, beach: Beach) -> None:
    """Execute the kernel-force command."""

    save_path = (
        args.save
        if args.save is not None
        else args.output_dir / "object_forces_kernel.csv"
    )
    records = beach.calc_object_forces_kernel(
        step=args.step,
        target_mesh_ids=_parse_mesh_ids(args.target_mesh_ids),
        theta=args.theta,
        leaf_max=args.leaf_max,
        order=args.order,
        library_path=args.library,
    )

    _write_csv(save_path, records)
    emit(
        args,
        {
            "saved": save_path,
            "objects": [
                {
                    "mesh_id": record.mesh_id,
                    "step": record.step,
                    "charge_C": record.total_charge_C,
                    "force_norm_N": float(np.linalg.norm(record.force_N)),
                    "torque_norm_Nm": float(np.linalg.norm(record.torque_Nm)),
                }
                for record in records
            ],
        },
    )


def _parse_mesh_ids(value: str | None) -> tuple[int, ...] | None:
    if value is None:
        return None
    mesh_ids = tuple(int(part.strip()) for part in value.split(",") if part.strip())
    return mesh_ids or None


def _write_csv(path: Path, records) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "mesh_id",
                "step",
                "total_charge_C",
                "center_x_m",
                "center_y_m",
                "center_z_m",
                "force_x_N",
                "force_y_N",
                "force_z_N",
                "torque_x_Nm",
                "torque_y_Nm",
                "torque_z_Nm",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "mesh_id": record.mesh_id,
                    "step": record.step,
                    "total_charge_C": record.total_charge_C,
                    "center_x_m": record.center_m[0],
                    "center_y_m": record.center_m[1],
                    "center_z_m": record.center_m[2],
                    "force_x_N": record.force_N[0],
                    "force_y_N": record.force_N[1],
                    "force_z_N": record.force_N[2],
                    "torque_x_Nm": record.torque_Nm[0],
                    "torque_y_Nm": record.torque_Nm[1],
                    "torque_z_Nm": record.torque_Nm[2],
                }
            )

