"""CLI for object-wise forces from the Fortran field kernel."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

import numpy as np

from beach import Beach, FieldKernelError

from ._shared import (
    configure_entry_parser,
    require_nonnegative_finite,
    require_positive_finite,
)

COMMAND_NAME = "kernel-forces"


def _configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="path to beach.toml (default: auto-detect near output_dir)",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=-1,
        help="history batch step used for charges (-1: latest history)",
    )
    parser.add_argument(
        "--softening",
        type=float,
        default=None,
        help="softening length [m] (default: read sim.softening near output_dir)",
    )
    parser.add_argument(
        "--target-mesh-ids",
        default=None,
        help="comma-separated mesh ids to evaluate (default: all objects)",
    )
    parser.add_argument(
        "--theta",
        type=float,
        default=None,
        help="FMM theta (default: read sim.tree_theta or 0.5)",
    )
    parser.add_argument(
        "--leaf-max",
        type=int,
        default=None,
        help="FMM leaf_max (default: read sim.tree_leaf_max or 16)",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=4,
        help="FMM Cartesian expansion order",
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=None,
        help="path to libbeach_field_kernel shared library",
    )
    parser.add_argument(
        "--save-csv",
        type=Path,
        default=None,
        help="output CSV path (default: <output_dir>/object_forces_kernel.csv)",
    )


def build_parser(*, prog: str | None = COMMAND_NAME) -> argparse.ArgumentParser:
    """Build the parser for the kernel-force CLI."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command under the unified ``beachx`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="write object-wise net forces from the Fortran field kernel",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def run(args: argparse.Namespace) -> None:
    """Execute the kernel-force command."""

    parser = args._parser
    require_nonnegative_finite(parser, args.softening, "--softening")
    require_positive_finite(parser, args.theta, "--theta")
    if args.leaf_max is not None and args.leaf_max <= 0:
        parser.error("--leaf-max must be > 0.")
    if args.order < 0:
        parser.error("--order must be >= 0.")

    output_dir = Path(args.output_dir)
    save_csv = (
        args.save_csv
        if args.save_csv is not None
        else output_dir / "object_forces_kernel.csv"
    )
    target_mesh_ids = _parse_mesh_ids(args.target_mesh_ids)

    beach = Beach(output_dir)
    try:
        beach.result
    except FileNotFoundError as exc:
        raise SystemExit(
            f'Fortran output files are missing under "{output_dir}". '
            "Expected at least summary.txt and charges.csv."
        ) from exc

    try:
        records = beach.calc_object_forces_kernel(
            step=args.step,
            target_mesh_ids=target_mesh_ids,
            softening=args.softening,
            theta=args.theta,
            leaf_max=args.leaf_max,
            order=args.order,
            config_path=args.config,
            library_path=args.library,
        )
    except (ValueError, FieldKernelError) as exc:
        raise SystemExit(str(exc)) from exc

    _write_csv(save_csv, records)
    print(f"saved_csv={save_csv}")
    for record in records:
        force_norm = float(np.linalg.norm(record.force_N))
        torque_norm = float(np.linalg.norm(record.torque_Nm))
        print(
            f"mesh_id={record.mesh_id} "
            f"charge_C={record.total_charge_C:.6e} "
            f"force_norm_N={force_norm:.6e} "
            f"torque_norm_Nm={torque_norm:.6e}"
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


def main(argv: Sequence[str] | None = None) -> None:
    """Run the kernel-force CLI entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
