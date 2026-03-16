"""CLI for Coulomb mobility analysis."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from beach import Beach

from ._shared import configure_entry_parser

COMMAND_NAME = "mobility"


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
        default=0.0,
        help="softening length [m] for Coulomb-force evaluation",
    )
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
    parser.add_argument(
        "--save-csv",
        type=Path,
        default=None,
        help="output CSV path (default: <output_dir>/mobility_summary.csv)",
    )


def build_parser(*, prog: str | None = COMMAND_NAME) -> argparse.ArgumentParser:
    """Build the argument parser for Coulomb mobility CLI."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command under the unified ``beachx`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="analyze object-wise Coulomb mobility indicators",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def _parse_kind_list(value: str | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    return tuple(part.strip().lower() for part in value.split(",") if part.strip())


def _fmt_scalar(value: float | None) -> str:
    if value is None:
        return "n/a"
    if value == float("inf"):
        return "inf"
    return f"{value:.3e}"


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


def run(args: argparse.Namespace) -> None:
    """Execute the Coulomb mobility command."""

    parser = args._parser
    if args.softening < 0.0:
        parser.error("--softening must be >= 0.")
    if args.adhesion_force_N < 0.0:
        parser.error("--adhesion-force-N must be >= 0.")

    output_dir = Path(args.output_dir)
    save_csv = (
        args.save_csv
        if args.save_csv is not None
        else output_dir / "mobility_summary.csv"
    )
    support_kinds = _parse_kind_list(args.support_kinds)
    target_kinds = _parse_kind_list(args.target_kinds)

    beach = Beach(output_dir)
    try:
        beach.result
    except FileNotFoundError as exc:
        raise SystemExit(
            f'Fortran output files are missing under "{output_dir}". '
            "Expected at least summary.txt and charges.csv."
        ) from exc

    try:
        analysis = beach.analyze_coulomb_mobility(
            step=args.step,
            softening=args.softening,
            config_path=args.config,
            gravity=args.gravity,
            support_normal=args.support_normal,
            support_kinds=support_kinds,
            target_kinds=target_kinds,
            density_kg_m3=args.density_kg_m3,
            mu_static=args.mu_static,
            mu_roll=args.mu_roll,
            adhesion_force_N=args.adhesion_force_N,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    _write_csv(save_csv, analysis)

    print(f"saved_csv={save_csv}")
    print(f"support_normal={analysis.support_normal_m.tolist()}")
    for record in analysis.records:
        force_norm = float((record.force_N @ record.force_N) ** 0.5)
        torque_norm = float((record.torque_Nm @ record.torque_Nm) ** 0.5)
        print(
            "object="
            f"{record.label} "
            f"mesh_id={record.mesh_id} "
            f"kind={record.kind} "
            f"force_norm_N={force_norm:.3e} "
            f"torque_norm_Nm={torque_norm:.3e} "
            f"lift_ratio={_fmt_scalar(record.lift_ratio)} "
            f"slide_ratio={_fmt_scalar(record.slide_ratio)} "
            f"roll_ratio={_fmt_scalar(record.roll_ratio)}"
        )
        if record.notes:
            print(f"notes[{record.label}]={'; '.join(record.notes)}")


def main(argv: Sequence[str] | None = None) -> None:
    """Run the Coulomb mobility CLI entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
