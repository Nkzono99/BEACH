"""CLI for generating close-packed sphere model configs."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from beach.closepack_config import (
    ClosePackSpec,
    build_closepack_config,
    load_base_config,
    render_closepack_toml,
    total_sphere_count,
)

from ._shared import configure_entry_parser

COMMAND_NAME = "close-pack"


def _configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "base_config",
        nargs="?",
        default=Path("beach.toml"),
        type=Path,
        help="base beach.toml used for sim/particles/output defaults (default: ./beach.toml)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("beach.close-pack.toml"),
        help="destination TOML path (default: ./beach.close-pack.toml)",
    )
    parser.add_argument("--layers", type=int, required=True, help="number of sphere layers")
    parser.add_argument("--radius", type=float, required=True, help="sphere radius [m]")
    parser.add_argument(
        "--cells-x",
        type=int,
        required=True,
        help="number of close-pack unit cells repeated in x",
    )
    parser.add_argument(
        "--cells-y",
        type=int,
        help="number of close-pack unit cells repeated in y (default: same as --cells-x)",
    )
    parser.add_argument(
        "--floor-z",
        type=float,
        default=0.02,
        help="floor plane z position [m]",
    )
    parser.add_argument(
        "--plane-nx",
        type=int,
        default=20,
        help="plane mesh subdivisions in x",
    )
    parser.add_argument(
        "--plane-ny",
        type=int,
        default=20,
        help="plane mesh subdivisions in y",
    )
    parser.add_argument(
        "--sphere-n-lon",
        type=int,
        default=14,
        help="sphere mesh longitude count",
    )
    parser.add_argument(
        "--sphere-n-lat",
        type=int,
        default=6,
        help="sphere mesh latitude count",
    )
    box_group = parser.add_mutually_exclusive_group()
    box_group.add_argument(
        "--box-height",
        type=float,
        help="explicit box height in z from sim.box_min[2] [m]",
    )
    box_group.add_argument(
        "--top-clearance",
        type=float,
        help="extra headroom above the top sphere [m]",
    )
    parser.add_argument(
        "--output-dir",
        help="override [output].dir in the generated TOML",
    )


def build_parser(*, prog: str | None = None) -> argparse.ArgumentParser:
    """Build the argument parser for close-packed model generation."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this generator under a parent ``model`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="generate a close-packed sphere-stack beach.toml",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def run(args: argparse.Namespace) -> None:
    """Execute close-packed model generation."""

    cells_y = args.cells_y if args.cells_y is not None else args.cells_x
    spec = ClosePackSpec(
        layers=args.layers,
        radius=args.radius,
        cells_x=args.cells_x,
        cells_y=cells_y,
        floor_z=args.floor_z,
        plane_nx=args.plane_nx,
        plane_ny=args.plane_ny,
        sphere_n_lon=args.sphere_n_lon,
        sphere_n_lat=args.sphere_n_lat,
        box_height=args.box_height,
        top_clearance=args.top_clearance,
    )

    try:
        base_config = load_base_config(args.base_config)
    except FileNotFoundError as exc:
        raise SystemExit(f"base config file not found: {args.base_config}") from exc
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    try:
        config = build_closepack_config(
            spec,
            base_config=base_config,
            output_dir=args.output_dir,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    text = render_closepack_toml(
        config,
        spec=spec,
        base_config_path=args.base_config,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")

    print(f"saved={args.output}")
    print(f"base_config={args.base_config}")
    print(f"sphere_count={total_sphere_count(spec)}")
    print(f"box_min={config['sim']['box_min']}")
    print(f"box_max={config['sim']['box_max']}")


def main(argv: Sequence[str] | None = None) -> None:
    """Run the close-packed model generator CLI."""

    args = build_parser(prog="beachx model close-pack").parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
