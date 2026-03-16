"""CLI for plotting object-wise Coulomb-force matrices."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from beach import Beach

from ._shared import configure_entry_parser

COMMAND_NAME = "coulomb"
LEGACY_COMMAND_NAME = "beach-plot-coulomb-force-matrix"


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
        "--component",
        choices=("x", "y", "z"),
        default="z",
        help="force component rendered in the matrix",
    )
    parser.add_argument(
        "--target-kinds",
        default=None,
        help="comma-separated template kinds used as targets (default: all objects)",
    )
    parser.add_argument(
        "--softening",
        type=float,
        default=0.0,
        help="softening length [m] for Coulomb-force evaluation",
    )
    parser.add_argument(
        "--cmap",
        default="coolwarm",
        help="matplotlib colormap name",
    )
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="output image path (default: <output_dir>/coulomb_force_<component>.png)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="saved image DPI",
    )
    parser.add_argument(
        "--no-annotate",
        action="store_true",
        help="do not draw numeric values in each matrix cell",
    )
    parser.add_argument("--show", action="store_true", help="display matplotlib window")


def build_parser(*, prog: str | None = LEGACY_COMMAND_NAME) -> argparse.ArgumentParser:
    """Build the argument parser for Coulomb-matrix plotting CLI."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command under the unified ``beachx`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="plot object-wise Coulomb force matrix",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def _parse_target_kinds(value: str | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    return tuple(part.strip().lower() for part in value.split(",") if part.strip())


def run(args: argparse.Namespace) -> None:
    """Execute the Coulomb-matrix plotting command."""

    parser = args._parser
    if args.dpi <= 0:
        parser.error("--dpi must be > 0.")
    if args.softening < 0.0:
        parser.error("--softening must be >= 0.")

    target_kinds = _parse_target_kinds(args.target_kinds)
    output_dir = Path(args.output_dir)
    save_path = (
        args.save
        if args.save is not None
        else output_dir / f"coulomb_force_{args.component}.png"
    )

    beach = Beach(output_dir)
    try:
        beach.result
    except FileNotFoundError as exc:
        raise SystemExit(
            f'Fortran output files are missing under "{output_dir}". '
            "Expected at least summary.txt and charges.csv."
        ) from exc

    try:
        fig, ax = beach.plot_coulomb_force_matrix(
            step=args.step,
            component=args.component,
            softening=args.softening,
            target_kinds=target_kinds,
            config_path=args.config,
            cmap=args.cmap,
            annotate=not args.no_annotate,
        )
        fig.savefig(save_path, dpi=args.dpi)
    except ModuleNotFoundError as exc:
        if exc.name is not None and exc.name.startswith("matplotlib"):
            raise SystemExit(
                "matplotlib is required for visualization. "
                "Install dependencies with `python -m pip install -e . --no-build-isolation`."
            ) from exc
        raise
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    metadata = getattr(ax, "_beach_coulomb_matrix", {})
    print(f"saved={save_path}")
    print(f"component={metadata.get('component', args.component)}")
    print(f"targets={list(metadata.get('target_labels', ())) }")
    print(f"sources={list(metadata.get('source_labels', ())) }")

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()
    else:
        fig.clf()


def main(argv: Sequence[str] | None = None) -> None:
    """Run the Coulomb-matrix plotting CLI entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
