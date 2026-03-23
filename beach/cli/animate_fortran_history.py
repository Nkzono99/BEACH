"""CLI for animating Fortran charge/potential history."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from beach import Beach

from ._shared import configure_entry_parser

COMMAND_NAME = "animate"
LEGACY_COMMAND_NAME = "beach-animate-history"


def _configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument(
        "--quantity",
        choices=("charge", "potential"),
        default="charge",
        help="value to visualize on mesh snapshots",
    )
    parser.add_argument(
        "--save-gif",
        type=Path,
        default=None,
        help="output GIF path (default: <output_dir>/<quantity>_history.gif)",
    )
    parser.add_argument("--fps", type=int, default=10, help="GIF frame rate")
    parser.add_argument(
        "--frame-stride",
        type=int,
        default=1,
        help="use every N-th history snapshot",
    )
    parser.add_argument(
        "--total-frames",
        type=int,
        default=None,
        help="target number of frames sampled evenly from history snapshots",
    )
    parser.add_argument(
        "--cmap",
        default=None,
        help="matplotlib colormap name (defaults: charge=coolwarm, potential=viridis)",
    )
    parser.add_argument(
        "--potential-softening",
        type=float,
        default=None,
        help="smoothing length [m] for potential reconstruction; default uses sim.softening when available",
    )
    parser.add_argument(
        "--potential-self-term",
        choices=("auto", "area-equivalent", "exclude", "softened-point"),
        default="auto",
        help="self-term treatment for potential reconstruction",
    )
    parser.add_argument(
        "--apply-periodic2-mesh",
        action="store_true",
        help=(
            "wrap each plotted triangle into the periodic2 cell from nearby "
            "beach.toml using its centroid"
        ),
    )


def build_parser(*, prog: str | None = LEGACY_COMMAND_NAME) -> argparse.ArgumentParser:
    """Build the argument parser for history animation CLI."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command under the unified ``beachx`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="animate charge or potential history",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def run(args: argparse.Namespace) -> None:
    """Execute the history-animation command."""

    parser = args._parser
    if args.total_frames is not None and args.frame_stride != 1:
        parser.error("--frame-stride and --total-frames cannot be used together.")

    output_dir = Path(args.output_dir)
    save_gif = (
        args.save_gif
        if args.save_gif is not None
        else output_dir / f"{args.quantity}_history.gif"
    )
    self_term = args.potential_self_term.replace("-", "_")
    beach = Beach(output_dir)
    try:
        result = beach.result
    except FileNotFoundError as exc:
        raise SystemExit(
            f'Fortran output files are missing under "{output_dir}". '
            "Expected at least summary.txt and charges.csv."
        ) from exc

    try:
        written = beach.animate_mesh(
            output_path=save_gif,
            quantity=args.quantity,
            fps=args.fps,
            frame_stride=args.frame_stride,
            total_frames=args.total_frames,
            cmap=args.cmap,
            softening=args.potential_softening,
            self_term=self_term,
            apply_periodic2_mesh=args.apply_periodic2_mesh,
        )
    except ModuleNotFoundError as exc:
        if exc.name is not None and (
            exc.name.startswith("matplotlib") or exc.name.startswith("PIL")
        ):
            raise SystemExit(
                "matplotlib and pillow are required for GIF animation. "
                "Install dependencies with `python -m pip install -e . --no-build-isolation`."
            ) from exc
        raise

    print(f"saved_gif={written}")
    print(f"quantity={args.quantity}")
    snapshot_count = (
        len(result.history)
        if result.history is not None and result.history.has_data
        else 0
    )
    if args.total_frames is None:
        rendered_frames = (
            (snapshot_count + args.frame_stride - 1) // args.frame_stride
            if snapshot_count > 0
            else 0
        )
    else:
        rendered_frames = min(snapshot_count, args.total_frames)
    print(f"snapshots={snapshot_count}")
    print(f"rendered_frames={rendered_frames}")


def main(argv: Sequence[str] | None = None) -> None:
    """Run the history-animation CLI entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
