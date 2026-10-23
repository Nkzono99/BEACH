"""CLI for animating Fortran charge/potential history."""

from __future__ import annotations

import argparse

from beach import Beach

from . import _options
from ._output import emit


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``animate`` options."""

    _options.add_output_dir(parser)
    parser.add_argument(
        "--quantity",
        choices=("charge", "potential"),
        default="charge",
        help="value to visualize on mesh snapshots",
    )
    _options.add_save(parser, default="<output_dir>/<quantity>_history.gif")
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
    _options.add_cmap(parser, "charge", "potential")
    _options.add_library(parser)
    _options.add_periodic2_mesh(parser, mesh_mode=False)
    _options.add_view(parser)


def run(args: argparse.Namespace, beach: Beach) -> None:
    """Execute the history-animation command."""

    save_path = (
        args.save
        if args.save is not None
        else args.output_dir / f"{args.quantity}_history.gif"
    )
    written = beach.animate_mesh(
        output_path=save_path,
        quantity=args.quantity,
        fps=args.fps,
        frame_stride=args.frame_stride,
        total_frames=args.total_frames,
        cmap=args.cmap,
        apply_periodic2_mesh=args.apply_periodic2_mesh,
        periodic2_repeat=args.periodic2_repeat,
        view_elev=args.view_elev,
        view_azim=args.view_azim,
    )

    result = beach.result
    snapshot_count = (
        len(result.history.charge)
        if result.history.charge is not None and result.history.charge.has_data
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
    emit(
        args,
        {
            "saved": written,
            "quantity": args.quantity,
            "snapshots": snapshot_count,
            "rendered_frames": rendered_frames,
        },
    )
