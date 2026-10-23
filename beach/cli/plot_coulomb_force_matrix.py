"""CLI for plotting object-wise Coulomb-force matrices."""

from __future__ import annotations

import argparse

from beach import Beach

from . import _options
from ._output import emit


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``coulomb`` options."""

    _options.add_output_dir(parser)
    _options.add_config(parser)
    _options.add_step(parser)
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
    _options.add_library(parser)
    _options.add_cmap(parser, "force")
    _options.add_save(parser, default="<output_dir>/coulomb_force_<component>.png")
    _options.add_dpi(parser)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="worker processes for Coulomb-force entries (default: 1)",
    )
    parser.add_argument(
        "--no-annotate",
        action="store_true",
        help="do not draw numeric values in each matrix cell",
    )
    _options.add_show(parser)
    _options.add_json(parser)


def _parse_target_kinds(value: str | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    return tuple(part.strip().lower() for part in value.split(",") if part.strip())


def run(args: argparse.Namespace, beach: Beach) -> None:
    """Execute the Coulomb-matrix plotting command."""

    save_path = (
        args.save
        if args.save is not None
        else args.output_dir / f"coulomb_force_{args.component}.png"
    )
    fig, ax = beach.plot_coulomb_force_matrix(
        step=args.step,
        component=args.component,
        target_kinds=_parse_target_kinds(args.target_kinds),
        cmap=args.cmap,
        annotate=not args.no_annotate,
        workers=args.workers,
    )
    fig.savefig(save_path, dpi=args.dpi)

    metadata = ax._beach_coulomb_matrix
    emit(
        args,
        {
            "saved": save_path,
            "component": metadata["component"],
            "step": metadata["step"],
            "targets": list(metadata["target_labels"]),
            "sources": list(metadata["source_labels"]),
            "force_matrix": [
                {"target": label, f"force_{args.component}_N": row}
                for label, row in zip(metadata["target_labels"], metadata["matrix"])
            ],
        },
    )

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()
    else:
        fig.clf()
