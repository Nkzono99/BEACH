"""CLI for plotting object-wise Coulomb-force matrices."""

from __future__ import annotations

import argparse

from beach import Beach

from . import _options


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
        library_path=args.library,
        cmap=args.cmap,
        annotate=not args.no_annotate,
        workers=args.workers,
    )
    fig.savefig(save_path, dpi=args.dpi)

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
