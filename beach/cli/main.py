"""Unified post-processing CLI for BEACH."""

from __future__ import annotations

import argparse
from typing import Sequence

from . import (
    analyze_coulomb_mobility,
    animate_fortran_history,
    estimate_fortran_workload,
    inspect_fortran_output,
    model,
    plot_coulomb_force_matrix,
    plot_fortran_potential_slices,
    plot_performance_profile,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the root parser for ``beachx``."""

    parser = argparse.ArgumentParser(
        prog="beachx",
        description="BEACH tools for post-processing, visualization, and model generation.",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        metavar="command",
        required=True,
    )
    inspect_fortran_output.add_subparser(subparsers)
    animate_fortran_history.add_subparser(subparsers)
    plot_coulomb_force_matrix.add_subparser(subparsers)
    analyze_coulomb_mobility.add_subparser(subparsers)
    plot_fortran_potential_slices.add_subparser(subparsers)
    estimate_fortran_workload.add_subparser(subparsers)
    plot_performance_profile.add_subparser(subparsers)
    model.add_subparser(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Run the unified ``beachx`` entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)
