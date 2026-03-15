"""Unified post-processing CLI for BEACH."""

from __future__ import annotations

import argparse
from typing import Sequence

from . import (
    animate_fortran_history,
    estimate_fortran_workload,
    inspect_fortran_output,
    plot_fortran_potential_slices,
    plot_performance_profile,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the root parser for ``beachx``."""

    parser = argparse.ArgumentParser(
        prog="beachx",
        description="BEACH post-processing tools for inspection and visualization.",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        metavar="command",
        required=True,
    )
    inspect_fortran_output.add_subparser(subparsers)
    animate_fortran_history.add_subparser(subparsers)
    plot_fortran_potential_slices.add_subparser(subparsers)
    estimate_fortran_workload.add_subparser(subparsers)
    plot_performance_profile.add_subparser(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Run the unified ``beachx`` entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)
