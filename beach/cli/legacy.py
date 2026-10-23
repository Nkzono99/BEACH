"""Legacy console-script entry points for BEACH post-processing CLIs."""

from __future__ import annotations

import sys
from typing import Sequence

from ._registry import main


def _forward(legacy_name: str, subcommand: str, argv: Sequence[str] | None) -> None:
    print(
        f"WARNING: `{legacy_name}` is deprecated; use `beachx {subcommand}` instead.",
        file=sys.stderr,
    )
    main([subcommand, *(sys.argv[1:] if argv is None else argv)])


def inspect_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-inspect`` alias."""

    _forward("beach-inspect", "inspect", argv)


def animate_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-animate-history`` alias."""

    _forward("beach-animate-history", "animate", argv)


def coulomb_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-plot-coulomb-force-matrix`` alias."""

    _forward("beach-plot-coulomb-force-matrix", "coulomb", argv)


def slices_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-plot-potential-slices`` alias."""

    _forward("beach-plot-potential-slices", "slices", argv)


def workload_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-estimate-workload`` alias."""

    _forward("beach-estimate-workload", "workload", argv)


def profile_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-plot-performance-profile`` alias."""

    _forward("beach-plot-performance-profile", "profile", argv)
