"""Legacy console-script entry points for BEACH post-processing CLIs."""

from __future__ import annotations

from typing import Sequence

from . import (
    animate_fortran_history,
    estimate_fortran_workload,
    inspect_fortran_output,
    plot_fortran_potential_slices,
    plot_performance_profile,
)
from ._shared import print_legacy_warning


def inspect_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-inspect`` alias."""

    print_legacy_warning("beach-inspect", "inspect")
    inspect_fortran_output.main(argv)


def animate_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-animate-history`` alias."""

    print_legacy_warning("beach-animate-history", "animate")
    animate_fortran_history.main(argv)


def slices_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-plot-potential-slices`` alias."""

    print_legacy_warning("beach-plot-potential-slices", "slices")
    plot_fortran_potential_slices.main(argv)


def workload_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-estimate-workload`` alias."""

    print_legacy_warning("beach-estimate-workload", "workload")
    estimate_fortran_workload.main(argv)


def profile_main(argv: Sequence[str] | None = None) -> None:
    """Run the deprecated ``beach-plot-performance-profile`` alias."""

    print_legacy_warning("beach-plot-performance-profile", "profile")
    plot_performance_profile.main(argv)
