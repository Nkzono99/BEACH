"""Option groups and argument types shared by ``beachx`` subcommands.

Each subcommand's ``configure(parser)`` adds the groups it needs from here, so an
option that several commands accept has one spelling, one default, and one help
text. Plotting defaults come from :mod:`beach.plot._defaults` through this module
only.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from beach.plot._defaults import (
    DEFAULT_CHARGE_CMAP,
    DEFAULT_FORCE_CMAP,
    DEFAULT_MESH_VIEW_AZIM,
    DEFAULT_MESH_VIEW_ELEV,
    DEFAULT_POTENTIAL_CMAP,
)

_CMAPS = {
    "charge": DEFAULT_CHARGE_CMAP,
    "potential": DEFAULT_POTENTIAL_CMAP,
    "force": DEFAULT_FORCE_CMAP,
}


def finite_float(value: str) -> float:
    """Parse a finite float."""

    parsed = float(value)
    if not math.isfinite(parsed):
        raise argparse.ArgumentTypeError("must be finite")
    return parsed


def positive_float(value: str) -> float:
    """Parse a finite float greater than zero."""

    parsed = finite_float(value)
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError("must be > 0")
    return parsed


def nonnegative_float(value: str) -> float:
    """Parse a finite float greater than or equal to zero."""

    parsed = finite_float(value)
    if parsed < 0.0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return parsed


def positive_int(value: str) -> int:
    """Parse an integer greater than zero."""

    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return parsed


def nonnegative_int(value: str) -> int:
    """Parse an integer greater than or equal to zero."""

    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return parsed


def charge_step(value: str) -> int | None:
    """Parse ``final`` (``charges.csv``) or a history batch index."""

    normalized = value.strip().lower()
    if normalized == "final":
        return None
    try:
        return int(normalized)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "step must be final or an integer batch."
        ) from exc


def add_output_dir(parser: argparse.ArgumentParser) -> None:
    """Add the run directory; the registry loads it before ``run`` is called."""

    parser.add_argument(
        "output_dir",
        nargs="?",
        type=Path,
        default=Path("outputs/latest"),
        help="Fortran output directory (default: %(default)s)",
    )
    parser.set_defaults(_opens_run=True)


def add_config(parser: argparse.ArgumentParser, *, required: bool = False) -> None:
    """Add ``--config``."""

    parser.add_argument(
        "--config",
        type=Path,
        required=required,
        default=None,
        help=(
            "beach.toml used by the simulation"
            if required
            else "path to beach.toml (default: auto-detect near output_dir)"
        ),
    )


def add_library(parser: argparse.ArgumentParser) -> None:
    """Add ``--library``."""

    parser.add_argument(
        "--library",
        type=Path,
        default=None,
        help="path to libbeach_field_kernel shared library",
    )


def add_step(parser: argparse.ArgumentParser, *, default: int | None = -1) -> None:
    """Add ``--step`` for the charge snapshot."""

    parser.add_argument(
        "--step",
        type=charge_step,
        default=default,
        metavar="final|BATCH",
        help=(
            "charge snapshot: final uses charges.csv, an integer selects that stored "
            "history batch, and -1 selects the last stored history row "
            f"(default: {'final' if default is None else default})"
        ),
    )


def add_periodic2_mesh(
    parser: argparse.ArgumentParser, *, mesh_mode: bool = True
) -> None:
    """Add ``--apply-periodic2-mesh``, ``--periodic2-mesh-mode``, and ``--periodic2-repeat``."""

    parser.add_argument(
        "--apply-periodic2-mesh",
        action="store_true",
        help=(
            "wrap each plotted triangle into the periodic2 cell from nearby "
            "beach.toml using its centroid"
        ),
    )
    if mesh_mode:
        parser.add_argument(
            "--periodic2-mesh-mode",
            choices=("centroid", "mesh"),
            default="centroid",
            help=(
                "periodic2 wrapping unit for mesh plots; 'mesh' keeps each mesh_id "
                "object intact"
            ),
        )
    parser.add_argument(
        "--periodic2-repeat",
        type=nonnegative_int,
        default=0,
        help=(
            "replicate the mesh over periodic images; n produces (2n+1)^2 copies "
            "(default: 0 = no replication)"
        ),
    )


def add_view(parser: argparse.ArgumentParser) -> None:
    """Add ``--view-elev`` and ``--view-azim`` for 3D mesh plots."""

    parser.add_argument(
        "--view-elev",
        type=finite_float,
        default=DEFAULT_MESH_VIEW_ELEV,
        help="elevation angle [deg] for 3D mesh plots (default: %(default)s)",
    )
    parser.add_argument(
        "--view-azim",
        type=finite_float,
        default=DEFAULT_MESH_VIEW_AZIM,
        help="azimuth angle [deg] for 3D mesh plots (default: %(default)s)",
    )


def add_axis_unit(parser: argparse.ArgumentParser) -> None:
    """Add ``--axis-unit`` for plot coordinates."""

    parser.add_argument(
        "--axis-unit",
        choices=("m", "um", "nm"),
        default="m",
        help="coordinate unit for plot axes (default: %(default)s)",
    )


def add_cmap(parser: argparse.ArgumentParser, *quantities: str) -> None:
    """Add ``--cmap``; several quantities keep each quantity's own default."""

    if len(quantities) == 1:
        parser.add_argument(
            "--cmap",
            default=_CMAPS[quantities[0]],
            help="matplotlib colormap name (default: %(default)s)",
        )
        return
    defaults = ", ".join(f"{name}={_CMAPS[name]}" for name in quantities)
    parser.add_argument(
        "--cmap",
        default=None,
        help=f"matplotlib colormap name (defaults: {defaults})",
    )


def add_save(parser: argparse.ArgumentParser, *, default: str | None = None) -> None:
    """Add ``--save PATH`` for a command that writes one file.

    ``default`` describes where the file goes when ``--save`` is omitted; without
    it the command writes nothing unless ``--save`` is given.
    """

    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        metavar="PATH",
        help=f"output file path (default: {default})" if default else "output file path",
    )


def add_save_dir(parser: argparse.ArgumentParser, *, default: str | None = None) -> None:
    """Add ``--save-dir DIR`` for a command that writes a fixed set of files."""

    parser.add_argument(
        "--save-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            f"directory for the output files (default: {default})"
            if default
            else "directory for the output files"
        ),
    )


def add_dpi(parser: argparse.ArgumentParser) -> None:
    """Add ``--dpi`` for saved figures."""

    parser.add_argument(
        "--dpi", type=positive_int, default=150, help="saved image DPI (default: %(default)s)"
    )


def add_show(
    parser: argparse.ArgumentParser, *, help: str = "display matplotlib window"
) -> None:
    """Add ``--show``."""

    parser.add_argument("--show", action="store_true", help=help)


def add_json(parser: argparse.ArgumentParser) -> None:
    """Add ``--json``."""

    parser.add_argument(
        "--json",
        action="store_true",
        help="print the report as one JSON object instead of key=value lines",
    )
