"""CLI for plotting simulation-box potential slices."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Sequence

from beach import Beach

from ._shared import configure_entry_parser

COMMAND_NAME = "slices"
LEGACY_COMMAND_NAME = "beach-plot-potential-slices"


def _configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="path to beach.toml (default: auto-detect near output_dir)",
    )
    parser.add_argument(
        "--grid-n",
        type=int,
        default=200,
        help="grid count for each slice axis (N x N)",
    )
    parser.add_argument(
        "--xy-z",
        type=float,
        default=None,
        help="z coordinate [m] for XY slice (default: box center)",
    )
    parser.add_argument(
        "--yz-x",
        type=float,
        default=None,
        help="x coordinate [m] for YZ slice (default: box center)",
    )
    parser.add_argument(
        "--xz-y",
        type=float,
        default=None,
        help="y coordinate [m] for XZ slice (default: box center)",
    )
    parser.add_argument(
        "--potential-softening",
        type=float,
        default=None,
        help="smoothing length [m] used for potential sampling; default uses sim.softening when available",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2048,
        help="number of sampling points per chunk",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="matplotlib colormap name",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="color scale lower bound [V] (default: auto from sampled data)",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="color scale upper bound [V] (default: auto from sampled data)",
    )
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="output image path (default: <output_dir>/potential_slices.png)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="saved image DPI",
    )
    parser.add_argument("--show", action="store_true", help="display matplotlib window")


def build_parser(*, prog: str | None = LEGACY_COMMAND_NAME) -> argparse.ArgumentParser:
    """Build the argument parser for potential-slice plotting CLI."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def _load_toml(path: Path) -> dict[str, Any]:
    try:
        import tomllib  # py311+

        with path.open("rb") as stream:
            return tomllib.load(stream)
    except ModuleNotFoundError:
        try:
            import tomli  # type: ignore

            with path.open("rb") as stream:
                return tomli.load(stream)
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "TOML parser is missing. Use Python 3.11+ or install tomli: "
                "`python -m pip install tomli`."
            ) from exc


def _coerce_vec3(value: object, *, name: str) -> list[float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise SystemExit(f"{name} must be an array of 3 numbers in beach.toml.")
    try:
        return [float(v) for v in value]
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"{name} must contain numeric values.") from exc


def _find_config_path(output_dir: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        if not explicit.exists():
            raise SystemExit(f'Config file is not found: "{explicit}".')
        return explicit

    candidates = (
        output_dir / "beach.toml",
        output_dir.parent / "beach.toml",
        output_dir.parent.parent / "beach.toml",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise SystemExit(
        "beach.toml is required to read sim.box_min/box_max. "
        "Pass --config <path/to/beach.toml>."
    )


def _resolve_periodic2_from_sim(
    sim: dict[str, Any],
    *,
    box_min: list[float],
    box_max: list[float],
) -> dict[str, object] | None:
    field_bc_mode = str(sim.get("field_bc_mode", "free")).strip().lower()
    if field_bc_mode != "periodic2":
        return None

    periodic_axes: list[int] = []
    for axis_idx, axis_name in enumerate(("x", "y", "z")):
        low_raw = str(sim.get(f"bc_{axis_name}_low", "open")).strip().lower()
        high_raw = str(sim.get(f"bc_{axis_name}_high", "open")).strip().lower()
        low = "open" if low_raw in {"open", "outflow", "escape"} else low_raw
        high = "open" if high_raw in {"open", "outflow", "escape"} else high_raw
        if (low == "periodic") != (high == "periodic"):
            raise SystemExit(
                "periodic2 requires bc_low(axis)=bc_high(axis)=periodic for periodic axes."
            )
        if low == "periodic":
            periodic_axes.append(axis_idx)

    if len(periodic_axes) != 2:
        raise SystemExit('sim.field_bc_mode="periodic2" requires exactly two periodic axes.')

    lengths = [box_max[axis] - box_min[axis] for axis in periodic_axes]
    if lengths[0] <= 0.0 or lengths[1] <= 0.0:
        raise SystemExit("periodic2 requires positive box length on periodic axes.")

    try:
        image_layers = int(sim.get("field_periodic_image_layers", 1))
        far_correction = str(sim.get("field_periodic_far_correction", "none")).strip().lower()
        ewald_alpha = float(sim.get("field_periodic_ewald_alpha", 0.0))
        ewald_layers = int(sim.get("field_periodic_ewald_layers", 4))
    except (TypeError, ValueError) as exc:
        raise SystemExit("invalid periodic2 potential settings in [sim].") from exc

    if far_correction in {"none", "m2l_root"}:
        far_correction = "m2l_root_trunc"
        ewald_layers = max(1, ewald_layers)

    return {
        "axes": tuple(periodic_axes),
        "lengths": tuple(lengths),
        "origins": tuple(box_min[axis] for axis in periodic_axes),
        "image_layers": image_layers,
        "far_correction": far_correction,
        "ewald_alpha": ewald_alpha,
        "ewald_layers": ewald_layers,
    }


def _load_sim_box(
    config_path: Path,
) -> tuple[list[float], list[float], bool, dict[str, object] | None]:
    config = _load_toml(config_path)
    sim = config.get("sim")
    if not isinstance(sim, dict):
        raise SystemExit("[sim] section is missing in beach.toml.")

    if "box_min" not in sim or "box_max" not in sim:
        raise SystemExit("sim.box_min and sim.box_max are required in beach.toml.")

    box_min = _coerce_vec3(sim["box_min"], name="sim.box_min")
    box_max = _coerce_vec3(sim["box_max"], name="sim.box_max")
    for axis, low, high in zip(("x", "y", "z"), box_min, box_max):
        if high <= low:
            raise SystemExit(
                f"sim.box_max[{axis}] must be greater than sim.box_min[{axis}]."
            )

    use_box = bool(sim.get("use_box", False))
    periodic2 = _resolve_periodic2_from_sim(
        sim,
        box_min=box_min,
        box_max=box_max,
    )
    return box_min, box_max, use_box, periodic2


def _default_slice_values(
    box_min: list[float],
    box_max: list[float],
    *,
    xy_z: float | None,
    yz_x: float | None,
    xz_y: float | None,
) -> tuple[float, float, float]:
    xy = 0.5 * (box_min[2] + box_max[2]) if xy_z is None else float(xy_z)
    yz = 0.5 * (box_min[0] + box_max[0]) if yz_x is None else float(yz_x)
    xz = 0.5 * (box_min[1] + box_max[1]) if xz_y is None else float(xz_y)
    return xy, yz, xz


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command under the unified ``beachx`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="plot potential slices inside the simulation box",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def run(args: argparse.Namespace) -> None:
    """Execute the potential-slice plotting command."""

    parser = args._parser
    if args.grid_n < 2:
        parser.error("--grid-n must be >= 2.")
    if args.chunk_size <= 0:
        parser.error("--chunk-size must be > 0.")
    if args.dpi <= 0:
        parser.error("--dpi must be > 0.")
    if (
        args.vmin is not None
        and args.vmax is not None
        and float(args.vmin) >= float(args.vmax)
    ):
        parser.error("--vmin must be smaller than --vmax.")

    output_dir = Path(args.output_dir)
    save_path = (
        args.save if args.save is not None else output_dir / "potential_slices.png"
    )

    beach = Beach(output_dir)
    try:
        beach.result
    except FileNotFoundError as exc:
        raise SystemExit(
            f'Fortran output files are missing under "{output_dir}". '
            "Expected at least summary.txt and charges.csv."
        ) from exc

    config_path = _find_config_path(output_dir, args.config)
    box_min, box_max, use_box, periodic2 = _load_sim_box(config_path)
    xy_z, yz_x, xz_y = _default_slice_values(
        box_min,
        box_max,
        xy_z=args.xy_z,
        yz_x=args.yz_x,
        xz_y=args.xz_y,
    )

    try:
        fig, _ = beach.plot_potential_slices(
            box_min=box_min,
            box_max=box_max,
            grid_n=args.grid_n,
            xy_z=xy_z,
            yz_x=yz_x,
            xz_y=xz_y,
            softening=args.potential_softening,
            chunk_size=args.chunk_size,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            periodic2=periodic2,
        )
    except ModuleNotFoundError as exc:
        if exc.name is not None and exc.name.startswith("matplotlib"):
            raise SystemExit(
                "matplotlib is required for visualization. "
                "Install dependencies with `python -m pip install -e . --no-build-isolation`."
            ) from exc
        raise
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=args.dpi)

    print(f"saved_potential_slices={save_path}")
    print(f"config={config_path}")
    print(f"grid_n={args.grid_n}")
    print(f"box_min={box_min} box_max={box_max}")
    print(f"periodic2={periodic2}")
    print(f"xy_z={xy_z:.6e} yz_x={yz_x:.6e} xz_y={xz_y:.6e}")
    print(f"vmin={args.vmin} vmax={args.vmax}")
    print(f"use_box={use_box}")

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()
    else:
        import matplotlib.pyplot as plt

        plt.close(fig)


def main(argv: Sequence[str] | None = None) -> None:
    """Run the potential-slice plotting CLI entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
