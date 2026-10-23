"""CLI for plotting simulation-box potential slices."""

from __future__ import annotations

import argparse
from pathlib import Path

from beach import Beach
from beach.analysis.periodic import periodic2_from_sim
from beach.config import ConfigError, load_config_file

from . import _options


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``slices`` options."""

    _options.add_output_dir(parser)
    _options.add_config(parser)
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
    _options.add_library(parser)
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2048,
        help="number of sampling points per chunk",
    )
    _options.add_cmap(parser, "potential")
    parser.add_argument(
        "--vmin",
        type=_options.finite_float,
        default=None,
        help="color scale lower bound [V] (default: auto from sampled data)",
    )
    parser.add_argument(
        "--vmax",
        type=_options.finite_float,
        default=None,
        help="color scale upper bound [V] (default: auto from sampled data)",
    )
    _options.add_save(parser, default="<output_dir>/potential_slices.png")
    _options.add_dpi(parser)
    _options.add_axis_unit(parser)
    _options.add_show(parser)


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
        "beach.toml is required to read domain.box_min/box_max. "
        "Pass --config <path/to/beach.toml>."
    )


def _load_sim_box(
    config_path: Path,
    *,
    allow_historical_root_oracle: bool = False,
) -> tuple[list[float], list[float], bool, dict[str, object] | None]:
    try:
        config = load_config_file(config_path)
    except ConfigError as exc:
        raise SystemExit(f"{config_path}: {exc}") from exc
    sim = config.get("sim")
    if not isinstance(sim, dict):
        raise SystemExit("[sim] section is missing in beach.toml.")
    domain = config.get("domain")
    if isinstance(domain, dict):
        box = domain
        prefix = "domain"
    else:
        box = sim
        prefix = "sim"
    if "box_min" not in box or "box_max" not in box:
        raise SystemExit(
            f"{prefix}.box_min and {prefix}.box_max are required in beach.toml."
        )

    box_min = _coerce_vec3(box["box_min"], name=f"{prefix}.box_min")
    box_max = _coerce_vec3(box["box_max"], name=f"{prefix}.box_max")
    for axis, low, high in zip(("x", "y", "z"), box_min, box_max):
        if high <= low:
            raise SystemExit(
                f"{prefix}.box_max[{axis}] must be greater than {prefix}.box_min[{axis}]."
            )

    use_box = isinstance(domain, dict) or bool(sim.get("use_box", False))
    periodic_sim = dict(sim)
    if isinstance(domain, dict):
        periodic_sim["box_min"] = box_min
        periodic_sim["box_max"] = box_max
        axes = domain.get("periodic_axes", [])
        for axis in ("x", "y", "z"):
            mode = "periodic" if axis in axes else "open"
            periodic_sim[f"bc_{axis}_low"] = mode
            periodic_sim[f"bc_{axis}_high"] = mode
    field_boundary = config.get("field_boundary")
    if isinstance(field_boundary, dict):
        periodic_sim["field_bc_mode"] = field_boundary.get("mode", "free")
    try:
        periodic2 = periodic2_from_sim(
            periodic_sim,
            allow_historical_root_oracle=allow_historical_root_oracle,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
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


def run(args: argparse.Namespace, beach: Beach) -> None:
    """Execute the potential-slice plotting command."""

    # The color limits are checked here because the plot checks them only after
    # sampling every slice through the field kernel.
    if (
        args.vmin is not None
        and args.vmax is not None
        and float(args.vmin) >= float(args.vmax)
    ):
        args._parser.error("--vmin must be smaller than --vmax.")

    output_dir = args.output_dir
    save_path = (
        args.save if args.save is not None else output_dir / "potential_slices.png"
    )

    config_path = _find_config_path(output_dir, args.config)
    box_min, box_max, use_box, periodic2 = _load_sim_box(
        config_path,
        allow_historical_root_oracle=args.config is None,
    )
    xy_z, yz_x, xz_y = _default_slice_values(
        box_min,
        box_max,
        xy_z=args.xy_z,
        yz_x=args.yz_x,
        xz_y=args.xz_y,
    )

    fig, _ = beach.plot_potential_slices(
        box_min=box_min,
        box_max=box_max,
        grid_n=args.grid_n,
        xy_z=xy_z,
        yz_x=yz_x,
        xz_y=xz_y,
        chunk_size=args.chunk_size,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        periodic2=periodic2,
        axis_unit=args.axis_unit,
        config_path=config_path,
        library_path=args.library,
    )

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

    import matplotlib.pyplot as plt

    if args.show:
        plt.show()
    else:
        plt.close(fig)

