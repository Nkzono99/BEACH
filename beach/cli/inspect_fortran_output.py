"""CLI for inspecting Fortran output files."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from beach import Beach, list_fortran_runs

from ._shared import (
    configure_entry_parser,
    require_finite,
)

COMMAND_NAME = "inspect"
LEGACY_COMMAND_NAME = "beach-inspect"


def _configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument(
        "--show",
        action="store_true",
        help="display all matplotlib plots, including potential plots",
    )
    parser.add_argument(
        "--save-bar", type=Path, default=None, help="save bar-chart figure path"
    )
    parser.add_argument(
        "--save-mesh", type=Path, default=None, help="save 3D mesh figure path"
    )
    parser.add_argument(
        "--save-potential-mesh",
        type=Path,
        default=None,
        help="save 3D electric-potential mesh figure path; may calculate potential",
    )
    parser.add_argument(
        "--recompute-potential",
        action="store_true",
        help=(
            "recompute potential plots and potential min/max through "
            "Beach.compute_potential instead of using only mesh_potential.csv; "
            "may be expensive"
        ),
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=None,
        help="path to libbeach_field_kernel shared library",
    )
    parser.add_argument(
        "--view-elev",
        type=float,
        default=24.0,
        help="elevation angle [deg] for 3D mesh plots",
    )
    parser.add_argument(
        "--view-azim",
        type=float,
        default=-58.0,
        help="azimuth angle [deg] for 3D mesh plots",
    )
    parser.add_argument(
        "--apply-periodic2-mesh",
        action="store_true",
        help=(
            "wrap each plotted triangle into the periodic2 cell from nearby "
            "beach.toml using its centroid"
        ),
    )
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
        type=int,
        default=0,
        help=(
            "replicate the mesh over periodic images; n produces (2n+1)^2 copies "
            "(default: 0 = no replication)"
        ),
    )
    parser.add_argument(
        "--axis-unit",
        choices=("m", "um", "nm"),
        default="m",
        help="coordinate unit for 3D mesh plot axes",
    )


def build_parser(*, prog: str | None = LEGACY_COMMAND_NAME) -> argparse.ArgumentParser:
    """Build the argument parser for output inspection CLI."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command under the unified ``beachx`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="inspect Fortran output files",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def run(args: argparse.Namespace) -> None:
    """Execute the output-inspection command."""

    parser = args._parser
    require_finite(parser, args.view_elev, "--view-elev")
    require_finite(parser, args.view_azim, "--view-azim")
    if args.periodic2_repeat < 0:
        parser.error("--periodic2-repeat must be >= 0.")

    reference_point = "species1_injection_center"

    beach = Beach(args.output_dir)
    try:
        result = beach.result
    except FileNotFoundError as exc:
        raise SystemExit(
            f'Fortran output files are missing under "{args.output_dir}". '
            "Expected at least summary.txt and charges.csv."
        ) from exc
    print(f"directory={result.directory}")
    print(f"mesh_nelem={result.mesh_nelem}")
    print(f"processed_particles={result.processed_particles}")
    print(f"absorbed={result.absorbed} escaped={result.escaped}")
    print(f"batches={result.batches} last_rel_change={result.last_rel_change:.6e}")
    if result.simulated_time_s is not None:
        print(f"simulated_time_s={result.simulated_time_s:.6e}")
    if (
        result.periodic2_max_nonzero_mode_potential_step_v is not None
        and result.periodic2_max_nonzero_mode_potential_step_v > 0.0
    ):
        print(
            "adaptive_nonzero_mode="
            f"limit_V:{result.periodic2_max_nonzero_mode_potential_step_v:.6e} "
            f"last_step_V:{result.adaptive_nonzero_mode_last_potential_step_v or 0.0:.6e} "
            f"last_duration_s:{result.adaptive_nonzero_mode_last_batch_duration_s or 0.0:.6e} "
            f"rejected_trials:{result.adaptive_nonzero_mode_rejected_trials} "
            f"omp_threads:{result.adaptive_nonzero_mode_omp_threads or 0}"
        )
    print(f"charge_sum={result.charges.sum():.6e}")
    if result.mesh_ids is not None:
        mesh_ids = sorted(int(v) for v in set(result.mesh_ids.tolist()))
        print(f"mesh_ids={mesh_ids}")
    if result.mesh_sources is not None:
        for mesh_id in sorted(result.mesh_sources):
            src = result.mesh_sources[mesh_id]
            print(
                "mesh_source="
                f"id:{src.mesh_id} source:{src.source_kind} "
                f"template:{src.template_kind} surface:{src.surface_model} "
                f"epsilon_r:{src.epsilon_r:g} elems:{src.elem_count}"
            )
    potential = None
    if args.recompute_potential:
        potential = beach.compute_potential(
            reference_point=reference_point,
            library_path=args.library,
        )
    elif result.mesh_potential_v is not None:
        potential = result.mesh_potential_v

    if potential is not None and potential.size > 0:
        print(f"potential_min={potential.min():.6e}")
        print(f"potential_max={potential.max():.6e}")
    if result.history is not None and result.history.has_data:
        snapshot_count = len(result.history)
        print(f"charge_history_shape=({result.mesh_nelem}, {snapshot_count})")
        print(f"batch_indices={result.history.batch_indices}")
        print(
            "processed_particles_by_batch="
            f"{result.history.processed_particles_by_batch}"
        )

    need_bar_plot = args.save_bar is not None or args.show
    need_mesh_plot = args.save_mesh is not None or args.show
    need_potential_mesh_plot = args.save_potential_mesh is not None or args.show

    try:
        if need_bar_plot:
            bar_fig, _ = beach.plot_bar()
            if args.save_bar is not None:
                bar_fig.savefig(args.save_bar, dpi=150)
                print(f"saved_bar={args.save_bar}")

        if need_mesh_plot:
            if result.triangles is not None:
                mesh_fig, _ = beach.plot_mesh(
                    view_elev=args.view_elev,
                    view_azim=args.view_azim,
                    apply_periodic2_mesh=args.apply_periodic2_mesh,
                    periodic2_mesh_mode=args.periodic2_mesh_mode,
                    periodic2_repeat=args.periodic2_repeat,
                    axis_unit=args.axis_unit,
                )
                if args.save_mesh is not None:
                    mesh_fig.savefig(args.save_mesh, dpi=150)
                    print(f"saved_mesh={args.save_mesh}")
            else:
                print("mesh_triangles.csv not found; mesh visualization is skipped")

        if need_potential_mesh_plot:
            if result.triangles is not None:
                potential_mesh_fig, _ = beach.plot_potential(
                    reference_point=reference_point,
                    view_elev=args.view_elev,
                    view_azim=args.view_azim,
                    apply_periodic2_mesh=args.apply_periodic2_mesh,
                    periodic2_mesh_mode=args.periodic2_mesh_mode,
                    periodic2_repeat=args.periodic2_repeat,
                    axis_unit=args.axis_unit,
                    library_path=args.library,
                )
                if args.save_potential_mesh is not None:
                    potential_mesh_fig.savefig(args.save_potential_mesh, dpi=150)
                    print(f"saved_potential_mesh={args.save_potential_mesh}")
            else:
                print(
                    "mesh_triangles.csv not found; potential mesh visualization is skipped"
                )
    except ModuleNotFoundError as exc:
        if exc.name is not None and exc.name.startswith("matplotlib"):
            raise SystemExit(
                "matplotlib is required for visualization. "
                "Install dependencies with `python -m pip install -e . --no-build-isolation`."
            ) from exc
        raise

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()

    runs = list_fortran_runs(Path(args.output_dir).parent)
    if runs:
        print("sibling_runs=")
        for run in runs:
            print(f"  - {run}")


def main(argv: Sequence[str] | None = None) -> None:
    """Run the output-inspection CLI entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
