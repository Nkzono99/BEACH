"""CLI for inspecting Fortran output files."""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace

from beach import Beach, compute_potential_mesh, list_fortran_runs

from . import _options
from ._output import emit


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``inspect`` options."""

    _options.add_output_dir(parser)
    _options.add_show(
        parser, help="display all matplotlib plots, including potential plots"
    )
    _options.add_save_dir(parser)
    parser.add_argument(
        "--recompute-potential",
        action="store_true",
        help=(
            "recompute the absolute potential for potential min/max with the "
            "native field kernel instead of reading mesh_potential.csv; "
            "may be expensive"
        ),
    )
    _options.add_library(parser)
    _options.add_view(parser)
    _options.add_periodic2_mesh(parser)
    _options.add_axis_unit(parser)
    _options.add_json(parser)
    parser.epilog = (
        "--save-dir writes charge_bar.png, charge_mesh.png, and potential_mesh.png "
        "(the figures --show displays); the potential mesh uses mesh_potential.csv "
        "when present and otherwise computes the potential."
    )


def run(args: argparse.Namespace, beach: Beach) -> None:
    """Execute the output-inspection command."""

    result = beach.result
    report: dict[str, object] = {
        "directory": result.directory,
        "mesh_nelem": result.mesh_nelem,
        "processed_particles": result.processed_particles,
        "absorbed": result.absorbed,
        "escaped": result.escaped,
        "batches": result.batches,
        "last_rel_change": result.last_rel_change,
    }
    if result.simulated_time_s is not None:
        report["simulated_time_s"] = result.simulated_time_s
    if (
        result.periodic2_max_nonzero_mode_potential_step_v is not None
        and result.periodic2_max_nonzero_mode_potential_step_v > 0.0
    ):
        report.update(
            {
                "adaptive_nonzero_mode_limit_V": (
                    result.periodic2_max_nonzero_mode_potential_step_v
                ),
                "adaptive_nonzero_mode_last_step_V": (
                    result.adaptive_nonzero_mode_last_potential_step_v or 0.0
                ),
                "adaptive_nonzero_mode_last_duration_s": (
                    result.adaptive_nonzero_mode_last_batch_duration_s or 0.0
                ),
                "adaptive_nonzero_mode_rejected_trials": (
                    result.adaptive_nonzero_mode_rejected_trials
                ),
                "adaptive_nonzero_mode_omp_threads": (
                    result.adaptive_nonzero_mode_omp_threads or 0
                ),
            }
        )
    report["charge_sum"] = float(result.charges.sum())
    if result.mesh_ids is not None:
        report["mesh_ids"] = sorted(int(v) for v in set(result.mesh_ids.tolist()))
    if result.mesh_sources is not None:
        report["mesh_sources"] = [
            {
                "mesh_id": src.mesh_id,
                "source": src.source_kind,
                "template": src.template_kind,
                "surface": src.surface_model,
                "epsilon_r": src.epsilon_r,
                "elems": src.elem_count,
            }
            for src in (result.mesh_sources[key] for key in sorted(result.mesh_sources))
        ]
    potential = None
    if args.recompute_potential:
        # Without saved values compute_potential_mesh evaluates the kernel.
        potential = compute_potential_mesh(
            replace(result, mesh_potential_v=None),
            library_path=args.library,
        )
    elif result.mesh_potential_v is not None:
        potential = result.mesh_potential_v
    if potential is not None and potential.size > 0:
        report["potential_min"] = float(potential.min())
        report["potential_max"] = float(potential.max())
    if result.history is not None and result.history.has_data:
        report["charge_history_shape"] = [result.mesh_nelem, len(result.history)]
        report["batch_indices"] = result.history.batch_indices
        report["processed_particles_by_batch"] = (
            result.history.processed_particles_by_batch
        )

    if args.save_dir is not None or args.show:
        report.update(_plot_figures(args, beach))

    runs = list_fortran_runs(args.output_dir.parent)
    if runs:
        report["sibling_runs"] = runs
    emit(args, report)

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


def _plot_figures(args: argparse.Namespace, beach: Beach) -> dict[str, object]:
    """Draw the inspect figures and save each one as soon as it is drawn."""

    saved: dict[str, object] = {}
    if args.save_dir is not None:
        args.save_dir.mkdir(parents=True, exist_ok=True)

    def keep(name: str, fig) -> None:
        if args.save_dir is not None:
            path = args.save_dir / f"{name}.png"
            fig.savefig(path, dpi=150)
            saved[f"saved_{name}"] = path

    keep("charge_bar", beach.plot_bar()[0])
    if beach.result.triangles is None:
        print("mesh_triangles.csv not found; mesh plots are skipped", file=sys.stderr)
        return saved
    mesh_options = {
        "view_elev": args.view_elev,
        "view_azim": args.view_azim,
        "apply_periodic2_mesh": args.apply_periodic2_mesh,
        "periodic2_mesh_mode": args.periodic2_mesh_mode,
        "periodic2_repeat": args.periodic2_repeat,
        "axis_unit": args.axis_unit,
    }
    keep("charge_mesh", beach.plot_mesh(**mesh_options)[0])
    keep(
        "potential_mesh",
        beach.plot_potential(**mesh_options, library_path=args.library)[0],
    )
    return saved
