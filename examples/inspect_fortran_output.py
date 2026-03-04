"""Inspect Fortran output files from Python."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from beach import (
    compute_potential_mesh,
    list_fortran_runs,
    load_fortran_result,
    plot_charge_mesh,
    plot_charges,
    plot_potential_mesh,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument("--show", action="store_true", help="display matplotlib window")
    parser.add_argument("--save-bar", type=Path, default=None, help="save bar-chart figure path")
    parser.add_argument("--save-mesh", type=Path, default=None, help="save 3D mesh figure path")
    parser.add_argument(
        "--save-potential-mesh",
        type=Path,
        default=None,
        help="save 3D electric-potential mesh figure path",
    )
    parser.add_argument(
        "--potential-softening",
        type=float,
        default=1.0e-6,
        help="softening length [m] used for potential reconstruction",
    )
    args = parser.parse_args()

    result = load_fortran_result(args.output_dir)
    print(f"directory={result.directory}")
    print(f"mesh_nelem={result.mesh_nelem}")
    print(f"processed_particles={result.processed_particles}")
    print(f"absorbed={result.absorbed} escaped={result.escaped}")
    print(f"batches={result.batches} last_rel_change={result.last_rel_change:.6e}")
    print(f"charge_sum={result.charges.sum():.6e}")
    if result.triangles is not None:
        potential = compute_potential_mesh(result, softening=args.potential_softening)
        print(f"potential_min={potential.min():.6e}")
        print(f"potential_max={potential.max():.6e}")
    if result.charge_history is not None:
        print(f"charge_history_shape={result.charge_history.shape}")
        print(f"batch_indices={result.batch_indices}")
        print(f"processed_particles_by_batch={result.processed_particles_by_batch}")

    need_bar_plot = args.save_bar is not None or args.show
    need_mesh_plot = args.save_mesh is not None or args.show
    need_potential_mesh_plot = args.save_potential_mesh is not None or args.show

    try:
        if need_bar_plot:
            bar_fig, _ = plot_charges(result)
            if args.save_bar is not None:
                bar_fig.savefig(args.save_bar, dpi=150)
                print(f"saved_bar={args.save_bar}")

        if need_mesh_plot:
            if result.triangles is not None:
                mesh_fig, _ = plot_charge_mesh(result)
                if args.save_mesh is not None:
                    mesh_fig.savefig(args.save_mesh, dpi=150)
                    print(f"saved_mesh={args.save_mesh}")
            else:
                print("mesh_triangles.csv not found; mesh visualization is skipped")

        if need_potential_mesh_plot:
            if result.triangles is not None:
                potential_mesh_fig, _ = plot_potential_mesh(
                    result,
                    softening=args.potential_softening,
                )
                if args.save_potential_mesh is not None:
                    potential_mesh_fig.savefig(args.save_potential_mesh, dpi=150)
                    print(f"saved_potential_mesh={args.save_potential_mesh}")
            else:
                print("mesh_triangles.csv not found; potential mesh visualization is skipped")
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


if __name__ == "__main__":
    main()
