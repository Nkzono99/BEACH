"""Inspect Fortran output files from Python."""

from __future__ import annotations

import argparse
from pathlib import Path

from bemtracer import list_fortran_runs, load_fortran_result, plot_charges


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument("--show", action="store_true", help="display matplotlib window")
    parser.add_argument("--save", type=Path, default=None, help="save figure path")
    args = parser.parse_args()

    result = load_fortran_result(args.output_dir)
    print(f"directory={result.directory}")
    print(f"mesh_nelem={result.mesh_nelem}")
    print(f"processed_particles={result.processed_particles}")
    print(f"absorbed={result.absorbed} escaped={result.escaped}")
    print(f"batches={result.batches} last_rel_change={result.last_rel_change:.6e}")
    print(f"charge_sum={result.charges.sum():.6e}")

    fig, _ = plot_charges(result)
    if args.save is not None:
        fig.savefig(args.save, dpi=150)
        print(f"saved={args.save}")
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
