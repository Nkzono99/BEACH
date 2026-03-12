"""Plot area-weighted mesh-source boxplots from Fortran outputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument(
        "--quantity",
        choices=("charge", "potential"),
        default="charge",
        help="boxplot quantity",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=-1,
        help="history step for charge snapshot (-1 means latest)",
    )
    parser.add_argument(
        "--softening",
        type=float,
        default=0.0,
        help="softening length [m] in potential mode",
    )
    parser.add_argument(
        "--self-term",
        choices=("area-equivalent", "exclude", "softened-point"),
        default="area-equivalent",
        help="potential self-term model in potential mode",
    )
    parser.add_argument(
        "--hide-fliers",
        action="store_true",
        help="hide outlier markers",
    )
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="path to save the figure",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="show matplotlib window",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    from beach import Beach

    args = build_parser().parse_args(argv)
    self_term = args.self_term.replace("-", "_")
    run = Beach(args.output_dir)
    fig, _ = run.plot_mesh_source_boxplot(
        quantity=args.quantity,
        step=args.step,
        softening=args.softening,
        self_term=self_term,
        showfliers=not args.hide_fliers,
    )

    if args.save is not None:
        fig.savefig(args.save, dpi=150)
        print(f"saved={args.save}")

    if args.show or args.save is None:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
