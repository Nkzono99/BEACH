"""Animate Fortran charge history on the mesh and save as GIF."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from beach import Beach


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="outputs/latest")
    parser.add_argument(
        "--quantity",
        choices=("charge", "potential"),
        default="charge",
        help="value to visualize on mesh snapshots",
    )
    parser.add_argument(
        "--save-gif",
        type=Path,
        default=None,
        help="output GIF path (default: <output_dir>/<quantity>_history.gif)",
    )
    parser.add_argument("--fps", type=int, default=10, help="GIF frame rate")
    parser.add_argument(
        "--frame-stride",
        type=int,
        default=1,
        help="use every N-th history snapshot",
    )
    parser.add_argument(
        "--cmap",
        default=None,
        help="matplotlib colormap name (defaults: charge=coolwarm, potential=viridis)",
    )
    parser.add_argument(
        "--potential-softening",
        type=float,
        default=0.0,
        help="smoothing length [m] for potential reconstruction",
    )
    parser.add_argument(
        "--potential-self-term",
        choices=("area-equivalent", "exclude", "softened-point"),
        default="area-equivalent",
        help="self-term treatment for potential reconstruction",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    save_gif = (
        args.save_gif
        if args.save_gif is not None
        else output_dir / f"{args.quantity}_history.gif"
    )
    self_term = args.potential_self_term.replace("-", "_")
    beach = Beach(output_dir)
    result = beach.result

    try:
        written = beach.animate_mesh(
            output_path=save_gif,
            quantity=args.quantity,
            fps=args.fps,
            frame_stride=args.frame_stride,
            cmap=args.cmap,
            softening=args.potential_softening,
            self_term=self_term,
        )
    except ModuleNotFoundError as exc:
        if exc.name is not None and (
            exc.name.startswith("matplotlib") or exc.name.startswith("PIL")
        ):
            raise SystemExit(
                "matplotlib and pillow are required for GIF animation. "
                "Install dependencies with `python -m pip install -e . --no-build-isolation`."
            ) from exc
        raise

    print(f"saved_gif={written}")
    print(f"quantity={args.quantity}")
    frame_count = (
        result.charge_history.shape[1] if result.charge_history is not None else 0
    )
    print(f"frames={frame_count}")


if __name__ == "__main__":
    main()
