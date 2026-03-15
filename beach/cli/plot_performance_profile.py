"""CLI for plotting BEACH performance profiles."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Sequence


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for performance-profile plotting CLI."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "profile_path",
        nargs="?",
        default="outputs/latest/performance_profile.csv",
        help="path to performance_profile.csv or its parent directory",
    )
    parser.add_argument(
        "--metric",
        choices=("rank_max_s", "rank_mean_s", "rank_min_s"),
        default="rank_max_s",
        help="time metric used for the primary bar chart",
    )
    parser.add_argument(
        "--secondary",
        choices=("auto", "share", "imbalance", "calls"),
        default="auto",
        help="secondary subplot quantity",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=0,
        help="show only top N regions after sorting (0 means all)",
    )
    parser.add_argument(
        "--include-zero",
        action="store_true",
        help="include zero-time rows in the plot",
    )
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        help="path to save the figure (default: profile path with .png suffix when not showing)",
    )
    parser.add_argument("--dpi", type=int, default=150, help="save DPI")
    parser.add_argument("--show", action="store_true", help="display matplotlib window")
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Run the performance-profile plotting CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)

    profile_path = _resolve_profile_path(args.profile_path)
    try:
        metadata, rows = load_performance_profile(profile_path)
    except FileNotFoundError as exc:
        raise SystemExit(
            f'Performance profile file is missing under "{profile_path}". '
            "Expected performance_profile.csv."
        ) from exc
    except ValueError as exc:
        raise SystemExit(f'Failed to parse "{profile_path}": {exc}') from exc

    try:
        fig, _ = plot_performance_profile(
            metadata,
            rows,
            metric=args.metric,
            secondary=args.secondary,
            top=args.top,
            include_zero=args.include_zero,
            profile_path=profile_path,
        )
    except ModuleNotFoundError as exc:
        if exc.name is not None and exc.name.startswith("matplotlib"):
            raise SystemExit(
                "matplotlib is required for visualization. "
                "Install dependencies with `python -m pip install -e . --no-build-isolation`."
            ) from exc
        raise

    save_path = args.save
    if save_path is None and not args.show:
        save_path = profile_path.with_suffix(".png")
    if save_path is not None:
        fig.savefig(save_path, dpi=args.dpi)
        print(f"saved={save_path}")

    _print_top_regions(rows, metric=args.metric, top=args.top)

    import matplotlib.pyplot as plt

    if args.show:
        plt.show()
    else:
        plt.close(fig)


def load_performance_profile(path: Path) -> tuple[dict[str, str], list[dict[str, Any]]]:
    """Load one ``performance_profile.csv`` file."""

    metadata: dict[str, str] = {}
    csv_lines: list[str] = []

    if not path.exists():
        raise FileNotFoundError(path)

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            comment = stripped[1:].strip()
            if "=" in comment:
                key, value = comment.split("=", 1)
                metadata[key.strip()] = value.strip()
            continue
        csv_lines.append(line)

    if not csv_lines:
        raise ValueError("CSV body is empty.")

    reader = csv.DictReader(csv_lines)
    rows: list[dict[str, Any]] = []
    for row in reader:
        if row is None:
            continue
        rows.append(
            {
                "region": _required_field(row, "region"),
                "calls_sum": float(_required_field(row, "calls_sum")),
                "calls_mean": float(_required_field(row, "calls_mean")),
                "rank_min_s": float(_required_field(row, "rank_min_s")),
                "rank_mean_s": float(_required_field(row, "rank_mean_s")),
                "rank_max_s": float(_required_field(row, "rank_max_s")),
                "imbalance_ratio": float(_required_field(row, "imbalance_ratio")),
            }
        )

    if not rows:
        raise ValueError("No profile rows were found.")
    return metadata, rows


def plot_performance_profile(
    metadata: dict[str, str],
    rows: list[dict[str, Any]],
    *,
    metric: str = "rank_max_s",
    secondary: str = "auto",
    top: int = 0,
    include_zero: bool = False,
    profile_path: Path | None = None,
):
    """Plot a performance profile as a two-panel horizontal bar chart."""

    import matplotlib.pyplot as plt

    rows_to_plot = list(rows)
    if not include_zero:
        rows_to_plot = [
            row
            for row in rows_to_plot
            if float(row[metric]) > 0.0 or float(row["calls_sum"]) > 0.0
        ]
    rows_to_plot.sort(key=lambda row: (float(row[metric]), row["region"]), reverse=True)
    if top > 0:
        rows_to_plot = rows_to_plot[:top]
    if not rows_to_plot:
        raise ValueError("No non-zero rows are available to plot.")

    world_size = int(metadata.get("mpi_world_size", "1"))
    secondary_mode = secondary
    if secondary_mode == "auto":
        secondary_mode = "imbalance" if world_size > 1 else "share"

    simulation_total = _simulation_total(rows, metric=metric)
    labels = [str(row["region"]) for row in reversed(rows_to_plot)]
    metric_values = [float(row[metric]) for row in reversed(rows_to_plot)]
    secondary_values, secondary_label = _secondary_series(
        rows_to_plot,
        metric=metric,
        mode=secondary_mode,
        simulation_total=simulation_total,
    )
    secondary_values = list(reversed(secondary_values))
    colors = [_region_color(name) for name in labels]

    height = max(4.5, 0.42 * len(labels) + 2.0)
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13.0, height),
        width_ratios=(1.8, 1.1),
        constrained_layout=True,
    )
    ax_time, ax_secondary = axes
    y_pos = list(range(len(labels)))

    ax_time.barh(y_pos, metric_values, color=colors, edgecolor="black", linewidth=0.4)
    ax_time.set_yticks(y_pos, labels=labels)
    ax_time.set_xlabel(f"{metric} [s]")
    ax_time.set_title("Time")
    ax_time.grid(axis="x", linestyle=":", alpha=0.35)

    max_metric = max(metric_values) if metric_values else 0.0
    time_pad = 0.04 * max_metric if max_metric > 0.0 else 1.0
    ax_time.set_xlim(0.0, max_metric + 3.0 * time_pad)
    for idx, row in enumerate(reversed(rows_to_plot)):
        ax_time.text(
            metric_values[idx] + time_pad,
            y_pos[idx],
            f"{metric_values[idx]:.3g}s  calls={float(row['calls_sum']):.0f}",
            va="center",
            fontsize=8.5,
        )

    ax_secondary.barh(
        y_pos,
        secondary_values,
        color=colors,
        edgecolor="black",
        linewidth=0.4,
    )
    ax_secondary.set_yticks(y_pos, labels=[])
    ax_secondary.set_xlabel(secondary_label)
    ax_secondary.set_title(secondary_label)
    ax_secondary.grid(axis="x", linestyle=":", alpha=0.35)

    max_secondary = max(secondary_values) if secondary_values else 0.0
    secondary_pad = 0.04 * max_secondary if max_secondary > 0.0 else 1.0
    ax_secondary.set_xlim(0.0, max_secondary + 3.0 * secondary_pad)
    for idx, value in enumerate(secondary_values):
        ax_secondary.text(
            value + secondary_pad,
            y_pos[idx],
            _secondary_label(value, secondary_mode),
            va="center",
            fontsize=8.5,
        )

    title_parts = [
        "BEACH Performance Profile",
        f"mpi={metadata.get('mpi_world_size', '?')}",
        f"omp={metadata.get('omp_max_threads', '?')}",
        f"detail={metadata.get('detail_enabled', '?')}",
    ]
    if profile_path is not None:
        title_parts.append(str(profile_path))
    fig.suptitle(" | ".join(title_parts), fontsize=12)
    return fig, axes


def _resolve_profile_path(value: str) -> Path:
    path = Path(value)
    if path.is_dir() or path.suffix == "":
        return path / "performance_profile.csv"
    return path


def _required_field(row: dict[str, str], key: str) -> str:
    value = row.get(key)
    if value is None or value.strip() == "":
        raise ValueError(f"Missing required column: {key}")
    return value.strip()


def _simulation_total(rows: list[dict[str, Any]], *, metric: str) -> float:
    for row in rows:
        if row["region"] == "simulation_total":
            return max(0.0, float(row[metric]))
    return max((float(row[metric]) for row in rows), default=0.0)


def _secondary_series(
    rows: list[dict[str, Any]],
    *,
    metric: str,
    mode: str,
    simulation_total: float,
) -> tuple[list[float], str]:
    if mode == "imbalance":
        return [float(row["imbalance_ratio"]) for row in rows], "imbalance_ratio"
    if mode == "calls":
        return [float(row["calls_sum"]) for row in rows], "calls_sum"
    scale = max(simulation_total, 1.0e-30)
    return [100.0 * float(row[metric]) / scale for row in rows], "% of simulation_total"


def _secondary_label(value: float, mode: str) -> str:
    if mode == "imbalance":
        return f"x{value:.2f}"
    if mode == "calls":
        return f"{value:.0f}"
    return f"{value:.1f}%"


def _region_color(region: str) -> str:
    if region == "program_total":
        return "#4c6ef5"
    if region == "simulation_total":
        return "#0f766e"
    if region.startswith("particle_"):
        return "#d97706"
    if region.startswith("write") or "history" in region:
        return "#9333ea"
    if "field" in region:
        return "#2563eb"
    if "mpi" in region:
        return "#dc2626"
    return "#475569"


def _print_top_regions(rows: list[dict[str, Any]], *, metric: str, top: int) -> None:
    sorted_rows = sorted(
        rows,
        key=lambda row: (float(row[metric]), row["region"]),
        reverse=True,
    )
    limit = top if top > 0 else min(8, len(sorted_rows))
    print(f"top_regions_by_{metric}=")
    for row in sorted_rows[:limit]:
        print(
            f"  {row['region']}: "
            f"{float(row[metric]):.6e}s "
            f"(calls={float(row['calls_sum']):.0f}, "
            f"imbalance={float(row['imbalance_ratio']):.3f})"
        )


if __name__ == "__main__":
    main()
