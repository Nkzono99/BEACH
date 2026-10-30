"""Area-weighted boxplots of charge or potential per mesh source."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Mapping

import numpy as np

from beach.analysis.context import RunContext
from beach.analysis.mesh import _triangle_areas
from beach.analysis.potential import compute_potential_mesh
from beach.io.types import RunResult, MeshSource

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def plot_mesh_source_boxplot(
    result: RunResult | object,
    *,
    quantity: str = "charge",
    step: int | None = -1,
    showfliers: bool = True,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
) -> tuple[Figure, Axes]:
    """Plot area-weighted boxplots per mesh source (heavy for potential not saved in mesh_potential.csv).

    Parameters
    ----------
    result : RunResult or Beach-like object
        Run result or object exposing ``result`` as ``RunResult``.
    quantity : {"charge", "potential"}, default "charge"
        Quantity used in boxplot values.
    step : int or None, default -1
        History batch step used for charge snapshot.
        ``None`` uses final charges from ``charges.csv``.
    showfliers : bool, default True
        Whether outlier markers are rendered.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.

    Raises
    ------
    ValueError
        If triangle geometry is unavailable or plot arguments are invalid.
    """

    import matplotlib.pyplot as plt

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    triangles = context.require_triangles()
    values = _boxplot_values_for_quantity(
        context,
        quantity=quantity,
        step=step,
        library_path=library_path,
    )
    mesh_ids = context.mesh_ids_or_default()
    areas = _triangle_areas(triangles)
    stats = _box_stats_by_mesh_source(
        values=values,
        weights=areas,
        mesh_ids=mesh_ids,
        mesh_sources=resolved.mesh.mesh_sources,
    )
    if len(stats) == 0:
        raise ValueError("no mesh elements are available for boxplot.")

    fig_width = max(7.0, 1.4 * len(stats))
    fig, ax = plt.subplots(figsize=(fig_width, 4.5))
    artists = ax.bxp(stats, showfliers=showfliers, patch_artist=True)
    for patch in artists["boxes"]:
        patch.set_facecolor("#9ecae1")
        patch.set_alpha(0.85)

    quantity_name, quantity_unit = _quantity_axis_labels(quantity)
    ax.set_xlabel("mesh source")
    ax.set_ylabel(quantity_unit)
    ax.set_title(f"Area-weighted {quantity_name} by mesh source")
    ax.tick_params(axis="x", rotation=24)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.45)
    ax._beach_box_stats = stats
    fig.tight_layout()
    return fig, ax


def _boxplot_values_for_quantity(
    context: RunContext,
    *,
    quantity: str,
    step: int | None,
    library_path: str | Path | None,
) -> np.ndarray:
    if quantity == "charge":
        return context.charges_at(step)
    if quantity == "potential":
        potential_result = replace(context.result, charges=context.charges_at(step))
        return compute_potential_mesh(
            RunContext(
                result=potential_result,
                requested_config_path=context.requested_config_path,
            ),
            library_path=library_path,
            source="auto" if step is None else "kernel",
        )
    raise ValueError("quantity must be one of {'charge', 'potential'}.")


def _box_stats_by_mesh_source(
    *,
    values: np.ndarray,
    weights: np.ndarray,
    mesh_ids: np.ndarray,
    mesh_sources: Mapping[int, MeshSource] | None,
) -> list[dict[str, object]]:
    unique_mesh_ids = np.unique(mesh_ids)
    stats: list[dict[str, object]] = []
    for raw_mesh_id in unique_mesh_ids:
        mesh_id = int(raw_mesh_id)
        mask = mesh_ids == mesh_id
        item = _weighted_box_stats(values[mask], weights[mask])
        item["label"] = _mesh_source_label(mesh_id, mesh_sources)
        stats.append(item)
    return stats


def _weighted_box_stats(values: np.ndarray, weights: np.ndarray) -> dict[str, object]:
    scalar_values = np.asarray(values, dtype=float)
    scalar_weights = np.asarray(weights, dtype=float)
    if scalar_values.size == 0:
        raise ValueError("weighted boxplot requires at least one value.")
    if scalar_values.shape != scalar_weights.shape:
        raise ValueError("values and weights must have the same shape.")
    if not np.all(np.isfinite(scalar_values)):
        raise ValueError("values must be finite.")
    if not np.all(np.isfinite(scalar_weights)):
        raise ValueError("weights must be finite.")
    if np.any(scalar_weights < 0.0):
        raise ValueError("weights must be >= 0.")

    positive_mask = scalar_weights > 0.0
    if not np.any(positive_mask):
        raise ValueError("at least one positive weight is required.")
    scalar_values = scalar_values[positive_mask]
    scalar_weights = scalar_weights[positive_mask]

    q1, med, q3 = _weighted_quantiles(scalar_values, scalar_weights, (0.25, 0.50, 0.75))
    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr
    upper_fence = q3 + 1.5 * iqr
    inside = (scalar_values >= lower_fence) & (scalar_values <= upper_fence)

    if np.any(inside):
        whislo = float(np.min(scalar_values[inside]))
        whishi = float(np.max(scalar_values[inside]))
    else:
        whislo = q1
        whishi = q3

    return {
        "whislo": whislo,
        "q1": q1,
        "med": med,
        "q3": q3,
        "whishi": whishi,
        "fliers": np.sort(scalar_values[~inside]),
    }


def _weighted_quantiles(
    values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...]
) -> tuple[float, ...]:
    """Lower weighted quantiles from one sort and one cumulative distribution."""

    cdf: np.ndarray | None = None
    sorted_values = values
    results: list[float] = []
    for q in quantiles:
        quantile = float(q)
        if quantile <= 0.0:
            results.append(float(np.min(values)))
            continue
        if quantile >= 1.0:
            results.append(float(np.max(values)))
            continue
        if cdf is None:
            order = np.argsort(values)
            sorted_values = values[order]
            sorted_weights = weights[order]
            total_weight = float(np.sum(sorted_weights))
            if total_weight <= 0.0:
                raise ValueError("weights must include at least one positive value.")
            cdf = np.cumsum(sorted_weights) / total_weight
        idx = int(np.searchsorted(cdf, quantile, side="left"))
        if idx >= sorted_values.size:
            idx = sorted_values.size - 1
        results.append(float(sorted_values[idx]))
    return tuple(results)


def _mesh_source_label(
    mesh_id: int, mesh_sources: Mapping[int, MeshSource] | None
) -> str:
    if mesh_sources is None or mesh_id not in mesh_sources:
        return f"id={mesh_id}"

    source = mesh_sources[mesh_id]
    source_kind = source.source_kind or "unknown"
    template_kind = source.template_kind or "n/a"
    surface_model = source.surface_model or "insulator"
    return f"id={mesh_id} ({source_kind}/{template_kind}/{surface_model}/eps={source.epsilon_r:g})"


def _quantity_axis_labels(quantity: str) -> tuple[str, str]:
    if quantity == "charge":
        return "charge", "charge [C]"
    if quantity == "potential":
        return "potential", "potential [V]"
    raise ValueError("quantity must be one of {'charge', 'potential'}.")
