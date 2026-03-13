"""Plotting helpers for Fortran results."""

from __future__ import annotations

from dataclasses import replace
from typing import Iterable, Mapping

import numpy as np

from .mesh import _plot_scalar_mesh, _surface_charge_density, _triangle_areas
from .potential import (
    _resolve_reference_point,
    compute_potential_mesh,
    compute_potential_slices,
)
from .selection import (
    _charges_for_step,
    _mesh_ids_or_default,
    _require_triangles,
    _resolve_result,
)
from .types import FortranRunResult, MeshSource


def plot_charges(result: FortranRunResult | object):
    """Plot per-element charges as a bar chart.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.
    """

    import matplotlib.pyplot as plt

    resolved = _resolve_result(result)
    fig, ax = plt.subplots(figsize=(8, 3))
    x = np.arange(resolved.charges.size)
    ax.bar(x, resolved.charges)
    ax.set_xlabel("element index")
    ax.set_ylabel("charge [C]")
    ax.set_title(f"Fortran charge distribution: {resolved.directory}")
    fig.tight_layout()
    return fig, ax


def plot_charge_mesh(result: FortranRunResult | object, *, cmap: str = "coolwarm"):
    """Plot a 3D mesh colored by surface charge density.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    cmap : str, default "coolwarm"
        Matplotlib colormap name.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.

    Raises
    ------
    ValueError
        If triangle geometry is unavailable.
    """

    resolved = _resolve_result(result)
    triangles = _require_triangles(resolved)
    sigma = _surface_charge_density(resolved.charges, triangles)
    return _plot_scalar_mesh(
        triangles,
        sigma,
        title=f"Surface charge density mesh: {resolved.directory}",
        colorbar_label="surface charge density [C/m^2]",
        cmap=cmap,
    )


def plot_potential_slices(
    result: FortranRunResult | object,
    *,
    box_min: Iterable[float],
    box_max: Iterable[float],
    grid_n: int = 200,
    xy_z: float | None = None,
    yz_x: float | None = None,
    xz_y: float | None = None,
    softening: float | None = None,
    chunk_size: int = 2048,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    periodic2: Mapping[str, object] | None = None,
    reference_point: Iterable[float] | str | None = "species1_injection_center",
):
    """Plot XY/YZ/XZ potential slices in one figure.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    box_min : iterable of float
        Lower simulation-box corner ``[x, y, z]`` in meters.
    box_max : iterable of float
        Upper simulation-box corner ``[x, y, z]`` in meters.
    grid_n : int, default 200
        Grid size per slice axis.
    xy_z : float or None, default None
        Z coordinate of XY slice. ``None`` uses box center.
    yz_x : float or None, default None
        X coordinate of YZ slice. ``None`` uses box center.
    xz_y : float or None, default None
        Y coordinate of XZ slice. ``None`` uses box center.
    softening : float or None, default None
        Softening length in meters.
    chunk_size : int, default 2048
        Sampling chunk size for potential computation.
    cmap : str, default "viridis"
        Matplotlib colormap name.
    vmin : float or None, default None
        Lower color limit. ``None`` selects automatically.
    vmax : float or None, default None
        Upper color limit. ``None`` selects automatically.
    periodic2 : mapping or None, default None
        Two-axis periodic setting. See
        :func:`beach.fortran_results.compute_potential_points`.
    reference_point : iterable of float, {"species1_injection_center"}, or None, default "species1_injection_center"
        基準電位を差し引く参照点。既定では species1 の注入面中心を使う。

    Returns
    -------
    tuple
        ``(figure, axes_array)`` from matplotlib.

    Raises
    ------
    ValueError
        If potential sampling or color-limit arguments are invalid.
    """

    import matplotlib.pyplot as plt

    resolved = _resolve_result(result)
    resolved_reference = _resolve_reference_point(resolved, reference_point)
    slices = compute_potential_slices(
        result,
        box_min=box_min,
        box_max=box_max,
        grid_n=grid_n,
        xy_z=xy_z,
        yz_x=yz_x,
        xz_y=xz_y,
        softening=softening,
        chunk_size=chunk_size,
        periodic2=periodic2,
        reference_point=reference_point,
    )

    ordered_planes = ("xy", "yz", "xz")
    data_min = min(float(np.min(slices[name].potential_V)) for name in ordered_planes)
    data_max = max(float(np.max(slices[name].potential_V)) for name in ordered_planes)
    global_min, global_max = _resolve_slice_color_limits(
        data_min=data_min,
        data_max=data_max,
        vmin=vmin,
        vmax=vmax,
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    mappable = None
    for ax, plane in zip(axes, ordered_planes):
        slc = slices[plane]
        mappable = ax.pcolormesh(
            slc.u_values_m,
            slc.v_values_m,
            slc.potential_V,
            shading="auto",
            cmap=cmap,
            vmin=global_min,
            vmax=global_max,
        )
        ax.set_xlabel(f"{slc.axis_u} [m]")
        ax.set_ylabel(f"{slc.axis_v} [m]")
        ax.set_title(
            f"{slc.plane.upper()} ({slc.fixed_axis}={slc.fixed_value_m:.3e} m)"
        )
        ax.set_aspect("equal", adjustable="box")

    if mappable is not None:
        label = "potential [V]" if resolved_reference is None else "potential difference [V]"
        fig.colorbar(mappable, ax=axes, shrink=0.9, label=label)
    return fig, axes


def plot_potential_mesh(
    result: FortranRunResult | object,
    *,
    softening: float | None = None,
    self_term: str = "auto",
    cmap: str = "viridis",
    periodic2: Mapping[str, object] | None = None,
    reference_point: Iterable[float] | str | None = "species1_injection_center",
):
    """Plot a 3D mesh colored by reconstructed electric potential.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    softening : float or None, default None
        Softening length in meters.
    self_term : {"auto", "area_equivalent", "exclude", "softened_point"}, default "auto"
        Self-interaction model for potential reconstruction.
    cmap : str, default "viridis"
        Matplotlib colormap name.
    periodic2 : mapping or None, default None
        Two-axis periodic setting. See
        :func:`beach.fortran_results.compute_potential_mesh`.
    reference_point : iterable of float, {"species1_injection_center"}, or None, default "species1_injection_center"
        基準電位を差し引く参照点。既定では species1 の注入面中心を使う。

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.

    Raises
    ------
    ValueError
        If potential reconstruction arguments are invalid.
    """

    resolved = _resolve_result(result)
    triangles = _require_triangles(resolved)
    resolved_reference = _resolve_reference_point(resolved, reference_point)
    phi = compute_potential_mesh(
        resolved,
        softening=softening,
        self_term=self_term,
        periodic2=periodic2,
        reference_point=reference_point,
    )
    title = "Electric potential mesh" if resolved_reference is None else "Electric potential difference mesh"
    colorbar_label = "potential [V]" if resolved_reference is None else "potential difference [V]"
    return _plot_scalar_mesh(
        triangles,
        phi,
        title=f"{title}: {resolved.directory}",
        colorbar_label=colorbar_label,
        cmap=cmap,
    )


def plot_mesh_source_boxplot(
    result: FortranRunResult | object,
    *,
    quantity: str = "charge",
    step: int | None = -1,
    softening: float | None = None,
    self_term: str = "auto",
    showfliers: bool = True,
):
    """Plot area-weighted boxplots per mesh source.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    quantity : {"charge", "potential"}, default "charge"
        Quantity used in boxplot values.
    step : int or None, default -1
        History batch step used for charge snapshot.
        ``None`` uses final charges from ``charges.csv``.
    softening : float or None, default None
        Softening length in meters (potential mode only).
    self_term : {"auto", "area_equivalent", "exclude", "softened_point"}, default "auto"
        Self-interaction model (potential mode only).
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

    resolved = _resolve_result(result)
    triangles = _require_triangles(resolved)
    values = _boxplot_values_for_quantity(
        resolved,
        quantity=quantity,
        step=step,
        softening=softening,
        self_term=self_term,
    )
    mesh_ids = _mesh_ids_or_default(resolved)
    areas = _triangle_areas(triangles)
    stats = _box_stats_by_mesh_source(
        values=values,
        weights=areas,
        mesh_ids=mesh_ids,
        mesh_sources=resolved.mesh_sources,
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


def _resolve_slice_color_limits(
    *,
    data_min: float,
    data_max: float,
    vmin: float | None,
    vmax: float | None,
) -> tuple[float, float]:
    if vmin is None and vmax is None and np.isclose(data_min, data_max):
        center = 0.5 * (data_min + data_max)
        width = max(abs(center) * 0.05, 1.0)
        return center - width, center + width

    lower = data_min if vmin is None else float(vmin)
    upper = data_max if vmax is None else float(vmax)
    if not np.isfinite(lower) or not np.isfinite(upper):
        raise ValueError("vmin/vmax must be finite values.")
    if lower >= upper:
        raise ValueError("vmin must be smaller than vmax.")
    return lower, upper


def _boxplot_values_for_quantity(
    result: FortranRunResult,
    *,
    quantity: str,
    step: int | None,
    softening: float | None,
    self_term: str,
) -> np.ndarray:
    if quantity == "charge":
        return _charges_for_step(result, step=step)
    if quantity == "potential":
        charges = _charges_for_step(result, step=step)
        potential_result = replace(result, charges=charges)
        return compute_potential_mesh(
            potential_result,
            softening=softening,
            self_term=self_term,
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

    q1 = _weighted_quantile(scalar_values, scalar_weights, 0.25)
    med = _weighted_quantile(scalar_values, scalar_weights, 0.50)
    q3 = _weighted_quantile(scalar_values, scalar_weights, 0.75)
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


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    quantile = float(q)
    if quantile <= 0.0:
        return float(np.min(values))
    if quantile >= 1.0:
        return float(np.max(values))

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
    return float(sorted_values[idx])


def _mesh_source_label(
    mesh_id: int, mesh_sources: Mapping[int, MeshSource] | None
) -> str:
    if mesh_sources is None or mesh_id not in mesh_sources:
        return f"id={mesh_id}"

    source = mesh_sources[mesh_id]
    source_kind = source.source_kind or "unknown"
    template_kind = source.template_kind or "n/a"
    return f"id={mesh_id} ({source_kind}/{template_kind})"


def _quantity_axis_labels(quantity: str) -> tuple[str, str]:
    if quantity == "charge":
        return "charge", "charge [C]"
    if quantity == "potential":
        return "potential", "potential [V]"
    raise ValueError("quantity must be one of {'charge', 'potential'}.")
