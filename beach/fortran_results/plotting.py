"""Plotting helpers for Fortran results."""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .coulomb import calc_coulomb
from .mesh import _plot_scalar_mesh, _surface_charge_density, _triangle_areas
from .potential import (
    _find_config_path_near_output,
    _load_toml,
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


def plot_coulomb_force_matrix(
    result: FortranRunResult | object,
    *,
    step: int | None = -1,
    component: str = "z",
    softening: float = 0.0,
    target_kinds: Iterable[str] | None = None,
    config_path: str | Path | None = None,
    cmap: str = "coolwarm",
    annotate: bool = True,
):
    """Plot an object-wise Coulomb-force matrix.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    step : int or None, default -1
        History batch step used for charge snapshot.
        ``None`` uses final charges from ``charges.csv``.
    component : {"x", "y", "z"}, default "z"
        Force component rendered in the matrix.
    softening : float, default 0.0
        Softening length in meters for Coulomb-force evaluation.
    target_kinds : iterable of str or None, default None
        Template kinds used as target objects. ``None`` means all objects.
    config_path : str, pathlib.Path, or None, default None
        Optional ``beach.toml`` path used for object labels/order.
        ``None`` auto-detects near ``result.directory`` when possible.
    cmap : str, default "coolwarm"
        Matplotlib colormap name.
    annotate : bool, default True
        Whether to draw scientific-notation values inside each cell.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.

    Raises
    ------
    ValueError
        If no suitable objects are available or arguments are invalid.
    """

    import matplotlib.pyplot as plt

    resolved = _resolve_result(result)
    component_idx, component_label = _force_component_info(component)
    object_specs = _resolve_coulomb_object_specs(resolved, config_path=config_path)

    resolved_target_kinds = _normalize_target_kinds(target_kinds)
    target_specs = [
        spec
        for spec in object_specs
        if resolved_target_kinds is None or spec["kind"] in resolved_target_kinds
    ]
    if len(target_specs) == 0:
        available_kinds = sorted({str(spec["kind"]) for spec in object_specs})
        raise ValueError(
            "target_kinds did not match any objects. "
            f"available={available_kinds}"
        )

    matrix = np.zeros((len(object_specs) + 1, len(target_specs)), dtype=float)
    source_labels = [str(spec["label"]) for spec in object_specs] + ["net"]
    target_labels = [str(spec["label"]) for spec in target_specs]

    for col_idx, target_spec in enumerate(target_specs):
        target_mesh_id = int(target_spec["mesh_id"])
        for row_idx, source_spec in enumerate(object_specs):
            source_mesh_id = int(source_spec["mesh_id"])
            if source_mesh_id == target_mesh_id:
                continue
            interaction = calc_coulomb(
                resolved,
                target=target_mesh_id,
                source=source_mesh_id,
                step=step,
                softening=softening,
            )
            matrix[row_idx, col_idx] = float(
                interaction.force_on_a_N[component_idx]
            )

        net_sources = [
            int(spec["mesh_id"])
            for spec in object_specs
            if int(spec["mesh_id"]) != target_mesh_id
        ]
        if net_sources:
            interaction = calc_coulomb(
                resolved,
                target=target_mesh_id,
                source=net_sources,
                step=step,
                softening=softening,
            )
            matrix[-1, col_idx] = float(interaction.force_on_a_N[component_idx])

    fig_width = max(6.0, 1.2 * len(target_specs) + 2.8)
    fig_height = max(4.8, 0.7 * len(source_labels) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    max_abs = float(np.max(np.abs(matrix)))
    if max_abs <= np.finfo(float).tiny:
        max_abs = 1.0

    image = ax.imshow(
        matrix,
        origin="lower",
        cmap=cmap,
        vmin=-max_abs,
        vmax=max_abs,
        aspect="auto",
    )
    ax.set_xticks(np.arange(len(target_labels)))
    ax.set_xticklabels(target_labels, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(source_labels)))
    ax.set_yticklabels(source_labels)
    ax.set_xlabel("target")
    ax.set_ylabel("source")
    ax.set_title(f"Coulomb force matrix ({component_label})")

    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(f"{component_label} [N]")

    if annotate:
        fontsize = max(5, min(8, int(12 - 0.35 * max(matrix.shape))))
        for row_idx in range(matrix.shape[0]):
            for col_idx in range(matrix.shape[1]):
                value = float(matrix[row_idx, col_idx])
                text_color = "white" if abs(value) > 0.45 * max_abs else "black"
                ax.text(
                    col_idx,
                    row_idx,
                    f"{value:.1e}",
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    color=text_color,
                )

    ax._beach_coulomb_matrix = {
        "matrix": matrix.copy(),
        "component": component_label,
        "target_labels": tuple(target_labels),
        "source_labels": tuple(source_labels),
        "target_mesh_ids": tuple(int(spec["mesh_id"]) for spec in target_specs),
        "source_mesh_ids": tuple(int(spec["mesh_id"]) for spec in object_specs),
        "step": step,
        "softening": float(softening),
    }
    fig.tight_layout()
    return fig, ax


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


def _force_component_info(component: str) -> tuple[int, str]:
    normalized = str(component).strip().lower()
    if normalized == "x":
        return 0, "Fx"
    if normalized == "y":
        return 1, "Fy"
    if normalized == "z":
        return 2, "Fz"
    raise ValueError("component must be one of {'x', 'y', 'z'}.")


def _normalize_target_kinds(
    target_kinds: Iterable[str] | None,
) -> set[str] | None:
    if target_kinds is None:
        return None

    normalized = {
        str(kind).strip().lower()
        for kind in target_kinds
        if str(kind).strip()
    }
    if len(normalized) == 0 or "all" in normalized:
        return None
    return normalized


def _resolve_coulomb_object_specs(
    result: FortranRunResult,
    *,
    config_path: str | Path | None,
) -> list[dict[str, object]]:
    mesh_ids = tuple(int(v) for v in np.unique(_mesh_ids_or_default(result)))
    config_specs = _object_specs_from_config(result, mesh_ids=mesh_ids, config_path=config_path)
    if config_specs is not None:
        return config_specs

    if result.mesh_sources is not None:
        kinds = [
            _mesh_kind_from_source(result.mesh_sources.get(mesh_id))
            for mesh_id in mesh_ids
        ]
        labels = _labels_from_kinds(kinds)
        return [
            {"mesh_id": mesh_id, "kind": kind, "label": label}
            for mesh_id, kind, label in zip(mesh_ids, kinds, labels)
        ]

    fallback_kinds = [f"mesh{mesh_id}" for mesh_id in mesh_ids]
    return [
        {"mesh_id": mesh_id, "kind": kind, "label": kind}
        for mesh_id, kind in zip(mesh_ids, fallback_kinds)
    ]


def _object_specs_from_config(
    result: FortranRunResult,
    *,
    mesh_ids: tuple[int, ...],
    config_path: str | Path | None,
) -> list[dict[str, object]] | None:
    path: Path | None
    if config_path is None:
        path = _find_config_path_near_output(result.directory)
    else:
        path = Path(config_path)
        if not path.exists():
            raise ValueError(f'config file is not found: "{path}".')

    if path is None:
        return None

    config = _load_toml(path)
    mesh_cfg = config.get("mesh")
    if not isinstance(mesh_cfg, Mapping):
        return None
    templates = mesh_cfg.get("templates")
    if not isinstance(templates, list):
        return None

    enabled_kinds: list[str] = []
    for template in templates:
        if not isinstance(template, Mapping):
            continue
        if not bool(template.get("enabled", True)):
            continue
        kind = str(template.get("kind", "mesh")).strip().lower() or "mesh"
        enabled_kinds.append(kind)

    if len(enabled_kinds) != len(mesh_ids):
        return None

    labels = _labels_from_kinds(enabled_kinds)
    return [
        {"mesh_id": mesh_id, "kind": kind, "label": label}
        for mesh_id, kind, label in zip(mesh_ids, enabled_kinds, labels)
    ]


def _labels_from_kinds(kinds: Iterable[str]) -> list[str]:
    normalized = [str(kind).strip().lower() or "mesh" for kind in kinds]
    totals = Counter(normalized)
    seen: Counter[str] = Counter()
    labels: list[str] = []
    for kind in normalized:
        seen[kind] += 1
        if totals[kind] == 1:
            labels.append(kind)
        else:
            labels.append(f"{kind}{seen[kind]}")
    return labels


def _mesh_kind_from_source(source: MeshSource | None) -> str:
    if source is None:
        return "mesh"
    template_kind = str(source.template_kind).strip().lower()
    if template_kind:
        return template_kind
    source_kind = str(source.source_kind).strip().lower()
    if source_kind:
        return source_kind
    return "mesh"


def _quantity_axis_labels(quantity: str) -> tuple[str, str]:
    if quantity == "charge":
        return "charge", "charge [C]"
    if quantity == "potential":
        return "potential", "potential [V]"
    raise ValueError("quantity must be one of {'charge', 'potential'}.")
