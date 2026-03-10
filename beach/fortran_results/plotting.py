"""Plotting helpers for Fortran results."""

from __future__ import annotations

from typing import Iterable

import numpy as np

from .mesh import _plot_scalar_mesh, _surface_charge_density
from .potential import compute_potential_mesh, compute_potential_slices
from .selection import _require_triangles, _resolve_result
from .types import FortranRunResult


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
    softening: float = 0.0,
    chunk_size: int = 2048,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
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
    softening : float, default 0.0
        Softening length in meters.
    chunk_size : int, default 2048
        Sampling chunk size for potential computation.
    cmap : str, default "viridis"
        Matplotlib colormap name.
    vmin : float or None, default None
        Lower color limit. ``None`` selects automatically.
    vmax : float or None, default None
        Upper color limit. ``None`` selects automatically.

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
        fig.colorbar(mappable, ax=axes, shrink=0.9, label="potential [V]")
    return fig, axes


def plot_potential_mesh(
    result: FortranRunResult | object,
    *,
    softening: float = 0.0,
    self_term: str = "area_equivalent",
    cmap: str = "viridis",
):
    """Plot a 3D mesh colored by reconstructed electric potential.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    softening : float, default 0.0
        Softening length in meters.
    self_term : {"area_equivalent", "exclude", "softened_point"}, default "area_equivalent"
        Self-interaction model for potential reconstruction.
    cmap : str, default "viridis"
        Matplotlib colormap name.

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
    phi = compute_potential_mesh(resolved, softening=softening, self_term=self_term)
    return _plot_scalar_mesh(
        triangles,
        phi,
        title=f"Electric potential mesh: {resolved.directory}",
        colorbar_label="potential [V]",
        cmap=cmap,
    )


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
