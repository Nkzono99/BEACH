"""XY/YZ/XZ potential-slice figures."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Mapping

import numpy as np

from beach.analysis.context import RunContext
from beach.analysis.potential import _resolve_reference_point, compute_potential_slices
from beach.io.types import RunResult

from ._defaults import DEFAULT_POTENTIAL_CMAP
from ._mesh3d import _resolve_axis_unit

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from numpy.typing import NDArray

    from beach.analysis.field_kernel import FieldKernel


def plot_potential_slices(
    result: RunResult | object,
    *,
    box_min: Iterable[float],
    box_max: Iterable[float],
    grid_n: int = 200,
    xy_z: float | None = None,
    yz_x: float | None = None,
    xz_y: float | None = None,
    chunk_size: int = 2048,
    cmap: str = DEFAULT_POTENTIAL_CMAP,
    vmin: float | None = None,
    vmax: float | None = None,
    periodic2: Mapping[str, object] | None = None,
    reference_point: Iterable[float] | str | None = None,
    axis_unit: str = "m",
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    kernel: FieldKernel | None = None,
) -> tuple[Figure, NDArray[np.object_]]:
    """Plot XY/YZ/XZ potential slices in one figure (heavy: native kernel).

    Parameters
    ----------
    result : RunResult or Beach-like object
        Run result or object exposing ``result`` as ``RunResult``.
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
        :func:`beach.compute_potential_points`.
    reference_point : iterable of float, {"species1_injection_center"}, or None, default None
        基準電位を差し引く参照点。``None`` は絶対電位、``"species1_injection_center"``
        は species1 の注入面中心を基準にした電位差を描く。
    axis_unit : {"m", "um", "nm"}, default "m"
        Coordinate unit used for slice axes.
    config_path, library_path : str, pathlib.Path, or None, default None
        Explicit ``beach.toml`` and native kernel library.
    kernel : FieldKernel or None, default None
        Open kernel built from this run, used as given and not closed; see
        :func:`beach.compute_potential_points`. ``None`` builds one for this call.

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

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    axis_label, axis_scale = _resolve_axis_unit(axis_unit)
    resolved_reference = _resolve_reference_point(
        resolved,
        reference_point,
        context=context,
    )
    slices = compute_potential_slices(
        context,
        box_min=box_min,
        box_max=box_max,
        grid_n=grid_n,
        xy_z=xy_z,
        yz_x=yz_x,
        xz_y=xz_y,
        chunk_size=chunk_size,
        periodic2=periodic2,
        reference_point=resolved_reference,
        config_path=config_path,
        library_path=library_path,
        kernel=kernel,
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
            slc.u_values_m * axis_scale,
            slc.v_values_m * axis_scale,
            slc.potential_V,
            shading="auto",
            cmap=cmap,
            vmin=global_min,
            vmax=global_max,
        )
        ax.set_xlabel(f"{slc.axis_u} [{axis_label}]")
        ax.set_ylabel(f"{slc.axis_v} [{axis_label}]")
        ax.set_title(
            f"{slc.plane.upper()} "
            f"({slc.fixed_axis}={slc.fixed_value_m * axis_scale:.3g} {axis_label})"
        )
        ax.set_aspect("equal", adjustable="box")

    if mappable is not None:
        label = "potential [V]" if resolved_reference is None else "potential difference [V]"
        fig.colorbar(mappable, ax=axes, shrink=0.9, label=label)
    return fig, axes


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
