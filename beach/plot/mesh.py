"""3D mesh figures colored by charge density or potential."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Literal, Mapping

import numpy as np

from beach.analysis.context import RunContext, resolve_result
from beach.analysis.mesh import (
    _maybe_apply_periodic2_mesh,
    _maybe_replicate_periodic2,
    _surface_charge_density,
)
from beach.analysis.potential import _resolve_reference_point, compute_potential_mesh
from beach.io.types import RunResult

from ._defaults import (
    DEFAULT_CHARGE_CMAP,
    DEFAULT_MESH_VIEW_AZIM,
    DEFAULT_MESH_VIEW_ELEV,
    DEFAULT_POTENTIAL_CMAP,
)
from ._mesh3d import _plot_scalar_mesh

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d import Axes3D

    from beach.analysis.field_kernel import FieldKernel


def plot_charges(result: RunResult | object) -> tuple[Figure, Axes]:
    """Plot per-element charges as a bar chart (light: saved values only).

    Parameters
    ----------
    result : RunResult or Beach-like object
        Run result or object exposing ``result`` as ``RunResult``.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.
    """

    import matplotlib.pyplot as plt

    resolved = resolve_result(result)
    fig, ax = plt.subplots(figsize=(8, 3))
    x = np.arange(resolved.charges.size)
    ax.bar(x, resolved.charges)
    ax.set_xlabel("element index")
    ax.set_ylabel("charge [C]")
    ax.set_title(f"Fortran charge distribution: {resolved.directory}")
    fig.tight_layout()
    return fig, ax


def plot_charge_mesh(
    result: RunResult | object,
    *,
    cmap: str = DEFAULT_CHARGE_CMAP,
    view_elev: float = DEFAULT_MESH_VIEW_ELEV,
    view_azim: float = DEFAULT_MESH_VIEW_AZIM,
    periodic2: Mapping[str, object] | None = None,
    apply_periodic2_mesh: bool = False,
    periodic2_mesh_mode: str = "centroid",
    periodic2_repeat: int = 0,
    axis_unit: str = "m",
) -> tuple[Figure, Axes3D]:
    """Plot a 3D mesh colored by surface charge density (light: saved values only).

    Parameters
    ----------
    result : RunResult or Beach-like object
        Run result or object exposing ``result`` as ``RunResult``.
    cmap : str, default "coolwarm"
        Matplotlib colormap name.
    view_elev : float, default 24.0
        Elevation angle in degrees passed to ``Axes3D.view_init``.
    view_azim : float, default -58.0
        Azimuth angle in degrees passed to ``Axes3D.view_init``.
    periodic2 : mapping or None, default None
        Two-axis periodic setting used when ``apply_periodic2_mesh=True``.
        ``None`` の場合は出力ディレクトリ近傍の ``beach.toml`` から自動判定する。
    apply_periodic2_mesh : bool, default False
        ``True`` の場合、triangle 重心を周期セルへ wrap する平行移動を各 face に
        適用して描画する。
    periodic2_mesh_mode : {"centroid", "mesh"}, default "centroid"
        周期セルへ wrap する単位。``"centroid"`` は face ごと、``"mesh"`` は
        mesh_id ごとに丸ごと移動して粒子形状を保つ。
    periodic2_repeat : int, default 0
        周期イメージの複製数。``n`` を指定すると各周期軸方向に ``-n`` から ``+n``
        までの ``(2n+1)^2`` コピーを生成して描画する。``0`` は複製なし。
    axis_unit : {"m", "um", "nm"}, default "m"
        Coordinate unit used for 3D axes.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.

    Raises
    ------
    ValueError
        If triangle geometry is unavailable.
    """

    context = RunContext.from_value(result)
    resolved = context.result
    triangles = context.require_triangles()
    triangles = _maybe_apply_periodic2_mesh(
        resolved,
        triangles,
        periodic2=periodic2,
        apply_periodic2_mesh=apply_periodic2_mesh,
        periodic2_mesh_mode=periodic2_mesh_mode,
        context=context,
    )
    sigma = _surface_charge_density(resolved.charges, triangles)
    triangles, sigma = _maybe_replicate_periodic2(
        resolved,
        triangles,
        sigma,
        periodic2=periodic2,
        periodic2_repeat=periodic2_repeat,
        context=context,
    )
    return _plot_scalar_mesh(
        triangles,
        sigma,
        title=f"Surface charge density mesh: {_short_directory_label(resolved.directory)}",
        colorbar_label="surface charge density [C/m^2]",
        cmap=cmap,
        view_elev=view_elev,
        view_azim=view_azim,
        axis_unit=axis_unit,
    )


def plot_potential_mesh(
    result: RunResult | object,
    *,
    cmap: str = DEFAULT_POTENTIAL_CMAP,
    view_elev: float = DEFAULT_MESH_VIEW_ELEV,
    view_azim: float = DEFAULT_MESH_VIEW_AZIM,
    periodic2: Mapping[str, object] | None = None,
    apply_periodic2_mesh: bool = False,
    periodic2_mesh_mode: str = "centroid",
    periodic2_repeat: int = 0,
    reference_point: Iterable[float] | str | None = None,
    axis_unit: str = "m",
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    source: Literal["auto", "saved", "kernel"] = "auto",
    kernel: FieldKernel | None = None,
) -> tuple[Figure, Axes3D]:
    """Plot a 3D mesh colored by electric potential (heavy unless saved values are used).

    Parameters
    ----------
    result : RunResult or Beach-like object
        Run result or object exposing ``result`` as ``RunResult``.
    cmap : str, default "viridis"
        Matplotlib colormap name.
    view_elev : float, default 24.0
        Elevation angle in degrees passed to ``Axes3D.view_init``.
    view_azim : float, default -58.0
        Azimuth angle in degrees passed to ``Axes3D.view_init``.
    periodic2 : mapping or None, default None
        Two-axis periodic setting. See
        :func:`beach.compute_potential_mesh`.
    apply_periodic2_mesh : bool, default False
        ``True`` の場合、triangle 重心を周期セルへ wrap する平行移動を各 face に
        適用して描画する。``periodic2`` が ``None`` なら近傍の ``beach.toml`` から
        自動判定する。
    periodic2_mesh_mode : {"centroid", "mesh"}, default "centroid"
        周期セルへ wrap する単位。``"mesh"`` は粒子を割らずに表示したい場合に使う。
    periodic2_repeat : int, default 0
        周期イメージの複製数。``n`` を指定すると各周期軸方向に ``-n`` から ``+n``
        までの ``(2n+1)^2`` コピーを生成して描画する。``0`` は複製なし。
    reference_point : iterable of float, {"species1_injection_center"}, or None, default None
        基準電位を差し引く参照点。``None`` は絶対電位、``"species1_injection_center"``
        は species1 の注入面中心を基準にした電位差を描く。
    axis_unit : {"m", "um", "nm"}, default "m"
        Coordinate unit used for 3D axes.
    config_path, library_path : str, pathlib.Path, or None, default None
        Explicit ``beach.toml`` and native kernel library.
    source : {"auto", "saved", "kernel"}, default "auto"
        Saved ``mesh_potential.csv`` or kernel evaluation; see
        :func:`beach.compute_potential_mesh`.
    kernel : FieldKernel or None, default None
        Open kernel built from this run that evaluates the potential with its own
        plan and charges and is not closed. ``periodic2`` then only wraps and
        replicates the drawn mesh. ``None`` builds a kernel when one is needed.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.

    Raises
    ------
    ValueError
        If potential reconstruction arguments are invalid.
    """

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    triangles = context.require_triangles()
    plot_triangles = _maybe_apply_periodic2_mesh(
        resolved,
        triangles,
        periodic2=periodic2,
        apply_periodic2_mesh=apply_periodic2_mesh,
        periodic2_mesh_mode=periodic2_mesh_mode,
        context=context,
    )
    if reference_point is not None and not isinstance(reference_point, str):
        reference_point = tuple(reference_point)
    phi = compute_potential_mesh(
        context,
        periodic2=periodic2,
        reference_point=reference_point,
        config_path=config_path,
        library_path=library_path,
        source=source,
        kernel=kernel,
    )
    resolved_reference = _resolve_reference_point(
        resolved,
        reference_point,
        context=context,
    )
    plot_triangles, phi = _maybe_replicate_periodic2(
        resolved,
        plot_triangles,
        phi,
        periodic2=periodic2,
        periodic2_repeat=periodic2_repeat,
        context=context,
    )
    title = "Electric potential mesh" if resolved_reference is None else "Electric potential difference mesh"
    colorbar_label = "potential [V]" if resolved_reference is None else "potential difference [V]"
    return _plot_scalar_mesh(
        plot_triangles,
        phi,
        title=f"{title}: {_short_directory_label(resolved.directory)}",
        colorbar_label=colorbar_label,
        cmap=cmap,
        view_elev=view_elev,
        view_azim=view_azim,
        axis_unit=axis_unit,
    )


def _short_directory_label(directory: Path) -> str:
    parts = Path(directory).parts
    if len(parts) >= 2:
        return str(Path(*parts[-2:]))
    return str(directory)
