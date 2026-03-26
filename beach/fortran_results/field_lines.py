"""Electric field computation and 3D field-line tracing."""

from __future__ import annotations

from typing import Iterable, Mapping

import numpy as np

from .constants import K_COULOMB
from .mesh import _triangle_centers
from .selection import _require_triangles, _resolve_result
from .types import FortranRunResult


# ---------------------------------------------------------------------------
# Electric field at arbitrary points
# ---------------------------------------------------------------------------


def compute_electric_field_points(
    result: FortranRunResult | object,
    points: np.ndarray,
    *,
    softening: float | None = None,
    chunk_size: int = 2048,
    periodic2: Mapping[str, object] | None = None,
) -> np.ndarray:
    """Compute the electric field vector at arbitrary 3D points.

    The field is computed directly from the surface charges via Coulomb's law:
    ``E(r) = K * sum_j  q_j * (r - r_j) / |r - r_j|^3``.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    points : numpy.ndarray
        Sampling points with shape ``(n_points, 3)`` in meters.
    softening : float or None, default None
        Softening length in meters. ``None`` は ``sim.softening`` を自動参照する。
    chunk_size : int, default 2048
        Number of points processed per chunk.
    periodic2 : mapping or None, default None
        Two-axis periodic setting. ``None`` の場合は出力ディレクトリ近傍の
        ``beach.toml`` から自動判定。

    Returns
    -------
    numpy.ndarray
        Electric field vectors in V/m with shape ``(n_points, 3)``.
    """

    from .potential import (
        _auto_periodic2_from_result,
        _coerce_periodic2,
        _resolve_softening,
    )

    resolved = _resolve_result(result)
    resolved_softening = _resolve_softening(resolved, softening)
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0.")

    sample_points = np.asarray(points, dtype=float)
    if sample_points.ndim != 2 or sample_points.shape[1] != 3:
        raise ValueError("points must have shape (n_points, 3).")
    if sample_points.shape[0] == 0:
        return np.empty((0, 3), dtype=float)

    triangles = _require_triangles(resolved)
    centers = _triangle_centers(triangles)
    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)

    if periodic_cfg is None:
        return _efield_points_free(
            sample_points,
            centers,
            resolved.charges,
            softening=resolved_softening,
            chunk_size=chunk_size,
        )
    else:
        return _efield_points_periodic2(
            sample_points,
            centers,
            resolved.charges,
            softening=resolved_softening,
            chunk_size=chunk_size,
            periodic2=periodic_cfg,
        )


def _efield_points_free(
    points: np.ndarray,
    centers: np.ndarray,
    charges: np.ndarray,
    *,
    softening: float,
    chunk_size: int,
) -> np.ndarray:
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny
    efield = np.empty((points.shape[0], 3), dtype=float)
    for start in range(0, points.shape[0], chunk_size):
        stop = min(start + chunk_size, points.shape[0])
        # delta[i, j, :] = points[i] - centers[j]
        delta = points[start:stop, None, :] - centers[None, :, :]
        dist2 = np.sum(delta * delta, axis=2) + eps2  # (chunk, n_src)
        inv_r3 = 1.0 / (
            np.maximum(dist2, min_dist2)
            * np.sqrt(np.maximum(dist2, min_dist2))
        )  # (chunk, n_src)
        # coeff[i, j] = K * q_j / |r_ij|^3
        coeff = K_COULOMB * charges[None, :] * inv_r3  # (chunk, n_src)
        efield[start:stop] = np.einsum("ij,ijk->ik", coeff, delta)
    return efield


def _efield_points_periodic2(
    points: np.ndarray,
    centers: np.ndarray,
    charges: np.ndarray,
    *,
    softening: float,
    chunk_size: int,
    periodic2: tuple,
) -> np.ndarray:
    from .potential import _wrap_periodic2_points

    axes, lengths, origins, nimg, _far_correction, _alpha, _ewald_layers = periodic2
    axis1, axis2 = axes
    l1, l2 = lengths
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny
    efield = np.zeros((points.shape[0], 3), dtype=float)
    wrapped_points = _wrap_periodic2_points(
        points, axes=axes, lengths=lengths, origins=origins,
    )

    for ix in range(-nimg, nimg + 1):
        for iy in range(-nimg, nimg + 1):
            shifted = centers.copy()
            shifted[:, axis1] += float(ix) * l1
            shifted[:, axis2] += float(iy) * l2
            for start in range(0, points.shape[0], chunk_size):
                stop = min(start + chunk_size, points.shape[0])
                delta = wrapped_points[start:stop, None, :] - shifted[None, :, :]
                dist2 = np.sum(delta * delta, axis=2) + eps2
                inv_r3 = 1.0 / (
                    np.maximum(dist2, min_dist2)
                    * np.sqrt(np.maximum(dist2, min_dist2))
                )
                coeff = K_COULOMB * charges[None, :] * inv_r3
                efield[start:stop] += np.einsum("ij,ijk->ik", coeff, delta)

    return efield


# ---------------------------------------------------------------------------
# Field-line tracing (RK4)
# ---------------------------------------------------------------------------


def trace_field_lines(
    result: FortranRunResult | object,
    seed_points: np.ndarray,
    *,
    ds: float | None = None,
    max_steps: int = 500,
    softening: float | None = None,
    periodic2: Mapping[str, object] | None = None,
    direction: str = "both",
    box_min: Iterable[float] | None = None,
    box_max: Iterable[float] | None = None,
) -> list[np.ndarray]:
    """Trace electric field lines from seed points using RK4 integration.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    seed_points : numpy.ndarray
        Starting points with shape ``(n_seeds, 3)`` in meters.
    ds : float or None, default None
        Integration step size in meters. ``None`` は三角形メッシュの代表長さ
        (平均辺長の 0.5 倍) から自動設定する。
    max_steps : int, default 500
        Maximum number of integration steps per direction.
    softening : float or None, default None
        Softening length in meters. ``None`` は ``sim.softening`` を自動参照する。
    periodic2 : mapping or None, default None
        Two-axis periodic setting. ``None`` で自動判定。
    direction : {"both", "forward", "backward"}, default "both"
        ``"forward"`` は電場方向、``"backward"`` は逆方向、``"both"`` は両方。
    box_min : iterable of float or None, default None
        Bounding box lower corner. 力線がこの範囲外に出たら打ち切る。
    box_max : iterable of float or None, default None
        Bounding box upper corner.

    Returns
    -------
    list of numpy.ndarray
        Each element has shape ``(n_points_i, 3)`` representing one field line.
    """

    from .potential import (
        _auto_periodic2_from_result,
        _coerce_periodic2,
        _resolve_softening,
    )

    resolved = _resolve_result(result)
    resolved_softening = _resolve_softening(resolved, softening)

    seeds = np.asarray(seed_points, dtype=float)
    if seeds.ndim == 1 and seeds.shape[0] == 3:
        seeds = seeds.reshape(1, 3)
    if seeds.ndim != 2 or seeds.shape[1] != 3:
        raise ValueError("seed_points must have shape (n_seeds, 3).")

    triangles = _require_triangles(resolved)
    centers = _triangle_centers(triangles)
    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)

    if ds is None:
        edges = np.concatenate([
            np.linalg.norm(triangles[:, 1] - triangles[:, 0], axis=1),
            np.linalg.norm(triangles[:, 2] - triangles[:, 1], axis=1),
            np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1),
        ])
        ds = float(np.mean(edges)) * 0.5

    bb_min = np.asarray(box_min, dtype=float) if box_min is not None else None
    bb_max = np.asarray(box_max, dtype=float) if box_max is not None else None

    lines: list[np.ndarray] = []
    for idx in range(seeds.shape[0]):
        seed = seeds[idx]
        parts: list[np.ndarray] = []
        if direction in ("forward", "both"):
            fwd = _rk4_trace(
                seed,
                centers,
                resolved.charges,
                softening=resolved_softening,
                periodic2=periodic_cfg,
                ds=ds,
                max_steps=max_steps,
                sign=1.0,
                bb_min=bb_min,
                bb_max=bb_max,
            )
            parts.append(fwd)
        if direction in ("backward", "both"):
            bwd = _rk4_trace(
                seed,
                centers,
                resolved.charges,
                softening=resolved_softening,
                periodic2=periodic_cfg,
                ds=ds,
                max_steps=max_steps,
                sign=-1.0,
                bb_min=bb_min,
                bb_max=bb_max,
            )
            parts.append(bwd[::-1])

        if direction == "both" and len(parts) == 2:
            # backward (reversed) + forward, skip duplicate seed
            line = np.concatenate([parts[1], parts[0][1:]], axis=0)
        elif parts:
            line = parts[0]
        else:
            line = seed.reshape(1, 3)
        lines.append(line)

    return lines


def _rk4_trace(
    seed: np.ndarray,
    centers: np.ndarray,
    charges: np.ndarray,
    *,
    softening: float,
    periodic2: tuple | None,
    ds: float,
    max_steps: int,
    sign: float,
    bb_min: np.ndarray | None,
    bb_max: np.ndarray | None,
) -> np.ndarray:
    """Trace a single field line using 4th-order Runge-Kutta."""

    points_list = [seed.copy()]
    pos = seed.copy()

    for _ in range(max_steps):
        k1 = _efield_single_point(
            pos, centers, charges, softening=softening, periodic2=periodic2,
        )
        norm1 = np.linalg.norm(k1)
        if norm1 < 1e-30:
            break
        k1 = k1 / norm1

        k2 = _efield_single_point(
            pos + 0.5 * ds * sign * k1, centers, charges,
            softening=softening, periodic2=periodic2,
        )
        norm2 = np.linalg.norm(k2)
        if norm2 < 1e-30:
            break
        k2 = k2 / norm2

        k3 = _efield_single_point(
            pos + 0.5 * ds * sign * k2, centers, charges,
            softening=softening, periodic2=periodic2,
        )
        norm3 = np.linalg.norm(k3)
        if norm3 < 1e-30:
            break
        k3 = k3 / norm3

        k4 = _efield_single_point(
            pos + ds * sign * k3, centers, charges,
            softening=softening, periodic2=periodic2,
        )
        norm4 = np.linalg.norm(k4)
        if norm4 < 1e-30:
            break
        k4 = k4 / norm4

        pos = pos + (ds * sign / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        points_list.append(pos.copy())

        if bb_min is not None and np.any(pos < bb_min):
            break
        if bb_max is not None and np.any(pos > bb_max):
            break

    return np.array(points_list)


def _efield_single_point(
    point: np.ndarray,
    centers: np.ndarray,
    charges: np.ndarray,
    *,
    softening: float,
    periodic2: tuple | None,
) -> np.ndarray:
    """Compute the electric field vector at a single 3D point (internal)."""

    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny

    if periodic2 is None:
        delta = point - centers  # (n_src, 3)
        dist2 = np.sum(delta * delta, axis=1) + eps2
        inv_r3 = 1.0 / (
            np.maximum(dist2, min_dist2)
            * np.sqrt(np.maximum(dist2, min_dist2))
        )
        coeff = K_COULOMB * charges * inv_r3
        return np.sum(coeff[:, None] * delta, axis=0)

    from .potential import _wrap_periodic2_points

    axes, lengths, origins, nimg, _fc, _alpha, _el = periodic2
    axis1, axis2 = axes
    l1, l2 = lengths
    pt = _wrap_periodic2_points(
        point.reshape(1, 3), axes=axes, lengths=lengths, origins=origins,
    )[0]

    efield = np.zeros(3, dtype=float)
    for ix in range(-nimg, nimg + 1):
        for iy in range(-nimg, nimg + 1):
            shifted = centers.copy()
            shifted[:, axis1] += float(ix) * l1
            shifted[:, axis2] += float(iy) * l2
            delta = pt - shifted
            dist2 = np.sum(delta * delta, axis=1) + eps2
            inv_r3 = 1.0 / (
                np.maximum(dist2, min_dist2)
                * np.sqrt(np.maximum(dist2, min_dist2))
            )
            coeff = K_COULOMB * charges * inv_r3
            efield += np.sum(coeff[:, None] * delta, axis=0)

    return efield


# ---------------------------------------------------------------------------
# 3D field-line plotting
# ---------------------------------------------------------------------------


def plot_field_lines_3d(
    result: FortranRunResult | object,
    seed_points: np.ndarray,
    *,
    ds: float | None = None,
    max_steps: int = 500,
    softening: float | None = None,
    periodic2: Mapping[str, object] | None = None,
    direction: str = "both",
    box_min: Iterable[float] | None = None,
    box_max: Iterable[float] | None = None,
    show_mesh: bool = True,
    mesh_alpha: float = 0.25,
    mesh_cmap: str = "coolwarm",
    line_color: str | None = None,
    line_cmap: str = "plasma",
    line_width: float = 1.2,
    view_elev: float = 24.0,
    view_azim: float = -58.0,
    title: str = "Electric field lines",
    figsize: tuple[float, float] = (9, 7),
):
    """Plot 3D electric field lines with optional surface mesh overlay.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result.
    seed_points : numpy.ndarray
        Starting points with shape ``(n_seeds, 3)``.
    ds : float or None
        Integration step size.
    max_steps : int
        Maximum steps per direction.
    softening : float or None
        Softening length.
    periodic2 : mapping or None
        Periodic boundary setting.
    direction : {"both", "forward", "backward"}
        Tracing direction.
    box_min, box_max : iterable of float or None
        Bounding box for line termination.
    show_mesh : bool, default True
        Whether to overlay the triangle mesh.
    mesh_alpha : float, default 0.25
        Mesh face transparency.
    mesh_cmap : str, default "coolwarm"
        Colormap for mesh surface charge density.
    line_color : str or None, default None
        Fixed color for field lines. ``None`` で各力線を ``line_cmap`` で
        色分けする。
    line_cmap : str, default "plasma"
        Colormap used when ``line_color`` is ``None``.
    line_width : float, default 1.2
        Line width.
    view_elev : float, default 24.0
        Elevation angle.
    view_azim : float, default -58.0
        Azimuth angle.
    title : str, default "Electric field lines"
        Plot title.
    figsize : tuple of float, default (9, 7)
        Figure size.

    Returns
    -------
    tuple
        ``(figure, axes)`` from matplotlib.
    """

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    from .mesh import (
        _configure_mesh_axes,
        _surface_charge_density,
    )

    resolved = _resolve_result(result)
    lines = trace_field_lines(
        resolved,
        seed_points,
        ds=ds,
        max_steps=max_steps,
        softening=softening,
        periodic2=periodic2,
        direction=direction,
        box_min=box_min,
        box_max=box_max,
    )

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    triangles = _require_triangles(resolved)

    # Optionally draw the mesh
    if show_mesh:
        charges = resolved.charges
        scd = _surface_charge_density(charges, triangles)
        max_abs = float(np.max(np.abs(scd)))
        if max_abs <= np.finfo(float).tiny:
            max_abs = 1.0
        norm = plt.Normalize(vmin=-max_abs, vmax=max_abs)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=mesh_cmap)
        facecolors = sm.to_rgba(scd)
        facecolors[:, 3] = mesh_alpha
        mesh_coll = Poly3DCollection(
            triangles,
            facecolors=facecolors,
            edgecolor=(0.0, 0.0, 0.0, 0.15),
            linewidth=0.2,
        )
        ax.add_collection3d(mesh_coll)

    # Draw field lines
    n_lines = len(lines)
    if line_color is not None:
        colors = [line_color] * n_lines
    else:
        cmap = plt.get_cmap(line_cmap)
        colors = [cmap(i / max(n_lines - 1, 1)) for i in range(n_lines)]

    for line, c in zip(lines, colors):
        if line.shape[0] < 2:
            continue
        ax.plot(line[:, 0], line[:, 1], line[:, 2], color=c, linewidth=line_width)
        # Arrow at midpoint
        mid = line.shape[0] // 2
        if mid > 0 and mid < line.shape[0] - 1:
            ax.quiver(
                line[mid, 0], line[mid, 1], line[mid, 2],
                line[mid + 1, 0] - line[mid, 0],
                line[mid + 1, 1] - line[mid, 1],
                line[mid + 1, 2] - line[mid, 2],
                color=c,
                arrow_length_ratio=0.4,
                linewidth=line_width * 1.5,
            )

    # Draw seed points
    seeds = np.asarray(seed_points, dtype=float)
    if seeds.ndim == 1:
        seeds = seeds.reshape(1, 3)
    ax.scatter(seeds[:, 0], seeds[:, 1], seeds[:, 2],
               c="red", s=20, zorder=5, depthshade=False)

    _configure_mesh_axes(ax, triangles, view_elev=view_elev, view_azim=view_azim)
    ax.set_title(title)
    fig.tight_layout()
    return fig, ax
