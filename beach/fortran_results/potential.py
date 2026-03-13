"""Potential computation helpers for mesh and sampled points."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .constants import K_COULOMB
from .mesh import _triangle_areas, _triangle_centers
from .selection import _require_triangles, _resolve_result
from .types import FortranRunResult, PotentialSlice2D


def compute_potential_mesh(
    result: FortranRunResult | object,
    *,
    softening: float | None = None,
    self_term: str = "auto",
    periodic2: Mapping[str, object] | None = None,
) -> np.ndarray:
    """Compute electric potential at each triangle centroid.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    softening : float or None, default None
        Softening length in meters for off-diagonal interactions. ``None`` は
        ``sim.softening`` を自動参照し、見つからなければ ``0.0`` を使う。
    self_term : {"auto", "area_equivalent", "exclude", "softened_point"}, default "auto"
        Self-interaction model at each centroid. ``"auto"`` は softening が
        正なら ``"softened_point"``, そうでなければ ``"area_equivalent"`` を選ぶ。
    periodic2 : mapping or None, default None
        Two-axis periodic setting for potential reconstruction. Use a mapping with
        keys ``axes`` (length-2, 0-based axis indices), ``lengths`` (length-2,
        positive box lengths), and optional ``origins`` (length-2 periodic-box
        origins) or ``box_min`` (length-3), ``image_layers`` (int, default 1),
        ``far_correction`` (``"none"`` or ``"ewald_like"``), ``ewald_alpha``
        (float, default 0.0 for auto), ``ewald_layers`` (int, default 4).
        If ``None``, ``result.directory`` 近傍の ``beach.toml`` を探索し、
        ``sim.field_bc_mode="periodic2"`` なら自動適用する。

    Returns
    -------
    numpy.ndarray
        Potential values in volts with shape ``(mesh_nelem,)``.

    Raises
    ------
    ValueError
        If softening is negative, triangles are unavailable, or ``self_term`` is invalid.
    """

    resolved = _resolve_result(result)
    resolved_softening = _resolve_softening(resolved, softening)
    self_term_key = _resolve_self_term(self_term, resolved_softening)

    triangles = _require_triangles(resolved)
    centers = _triangle_centers(triangles)
    self_coeff = _self_potential_coefficients(
        triangles, self_term=self_term_key, softening=resolved_softening
    )
    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)
    if periodic_cfg is None:
        offdiag_kernel = _potential_offdiag_kernel(centers, softening=resolved_softening)
        potential = offdiag_kernel @ resolved.charges + self_coeff * resolved.charges
    else:
        potential = _compute_potential_mesh_periodic2(
            centers,
            resolved.charges,
            softening=resolved_softening,
            self_coeff=self_coeff,
            periodic2=periodic_cfg,
        )
    return K_COULOMB * potential


def compute_potential_points(
    result: FortranRunResult | object,
    points: np.ndarray,
    *,
    softening: float | None = None,
    chunk_size: int = 2048,
    periodic2: Mapping[str, object] | None = None,
) -> np.ndarray:
    """Compute electric potential at arbitrary 3D sampling points.

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
        Two-axis periodic setting for potential sampling. Use a mapping with
        keys ``axes`` (length-2, 0-based axis indices), ``lengths`` (length-2,
        positive box lengths), and optional ``origins`` (length-2 periodic-box
        origins) or ``box_min`` (length-3), ``image_layers`` (int, default 1),
        ``far_correction`` (``"none"`` or ``"ewald_like"``), ``ewald_alpha``
        (float, default 0.0 for auto), ``ewald_layers`` (int, default 4).
        If ``None``, ``result.directory`` 近傍の ``beach.toml`` から自動判定する。

    Returns
    -------
    numpy.ndarray
        Potential values in volts with shape ``(n_points,)``.

    Raises
    ------
    ValueError
        If argument shapes/ranges are invalid.
    """

    resolved = _resolve_result(result)
    resolved_softening = _resolve_softening(resolved, softening)
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0.")

    sample_points = np.asarray(points, dtype=float)
    if sample_points.ndim != 2 or sample_points.shape[1] != 3:
        raise ValueError("points must have shape (n_points, 3).")
    if sample_points.shape[0] == 0:
        return np.empty(0, dtype=float)

    triangles = _require_triangles(resolved)
    centers = _triangle_centers(triangles)
    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)
    if periodic_cfg is None:
        potential = _compute_potential_points_free(
            sample_points,
            centers,
            resolved.charges,
            softening=resolved_softening,
            chunk_size=chunk_size,
        )
    else:
        potential = _compute_potential_points_periodic2(
            sample_points,
            centers,
            resolved.charges,
            softening=resolved_softening,
            chunk_size=chunk_size,
            periodic2=periodic_cfg,
        )

    return K_COULOMB * potential


def compute_potential_slices(
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
    periodic2: Mapping[str, object] | None = None,
) -> dict[str, PotentialSlice2D]:
    """Sample potential on XY/YZ/XZ slices inside a simulation box.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    box_min : iterable of float
        Lower simulation-box corner ``[x, y, z]`` in meters.
    box_max : iterable of float
        Upper simulation-box corner ``[x, y, z]`` in meters.
    grid_n : int, default 200
        Grid size per slice axis (``grid_n x grid_n``).
    xy_z : float or None, default None
        Z coordinate for the XY slice. ``None`` uses box center.
    yz_x : float or None, default None
        X coordinate for the YZ slice. ``None`` uses box center.
    xz_y : float or None, default None
        Y coordinate for the XZ slice. ``None`` uses box center.
    softening : float or None, default None
        Softening length in meters. ``None`` は ``sim.softening`` を自動参照する。
    chunk_size : int, default 2048
        Number of sampling points processed per chunk.
    periodic2 : mapping or None, default None
        Two-axis periodic setting for potential sampling. See
        :func:`compute_potential_points` for supported keys. ``None`` 時は
        ``result.directory`` 近傍の ``beach.toml`` から自動判定する。

    Returns
    -------
    dict[str, PotentialSlice2D]
        Slice data keyed by ``"xy"``, ``"yz"``, and ``"xz"``.

    Raises
    ------
    ValueError
        If box bounds, grid settings, or slice coordinates are invalid.
    """

    if grid_n < 2:
        raise ValueError("grid_n must be >= 2.")

    resolved = _resolve_result(result)
    resolved_softening = _resolve_softening(resolved, softening)
    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)
    min_corner, max_corner = _coerce_box_bounds(box_min, box_max)
    x = np.linspace(min_corner[0], max_corner[0], grid_n, dtype=float)
    y = np.linspace(min_corner[1], max_corner[1], grid_n, dtype=float)
    z = np.linspace(min_corner[2], max_corner[2], grid_n, dtype=float)

    z_xy = _coerce_slice_coordinate(
        xy_z,
        lower=min_corner[2],
        upper=max_corner[2],
        name="xy_z",
    )
    x_yz = _coerce_slice_coordinate(
        yz_x,
        lower=min_corner[0],
        upper=max_corner[0],
        name="yz_x",
    )
    y_xz = _coerce_slice_coordinate(
        xz_y,
        lower=min_corner[1],
        upper=max_corner[1],
        name="xz_y",
    )

    xx, yy = np.meshgrid(x, y, indexing="xy")
    points_xy = np.column_stack((xx.reshape(-1), yy.reshape(-1), np.full(xx.size, z_xy)))
    potential_xy = compute_potential_points(
        resolved,
        points_xy,
        softening=resolved_softening,
        chunk_size=chunk_size,
        periodic2=periodic_cfg,
    ).reshape(grid_n, grid_n)

    yy2, zz = np.meshgrid(y, z, indexing="xy")
    points_yz = np.column_stack((np.full(yy2.size, x_yz), yy2.reshape(-1), zz.reshape(-1)))
    potential_yz = compute_potential_points(
        resolved,
        points_yz,
        softening=resolved_softening,
        chunk_size=chunk_size,
        periodic2=periodic_cfg,
    ).reshape(grid_n, grid_n)

    xx2, zz2 = np.meshgrid(x, z, indexing="xy")
    points_xz = np.column_stack((xx2.reshape(-1), np.full(xx2.size, y_xz), zz2.reshape(-1)))
    potential_xz = compute_potential_points(
        resolved,
        points_xz,
        softening=resolved_softening,
        chunk_size=chunk_size,
        periodic2=periodic_cfg,
    ).reshape(grid_n, grid_n)

    return {
        "xy": PotentialSlice2D(
            plane="xy",
            axis_u="x",
            axis_v="y",
            fixed_axis="z",
            fixed_value_m=float(z_xy),
            u_values_m=x,
            v_values_m=y,
            potential_V=potential_xy,
        ),
        "yz": PotentialSlice2D(
            plane="yz",
            axis_u="y",
            axis_v="z",
            fixed_axis="x",
            fixed_value_m=float(x_yz),
            u_values_m=y,
            v_values_m=z,
            potential_V=potential_yz,
        ),
        "xz": PotentialSlice2D(
            plane="xz",
            axis_u="x",
            axis_v="z",
            fixed_axis="y",
            fixed_value_m=float(y_xz),
            u_values_m=x,
            v_values_m=z,
            potential_V=potential_xz,
        ),
    }


def _coerce_box_bounds(
    box_min: Iterable[float], box_max: Iterable[float]
) -> tuple[np.ndarray, np.ndarray]:
    min_corner = np.asarray(list(box_min), dtype=float)
    max_corner = np.asarray(list(box_max), dtype=float)
    if min_corner.shape != (3,) or max_corner.shape != (3,):
        raise ValueError("box_min and box_max must contain exactly 3 values.")
    if np.any(max_corner <= min_corner):
        raise ValueError("box_max must be greater than box_min on all axes.")
    return min_corner, max_corner


def _coerce_slice_coordinate(
    value: float | None, *, lower: float, upper: float, name: str
) -> float:
    if value is None:
        return 0.5 * (lower + upper)
    coordinate = float(value)
    if coordinate < lower or coordinate > upper:
        raise ValueError(
            f"{name} must be within [{lower:.6e}, {upper:.6e}] but got {coordinate:.6e}."
        )
    return coordinate


def _potential_offdiag_kernel(centers: np.ndarray, *, softening: float) -> np.ndarray:
    delta = centers[:, None, :] - centers[None, :, :]
    dist2 = np.sum(delta * delta, axis=2) + softening * softening
    min_dist2 = np.finfo(float).tiny
    inv_r = 1.0 / np.sqrt(np.maximum(dist2, min_dist2))
    np.fill_diagonal(inv_r, 0.0)
    return inv_r


def _potential_history(
    charges_history: np.ndarray,
    triangles: np.ndarray,
    *,
    softening: float,
    self_term: str,
    periodic2: tuple[
        tuple[int, int],
        tuple[float, float],
        tuple[float, float],
        int,
        str,
        float,
        int,
    ]
    | None = None,
) -> np.ndarray:
    if softening < 0.0:
        raise ValueError("softening must be >= 0.")

    centers = _triangle_centers(triangles)
    self_coeff = _self_potential_coefficients(
        triangles,
        self_term=self_term,
        softening=softening,
    )
    if periodic2 is None:
        offdiag_kernel = _potential_offdiag_kernel(centers, softening=softening)
        potential = offdiag_kernel @ charges_history + self_coeff[:, None] * charges_history
    else:
        potential = np.empty_like(charges_history, dtype=float)
        for frame_idx in range(charges_history.shape[1]):
            potential[:, frame_idx] = _compute_potential_mesh_periodic2(
                centers,
                charges_history[:, frame_idx],
                softening=softening,
                self_coeff=self_coeff,
                periodic2=periodic2,
            )
    return K_COULOMB * potential


def _compute_potential_points_free(
    points: np.ndarray,
    centers: np.ndarray,
    charges: np.ndarray,
    *,
    softening: float,
    chunk_size: int,
) -> np.ndarray:
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny
    potential = np.empty(points.shape[0], dtype=float)
    for start in range(0, points.shape[0], chunk_size):
        stop = min(start + chunk_size, points.shape[0])
        delta = points[start:stop, None, :] - centers[None, :, :]
        dist2 = np.sum(delta * delta, axis=2) + eps2
        inv_r = 1.0 / np.sqrt(np.maximum(dist2, min_dist2))
        potential[start:stop] = inv_r @ charges
    return potential


def _compute_potential_points_periodic2(
    points: np.ndarray,
    centers: np.ndarray,
    charges: np.ndarray,
    *,
    softening: float,
    chunk_size: int,
    periodic2: tuple[
        tuple[int, int],
        tuple[float, float],
        tuple[float, float],
        int,
        str,
        float,
        int,
    ],
) -> np.ndarray:
    axes, lengths, origins, nimg, far_correction, alpha, ewald_layers = periodic2
    axis1, axis2 = axes
    l1, l2 = lengths
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny
    potential = np.zeros(points.shape[0], dtype=float)
    wrapped_points = _wrap_periodic2_points(
        points,
        axes=axes,
        lengths=lengths,
        origins=origins,
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
                inv_r = 1.0 / np.sqrt(np.maximum(dist2, min_dist2))
                potential[start:stop] += inv_r @ charges

    if far_correction == "ewald_like":
        img_outer = nimg + ewald_layers
        for ix in range(-img_outer, img_outer + 1):
            for iy in range(-img_outer, img_outer + 1):
                if abs(ix) <= nimg and abs(iy) <= nimg:
                    continue
                shifted = centers.copy()
                shifted[:, axis1] += float(ix) * l1
                shifted[:, axis2] += float(iy) * l2
                for start in range(0, points.shape[0], chunk_size):
                    stop = min(start + chunk_size, points.shape[0])
                    delta = wrapped_points[start:stop, None, :] - shifted[None, :, :]
                    dist2 = np.sum(delta * delta, axis=2) + eps2
                    r = np.sqrt(np.maximum(dist2, min_dist2))
                    kernel = np.erfc(alpha * r) / r
                    potential[start:stop] += kernel @ charges

    return potential


def _compute_potential_mesh_periodic2(
    centers: np.ndarray,
    charges: np.ndarray,
    *,
    softening: float,
    self_coeff: np.ndarray,
    periodic2: tuple[
        tuple[int, int],
        tuple[float, float],
        tuple[float, float],
        int,
        str,
        float,
        int,
    ],
) -> np.ndarray:
    axes, lengths, origins, nimg, far_correction, alpha, ewald_layers = periodic2
    axis1, axis2 = axes
    l1, l2 = lengths
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny
    target_centers = _wrap_periodic2_points(
        centers,
        axes=axes,
        lengths=lengths,
        origins=origins,
    )
    potential = self_coeff * charges

    # Central-cell off-diagonal interaction.
    delta0 = target_centers[:, None, :] - centers[None, :, :]
    dist20 = np.sum(delta0 * delta0, axis=2) + eps2
    inv_r0 = 1.0 / np.sqrt(np.maximum(dist20, min_dist2))
    np.fill_diagonal(inv_r0, 0.0)
    potential += inv_r0 @ charges

    for ix in range(-nimg, nimg + 1):
        for iy in range(-nimg, nimg + 1):
            if ix == 0 and iy == 0:
                continue
            shifted = centers.copy()
            shifted[:, axis1] += float(ix) * l1
            shifted[:, axis2] += float(iy) * l2
            delta = target_centers[:, None, :] - shifted[None, :, :]
            dist2 = np.sum(delta * delta, axis=2) + eps2
            inv_r = 1.0 / np.sqrt(np.maximum(dist2, min_dist2))
            potential += inv_r @ charges

    if far_correction == "ewald_like":
        img_outer = nimg + ewald_layers
        for ix in range(-img_outer, img_outer + 1):
            for iy in range(-img_outer, img_outer + 1):
                if abs(ix) <= nimg and abs(iy) <= nimg:
                    continue
                shifted = centers.copy()
                shifted[:, axis1] += float(ix) * l1
                shifted[:, axis2] += float(iy) * l2
                delta = target_centers[:, None, :] - shifted[None, :, :]
                dist2 = np.sum(delta * delta, axis=2) + eps2
                r = np.sqrt(np.maximum(dist2, min_dist2))
                kernel = np.erfc(alpha * r) / r
                potential += kernel @ charges

    return potential


def _coerce_periodic2(
    periodic2: Mapping[str, object] | None,
) -> tuple[
    tuple[int, int],
    tuple[float, float],
    tuple[float, float],
    int,
    str,
    float,
    int,
] | None:
    if periodic2 is None:
        return None
    if not isinstance(periodic2, Mapping):
        raise ValueError("periodic2 must be a mapping or None.")

    if "axes" not in periodic2 or "lengths" not in periodic2:
        raise ValueError('periodic2 requires "axes" and "lengths".')
    axes_obj = periodic2["axes"]
    lengths_obj = periodic2["lengths"]
    if not isinstance(axes_obj, (list, tuple)) or len(axes_obj) != 2:
        raise ValueError("periodic2.axes must be a length-2 sequence.")
    if not isinstance(lengths_obj, (list, tuple)) or len(lengths_obj) != 2:
        raise ValueError("periodic2.lengths must be a length-2 sequence.")

    axes = (int(axes_obj[0]), int(axes_obj[1]))
    if axes[0] == axes[1] or any(axis < 0 or axis > 2 for axis in axes):
        raise ValueError("periodic2.axes must contain two distinct axis indices in {0,1,2}.")

    lengths = (float(lengths_obj[0]), float(lengths_obj[1]))
    if lengths[0] <= 0.0 or lengths[1] <= 0.0:
        raise ValueError("periodic2.lengths must be positive.")

    if "origins" in periodic2:
        origins_obj = periodic2["origins"]
        if not isinstance(origins_obj, (list, tuple)) or len(origins_obj) != 2:
            raise ValueError("periodic2.origins must be a length-2 sequence.")
        origins = (float(origins_obj[0]), float(origins_obj[1]))
    elif "box_min" in periodic2:
        box_min = _coerce_vec3(periodic2["box_min"], name="periodic2.box_min")
        origins = (box_min[axes[0]], box_min[axes[1]])
    else:
        origins = (0.0, 0.0)

    nimg = int(periodic2.get("image_layers", 1))
    if nimg < 0:
        raise ValueError("periodic2.image_layers must be >= 0.")

    far_correction = str(periodic2.get("far_correction", "none")).strip().lower()
    if far_correction not in {"none", "ewald_like"}:
        raise ValueError('periodic2.far_correction must be "none" or "ewald_like".')

    alpha = float(periodic2.get("ewald_alpha", 0.0))
    if (not math.isfinite(alpha)) or alpha < 0.0:
        raise ValueError("periodic2.ewald_alpha must be finite and >= 0.")
    if alpha <= 0.0 and far_correction == "ewald_like":
        alpha = 1.2 / (float(nimg + 1) * min(lengths))

    ewald_layers = int(periodic2.get("ewald_layers", 4))
    if ewald_layers < 0:
        raise ValueError("periodic2.ewald_layers must be >= 0.")
    if far_correction == "ewald_like" and ewald_layers < 1:
        raise ValueError("periodic2.ewald_layers must be >= 1 for ewald_like.")

    return axes, lengths, origins, nimg, far_correction, alpha, ewald_layers


def _auto_periodic2_from_result(
    resolved: FortranRunResult,
) -> tuple[
    tuple[int, int],
    tuple[float, float],
    tuple[float, float],
    int,
    str,
    float,
    int,
] | None:
    sim = _load_sim_near_output(resolved.directory)
    if sim is None:
        return None
    periodic2 = _periodic2_from_sim(sim)
    return _coerce_periodic2(periodic2)


def _resolve_softening(resolved: FortranRunResult, softening: float | None) -> float:
    if softening is None:
        sim = _load_sim_near_output(resolved.directory)
        raw = 0.0 if sim is None else sim.get("softening", 0.0)
    else:
        raw = softening

    value = float(raw)
    if (not math.isfinite(value)) or value < 0.0:
        raise ValueError("softening must be finite and >= 0.")
    return value


def _resolve_self_term(self_term: str, softening: float) -> str:
    normalized = str(self_term).strip().lower().replace("-", "_")
    if normalized == "auto":
        return "softened_point" if softening > 0.0 else "area_equivalent"
    if normalized in {"area_equivalent", "exclude", "softened_point"}:
        return normalized
    raise ValueError(
        "self_term must be one of {'auto', 'area_equivalent', 'exclude', 'softened_point'}."
    )


def _load_sim_near_output(output_dir: Path) -> Mapping[str, object] | None:
    config_path = _find_config_path_near_output(output_dir)
    if config_path is None:
        return None
    config = _load_toml(config_path)
    sim = config.get("sim")
    return sim if isinstance(sim, Mapping) else None


def _wrap_periodic2_points(
    points: np.ndarray,
    *,
    axes: tuple[int, int],
    lengths: tuple[float, float],
    origins: tuple[float, float],
) -> np.ndarray:
    wrapped = np.asarray(points, dtype=float).copy()
    wrapped[:, axes[0]] = origins[0] + np.mod(
        wrapped[:, axes[0]] - origins[0], lengths[0]
    )
    wrapped[:, axes[1]] = origins[1] + np.mod(
        wrapped[:, axes[1]] - origins[1], lengths[1]
    )
    return wrapped


def _find_config_path_near_output(output_dir: Path) -> Path | None:
    candidates = (
        output_dir / "beach.toml",
        output_dir.parent / "beach.toml",
        output_dir.parent.parent / "beach.toml",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _load_toml(path: Path) -> dict[str, object]:
    try:
        import tomllib  # py311+

        with path.open("rb") as stream:
            return tomllib.load(stream)
    except ModuleNotFoundError:
        try:
            import tomli  # type: ignore

            with path.open("rb") as stream:
                return tomli.load(stream)
        except ModuleNotFoundError as exc:
            raise ValueError(
                "TOML parser is missing. Use Python 3.11+ or install tomli: "
                "`python -m pip install tomli`."
            ) from exc


def _periodic2_from_sim(sim: Mapping[str, object]) -> dict[str, object] | None:
    field_bc_mode = str(sim.get("field_bc_mode", "free")).strip().lower()
    if field_bc_mode != "periodic2":
        return None

    if "box_min" not in sim or "box_max" not in sim:
        raise ValueError('periodic2 potential requires "sim.box_min" and "sim.box_max".')
    box_min = _coerce_vec3(sim["box_min"], name="sim.box_min")
    box_max = _coerce_vec3(sim["box_max"], name="sim.box_max")

    periodic_axes: list[int] = []
    for axis_idx, axis_name in enumerate(("x", "y", "z")):
        low = _canonical_boundary_mode(sim.get(f"bc_{axis_name}_low", "open"))
        high = _canonical_boundary_mode(sim.get(f"bc_{axis_name}_high", "open"))
        if (low == "periodic") != (high == "periodic"):
            raise ValueError("periodic2 requires bc_low(axis)=bc_high(axis)=periodic.")
        if low == "periodic":
            periodic_axes.append(axis_idx)

    if len(periodic_axes) != 2:
        raise ValueError('sim.field_bc_mode="periodic2" requires exactly two periodic axes.')

    lengths = [box_max[axis] - box_min[axis] for axis in periodic_axes]
    if lengths[0] <= 0.0 or lengths[1] <= 0.0:
        raise ValueError("periodic2 requires positive box length on periodic axes.")

    return {
        "axes": tuple(periodic_axes),
        "lengths": tuple(lengths),
        "origins": tuple(box_min[axis] for axis in periodic_axes),
        "image_layers": int(sim.get("field_periodic_image_layers", 1)),
        "far_correction": str(sim.get("field_periodic_far_correction", "none")),
        "ewald_alpha": float(sim.get("field_periodic_ewald_alpha", 0.0)),
        "ewald_layers": int(sim.get("field_periodic_ewald_layers", 4)),
    }


def _coerce_vec3(value: object, *, name: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain exactly 3 values.")
    try:
        return float(value[0]), float(value[1]), float(value[2])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain numeric values.") from exc


def _canonical_boundary_mode(value: object) -> str:
    mode = str(value).strip().lower()
    if mode in {"open", "outflow", "escape"}:
        return "open"
    return mode


def _self_potential_coefficients(
    triangles: np.ndarray, *, self_term: str, softening: float
) -> np.ndarray:
    if self_term == "exclude":
        return np.zeros(triangles.shape[0], dtype=float)

    if self_term == "softened_point":
        if softening <= 0.0:
            raise ValueError("softening must be > 0 when self_term='softened_point'.")
        return np.full(triangles.shape[0], 1.0 / softening, dtype=float)

    if self_term == "area_equivalent":
        areas = _triangle_areas(triangles)
        min_area = np.finfo(float).tiny
        safe_areas = np.maximum(areas, min_area)
        return 2.0 * np.sqrt(np.pi) / np.sqrt(safe_areas)

    raise ValueError(
        "self_term must be one of {'area_equivalent', 'exclude', 'softened_point'}."
    )
