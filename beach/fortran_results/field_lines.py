"""Triangle-panel electric-field evaluation and field-line tracing."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .context import RunContext
from .kernel import (
    FieldKernel,
    field_kernel_options_from_result,
    _require_total_field_config,
    _require_total_field_reconstruction,
)
from .selection import _require_triangle_source_model, _require_triangles
from .types import FortranRunResult


def compute_electric_field_points(
    result: FortranRunResult | object,
    points: np.ndarray,
    *,
    chunk_size: int = 2048,
    periodic2: Mapping[str, object] | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
) -> np.ndarray:
    """Evaluate the native triangle-P0 electric field at arbitrary points."""

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    _require_triangle_source_model(resolved)
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0.")
    sample_points = _points(points)
    if sample_points.shape[0] == 0:
        return np.empty((0, 3), dtype=float)

    _require_total_field_config(
        context,
        operation="electric-field recomputation",
    )
    options = field_kernel_options_from_result(
        context,
        periodic2=periodic2,
        config_path=config_path,
    )
    _require_total_field_reconstruction(
        context,
        options,
        operation="electric-field recomputation",
    )
    field = np.empty_like(sample_points)
    with FieldKernel.from_result(
        context,
        periodic2=periodic2,
        config_path=config_path,
        library_path=library_path,
    ) as kernel:
        for start in range(0, sample_points.shape[0], chunk_size):
            stop = min(start + chunk_size, sample_points.shape[0])
            field[start:stop] = kernel.eval_e(sample_points[start:stop])
    return field


def trace_field_lines(
    result: FortranRunResult | object,
    seed_points: np.ndarray,
    *,
    ds: float | None = None,
    max_steps: int = 500,
    periodic2: Mapping[str, object] | None = None,
    direction: str = "both",
    box_min: Iterable[float] | None = None,
    box_max: Iterable[float] | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
) -> list[np.ndarray]:
    """Trace field lines using the same triangle-P0 kernel as the simulator."""

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    _require_triangle_source_model(resolved)
    seeds = np.asarray(seed_points, dtype=float)
    if seeds.ndim == 1 and seeds.shape == (3,):
        seeds = seeds.reshape(1, 3)
    if seeds.ndim != 2 or seeds.shape[1] != 3:
        raise ValueError("seed_points must have shape (n_seeds, 3).")
    if not np.all(np.isfinite(seeds)):
        raise ValueError("seed_points must contain finite values.")
    if isinstance(max_steps, bool) or int(max_steps) != max_steps or max_steps < 0:
        raise ValueError("max_steps must be a non-negative integer.")
    direction_key = str(direction).strip().lower()
    if direction_key not in {"both", "forward", "backward"}:
        raise ValueError("direction must be one of {'both', 'forward', 'backward'}.")

    triangles = _require_triangles(resolved)
    if ds is None:
        edges = np.concatenate(
            [
                np.linalg.norm(triangles[:, 1] - triangles[:, 0], axis=1),
                np.linalg.norm(triangles[:, 2] - triangles[:, 1], axis=1),
                np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1),
            ]
        )
        ds_value = float(np.mean(edges)) * 0.5
    else:
        ds_value = float(ds)
    if not np.isfinite(ds_value) or ds_value <= 0.0:
        raise ValueError("ds must be finite and positive.")

    bb_min = _optional_vec3(box_min, "box_min")
    bb_max = _optional_vec3(box_max, "box_max")
    if (bb_min is None) != (bb_max is None):
        raise ValueError("box_min and box_max must be specified together.")
    if bb_min is not None and bb_max is not None and np.any(bb_max <= bb_min):
        raise ValueError("box_max must be greater than box_min on every axis.")

    _require_total_field_config(
        context,
        operation="field-line tracing",
    )
    options = field_kernel_options_from_result(
        context,
        periodic2=periodic2,
        config_path=config_path,
    )
    _require_total_field_reconstruction(
        context,
        options,
        operation="field-line tracing",
    )
    lines: list[np.ndarray] = []
    with FieldKernel.from_result(
        context,
        periodic2=periodic2,
        config_path=config_path,
        library_path=library_path,
    ) as kernel:
        for seed in seeds:
            forward: np.ndarray | None = None
            backward: np.ndarray | None = None
            if direction_key in {"forward", "both"}:
                forward = _rk4_trace(
                    kernel,
                    seed,
                    ds=ds_value,
                    max_steps=int(max_steps),
                    sign=1.0,
                    bb_min=bb_min,
                    bb_max=bb_max,
                )
            if direction_key in {"backward", "both"}:
                backward = _rk4_trace(
                    kernel,
                    seed,
                    ds=ds_value,
                    max_steps=int(max_steps),
                    sign=-1.0,
                    bb_min=bb_min,
                    bb_max=bb_max,
                )
            if forward is not None and backward is not None:
                lines.append(np.concatenate((backward[::-1], forward[1:]), axis=0))
            elif forward is not None:
                lines.append(forward)
            elif backward is not None:
                lines.append(backward[::-1])
    return lines


def _rk4_trace(
    kernel: FieldKernel,
    seed: np.ndarray,
    *,
    ds: float,
    max_steps: int,
    sign: float,
    bb_min: np.ndarray | None,
    bb_max: np.ndarray | None,
) -> np.ndarray:
    points = [np.asarray(seed, dtype=float).copy()]
    position = points[0].copy()
    for _ in range(max_steps):
        k1 = _unit_field(kernel, position)
        if k1 is None:
            break
        k2 = _unit_field(kernel, position + 0.5 * ds * sign * k1)
        if k2 is None:
            break
        k3 = _unit_field(kernel, position + 0.5 * ds * sign * k2)
        if k3 is None:
            break
        k4 = _unit_field(kernel, position + ds * sign * k3)
        if k4 is None:
            break
        position = position + (ds * sign / 6.0) * (
            k1 + 2.0 * k2 + 2.0 * k3 + k4
        )
        points.append(position.copy())
        if bb_min is not None and np.any(position < bb_min):
            break
        if bb_max is not None and np.any(position > bb_max):
            break
    return np.asarray(points)


def _unit_field(kernel: FieldKernel, point: np.ndarray) -> np.ndarray | None:
    field = kernel.eval_e(np.asarray(point, dtype=float).reshape(1, 3))[0]
    magnitude = float(np.linalg.norm(field))
    if not np.isfinite(magnitude) or magnitude < 1.0e-30:
        return None
    return field / magnitude


def plot_field_lines_3d(
    result: FortranRunResult | object,
    seed_points: np.ndarray,
    *,
    ds: float | None = None,
    max_steps: int = 500,
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
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
):
    """Plot native triangle-P0 field lines and an optional surface mesh."""

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    from .mesh import _configure_mesh_axes, _surface_charge_density

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    lines = trace_field_lines(
        context,
        seed_points,
        ds=ds,
        max_steps=max_steps,
        periodic2=periodic2,
        direction=direction,
        box_min=box_min,
        box_max=box_max,
        config_path=config_path,
        library_path=library_path,
    )
    triangles = _require_triangles(resolved)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    if show_mesh:
        density = _surface_charge_density(resolved.charges, triangles)
        max_abs = max(float(np.max(np.abs(density))), np.finfo(float).tiny)
        norm = plt.Normalize(vmin=-max_abs, vmax=max_abs)
        mapper = plt.cm.ScalarMappable(norm=norm, cmap=mesh_cmap)
        facecolors = mapper.to_rgba(density)
        facecolors[:, 3] = mesh_alpha
        ax.add_collection3d(
            Poly3DCollection(
                triangles,
                facecolors=facecolors,
                edgecolor=(0.0, 0.0, 0.0, 0.15),
                linewidth=0.2,
            )
        )

    if line_color is None:
        cmap = plt.get_cmap(line_cmap)
        colors = [
            cmap(index / max(len(lines) - 1, 1))
            for index in range(len(lines))
        ]
    else:
        colors = [line_color] * len(lines)
    for line, color in zip(lines, colors):
        if line.shape[0] < 2:
            continue
        ax.plot(
            line[:, 0],
            line[:, 1],
            line[:, 2],
            color=color,
            linewidth=line_width,
        )
        midpoint = line.shape[0] // 2
        if 0 < midpoint < line.shape[0] - 1:
            delta = line[midpoint + 1] - line[midpoint]
            ax.quiver(
                *line[midpoint],
                *delta,
                color=color,
                arrow_length_ratio=0.4,
                linewidth=line_width * 1.5,
            )

    seeds = np.asarray(seed_points, dtype=float)
    if seeds.ndim == 1:
        seeds = seeds.reshape(1, 3)
    ax.scatter(
        seeds[:, 0],
        seeds[:, 1],
        seeds[:, 2],
        c="red",
        s=20,
        zorder=5,
        depthshade=False,
    )
    _configure_mesh_axes(
        ax,
        triangles,
        view_elev=view_elev,
        view_azim=view_azim,
    )
    ax.set_title(title)
    fig.tight_layout()
    return fig, ax


def _points(value: np.ndarray) -> np.ndarray:
    points = np.asarray(value, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n_points, 3).")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must contain finite values.")
    return points


def _optional_vec3(
    value: Iterable[float] | None,
    name: str,
) -> np.ndarray | None:
    if value is None:
        return None
    vector = np.asarray(list(value), dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain exactly three finite values.")
    return vector
