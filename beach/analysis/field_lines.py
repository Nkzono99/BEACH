"""Triangle-panel electric-field evaluation and field-line tracing."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from beach.io._numeric import _points
from beach.io.types import FortranRunResult
from beach.native.kernel import FieldKernel

from .context import RunContext
from .field_kernel import (
    _require_total_field_config,
    _require_total_field_reconstruction,
    field_kernel_from_result,
    field_kernel_options_from_result,
)


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
    context.require_triangle_source_model()
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
    with field_kernel_from_result(
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
    context.require_triangle_source_model()
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

    triangles = context.require_triangles()
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
    with field_kernel_from_result(
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
