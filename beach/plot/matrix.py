"""Object-wise Coulomb-force matrix figures."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from beach.analysis.context import RunContext
from beach.analysis.coulomb import _coulomb_force_matrix
from beach.analysis.objects import normalize_kind_filter, resolve_object_specs
from beach.io.types import FortranRunResult

from ._defaults import DEFAULT_FORCE_CMAP


def plot_coulomb_force_matrix(
    result: FortranRunResult | object,
    *,
    step: int | None = -1,
    component: str = "z",
    target_kinds: Iterable[str] | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    cmap: str = DEFAULT_FORCE_CMAP,
    annotate: bool = True,
    workers: int = 1,
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
    target_kinds : iterable of str or None, default None
        Template kinds used as target objects. ``None`` means all objects.
    config_path : str, pathlib.Path, or None, default None
        Optional ``beach.toml`` path used for object labels/order.
        ``None`` auto-detects near ``result.directory`` when possible.
    cmap : str, default "coolwarm"
        Matplotlib colormap name.
    annotate : bool, default True
        Whether to draw scientific-notation values inside each cell.
    workers : int, default 1
        Number of worker processes used for source/target force entries.

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

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    if workers < 1:
        raise ValueError("workers must be >= 1.")
    component_idx, component_label = _force_component_info(component)
    object_specs = resolve_object_specs(
        resolved,
        config_path=config_path,
        context=context,
    )
    resolved_target_kinds = normalize_kind_filter(target_kinds)
    target_specs = [
        spec
        for spec in object_specs
        if resolved_target_kinds is None or spec.kind in resolved_target_kinds
    ]
    if len(target_specs) == 0:
        available_kinds = sorted({spec.kind for spec in object_specs})
        raise ValueError(
            "target_kinds did not match any objects. "
            f"available={available_kinds}"
        )

    source_labels = [spec.label for spec in object_specs] + ["net"]
    target_labels = [spec.label for spec in target_specs]
    matrix = _coulomb_force_matrix(
        context,
        source_mesh_ids=[int(spec.mesh_id) for spec in object_specs],
        target_mesh_ids=[int(spec.mesh_id) for spec in target_specs],
        component_idx=component_idx,
        step=step,
        library_path=library_path,
        workers=workers,
    )

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
        "target_mesh_ids": tuple(int(spec.mesh_id) for spec in target_specs),
        "source_mesh_ids": tuple(int(spec.mesh_id) for spec in object_specs),
        "step": step,
        "workers": int(workers),
    }
    fig.tight_layout()
    return fig, ax


def _force_component_info(component: str) -> tuple[int, str]:
    normalized = str(component).strip().lower()
    if normalized == "x":
        return 0, "Fx"
    if normalized == "y":
        return 1, "Fy"
    if normalized == "z":
        return 2, "Fz"
    raise ValueError("component must be one of {'x', 'y', 'z'}.")
