"""Triangle-panel force and torque computation utilities."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Iterable, Literal, Mapping

import numpy as np

from .context import RunContext
from .kernel import (
    FieldKernel,
    _options_from_result,
    _require_total_field_config,
    _require_total_field_reconstruction,
)
from .mesh import _triangle_centers
from .panel_quadrature import _quadrature_order, panel_target_quadrature
from .selection import (
    _coerce_group_selection,
    _require_triangle_source_model,
)
from .types import CoulombInteraction, FortranRunResult, MeshSelection


def calc_coulomb(
    result: FortranRunResult | object,
    target: int | MeshSelection | Iterable[int | MeshSelection],
    source: int | MeshSelection | Iterable[int | MeshSelection],
    *,
    step: int | None = -1,
    torque_origin: Literal[
        "target_center",
        "source_center",
        "origin",
        "group_a_center",
        "group_b_center",
    ] = "target_center",
    periodic2: Mapping[str, object] | None = None,
    quadrature_order: int = 7,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
) -> CoulombInteraction:
    """Compute force/torque on one panel group from another panel group."""

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    _require_triangle_source_model(resolved)
    sel_target = _coerce_group_selection(resolved, target, step=step)
    sel_source = _coerce_group_selection(resolved, source, step=step)
    if sel_target.elem_indices.size == 0:
        raise ValueError("target does not contain any mesh elements.")
    if sel_source.elem_indices.size == 0:
        raise ValueError("source does not contain any mesh elements.")
    if np.intersect1d(sel_target.elem_indices, sel_source.elem_indices).size:
        raise ValueError("target and source mesh selections must be disjoint.")

    origin_key = str(torque_origin).strip().lower()
    if origin_key == "group_a_center":
        origin_key = "target_center"
    elif origin_key == "group_b_center":
        origin_key = "source_center"
    if origin_key == "target_center":
        origin = _triangle_centers(sel_target.triangles).mean(axis=0)
    elif origin_key == "source_center":
        origin = _triangle_centers(sel_source.triangles).mean(axis=0)
    elif origin_key == "origin":
        origin = np.zeros(3, dtype=float)
    else:
        raise ValueError(
            "torque_origin must be one of "
            "{'target_center', 'source_center', 'origin'}."
        )

    order = _quadrature_order(quadrature_order)
    target_points, target_weights, _ = panel_target_quadrature(
        sel_target.triangles,
        sel_target.charges,
        order,
    )
    _require_total_field_config(
        context,
        operation="calc_coulomb",
    )
    options = _options_from_result(
        resolved,
        periodic2=periodic2,
        theta=None,
        leaf_max=None,
        order=4,
        config_path=config_path,
        context=context,
    )
    _require_total_field_reconstruction(
        context,
        options,
        operation="calc_coulomb",
    )
    options = replace(options, external_e0=(0.0, 0.0, 0.0))
    with FieldKernel(
        sel_source.triangles,
        sel_source.charges,
        options=options,
        library_path=library_path,
    ) as kernel:
        force_target, torque_target = kernel.force_on_charges(
            target_points,
            target_weights,
            origin=origin,
        )

    force_source = -force_target
    torque_source = -torque_target
    target_count = float(sel_target.elem_indices.size)
    return CoulombInteraction(
        group_a_mesh_ids=sel_target.mesh_ids,
        group_b_mesh_ids=sel_source.mesh_ids,
        step=sel_target.step,
        torque_origin_m=origin,
        force_on_a_N=force_target,
        force_on_b_N=force_source,
        torque_on_a_Nm=torque_target,
        torque_on_b_Nm=torque_source,
        mean_force_on_a_per_element_N=force_target / target_count,
        mean_torque_on_a_per_element_Nm=torque_target / target_count,
    )
