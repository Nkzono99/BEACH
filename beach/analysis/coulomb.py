"""Triangle-panel force and torque computation utilities."""

from __future__ import annotations

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Literal, Mapping, Sequence

import numpy as np

from beach.io.types import CoulombInteraction, FortranRunResult, MeshSelection
from beach.native.kernel import FieldKernel

from ._selection import _coerce_group_selection
from .context import RunContext
from .field_kernel import (
    _options_from_result,
    _require_total_field_config,
    _require_total_field_reconstruction,
)
from .mesh import _triangle_centers
from .panel_quadrature import _quadrature_order, panel_target_quadrature
from .periodic import auto_periodic2_from_result


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
    context.require_triangle_source_model()
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
        order=None,
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


_COULOMB_MATRIX_WORKER_CONTEXT: dict[str, object] = {}


def _coulomb_force_matrix(
    context: RunContext,
    *,
    source_mesh_ids: Sequence[int],
    target_mesh_ids: Sequence[int],
    component_idx: int,
    step: int | None,
    config_path: str | Path | None,
    library_path: str | Path | None,
    workers: int,
) -> np.ndarray:
    """Return one force component per source row and target column.

    The extra last row holds the net force from every other source object.
    """

    resolved = context.result
    periodic2 = auto_periodic2_from_result(resolved, context=context)
    matrix = np.zeros((len(source_mesh_ids) + 1, len(target_mesh_ids)), dtype=float)

    tasks: list[tuple[int, int, int, list[int]]] = []
    for col_idx, target_mesh_id in enumerate(target_mesh_ids):
        for row_idx, source_mesh_id in enumerate(source_mesh_ids):
            if source_mesh_id == target_mesh_id:
                continue
            tasks.append((row_idx, col_idx, target_mesh_id, [source_mesh_id]))

        net_sources = [
            mesh_id for mesh_id in source_mesh_ids if mesh_id != target_mesh_id
        ]
        if net_sources:
            tasks.append((len(source_mesh_ids), col_idx, target_mesh_id, net_sources))

    def _compute_entry_serial(
        task: tuple[int, int, int, list[int]],
    ) -> tuple[int, int, float]:
        row_idx, col_idx, target_mesh_id, task_source_ids = task
        interaction = calc_coulomb(
            resolved,
            target=target_mesh_id,
            source=task_source_ids,
            step=step,
            periodic2=periodic2,
            config_path=config_path,
            library_path=library_path,
        )
        return row_idx, col_idx, float(interaction.force_on_a_N[component_idx])

    if workers == 1 or len(tasks) <= 1:
        entries = map(_compute_entry_serial, tasks)
        for row_idx, col_idx, value in entries:
            matrix[row_idx, col_idx] = value
    else:
        max_workers = min(int(workers), len(tasks))
        try:
            mp_context = mp.get_context("fork")
        except ValueError:
            mp_context = None
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=mp_context,
            initializer=_init_coulomb_matrix_worker,
            initargs=(
                resolved,
                step,
                component_idx,
                periodic2,
                config_path,
                library_path,
            ),
        ) as executor:
            for row_idx, col_idx, value in executor.map(
                _compute_coulomb_matrix_entry_worker,
                tasks,
            ):
                matrix[row_idx, col_idx] = value
    return matrix


def _init_coulomb_matrix_worker(
    resolved: FortranRunResult,
    step: int | None,
    component_idx: int,
    periodic2: Mapping[str, object] | None,
    config_path: str | Path | None,
    library_path: str | Path | None,
) -> None:
    _COULOMB_MATRIX_WORKER_CONTEXT.clear()
    _COULOMB_MATRIX_WORKER_CONTEXT.update(
        {
            "resolved": resolved,
            "step": step,
            "component_idx": component_idx,
            "periodic2": periodic2,
            "config_path": config_path,
            "library_path": library_path,
        }
    )


def _compute_coulomb_matrix_entry_worker(
    task: tuple[int, int, int, list[int]],
) -> tuple[int, int, float]:
    row_idx, col_idx, target_mesh_id, source_mesh_ids = task
    interaction = calc_coulomb(
        _COULOMB_MATRIX_WORKER_CONTEXT["resolved"],
        target=target_mesh_id,
        source=source_mesh_ids,
        step=_COULOMB_MATRIX_WORKER_CONTEXT["step"],
        periodic2=_COULOMB_MATRIX_WORKER_CONTEXT["periodic2"],
        config_path=_COULOMB_MATRIX_WORKER_CONTEXT["config_path"],
        library_path=_COULOMB_MATRIX_WORKER_CONTEXT["library_path"],
    )
    component_idx = int(_COULOMB_MATRIX_WORKER_CONTEXT["component_idx"])
    return row_idx, col_idx, float(interaction.force_on_a_N[component_idx])
