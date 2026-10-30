"""Mesh-group selection for force and mesh-access helpers."""

from __future__ import annotations

from typing import Iterable

import numpy as np

from beach.io.types import RunResult, MeshSelection

from .context import RunContext


def _coerce_group_selection(
    result: RunResult,
    group: int | MeshSelection | Iterable[int | MeshSelection],
    *,
    step: int | None,
) -> MeshSelection:
    mesh_ids: list[int] = []
    mesh_steps: list[int] = []

    def _add_item(item: int | MeshSelection) -> None:
        if isinstance(item, MeshSelection):
            if item.directory != result.directory:
                raise ValueError("mesh selections must come from the same output directory.")
            mesh_ids.extend(item.mesh_ids)
            if item.step is not None:
                mesh_steps.append(item.step)
            return
        mesh_ids.append(int(item))

    if isinstance(group, MeshSelection) or isinstance(group, (int, np.integer)):
        _add_item(group)
    else:
        for item in group:
            _add_item(item)

    if not mesh_ids:
        raise ValueError("mesh group is empty.")

    requested_step = step
    if requested_step is None and mesh_steps:
        unique_steps = sorted(set(mesh_steps))
        if len(unique_steps) > 1:
            raise ValueError("mesh selections in one group must share the same step.")
        requested_step = unique_steps[0]

    return _build_mesh_selection(result, tuple(mesh_ids), step=requested_step)


def _mesh_element_indices(
    context: RunContext, mesh_ids: tuple[int, ...]
) -> tuple[tuple[int, ...], np.ndarray]:
    """Validate mesh ids and return them de-duplicated with their element indices."""

    unique_mesh_ids = tuple(dict.fromkeys(int(mid) for mid in mesh_ids))
    all_mesh_ids = context.mesh_ids_or_default()
    available = set(int(mid) for mid in np.unique(all_mesh_ids))
    missing = [mid for mid in unique_mesh_ids if mid not in available]
    if missing:
        raise ValueError(f"unknown mesh id(s): {missing}. available={sorted(available)}")
    mask = np.isin(all_mesh_ids, np.asarray(unique_mesh_ids, dtype=np.int64))
    return unique_mesh_ids, np.flatnonzero(mask)


def _build_mesh_selection(
    result: RunResult, mesh_ids: tuple[int, ...], *, step: int | None
) -> MeshSelection:
    context = RunContext.from_value(result)
    unique_mesh_ids, elem_indices = _mesh_element_indices(context, mesh_ids)
    triangles = context.require_triangles()
    charges = context.charges_at(step)

    return MeshSelection(
        directory=result.directory,
        mesh_ids=unique_mesh_ids,
        elem_indices=elem_indices,
        triangles=triangles[elem_indices],
        charges=charges[elem_indices],
        step=step,
    )
