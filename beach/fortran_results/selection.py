"""Result resolution and mesh-selection helpers."""

from __future__ import annotations

from typing import Iterable

import numpy as np

from .history import FortranChargeHistory
from .types import FortranRunResult, MeshSelection


def _resolve_result(result: FortranRunResult | object) -> FortranRunResult:
    if isinstance(result, FortranRunResult):
        return result

    candidate = getattr(result, "result", None)
    if isinstance(candidate, FortranRunResult):
        return candidate

    raise TypeError("result must be FortranRunResult or Beach.")


def _coerce_group_selection(
    result: FortranRunResult,
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


def _build_mesh_selection(
    result: FortranRunResult, mesh_ids: tuple[int, ...], *, step: int | None
) -> MeshSelection:
    unique_mesh_ids = tuple(dict.fromkeys(int(mid) for mid in mesh_ids))
    all_mesh_ids = _mesh_ids_or_default(result)
    available = set(int(mid) for mid in np.unique(all_mesh_ids))
    missing = [mid for mid in unique_mesh_ids if mid not in available]
    if missing:
        raise ValueError(f"unknown mesh id(s): {missing}. available={sorted(available)}")

    triangles = _require_triangles(result)
    charges = _charges_for_step(result, step=step)
    mask = np.isin(all_mesh_ids, np.asarray(unique_mesh_ids, dtype=np.int64))
    elem_indices = np.flatnonzero(mask)

    return MeshSelection(
        directory=result.directory,
        mesh_ids=unique_mesh_ids,
        elem_indices=elem_indices,
        triangles=triangles[elem_indices],
        charges=charges[elem_indices],
        step=step,
    )


def _charges_for_step(result: FortranRunResult, *, step: int | None) -> np.ndarray:
    if step is None:
        return result.charges

    history = result.history
    if history is None or not history.has_data:
        if step == -1:
            return result.charges
        raise ValueError(
            "charge_history.csv is required when step is specified and must not be empty."
        )
    if step == -1:
        return history.get_step(-1)
    return history.get_step(step)


def _mesh_ids_or_default(result: FortranRunResult) -> np.ndarray:
    if result.mesh_ids is None or result.mesh_ids.size != result.mesh_nelem:
        return np.ones(result.mesh_nelem, dtype=np.int64)
    return result.mesh_ids.astype(np.int64, copy=False)


def _require_triangles(result: FortranRunResult) -> np.ndarray:
    if result.triangles is None:
        raise ValueError(
            "mesh_triangles.csv is not found. Re-run Fortran with latest output format."
        )
    return result.triangles


def _require_history(result: FortranRunResult) -> FortranChargeHistory:
    history = result.history
    if history is None or not history.has_data:
        raise ValueError(
            "charge_history.csv is not found or empty. Enable history output and rerun."
        )
    return history
