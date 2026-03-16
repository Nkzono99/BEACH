"""Object metadata helpers for mesh-group analyses."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .potential import _find_config_path_near_output, _load_toml
from .selection import _mesh_ids_or_default
from .types import FortranRunResult, MeshSource


@dataclass(frozen=True)
class ResolvedObjectSpec:
    """Resolved metadata for one object-level mesh group."""

    mesh_id: int
    kind: str
    label: str
    template: Mapping[str, object] | None = None


def normalize_kind_filter(kinds: Iterable[str] | None) -> set[str] | None:
    """Normalize an optional iterable of kind names."""

    if kinds is None:
        return None

    normalized = {
        str(kind).strip().lower()
        for kind in kinds
        if str(kind).strip()
    }
    if len(normalized) == 0 or "all" in normalized:
        return None
    return normalized


def resolve_object_specs(
    result: FortranRunResult,
    *,
    config_path: str | Path | None,
) -> list[ResolvedObjectSpec]:
    """Resolve object labels, kinds, and optional template definitions."""

    mesh_ids = tuple(int(v) for v in np.unique(_mesh_ids_or_default(result)))
    config_specs = _object_specs_from_config(
        result,
        mesh_ids=mesh_ids,
        config_path=config_path,
    )
    if config_specs is not None:
        return config_specs

    if result.mesh_sources is not None:
        kinds = [
            _mesh_kind_from_source(result.mesh_sources.get(mesh_id))
            for mesh_id in mesh_ids
        ]
        labels = _labels_from_kinds(kinds)
        return [
            ResolvedObjectSpec(mesh_id=mesh_id, kind=kind, label=label)
            for mesh_id, kind, label in zip(mesh_ids, kinds, labels)
        ]

    fallback_kinds = [f"mesh{mesh_id}" for mesh_id in mesh_ids]
    return [
        ResolvedObjectSpec(mesh_id=mesh_id, kind=kind, label=kind)
        for mesh_id, kind in zip(mesh_ids, fallback_kinds)
    ]


def _object_specs_from_config(
    result: FortranRunResult,
    *,
    mesh_ids: tuple[int, ...],
    config_path: str | Path | None,
) -> list[ResolvedObjectSpec] | None:
    path: Path | None
    if config_path is None:
        path = _find_config_path_near_output(result.directory)
    else:
        path = Path(config_path)
        if not path.exists():
            raise ValueError(f'config file is not found: "{path}".')

    if path is None:
        return None

    config = _load_toml(path)
    mesh_cfg = config.get("mesh")
    if not isinstance(mesh_cfg, Mapping):
        return None
    templates = mesh_cfg.get("templates")
    if not isinstance(templates, list):
        return None

    enabled_templates: list[Mapping[str, object]] = []
    enabled_kinds: list[str] = []
    for template in templates:
        if not isinstance(template, Mapping):
            continue
        if not bool(template.get("enabled", True)):
            continue
        enabled_templates.append(template)
        kind = str(template.get("kind", "mesh")).strip().lower() or "mesh"
        enabled_kinds.append(kind)

    if len(enabled_kinds) != len(mesh_ids):
        return None

    labels = _labels_from_kinds(enabled_kinds)
    return [
        ResolvedObjectSpec(
            mesh_id=mesh_id,
            kind=kind,
            label=label,
            template=template,
        )
        for mesh_id, kind, label, template in zip(
            mesh_ids,
            enabled_kinds,
            labels,
            enabled_templates,
        )
    ]


def _labels_from_kinds(kinds: Iterable[str]) -> list[str]:
    normalized = [str(kind).strip().lower() or "mesh" for kind in kinds]
    totals = Counter(normalized)
    seen: Counter[str] = Counter()
    labels: list[str] = []
    for kind in normalized:
        seen[kind] += 1
        if totals[kind] == 1:
            labels.append(kind)
        else:
            labels.append(f"{kind}{seen[kind]}")
    return labels


def _mesh_kind_from_source(source: MeshSource | None) -> str:
    if source is None:
        return "mesh"
    template_kind = str(source.template_kind).strip().lower()
    if template_kind:
        return template_kind
    source_kind = str(source.source_kind).strip().lower()
    if source_kind:
        return source_kind
    return "mesh"
