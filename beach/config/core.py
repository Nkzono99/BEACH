"""Preset-based BEACH case rendering and validation."""

from __future__ import annotations

import copy
import os
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

from ._toml import load_toml_file, render_toml_document

CASE_FILENAME = "case.toml"
RENDERED_FILENAME = "beach.toml"
CASE_SCHEMA_VERSION = 1
SCHEMA_BASE_URL = "https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas"
BEACH_SCHEMA_URL = f"{SCHEMA_BASE_URL}/beach.schema.json"
CASE_SCHEMA_URL = f"{SCHEMA_BASE_URL}/beach.case.schema.json"
PRESET_SCHEMA_URL = f"{SCHEMA_BASE_URL}/beach.preset.schema.json"
DEFAULT_PRESET_NAMES = (
    "sim/periodic2_fmm",
    "species/solarwind_electron",
    "species/solarwind_ion",
    "mesh/plane_basic",
    "output/standard",
)
TOP_LEVEL_RENDER_ORDER = ("sim", "particles", "mesh", "output")
TOP_LEVEL_CASE_ORDER = ("schema_version", "title", "use_presets", "override")
SCHEMA_DIRECTIVE = (
    "#:schema "
    f"{BEACH_SCHEMA_URL}"
)
CASE_SCHEMA_DIRECTIVE = "#:schema " f"{CASE_SCHEMA_URL}"
PRESET_SCHEMA_DIRECTIVE = "#:schema " f"{PRESET_SCHEMA_URL}"
_FRAGMENT_TOP_LEVEL_KEYS = frozenset({"sim", "particles", "mesh", "output"})
_PRESET_FORBIDDEN_KEYS = frozenset(
    {"schema_version", "title", "use_presets", "override", "base_case"}
)
_FACE_SOURCE_MODES = frozenset({"reservoir_face", "photo_raycast"})


class ConfigError(ValueError):
    """Base error for preset-based BEACH config handling."""


class CaseSpecError(ConfigError):
    """Raised when ``case.toml`` itself violates the case-layer contract."""


class PresetResolutionError(ConfigError):
    """Raised when a referenced preset cannot be resolved or parsed."""


class MergeError(ConfigError):
    """Raised when preset fragments cannot be combined safely."""


class RenderValidationError(ConfigError):
    """Raised when the rendered ``beach.toml`` violates known constraints."""


@dataclass(frozen=True, slots=True)
class CaseDocument:
    """Parsed and validated ``case.toml`` payload."""

    schema_version: int
    use_presets: tuple[str, ...]
    override: dict[str, Any]
    title: str | None = None
    source_path: Path | None = None

    def to_dict(self) -> dict[str, Any]:
        """Materialize a mutable TOML document."""

        document: dict[str, Any] = {
            "schema_version": self.schema_version,
            "use_presets": list(self.use_presets),
        }
        if self.title is not None:
            document["title"] = self.title
        if self.override:
            document["override"] = copy.deepcopy(self.override)
        return document


@dataclass(frozen=True, slots=True)
class ResolvedPreset:
    """One resolved preset fragment and its origin."""

    name: str
    scope: str
    source_label: str
    data: dict[str, Any]
    shadowed_sources: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RenderResult:
    """Fully rendered config plus metadata about how it was built."""

    case: CaseDocument
    config: dict[str, Any]
    presets: tuple[ResolvedPreset, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ListedPreset:
    """One visible preset after local/user/built-in precedence resolution."""

    name: str
    scope: str
    source_path: Path
    shadowed_sources: tuple[str, ...] = ()


def beachx_home() -> Path:
    """Return the BEACHX state directory."""

    override = os.environ.get("BEACHX_HOME")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".config" / "beachx"


def saved_cases_dir() -> Path:
    """Return the directory for saved case templates."""

    return beachx_home() / "cases"


def user_presets_dir() -> Path:
    """Return the user-level preset directory."""

    return beachx_home() / "presets"


def builtin_presets_dir() -> Path:
    """Return the bundled preset directory shipped with the Python package."""

    return Path(__file__).resolve().parent / "presets"


def project_presets_dir(base_dir: str | Path) -> Path:
    """Return the project-local preset directory for one case directory."""

    return Path(base_dir) / ".beachx" / "presets"


def default_case_document(
    *,
    use_presets: Sequence[str] | None = None,
    title: str | None = None,
) -> CaseDocument:
    """Build the default editable ``case.toml`` payload."""

    raw_presets = tuple(use_presets) if use_presets is not None else DEFAULT_PRESET_NAMES
    if len(raw_presets) == 0:
        raise CaseSpecError("case spec error: use_presets must not be empty.")
    selected = tuple(normalize_preset_name(name) for name in raw_presets)
    return CaseDocument(
        schema_version=CASE_SCHEMA_VERSION,
        title=title,
        use_presets=selected,
        override={},
    )


def load_case_document(path: str | Path) -> CaseDocument:
    """Load and validate one ``case.toml`` document."""

    source_path = Path(path)
    raw = load_toml_file(source_path)
    allowed_top_level = {"schema_version", "title", "use_presets", "override"}
    unknown_keys = [key for key in raw if key not in allowed_top_level]
    if unknown_keys:
        raise CaseSpecError(
            "case spec error: unknown top-level key(s): "
            + ", ".join(sorted(unknown_keys))
        )

    schema_version = raw.get("schema_version", CASE_SCHEMA_VERSION)
    if not isinstance(schema_version, int):
        raise CaseSpecError("case spec error: schema_version must be an integer.")
    if schema_version != CASE_SCHEMA_VERSION:
        raise CaseSpecError(
            f"case spec error: unsupported schema_version={schema_version}. "
            f"Supported value is {CASE_SCHEMA_VERSION}."
        )

    title = raw.get("title")
    if title is not None and not isinstance(title, str):
        raise CaseSpecError("case spec error: title must be a string.")

    use_presets_raw = raw.get("use_presets")
    if not isinstance(use_presets_raw, list) or len(use_presets_raw) == 0:
        raise CaseSpecError(
            "case spec error: use_presets must be a non-empty array of strings."
        )
    use_presets: list[str] = []
    for value in use_presets_raw:
        if not isinstance(value, str):
            raise CaseSpecError(
                "case spec error: use_presets must be a non-empty array of strings."
            )
        use_presets.append(normalize_preset_name(value))

    override = raw.get("override", {})
    if override is None:
        override = {}
    if not isinstance(override, Mapping):
        raise CaseSpecError("case spec error: override must be a table.")
    override_document = copy.deepcopy(dict(override))
    _validate_fragment_structure(
        override_document,
        context="case override",
        allow_meta_keys=False,
    )
    _validate_high_level_fragment(override_document, context="case override")

    return CaseDocument(
        schema_version=schema_version,
        title=title,
        use_presets=tuple(use_presets),
        override=override_document,
        source_path=source_path,
    )


def render_case_file(path: str | Path) -> RenderResult:
    """Render one ``case.toml`` path into the final ``beach.toml`` document."""

    return render_case_document(load_case_document(path))


def render_case_document(case_document: CaseDocument) -> RenderResult:
    """Resolve presets, merge them, and validate the final config."""

    case_dir = (
        case_document.source_path.parent
        if case_document.source_path is not None
        else Path.cwd()
    )
    config: dict[str, Any] = {}
    presets: list[ResolvedPreset] = []
    warnings: list[str] = []
    preset_group_names: set[str] = set()

    for preset_name in case_document.use_presets:
        preset = resolve_preset(preset_name, case_dir=case_dir)
        duplicate_groups = preset_group_names.intersection(_fragment_group_names(preset.data))
        if duplicate_groups:
            names = ", ".join(sorted(duplicate_groups))
            raise MergeError(
                "merge error: mesh.groups defined by multiple presets: "
                f"{names}. Use [override.mesh.groups.<name>] for overrides."
            )
        preset_group_names.update(_fragment_group_names(preset.data))
        presets.append(preset)
        warnings.extend(_shadow_warnings(preset))
        config = merge_fragments(config, preset.data)

    config = merge_fragments(config, case_document.override)
    config = resolve_high_level_config(config)
    validate_rendered_config(config)

    return RenderResult(
        case=case_document,
        config=config,
        presets=tuple(presets),
        warnings=tuple(warnings),
    )


def resolve_preset(name: str, *, case_dir: str | Path) -> ResolvedPreset:
    """Resolve one preset name from project-local, user, or built-in storage."""

    normalized_name = normalize_preset_name(name)
    relative_path = preset_name_to_path(normalized_name)

    project_candidate = project_presets_dir(case_dir) / relative_path
    user_candidate = user_presets_dir() / relative_path
    builtin_candidate = builtin_presets_dir() / relative_path

    candidates = (
        ("project-local", project_candidate),
        ("user", user_candidate),
        ("built-in", builtin_candidate),
    )
    existing = [(scope, path) for scope, path in candidates if path.exists()]
    if not existing:
        raise PresetResolutionError(
            f'preset resolution error: preset "{normalized_name}" was not found.'
        )

    scope, resolved_path = existing[0]
    try:
        data = load_toml_file(resolved_path)
    except ValueError as exc:
        raise PresetResolutionError(
            f'preset resolution error: failed to parse "{normalized_name}" '
            f"from {resolved_path}: {exc}"
        ) from exc

    _validate_fragment_structure(
        data,
        context=f'preset "{normalized_name}"',
        allow_meta_keys=False,
    )
    _validate_high_level_fragment(data, context=f'preset "{normalized_name}"')

    shadowed_sources = tuple(str(path) for _, path in existing[1:])
    return ResolvedPreset(
        name=normalized_name,
        scope=scope,
        source_label=str(resolved_path),
        data=data,
        shadowed_sources=shadowed_sources,
    )


def merge_fragments(
    base: Mapping[str, Any] | None,
    overlay: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Merge two config fragments using BEACH config semantics."""

    merged = copy.deepcopy(dict(base) if base is not None else {})
    if overlay is None:
        return merged
    _merge_into(merged, dict(overlay), path=())
    return merged


def resolve_high_level_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve high-level spatial notation into plain beach.toml values."""

    resolved = copy.deepcopy(dict(config))
    sim = resolved.get("sim")
    if isinstance(sim, Mapping):
        resolved["sim"] = _resolve_sim_high_level(dict(sim))
    box_min, box_max = _resolved_box_bounds(resolved)

    particles = resolved.get("particles")
    if isinstance(particles, Mapping):
        species = particles.get("species")
        if isinstance(species, list):
            resolved_species = [
                _resolve_species_high_level(dict(item), box_min=box_min, box_max=box_max)
                for item in species
            ]
            particles_dict = dict(particles)
            particles_dict["species"] = resolved_species
            resolved["particles"] = particles_dict

    mesh = resolved.get("mesh")
    if isinstance(mesh, Mapping):
        resolved["mesh"] = _resolve_mesh_high_level(dict(mesh), box_min=box_min, box_max=box_max)

    return resolved


def _resolve_sim_high_level(sim: dict[str, Any]) -> dict[str, Any]:
    box_origin = sim.pop("box_origin", None)
    box_size = sim.pop("box_size", None)
    if (box_origin is None) != (box_size is None):
        raise ConfigError(
            "high-level config error: sim.box_origin and sim.box_size must be specified together."
        )
    if box_origin is not None and box_size is not None:
        origin = _coerce_numeric_sequence(box_origin, length=3, name="sim.box_origin")
        size = _coerce_numeric_sequence(box_size, length=3, name="sim.box_size")
        if any(component <= 0.0 for component in size):
            raise ConfigError("high-level config error: sim.box_size components must be > 0.")
        sim["box_min"] = origin
        sim["box_max"] = [origin[i] + size[i] for i in range(3)]
    return sim


def _resolve_species_high_level(
    species: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    mode = species.pop("inject_region_mode", None)
    uv_low = species.pop("uv_low", None)
    uv_high = species.pop("uv_high", None)
    source_mode = species.get("source_mode", "volume_seed")
    if mode is None:
        if uv_low is not None or uv_high is not None:
            raise ConfigError(
                "high-level config error: uv_low/uv_high require inject_region_mode=\"face_fraction\"."
            )
        return species
    if not isinstance(source_mode, str) or source_mode not in _FACE_SOURCE_MODES:
        raise ConfigError(
            "high-level config error: inject_region_mode is only supported for "
            'source_mode="reservoir_face" or "photo_raycast".'
        )
    if not isinstance(mode, str):
        raise ConfigError("high-level config error: inject_region_mode must be a string.")
    if mode == "absolute":
        if uv_low is not None or uv_high is not None:
            raise ConfigError(
                "high-level config error: inject_region_mode=\"absolute\" cannot use uv_low/uv_high."
            )
        return species
    if mode != "face_fraction":
        raise ConfigError(
            f"high-level config error: unsupported inject_region_mode={mode!r}."
        )
    if "pos_low" in species or "pos_high" in species:
        raise ConfigError(
            "high-level config error: face_fraction injection cannot be combined with pos_low/pos_high."
        )
    if box_min is None or box_max is None:
        raise ConfigError(
            "high-level config error: inject_region_mode=\"face_fraction\" requires sim.box_min/box_max."
        )
    if uv_low is None or uv_high is None:
        raise ConfigError(
            "high-level config error: face_fraction injection requires uv_low and uv_high."
        )
    uv_low_vec = _coerce_numeric_sequence(uv_low, length=2, name="uv_low")
    uv_high_vec = _coerce_numeric_sequence(uv_high, length=2, name="uv_high")
    if any(value < 0.0 or value > 1.0 for value in [*uv_low_vec, *uv_high_vec]):
        raise ConfigError("high-level config error: uv_low/uv_high must be inside [0, 1].")
    if any(uv_low_vec[i] > uv_high_vec[i] for i in range(2)):
        raise ConfigError("high-level config error: uv_low must be <= uv_high component-wise.")
    inject_face = species.get("inject_face")
    if not isinstance(inject_face, str) or inject_face == "":
        raise ConfigError(
            "high-level config error: face_fraction injection requires inject_face."
        )
    pos_low, pos_high = _resolve_face_fraction_region(
        inject_face=inject_face,
        uv_low=uv_low_vec,
        uv_high=uv_high_vec,
        box_min=box_min,
        box_max=box_max,
    )
    species["pos_low"] = pos_low
    species["pos_high"] = pos_high
    return species


def _resolve_mesh_high_level(
    mesh: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    groups_raw = mesh.pop("groups", {})
    group_definitions = _normalize_group_definitions(groups_raw)
    templates = mesh.get("templates")
    if not isinstance(templates, list):
        return mesh

    resolved_templates: list[dict[str, Any]] = []
    for template in templates:
        if not isinstance(template, Mapping):
            raise ConfigError("high-level config error: mesh.templates entries must be tables.")
        template_dict = dict(template)
        if "group" in template_dict:
            resolved_templates.append(
                _resolve_grouped_template(
                    template_dict,
                    group_definitions=group_definitions,
                    box_min=box_min,
                    box_max=box_max,
                )
            )
        else:
            resolved_templates.append(
                _resolve_direct_template(
                    template_dict,
                    box_min=box_min,
                    box_max=box_max,
                )
            )
    mesh["templates"] = resolved_templates
    return mesh


def _resolve_direct_template(
    template: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    placement_mode = template.pop("placement_mode", "absolute")
    anchor = template.pop("anchor", None)
    offset = template.pop("offset", None)
    offset_frac = template.pop("offset_frac", None)

    if placement_mode == "absolute":
        if anchor is not None or offset is not None or offset_frac is not None:
            raise ConfigError(
                "high-level config error: placement_mode=\"absolute\" cannot use anchor/offset/offset_frac."
            )
        if "center" not in template:
            raise ConfigError("high-level config error: direct mesh template requires center.")
    elif placement_mode == "box_anchor":
        if "center" in template:
            raise ConfigError(
                "high-level config error: placement_mode=\"box_anchor\" cannot be combined with center."
            )
        template["center"] = _resolve_anchor_position(
            anchor=anchor,
            offset=offset,
            offset_frac=offset_frac,
            box_min=box_min,
            box_max=box_max,
        )
    else:
        raise ConfigError(
            f"high-level config error: unsupported placement_mode={placement_mode!r}."
        )

    _resolve_template_size_high_level(template, box_min=box_min, box_max=box_max)
    return template


def _resolve_grouped_template(
    template: dict[str, Any],
    *,
    group_definitions: Mapping[str, Mapping[str, Any]],
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    group_name = template.pop("group")
    if not isinstance(group_name, str) or group_name == "":
        raise ConfigError("high-level config error: mesh template group must be a non-empty string.")
    if group_name not in group_definitions:
        raise ConfigError(
            f'high-level config error: mesh template references undefined group "{group_name}".'
        )
    for forbidden_key in (
        "center",
        "placement_mode",
        "anchor",
        "offset",
        "offset_frac",
        "size_mode",
        "size_frac",
    ):
        if forbidden_key in template:
            raise ConfigError(
                f'high-level config error: grouped mesh template "{group_name}" cannot define {forbidden_key}.'
            )
    center_local = template.pop("center_local", None)
    if center_local is None:
        raise ConfigError(
            f'high-level config error: grouped mesh template "{group_name}" requires center_local.'
        )
    center_local_vec = _coerce_numeric_sequence(center_local, length=3, name="center_local")
    group_origin = _resolve_group_origin(
        dict(group_definitions[group_name]),
        box_min=box_min,
        box_max=box_max,
    )
    group_scale = _resolve_group_scale(
        dict(group_definitions[group_name]),
        box_min=box_min,
        box_max=box_max,
    )
    template["center"] = [
        group_origin[i] + group_scale * center_local_vec[i] for i in range(3)
    ]
    _scale_template_lengths(template, group_scale)
    return template


def _resolve_template_size_high_level(
    template: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> None:
    size_mode = template.pop("size_mode", "absolute")
    size_frac = template.pop("size_frac", None)
    if size_mode == "absolute":
        if size_frac is not None:
            raise ConfigError(
                "high-level config error: size_frac requires size_mode=\"box_fraction\"."
            )
        return
    if size_mode != "box_fraction":
        raise ConfigError(f"high-level config error: unsupported size_mode={size_mode!r}.")
    if size_frac is None:
        raise ConfigError(
            "high-level config error: size_mode=\"box_fraction\" requires size_frac."
        )
    box_size = _require_box_size(box_min=box_min, box_max=box_max, context="size_frac")
    kind = template.get("kind", "plane")
    if not isinstance(kind, str):
        raise ConfigError("high-level config error: mesh template kind must be a string.")
    if kind in {"plane", "plane_hole", "plate_hole"}:
        frac = _coerce_numeric_sequence(size_frac, length=2, name="size_frac")
        template["size_x"] = frac[0] * box_size[0]
        template["size_y"] = frac[1] * box_size[1]
        return
    if kind == "box":
        frac = _coerce_numeric_sequence(size_frac, length=3, name="size_frac")
        template["size"] = [frac[i] * box_size[i] for i in range(3)]
        return
    if kind == "sphere":
        frac = _coerce_numeric_scalar(size_frac, name="size_frac")
        template["radius"] = frac * min(box_size)
        return
    if kind == "cylinder":
        frac = _coerce_numeric_sequence(size_frac, length=2, name="size_frac")
        template["radius"] = frac[0] * min(box_size[0], box_size[1])
        template["height"] = frac[1] * box_size[2]
        return
    raise ConfigError(
        f"high-level config error: size_mode=\"box_fraction\" is not supported for kind={kind!r}."
    )


def _resolve_face_fraction_region(
    *,
    inject_face: str,
    uv_low: Sequence[float],
    uv_high: Sequence[float],
    box_min: Sequence[float],
    box_max: Sequence[float],
) -> tuple[list[float], list[float]]:
    pos_low = [float(component) for component in box_min]
    pos_high = [float(component) for component in box_max]
    if inject_face == "x_low":
        boundary = float(box_min[0])
        ranges = ((1, box_min[1], box_max[1]), (2, box_min[2], box_max[2]))
        axis = 0
    elif inject_face == "x_high":
        boundary = float(box_max[0])
        ranges = ((1, box_min[1], box_max[1]), (2, box_min[2], box_max[2]))
        axis = 0
    elif inject_face == "y_low":
        boundary = float(box_min[1])
        ranges = ((0, box_min[0], box_max[0]), (2, box_min[2], box_max[2]))
        axis = 1
    elif inject_face == "y_high":
        boundary = float(box_max[1])
        ranges = ((0, box_min[0], box_max[0]), (2, box_min[2], box_max[2]))
        axis = 1
    elif inject_face == "z_low":
        boundary = float(box_min[2])
        ranges = ((0, box_min[0], box_max[0]), (1, box_min[1], box_max[1]))
        axis = 2
    elif inject_face == "z_high":
        boundary = float(box_max[2])
        ranges = ((0, box_min[0], box_max[0]), (1, box_min[1], box_max[1]))
        axis = 2
    else:
        raise ConfigError(f"high-level config error: invalid inject_face={inject_face!r}.")

    pos_low[axis] = boundary
    pos_high[axis] = boundary
    for uv_index, (coordinate_index, low_bound, high_bound) in enumerate(ranges):
        span = float(high_bound) - float(low_bound)
        pos_low[coordinate_index] = float(low_bound) + uv_low[uv_index] * span
        pos_high[coordinate_index] = float(low_bound) + uv_high[uv_index] * span
    return pos_low, pos_high


def _normalize_group_definitions(groups_raw: object) -> dict[str, dict[str, Any]]:
    if groups_raw in ({}, None):
        return {}
    if not isinstance(groups_raw, Mapping):
        raise ConfigError("high-level config error: mesh.groups must be a table of group definitions.")
    normalized: dict[str, dict[str, Any]] = {}
    for name, value in groups_raw.items():
        if not isinstance(value, Mapping):
            raise ConfigError(
                f'high-level config error: mesh.groups.{name} must be a table.'
            )
        normalized[name] = dict(value)
    return normalized


def _resolve_group_origin(
    group: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> list[float]:
    placement_mode = group.pop("placement_mode", "absolute")
    anchor = group.pop("anchor", None)
    offset = group.pop("offset", None)
    offset_frac = group.pop("offset_frac", None)
    group.pop("scale", None)
    group.pop("scale_from", None)
    group.pop("scale_factor", None)
    if placement_mode == "absolute":
        if anchor is not None:
            raise ConfigError(
                "high-level config error: group placement_mode=\"absolute\" cannot use anchor."
            )
        base = [0.0, 0.0, 0.0]
        offset_vec = _resolve_offset_vector(
            offset=offset,
            offset_frac=offset_frac,
            box_min=box_min,
            box_max=box_max,
            context="group absolute placement",
            allow_zero=True,
        )
        return [base[i] + offset_vec[i] for i in range(3)]
    if placement_mode != "box_anchor":
        raise ConfigError(
            f"high-level config error: unsupported group placement_mode={placement_mode!r}."
        )
    return _resolve_anchor_position(
        anchor=anchor,
        offset=offset,
        offset_frac=offset_frac,
        box_min=box_min,
        box_max=box_max,
    )


def _resolve_group_scale(
    group: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> float:
    scale = group.get("scale")
    scale_from = group.get("scale_from")
    scale_factor = group.get("scale_factor")
    if scale is not None and (scale_from is not None or scale_factor is not None):
        raise ConfigError(
            "high-level config error: group scale cannot combine scale with scale_from/scale_factor."
        )
    if scale is not None:
        scale_value = _coerce_numeric_scalar(scale, name="scale")
        if scale_value <= 0.0:
            raise ConfigError("high-level config error: group scale must be > 0.")
        return scale_value
    if scale_from is None and scale_factor is None:
        return 1.0
    if scale_from is None or scale_factor is None:
        raise ConfigError(
            "high-level config error: group scale_from and scale_factor must be specified together."
        )
    if not isinstance(scale_from, str):
        raise ConfigError("high-level config error: scale_from must be a string.")
    factor = _coerce_numeric_scalar(scale_factor, name="scale_factor")
    if factor <= 0.0:
        raise ConfigError("high-level config error: scale_factor must be > 0.")
    reference = _resolve_scale_reference(scale_from, box_min=box_min, box_max=box_max)
    return factor * reference


def _resolve_scale_reference(
    scale_from: str,
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> float:
    box_size = _require_box_size(box_min=box_min, box_max=box_max, context="scale_from")
    refs = {
        "box_x": box_size[0],
        "box_y": box_size[1],
        "box_z": box_size[2],
        "box_min_xy": min(box_size[0], box_size[1]),
        "box_max_xy": max(box_size[0], box_size[1]),
        "box_min_xyz": min(box_size),
        "box_max_xyz": max(box_size),
    }
    if scale_from not in refs:
        raise ConfigError(
            f"high-level config error: unsupported scale_from={scale_from!r}."
        )
    return refs[scale_from]


def _resolve_anchor_position(
    *,
    anchor: object,
    offset: object,
    offset_frac: object,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> list[float]:
    if not isinstance(anchor, str) or anchor == "":
        raise ConfigError("high-level config error: box_anchor placement requires anchor.")
    base = _resolve_anchor(anchor, box_min=box_min, box_max=box_max)
    offset_vec = _resolve_offset_vector(
        offset=offset,
        offset_frac=offset_frac,
        box_min=box_min,
        box_max=box_max,
        context=f"anchor {anchor}",
        allow_zero=True,
    )
    return [base[i] + offset_vec[i] for i in range(3)]


def _resolve_anchor(
    anchor: str,
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> list[float]:
    if box_min is None or box_max is None:
        raise ConfigError(
            "high-level config error: box_anchor placement requires sim.box_min/box_max."
        )
    center = [(box_min[i] + box_max[i]) * 0.5 for i in range(3)]
    anchors = {
        "box_center": center,
        "x_low_face_center": [box_min[0], center[1], center[2]],
        "x_high_face_center": [box_max[0], center[1], center[2]],
        "y_low_face_center": [center[0], box_min[1], center[2]],
        "y_high_face_center": [center[0], box_max[1], center[2]],
        "z_low_face_center": [center[0], center[1], box_min[2]],
        "z_high_face_center": [center[0], center[1], box_max[2]],
    }
    if anchor not in anchors:
        raise ConfigError(f"high-level config error: unsupported anchor={anchor!r}.")
    return list(anchors[anchor])


def _resolve_offset_vector(
    *,
    offset: object,
    offset_frac: object,
    box_min: list[float] | None,
    box_max: list[float] | None,
    context: str,
    allow_zero: bool,
) -> list[float]:
    if offset is not None and offset_frac is not None:
        raise ConfigError(
            f"high-level config error: {context} cannot combine offset and offset_frac."
        )
    if offset is None and offset_frac is None:
        return [0.0, 0.0, 0.0] if allow_zero else []
    if offset is not None:
        return _coerce_numeric_sequence(offset, length=3, name="offset")
    frac = _coerce_numeric_sequence(offset_frac, length=3, name="offset_frac")
    box_size = _require_box_size(box_min=box_min, box_max=box_max, context="offset_frac")
    return [frac[i] * box_size[i] for i in range(3)]


def _scale_template_lengths(template: dict[str, Any], scale: float) -> None:
    for key in ("size_x", "size_y", "radius", "inner_radius", "height"):
        if key in template:
            template[key] = scale * float(template[key])
    if "size" in template:
        template["size"] = [scale * float(value) for value in template["size"]]


def _resolved_box_bounds(config: Mapping[str, Any]) -> tuple[list[float] | None, list[float] | None]:
    sim = config.get("sim")
    if not isinstance(sim, Mapping):
        return None, None
    box_min = sim.get("box_min")
    box_max = sim.get("box_max")
    if box_min is None or box_max is None:
        return None, None
    return (
        _coerce_numeric_sequence(box_min, length=3, name="box_min"),
        _coerce_numeric_sequence(box_max, length=3, name="box_max"),
    )


def _require_box_size(
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
    context: str,
) -> list[float]:
    if box_min is None or box_max is None:
        raise ConfigError(
            f"high-level config error: {context} requires sim.box_min/box_max."
        )
    box_size = [box_max[i] - box_min[i] for i in range(3)]
    if any(component <= 0.0 for component in box_size):
        raise ConfigError(
            f"high-level config error: {context} requires positive box dimensions."
        )
    return box_size


def _coerce_numeric_sequence(value: object, *, length: int, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != length:
        raise ConfigError(f"high-level config error: {name} must be a {length}-element numeric array.")
    return [float(value[index]) for index in range(length)]


def _coerce_numeric_scalar(value: object, *, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigError(f"high-level config error: {name} must be numeric.")
    return float(value)


def validate_rendered_config(config: Mapping[str, Any]) -> None:
    """Validate the merged final config against known BEACH constraints."""

    final_config = copy.deepcopy(dict(config))
    _validate_fragment_structure(
        final_config,
        context="rendered config",
        allow_meta_keys=False,
    )

    for key in TOP_LEVEL_RENDER_ORDER:
        if key not in final_config:
            raise RenderValidationError(
                f"BEACH constraint error: rendered config is missing top-level [{key}] table."
            )

    sim = _require_table(final_config, "sim", context="rendered config")
    particles = _require_table(final_config, "particles", context="rendered config")
    mesh = _require_table(final_config, "mesh", context="rendered config")
    _require_table(final_config, "output", context="rendered config")

    species = particles.get("species")
    if not isinstance(species, list) or len(species) == 0 or not all(
        isinstance(item, Mapping) for item in species
    ):
        raise RenderValidationError(
            "BEACH constraint error: particles.species must be a non-empty array of tables."
        )

    resolved_batch_duration = _resolve_batch_duration(sim)
    use_box = bool(sim.get("use_box", False))

    field_bc_mode = sim.get("field_bc_mode", "free")
    field_solver = sim.get("field_solver", "auto")
    if field_bc_mode == "periodic2":
        if field_solver != "fmm":
            raise RenderValidationError(
                'BEACH constraint error: field_bc_mode="periodic2" requires field_solver="fmm".'
            )
        if not use_box:
            raise RenderValidationError(
                'BEACH constraint error: field_bc_mode="periodic2" requires sim.use_box=true.'
            )
        periodic_axes = _periodic_axis_count(sim)
        if periodic_axes != 2:
            raise RenderValidationError(
                "BEACH constraint error: field_bc_mode=\"periodic2\" requires exactly "
                "two periodic axes."
            )

    box_min = _maybe_vec3(sim.get("box_min"), name="sim.box_min")
    box_max = _maybe_vec3(sim.get("box_max"), name="sim.box_max")

    uses_face_sources = False
    has_volume_seed = False
    total_npcls_per_step = 0
    for index, item in enumerate(species, start=1):
        species_table = dict(item)
        source_mode = species_table.get("source_mode", "volume_seed")
        if not isinstance(source_mode, str):
            raise RenderValidationError(
                f"BEACH constraint error: particles.species[{index}] source_mode must be a string."
            )
        if "temperature_k" in species_table and "temperature_ev" in species_table:
            raise RenderValidationError(
                f"BEACH constraint error: particles.species[{index}] cannot define both "
                "temperature_k and temperature_ev."
            )

        if source_mode == "volume_seed":
            has_volume_seed = True
            npcls_per_step = species_table.get("npcls_per_step", 0)
            if not isinstance(npcls_per_step, int):
                raise RenderValidationError(
                    f"BEACH constraint error: particles.species[{index}].npcls_per_step "
                    "must be an integer."
                )
            total_npcls_per_step += npcls_per_step
            if "target_macro_particles_per_batch" in species_table:
                raise RenderValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="volume_seed" and cannot define '
                    "target_macro_particles_per_batch."
                )
            continue

        if source_mode == "reservoir_face":
            uses_face_sources = True
            _validate_face_source_common(
                species_table,
                index=index,
                source_mode=source_mode,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
                box_min=box_min,
                box_max=box_max,
            )
            if (
                "number_density_cm3" not in species_table
                and "number_density_m3" not in species_table
            ):
                raise RenderValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="reservoir_face" and requires number_density_cm3 '
                    "or number_density_m3."
                )
            if (
                "w_particle" in species_table
                and "target_macro_particles_per_batch" in species_table
            ):
                raise RenderValidationError(
                    f"BEACH constraint error: particles.species[{index}] cannot define both "
                    "w_particle and target_macro_particles_per_batch."
                )
            continue

        if source_mode == "photo_raycast":
            uses_face_sources = True
            _validate_face_source_common(
                species_table,
                index=index,
                source_mode=source_mode,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
                box_min=box_min,
                box_max=box_max,
            )
            current_density = species_table.get("emit_current_density_a_m2", 0.0)
            rays_per_batch = species_table.get("rays_per_batch", 0)
            if not isinstance(current_density, (int, float)) or float(current_density) <= 0.0:
                raise RenderValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="photo_raycast" and requires emit_current_density_a_m2 > 0.'
                )
            if not isinstance(rays_per_batch, int) or rays_per_batch <= 0:
                raise RenderValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="photo_raycast" and requires rays_per_batch > 0.'
                )
            forbidden = (
                "npcls_per_step",
                "number_density_cm3",
                "number_density_m3",
                "w_particle",
                "target_macro_particles_per_batch",
            )
            for key in forbidden:
                if key in species_table:
                    raise RenderValidationError(
                        f"BEACH constraint error: particles.species[{index}] uses "
                        f'source_mode="photo_raycast" and cannot define {key}.'
                    )
            continue

        raise RenderValidationError(
            f"BEACH constraint error: particles.species[{index}] has unsupported "
            f"source_mode={source_mode!r}."
        )

    if has_volume_seed and not uses_face_sources and total_npcls_per_step < 1:
        raise RenderValidationError(
            "BEACH constraint error: volume_seed species require total npcls_per_step >= 1."
        )

    _validate_rendered_mesh(mesh)


def render_beach_toml(
    config: Mapping[str, Any],
    *,
    source_case: str | Path | None = None,
) -> str:
    """Render one validated final config to ``beach.toml`` text."""

    header_comments = [SCHEMA_DIRECTIVE, "# Generated by beachx config render"]
    if source_case is not None:
        header_comments.append(f"# source_case={source_case}")
    return render_toml_document(
        config,
        header_comments=header_comments,
        top_level_order=TOP_LEVEL_RENDER_ORDER,
    )


def render_case_toml(case_document: CaseDocument) -> str:
    """Render an editable ``case.toml`` file."""

    header_comments = [
        CASE_SCHEMA_DIRECTIVE,
        "# Generated by beachx config init",
        "# Edit this file and run `beachx config render` to create beach.toml.",
    ]
    return render_toml_document(
        case_document.to_dict(),
        header_comments=header_comments,
        top_level_order=TOP_LEVEL_CASE_ORDER,
    )


def list_saved_cases() -> list[str]:
    """Return saved case names without the ``.toml`` suffix."""

    root = saved_cases_dir()
    if not root.exists():
        return []
    return sorted(path.stem for path in root.glob("*.toml") if path.is_file())


def load_saved_case_path(name: str) -> Path:
    """Resolve one saved case template path."""

    normalized = normalize_saved_case_name(name)
    path = saved_cases_dir() / f"{normalized}.toml"
    if not path.exists():
        raise FileNotFoundError(f"saved case not found: {normalized}")
    return path


def save_case_file(
    source_path: str | Path,
    name: str,
    *,
    force: bool = False,
) -> Path:
    """Copy one current case file into the saved-case store."""

    source = Path(source_path)
    normalized = normalize_saved_case_name(name)
    destination = saved_cases_dir() / f"{normalized}.toml"
    if destination.exists() and not force:
        raise FileExistsError(f"saved case already exists: {normalized}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return destination


def preset_category(name: str) -> str:
    """Return the category prefix for one preset name."""

    normalized = normalize_preset_name(name)
    return PurePosixPath(normalized).parts[0]


def preset_target_path(
    name: str,
    *,
    scope: str,
    case_dir: str | Path | None = None,
) -> Path:
    """Return the writable location for one user or project-local preset."""

    normalized = normalize_preset_name(name)
    relative_path = preset_name_to_path(normalized)
    if scope == "project-local":
        base_dir = Path.cwd() if case_dir is None else Path(case_dir)
        return project_presets_dir(base_dir) / relative_path
    if scope == "user":
        return user_presets_dir() / relative_path
    if scope == "built-in":
        return builtin_presets_dir() / relative_path
    raise ValueError(f"unsupported preset scope: {scope}")


def default_preset_fragment(name: str) -> dict[str, Any]:
    """Return a starter fragment for one new preset name."""

    normalized = normalize_preset_name(name)
    category = PurePosixPath(normalized).parts[0]
    if category == "sim":
        return {"sim": {}}
    if category == "species":
        return {"particles": {"species": [{}]}}
    if category == "mesh":
        return {"mesh": {}}
    if category == "output":
        return {"output": {}}
    raise CaseSpecError(
        f"case spec error: unsupported preset category {category!r}. "
        "Use sim/, species/, mesh/, or output/."
    )


def render_preset_toml(
    fragment: Mapping[str, Any],
    *,
    preset_name: str | None = None,
    source_label: str | Path | None = None,
) -> str:
    """Render one preset fragment as TOML text."""

    header_comments = [PRESET_SCHEMA_DIRECTIVE, "# BEACH preset fragment"]
    if preset_name is not None:
        header_comments.append(f"# preset_name={preset_name}")
    if source_label is not None:
        header_comments.append(f"# source={source_label}")
    return render_toml_document(
        fragment,
        header_comments=header_comments,
        top_level_order=TOP_LEVEL_RENDER_ORDER,
    )


def save_preset_fragment(
    name: str,
    fragment: Mapping[str, Any],
    *,
    scope: str,
    case_dir: str | Path | None = None,
    force: bool = False,
) -> Path:
    """Save one preset fragment into project-local or user storage."""

    materialized = copy.deepcopy(dict(fragment))
    _validate_fragment_structure(
        materialized,
        context=f'preset "{name}"',
        allow_meta_keys=False,
    )
    _validate_high_level_fragment(materialized, context=f'preset "{name}"')
    destination = preset_target_path(name, scope=scope, case_dir=case_dir)
    if destination.exists() and not force:
        raise FileExistsError(f"preset already exists: {name}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        render_preset_toml(materialized, preset_name=name, source_label=destination),
        encoding="utf-8",
    )
    return destination


def list_presets(*, case_dir: str | Path) -> list[ListedPreset]:
    """List visible presets with local/user/built-in precedence applied."""

    grouped: dict[str, list[tuple[str, Path]]] = {}
    search_roots = (
        ("project-local", project_presets_dir(case_dir)),
        ("user", user_presets_dir()),
        ("built-in", builtin_presets_dir()),
    )
    for scope, root in search_roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*.toml")):
            if not path.is_file():
                continue
            relative = path.relative_to(root)
            name = str(PurePosixPath(relative.as_posix()).with_suffix(""))
            grouped.setdefault(name, []).append((scope, path))

    listed: list[ListedPreset] = []
    for name in sorted(grouped):
        candidates = grouped[name]
        scope, path = candidates[0]
        listed.append(
            ListedPreset(
                name=name,
                scope=scope,
                source_path=path,
                shadowed_sources=tuple(str(candidate[1]) for candidate in candidates[1:]),
            )
        )
    return listed


def extract_preset_fragment(
    config: Mapping[str, Any],
    *,
    section: str,
    index: int | None = None,
) -> dict[str, Any]:
    """Extract one reusable preset fragment from a rendered config."""

    normalized = _normalize_preset_section(section)
    materialized = copy.deepcopy(dict(config))

    if normalized == "sim":
        sim = materialized.get("sim")
        if not isinstance(sim, Mapping):
            raise ConfigError('preset save error: rendered config does not contain [sim].')
        return {"sim": dict(sim)}

    if normalized == "mesh":
        mesh = materialized.get("mesh")
        if not isinstance(mesh, Mapping):
            raise ConfigError('preset save error: rendered config does not contain [mesh].')
        return {"mesh": dict(mesh)}

    if normalized == "output":
        output = materialized.get("output")
        if not isinstance(output, Mapping):
            raise ConfigError('preset save error: rendered config does not contain [output].')
        return {"output": dict(output)}

    if normalized == "species":
        species_list = _extract_array_of_tables(
            materialized,
            table_key="particles",
            array_key="species",
            section_label="particles.species",
        )
        selected = _select_indexed_item(
            species_list,
            index=index,
            section_label="particles.species",
        )
        return {"particles": {"species": [selected]}}

    templates = _extract_array_of_tables(
        materialized,
        table_key="mesh",
        array_key="templates",
        section_label="mesh.templates",
    )
    selected = _select_indexed_item(
        templates,
        index=index,
        section_label="mesh.templates",
    )
    return {"mesh": {"templates": [selected]}}


def semantic_diff(left: Any, right: Any) -> list[str]:
    """Return human-readable semantic differences between two TOML payloads."""

    lines: list[str] = []
    _append_semantic_diff(lines, (), left, right)
    return lines


def normalize_preset_name(name: str) -> str:
    """Validate and normalize one preset reference name."""

    normalized = name.strip().replace("\\", "/")
    if normalized == "":
        raise CaseSpecError("case spec error: preset names cannot be empty.")
    if normalized.startswith("/") or normalized.endswith(".toml"):
        raise CaseSpecError(
            f"case spec error: invalid preset name {name!r}. Use names like sim/periodic2_fmm."
        )
    path = PurePosixPath(normalized)
    if any(part in ("", ".", "..") for part in path.parts):
        raise CaseSpecError(
            f"case spec error: invalid preset name {name!r}. Use names like sim/periodic2_fmm."
        )
    return str(path)


def normalize_preset_section(section: str) -> str:
    """Validate and normalize one preset extraction section selector."""

    return _normalize_preset_section(section)


def preset_name_to_path(name: str) -> Path:
    """Convert one normalized preset name to a relative TOML path."""

    return Path(*PurePosixPath(name).parts).with_suffix(".toml")


def normalize_saved_case_name(name: str) -> str:
    """Validate one saved-case handle."""

    normalized = name.strip()
    if normalized == "":
        raise CaseSpecError("case spec error: saved case name cannot be empty.")
    invalid = {"/", "\\", ".."}
    if any(token in normalized for token in invalid):
        raise CaseSpecError(
            f"case spec error: invalid saved case name {name!r}. "
            "Use a simple name like cavity-base."
        )
    return normalized


def _fragment_group_names(fragment: Mapping[str, Any]) -> set[str]:
    mesh = fragment.get("mesh")
    if not isinstance(mesh, Mapping):
        return set()
    groups = mesh.get("groups")
    if not isinstance(groups, Mapping):
        return set()
    return {str(name) for name in groups}


def _merge_into(base: dict[str, Any], overlay: dict[str, Any], *, path: tuple[str, ...]) -> None:
    for key, overlay_value in overlay.items():
        current_path = (*path, key)
        if key not in base:
            base[key] = copy.deepcopy(overlay_value)
            continue

        base_value = base[key]
        if isinstance(base_value, Mapping):
            if not isinstance(overlay_value, Mapping):
                raise MergeError(
                    f"merge error: type conflict at {_format_path(current_path)}. "
                    "Cannot replace a table with a non-table value."
                )
            _merge_into(base_value, dict(overlay_value), path=current_path)
            continue

        if _is_array_of_tables(base_value) or _is_array_of_tables(overlay_value):
            if not _is_array_of_tables(base_value) or not _is_array_of_tables(overlay_value):
                raise MergeError(
                    f"merge error: type conflict at {_format_path(current_path)}. "
                    "Array-of-tables values must merge with the same kind."
                )
            base_value.extend(copy.deepcopy(list(overlay_value)))
            continue

        if isinstance(overlay_value, Mapping):
            raise MergeError(
                f"merge error: type conflict at {_format_path(current_path)}. "
                "Cannot replace a non-table value with a table."
            )

        base[key] = copy.deepcopy(overlay_value)


def _validate_fragment_structure(
    document: Mapping[str, Any],
    *,
    context: str,
    allow_meta_keys: bool,
) -> None:
    unknown_keys = [key for key in document if key not in _FRAGMENT_TOP_LEVEL_KEYS]
    if not allow_meta_keys:
        forbidden = [key for key in document if key in _PRESET_FORBIDDEN_KEYS]
        if forbidden:
            raise CaseSpecError(
                f"{context} error: reserved top-level key(s) are not allowed: "
                + ", ".join(sorted(forbidden))
            )
    if unknown_keys:
        raise CaseSpecError(
            f"{context} error: unsupported top-level key(s): " + ", ".join(sorted(unknown_keys))
        )

    for key in _FRAGMENT_TOP_LEVEL_KEYS.intersection(document):
        value = document[key]
        if not isinstance(value, Mapping):
            raise CaseSpecError(f"{context} error: top-level key {key!r} must be a table.")

    particles = document.get("particles")
    if isinstance(particles, Mapping) and "species" in particles:
        species = particles["species"]
        if not isinstance(species, list) or not all(isinstance(item, Mapping) for item in species):
            raise CaseSpecError(
                f"{context} error: particles.species must be an array of tables."
            )

    mesh = document.get("mesh")
    if isinstance(mesh, Mapping) and "templates" in mesh:
        templates = mesh["templates"]
        if not isinstance(templates, list) or not all(
            isinstance(item, Mapping) for item in templates
        ):
            raise CaseSpecError(
                f"{context} error: mesh.templates must be an array of tables."
            )


def _validate_high_level_fragment(
    document: Mapping[str, Any],
    *,
    context: str,
) -> None:
    sim = document.get("sim")
    if isinstance(sim, Mapping):
        if "box_origin" in sim and "box_min" in sim:
            raise CaseSpecError(
                f"{context} error: sim.box_origin and sim.box_min cannot be specified "
                "in the same fragment."
            )
        if "box_size" in sim and "box_max" in sim:
            raise CaseSpecError(
                f"{context} error: sim.box_size and sim.box_max cannot be specified "
                "in the same fragment."
            )

    particles = document.get("particles")
    species = particles.get("species") if isinstance(particles, Mapping) else None
    if isinstance(species, list):
        for index, item in enumerate(species, start=1):
            if not isinstance(item, Mapping):
                continue
            if not any(key in item for key in ("inject_region_mode", "uv_low", "uv_high")):
                continue
            source_mode = item.get("source_mode", "volume_seed")
            if not isinstance(source_mode, str) or source_mode not in _FACE_SOURCE_MODES:
                raise CaseSpecError(
                    f"{context} error: particles.species[{index}] uses inject_region_mode/uv_* "
                    'but source_mode must be "reservoir_face" or "photo_raycast".'
                )


def _resolve_batch_duration(sim: Mapping[str, Any]) -> float:
    dt = float(sim.get("dt", 0.0))
    has_batch_duration = "batch_duration" in sim
    has_batch_duration_step = "batch_duration_step" in sim
    if has_batch_duration and has_batch_duration_step:
        raise RenderValidationError(
            "BEACH constraint error: sim.batch_duration and sim.batch_duration_step "
            "cannot be specified together."
        )
    if has_batch_duration_step:
        return dt * float(sim["batch_duration_step"])
    return float(sim.get("batch_duration", 0.0))


def _periodic_axis_count(sim: Mapping[str, Any]) -> int:
    count = 0
    for axis in ("x", "y", "z"):
        low = sim.get(f"bc_{axis}_low", "open")
        high = sim.get(f"bc_{axis}_high", "open")
        if low == "periodic" and high == "periodic":
            count += 1
        elif low == "periodic" or high == "periodic":
            raise RenderValidationError(
                f"BEACH constraint error: bc_{axis}_low/high must both be periodic or both non-periodic."
            )
    return count


def _validate_face_source_common(
    species_table: Mapping[str, Any],
    *,
    index: int,
    source_mode: str,
    use_box: bool,
    batch_duration: float,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> None:
    if not use_box:
        raise RenderValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires sim.use_box=true.'
        )
    if batch_duration <= 0.0:
        raise RenderValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires batch_duration > 0.'
        )
    inject_face = species_table.get("inject_face")
    if not isinstance(inject_face, str) or inject_face == "":
        raise RenderValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires inject_face.'
        )
    if box_min is not None and box_max is not None:
        _validate_face_bounds(
            species_table,
            index=index,
            inject_face=inject_face,
            box_min=box_min,
            box_max=box_max,
        )


def _validate_face_bounds(
    species_table: Mapping[str, Any],
    *,
    index: int,
    inject_face: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
) -> None:
    pos_low = _maybe_vec3(
        species_table.get("pos_low"),
        name=f"particles.species[{index}].pos_low",
    )
    pos_high = _maybe_vec3(
        species_table.get("pos_high"),
        name=f"particles.species[{index}].pos_high",
    )
    if pos_low is None or pos_high is None:
        raise RenderValidationError(
            f"BEACH constraint error: particles.species[{index}] must define pos_low "
            "and pos_high on the inject_face."
        )

    axis_by_face = {
        "x_low": (0, float(box_min[0])),
        "x_high": (0, float(box_max[0])),
        "y_low": (1, float(box_min[1])),
        "y_high": (1, float(box_max[1])),
        "z_low": (2, float(box_min[2])),
        "z_high": (2, float(box_max[2])),
    }
    if inject_face not in axis_by_face:
        raise RenderValidationError(
            f"BEACH constraint error: particles.species[{index}] has invalid inject_face={inject_face!r}."
        )
    axis, boundary = axis_by_face[inject_face]
    if pos_low[axis] != boundary or pos_high[axis] != boundary:
        raise RenderValidationError(
            f"BEACH constraint error: particles.species[{index}] pos_low/pos_high must "
            f"lie on inject_face={inject_face!r}."
        )
    for other_axis in range(3):
        if other_axis == axis:
            continue
        if pos_low[other_axis] > pos_high[other_axis]:
            raise RenderValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_low must be <= "
                "pos_high along the inject-face coordinates."
            )
        low_bound = float(box_min[other_axis])
        high_bound = float(box_max[other_axis])
        if not (low_bound <= pos_low[other_axis] <= high_bound):
            raise RenderValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_low is outside the box."
            )
        if not (low_bound <= pos_high[other_axis] <= high_bound):
            raise RenderValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_high is outside the box."
            )


def _require_table(
    document: Mapping[str, Any],
    key: str,
    *,
    context: str,
) -> dict[str, Any]:
    value = document.get(key)
    if not isinstance(value, Mapping):
        raise RenderValidationError(
            f"BEACH constraint error: {context} requires [{key}] to be a table."
        )
    return dict(value)


def _validate_rendered_mesh(mesh: Mapping[str, Any]) -> None:
    templates = mesh.get("templates")
    if templates is None:
        return
    if not isinstance(templates, list) or not all(isinstance(item, Mapping) for item in templates):
        raise RenderValidationError(
            "BEACH constraint error: mesh.templates must be an array of tables."
        )
    for index, item in enumerate(templates, start=1):
        _validate_rendered_template(dict(item), index=index)


def _validate_rendered_template(template: Mapping[str, Any], *, index: int) -> None:
    kind_value = template.get("kind", "plane")
    if not isinstance(kind_value, str):
        raise RenderValidationError(
            f"BEACH constraint error: mesh.templates[{index}].kind must be a string."
        )
    kind = kind_value.strip().lower() or "plane"

    if "center" in template:
        _maybe_vec3(template.get("center"), name=f"mesh.templates[{index}].center")

    if kind == "plane":
        _positive_template_scalar(template, index=index, key="size_x", default=1.0)
        _positive_template_scalar(template, index=index, key="size_y", default=1.0)
        return

    if kind in {"plate_hole", "plane_hole"}:
        size_x = _positive_template_scalar(template, index=index, key="size_x", default=1.0)
        size_y = _positive_template_scalar(template, index=index, key="size_y", default=1.0)
        radius = _positive_template_scalar(template, index=index, key="radius", default=0.2)
        if radius >= 0.5 * min(size_x, size_y):
            raise RenderValidationError(
                f"BEACH constraint error: mesh.templates[{index}] radius must be smaller "
                "than half of min(size_x, size_y)."
            )
        return

    if kind == "disk":
        _positive_template_scalar(template, index=index, key="radius", default=0.5)
        return

    if kind == "annulus":
        radius = _positive_template_scalar(template, index=index, key="radius", default=0.5)
        inner_radius = _nonnegative_template_scalar(
            template,
            index=index,
            key="inner_radius",
            default=0.25,
        )
        if inner_radius >= radius:
            raise RenderValidationError(
                f"BEACH constraint error: mesh.templates[{index}].inner_radius must be "
                "smaller than radius."
            )
        return

    if kind == "box":
        size = _maybe_vec3(
            template.get("size", [1.0, 1.0, 1.0]),
            name=f"mesh.templates[{index}].size",
        )
        if size is None:
            raise RenderValidationError(
                f"BEACH constraint error: mesh.templates[{index}].size must be a 3-element array."
            )
        if any(component <= 0.0 for component in size):
            raise RenderValidationError(
                f"BEACH constraint error: mesh.templates[{index}].size must be positive on all axes."
            )
        return

    if kind == "cylinder":
        _positive_template_scalar(template, index=index, key="radius", default=0.5)
        _positive_template_scalar(template, index=index, key="height", default=1.0)
        return

    if kind == "sphere":
        _positive_template_scalar(template, index=index, key="radius", default=0.5)
        return

    raise RenderValidationError(
        f"BEACH constraint error: mesh.templates[{index}] has unsupported kind={kind_value!r}."
    )


def _positive_template_scalar(
    template: Mapping[str, Any],
    *,
    index: int,
    key: str,
    default: float,
) -> float:
    value = _template_scalar(template, index=index, key=key, default=default)
    if value <= 0.0:
        raise RenderValidationError(
            f"BEACH constraint error: mesh.templates[{index}].{key} must be > 0."
        )
    return value


def _nonnegative_template_scalar(
    template: Mapping[str, Any],
    *,
    index: int,
    key: str,
    default: float,
) -> float:
    value = _template_scalar(template, index=index, key=key, default=default)
    if value < 0.0:
        raise RenderValidationError(
            f"BEACH constraint error: mesh.templates[{index}].{key} must be >= 0."
        )
    return value


def _template_scalar(
    template: Mapping[str, Any],
    *,
    index: int,
    key: str,
    default: float,
) -> float:
    raw = template.get(key, default)
    if not isinstance(raw, (int, float)) or isinstance(raw, bool):
        raise RenderValidationError(
            f"BEACH constraint error: mesh.templates[{index}].{key} must be numeric."
        )
    return float(raw)


def _maybe_vec3(value: object, *, name: str) -> list[float] | None:
    if value is None:
        return None
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        raise RenderValidationError(f"BEACH constraint error: {name} must be a 3-element array.")
    return [float(value[0]), float(value[1]), float(value[2])]


def _normalize_preset_section(section: str) -> str:
    normalized = section.strip().lower()
    aliases = {
        "sim": "sim",
        "mesh": "mesh",
        "output": "output",
        "species": "species",
        "particles.species": "species",
        "particle-species": "species",
        "template": "mesh.templates",
        "templates": "mesh.templates",
        "mesh.templates": "mesh.templates",
        "mesh-template": "mesh.templates",
    }
    if normalized not in aliases:
        raise CaseSpecError(
            "case spec error: unsupported --section. Use sim, mesh, output, "
            "species, or mesh.templates."
        )
    return aliases[normalized]


def _extract_array_of_tables(
    document: Mapping[str, Any],
    *,
    table_key: str,
    array_key: str,
    section_label: str,
) -> list[dict[str, Any]]:
    table = document.get(table_key)
    if not isinstance(table, Mapping):
        raise ConfigError(f"preset save error: rendered config does not contain [{table_key}].")
    value = table.get(array_key)
    if not isinstance(value, list) or not all(isinstance(item, Mapping) for item in value):
        raise ConfigError(f"preset save error: rendered config does not contain {section_label}.")
    return [dict(item) for item in value]


def _select_indexed_item(
    items: Sequence[Mapping[str, Any]],
    *,
    index: int | None,
    section_label: str,
) -> dict[str, Any]:
    if index is None:
        raise ConfigError(
            f"preset save error: --section {section_label} requires --index <1-based index>."
        )
    if index < 1 or index > len(items):
        raise ConfigError(
            f"preset save error: --index {index} is out of range for {section_label} "
            f"(count={len(items)})."
        )
    return copy.deepcopy(dict(items[index - 1]))


def _shadow_warnings(preset: ResolvedPreset) -> list[str]:
    if not preset.shadowed_sources:
        return []
    return [
        f'warning: preset "{preset.name}" resolved to {preset.source_label} '
        f"and shadows {source}"
        for source in preset.shadowed_sources
    ]


def _append_semantic_diff(
    lines: list[str],
    path: tuple[str | int, ...],
    left: Any,
    right: Any,
) -> None:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        ordered_keys: list[str] = list(left.keys())
        for key in right.keys():
            if key not in left:
                ordered_keys.append(key)
        for key in ordered_keys:
            in_left = key in left
            in_right = key in right
            next_path = (*path, key)
            if in_left and in_right:
                _append_semantic_diff(lines, next_path, left[key], right[key])
            elif in_left:
                lines.append(f"- {_format_path(next_path)} = {_summarize_value(left[key])}")
            else:
                lines.append(f"+ {_format_path(next_path)} = {_summarize_value(right[key])}")
        return

    if _is_array_of_tables(left) and _is_array_of_tables(right):
        left_items = list(left)
        right_items = list(right)
        limit = max(len(left_items), len(right_items))
        for index in range(limit):
            next_path = (*path, index)
            if index < len(left_items) and index < len(right_items):
                _append_semantic_diff(lines, next_path, left_items[index], right_items[index])
            elif index < len(left_items):
                lines.append(f"- {_format_path(next_path)} = {_summarize_value(left_items[index])}")
            else:
                lines.append(f"+ {_format_path(next_path)} = {_summarize_value(right_items[index])}")
        return

    if left != right:
        lines.append(
            f"{_format_path(path)}: {_summarize_value(left)} -> {_summarize_value(right)}"
        )


def _format_path(path: tuple[str | int, ...]) -> str:
    if not path:
        return "<root>"
    fragments: list[str] = []
    for part in path:
        if isinstance(part, int):
            fragments[-1] = f"{fragments[-1]}[{part}]"
        else:
            fragments.append(part)
    return ".".join(fragments)


def _summarize_value(value: Any) -> str:
    if isinstance(value, str):
        return repr(value)
    return str(value)


def _is_array_of_tables(value: object) -> bool:
    return isinstance(value, list) and len(value) > 0 and all(
        isinstance(item, Mapping) for item in value
    )
