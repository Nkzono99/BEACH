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
    "https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json"
)
_FRAGMENT_TOP_LEVEL_KEYS = frozenset({"sim", "particles", "mesh", "output"})
_PRESET_FORBIDDEN_KEYS = frozenset(
    {"schema_version", "title", "use_presets", "override", "base_case"}
)


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

    for preset_name in case_document.use_presets:
        preset = resolve_preset(preset_name, case_dir=case_dir)
        presets.append(preset)
        warnings.extend(_shadow_warnings(preset))
        config = merge_fragments(config, preset.data)

    config = merge_fragments(config, case_document.override)
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
    _require_table(final_config, "mesh", context="rendered config")
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

    header_comments = ["# BEACH preset fragment"]
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


def _maybe_vec3(value: object, *, name: str) -> list[float] | None:
    if value is None:
        return None
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        raise RenderValidationError(f"BEACH constraint error: {name} must be a 3-element array.")
    return [float(value[0]), float(value[1]), float(value[2])]


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
