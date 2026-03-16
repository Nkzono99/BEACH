"""Helpers for generating BEACH close-packed sphere configs."""

from __future__ import annotations

import copy
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SQRT2 = math.sqrt(2.0)
DEFAULT_FLOOR_Z = 0.02
DEFAULT_TOP_CLEARANCE_FACTOR = 30.0

DEFAULT_BASE_CONFIG: dict[str, Any] = {
    "sim": {
        "dt": 2.0e-8,
        "batch_count": 100,
        "max_step": 100000,
        "tol_rel": 1.0e-8,
        "softening": 1.0e-6,
        "b0": [0.0, 0.0, 0.0],
        "use_box": True,
        "box_min": [0.0, 0.0, 0.0],
        "box_max": [1.0, 1.0, 10.0],
        "bc_x_low": "periodic",
        "bc_x_high": "periodic",
        "bc_y_low": "periodic",
        "bc_y_high": "periodic",
        "bc_z_low": "open",
        "bc_z_high": "open",
        "rng_seed": 12345,
        "field_solver": "treecode",
        "batch_duration_step": 60000,
    },
    "particles": {
        "species": [
            {
                "source_mode": "reservoir_face",
                "number_density_cm3": 5.0,
                "temperature_ev": 10.0,
                "q_particle": -1.602176634e-19,
                "m_particle": 9.10938356e-31,
                "target_macro_particles_per_batch": 5000,
                "inject_face": "z_high",
                "pos_low": [0.0, 0.0, 10.0],
                "pos_high": [1.0, 1.0, 10.0],
                "drift_velocity": [0.0, 0.0, -4.0e5],
            },
            {
                "source_mode": "reservoir_face",
                "number_density_cm3": 5.0,
                "temperature_ev": 10.0,
                "q_particle": 1.602176634e-19,
                "m_particle": 1.672482821616e-27,
                "target_macro_particles_per_batch": -1,
                "inject_face": "z_high",
                "pos_low": [0.0, 0.0, 10.0],
                "pos_high": [1.0, 1.0, 10.0],
                "drift_velocity": [0.0, 0.0, -4.0e5],
            },
        ]
    },
    "mesh": {
        "mode": "template",
    },
    "output": {
        "write_files": True,
        "dir": "./outputs/latest",
        "history_stride": 10,
    },
}

_TOP_LEVEL_SECTION_ORDER = ("sim", "particles", "mesh", "output")


@dataclass(frozen=True, slots=True)
class ClosePackSpec:
    """Geometry and meshing inputs for the close-packed sphere generator."""

    layers: int
    radius: float
    cells_x: int
    cells_y: int
    floor_z: float = DEFAULT_FLOOR_Z
    plane_nx: int = 20
    plane_ny: int = 20
    sphere_n_lon: int = 14
    sphere_n_lat: int = 6
    box_height: float | None = None
    top_clearance: float | None = None

    def validate(self) -> None:
        if self.layers < 1:
            raise ValueError("layers must be >= 1.")
        if self.radius <= 0.0:
            raise ValueError("radius must be > 0.")
        if self.cells_x < 1 or self.cells_y < 1:
            raise ValueError("cells_x and cells_y must be >= 1.")
        if self.floor_z < 0.0:
            raise ValueError("floor_z must be >= 0.")
        if self.plane_nx < 1 or self.plane_ny < 1:
            raise ValueError("plane_nx and plane_ny must be >= 1.")
        if self.sphere_n_lon < 3:
            raise ValueError("sphere_n_lon must be >= 3.")
        if self.sphere_n_lat < 2:
            raise ValueError("sphere_n_lat must be >= 2.")
        if self.box_height is not None and self.box_height <= 0.0:
            raise ValueError("box_height must be > 0 when provided.")
        if self.top_clearance is not None and self.top_clearance <= 0.0:
            raise ValueError("top_clearance must be > 0 when provided.")
        if self.box_height is not None and self.top_clearance is not None:
            raise ValueError("box_height and top_clearance cannot be used together.")


def default_base_config() -> dict[str, Any]:
    """Return a deep copy of the bundled base config."""

    return copy.deepcopy(DEFAULT_BASE_CONFIG)


def load_base_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load a base config from TOML, or use the bundled defaults."""

    if path is None:
        return default_base_config()
    config = _load_toml(Path(path))
    if not isinstance(config, dict):
        raise ValueError("base config must decode to a TOML table.")
    return config


def unit_cell_pitch(radius: float) -> float:
    """Return the lateral period for the A/B square close-packing cell."""

    if radius <= 0.0:
        raise ValueError("radius must be > 0.")
    return radius * (2.0 + SQRT2)


def layer_spacing(radius: float) -> float:
    """Return the vertical separation between neighboring close-packed layers."""

    if radius <= 0.0:
        raise ValueError("radius must be > 0.")
    return radius * SQRT2


def total_sphere_count(spec: ClosePackSpec) -> int:
    """Return the number of sphere templates produced by a spec."""

    spec.validate()
    return 2 * spec.layers * spec.cells_x * spec.cells_y


def geometry_top_z(spec: ClosePackSpec) -> float:
    """Return the highest z coordinate touched by the stacked spheres."""

    spec.validate()
    return (
        spec.floor_z
        + 2.0 * spec.radius
        + (spec.layers - 1) * layer_spacing(spec.radius)
    )


def generate_closepack_sphere_templates(
    spec: ClosePackSpec,
    *,
    box_min: Sequence[float] = (0.0, 0.0, 0.0),
) -> list[dict[str, Any]]:
    """Generate sphere templates matching the current `exp_bench_tree` pattern."""

    spec.validate()
    box_min_xyz = _coerce_vec3(box_min, name="box_min")
    pitch = unit_cell_pitch(spec.radius)
    basis_a = (
        (spec.radius, spec.radius),
        (pitch - spec.radius, pitch - spec.radius),
    )
    basis_b = (
        (pitch - spec.radius, spec.radius),
        (spec.radius, pitch - spec.radius),
    )
    z0 = spec.floor_z + spec.radius
    dz = layer_spacing(spec.radius)
    templates: list[dict[str, Any]] = []
    for layer_index in range(spec.layers):
        basis = basis_a if layer_index % 2 == 0 else basis_b
        center_z = z0 + layer_index * dz
        for iy in range(spec.cells_y):
            y_shift = box_min_xyz[1] + iy * pitch
            for ix in range(spec.cells_x):
                x_shift = box_min_xyz[0] + ix * pitch
                for basis_x, basis_y in basis:
                    templates.append(
                        {
                            "kind": "sphere",
                            "enabled": True,
                            "radius": spec.radius,
                            "center": [x_shift + basis_x, y_shift + basis_y, center_z],
                            "n_lon": spec.sphere_n_lon,
                            "n_lat": spec.sphere_n_lat,
                        }
                    )
    return templates


def build_closepack_config(
    spec: ClosePackSpec,
    *,
    base_config: Mapping[str, Any] | None = None,
    output_dir: str | None = None,
) -> dict[str, Any]:
    """Build a BEACH config with a close-packed floor + sphere stack."""

    spec.validate()
    config = copy.deepcopy(dict(base_config) if base_config is not None else default_base_config())
    sim = _ensure_table(config, "sim")
    particles = _ensure_table(config, "particles")
    mesh = _ensure_table(config, "mesh")
    output = _ensure_table(config, "output")

    box_min = _coerce_vec3(sim.get("box_min", [0.0, 0.0, 0.0]), name="sim.box_min")
    box_max_existing = _coerce_vec3(
        sim.get("box_max", [0.0, 0.0, 0.0]),
        name="sim.box_max",
    )

    width_x = spec.cells_x * unit_cell_pitch(spec.radius)
    width_y = spec.cells_y * unit_cell_pitch(spec.radius)
    required_top_z = geometry_top_z(spec)
    box_max_z = _resolve_box_top_z(
        spec=spec,
        box_min_z=box_min[2],
        existing_box_top_z=box_max_existing[2],
        required_top_z=required_top_z,
    )
    box_max = [box_min[0] + width_x, box_min[1] + width_y, box_max_z]

    if spec.floor_z <= box_min[2]:
        raise ValueError("floor_z must be above sim.box_min[2].")
    if box_max_z <= required_top_z:
        raise ValueError("The generated box is not tall enough for the sphere stack.")

    sim["use_box"] = True
    sim["box_min"] = box_min
    sim["box_max"] = box_max

    species_list = particles.get("species", [])
    if not isinstance(species_list, list):
        raise ValueError("particles.species must be an array of tables.")
    for species in species_list:
        if not isinstance(species, dict):
            raise ValueError("particles.species entries must be TOML tables.")
        _fit_species_face_to_box(species, box_min=box_min, box_max=box_max)

    mesh.pop("obj_path", None)
    mesh["mode"] = "template"
    mesh["templates"] = [
        _build_floor_plane_template(spec=spec, box_min=box_min, width_x=width_x, width_y=width_y),
        *generate_closepack_sphere_templates(spec, box_min=box_min),
    ]

    if output_dir is not None:
        output["dir"] = output_dir

    return config


def render_closepack_toml(
    config: Mapping[str, Any],
    *,
    spec: ClosePackSpec | None = None,
    base_config_path: str | Path | None = None,
) -> str:
    """Render a generated config as TOML text."""

    lines: list[str] = ["# Generated by examples/generate_closepack_config.py"]
    if spec is not None:
        lines.append(
            "# closepack "
            f"layers={spec.layers} radius={spec.radius} "
            f"cells_x={spec.cells_x} cells_y={spec.cells_y} "
            f"spheres={total_sphere_count(spec)}"
        )
        lines.append(
            "# lattice "
            f"pitch={unit_cell_pitch(spec.radius)} "
            f"layer_spacing={layer_spacing(spec.radius)}"
        )
    if base_config_path is not None:
        lines.append(f"# base_config={base_config_path}")
    lines.append("")
    lines.extend(_render_toml_tables(config))
    return "\n".join(lines).rstrip() + "\n"


def _resolve_box_top_z(
    *,
    spec: ClosePackSpec,
    box_min_z: float,
    existing_box_top_z: float,
    required_top_z: float,
) -> float:
    if spec.box_height is not None:
        return box_min_z + spec.box_height
    clearance = (
        spec.top_clearance
        if spec.top_clearance is not None
        else DEFAULT_TOP_CLEARANCE_FACTOR * spec.radius
    )
    resolved = required_top_z + clearance
    if spec.top_clearance is None:
        resolved = max(resolved, existing_box_top_z)
    return resolved


def _build_floor_plane_template(
    *,
    spec: ClosePackSpec,
    box_min: Sequence[float],
    width_x: float,
    width_y: float,
) -> dict[str, Any]:
    return {
        "kind": "plane",
        "enabled": True,
        "size_x": width_x,
        "size_y": width_y,
        "nx": spec.plane_nx,
        "ny": spec.plane_ny,
        "center": [
            float(box_min[0]) + 0.5 * width_x,
            float(box_min[1]) + 0.5 * width_y,
            spec.floor_z,
        ],
    }


def _fit_species_face_to_box(
    species: dict[str, Any],
    *,
    box_min: Sequence[float],
    box_max: Sequence[float],
) -> None:
    inject_face = species.get("inject_face")
    if not isinstance(inject_face, str) or inject_face == "":
        return
    axis, boundary = _axis_and_boundary_for_face(inject_face, box_min, box_max)
    pos_low = [float(v) for v in box_min]
    pos_high = [float(v) for v in box_max]
    pos_low[axis] = boundary
    pos_high[axis] = boundary
    species["pos_low"] = pos_low
    species["pos_high"] = pos_high


def _axis_and_boundary_for_face(
    inject_face: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
) -> tuple[int, float]:
    if inject_face == "x_low":
        return 0, float(box_min[0])
    if inject_face == "x_high":
        return 0, float(box_max[0])
    if inject_face == "y_low":
        return 1, float(box_min[1])
    if inject_face == "y_high":
        return 1, float(box_max[1])
    if inject_face == "z_low":
        return 2, float(box_min[2])
    if inject_face == "z_high":
        return 2, float(box_max[2])
    raise ValueError(f"unknown inject_face: {inject_face}")


def _ensure_table(config: dict[str, Any], key: str) -> dict[str, Any]:
    value = config.get(key)
    if value is None:
        table: dict[str, Any] = {}
        config[key] = table
        return table
    if not isinstance(value, dict):
        raise ValueError(f"{key} must be a TOML table.")
    return value


def _coerce_vec3(value: object, *, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a 3-element array.")
    if len(value) != 3:
        raise ValueError(f"{name} must be a 3-element array.")
    return [float(value[0]), float(value[1]), float(value[2])]


def _load_toml(path: Path) -> dict[str, Any]:
    try:
        import tomllib  # py311+

        with path.open("rb") as stream:
            return tomllib.load(stream)
    except ModuleNotFoundError:
        try:
            import tomli  # type: ignore

            with path.open("rb") as stream:
                return tomli.load(stream)
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "TOML parser is missing. Use Python 3.11+ or install tomli: "
                "`python -m pip install tomli`."
            ) from exc


def _render_toml_tables(config: Mapping[str, Any]) -> list[str]:
    lines: list[str] = []
    remaining_keys = [key for key in config if key not in _TOP_LEVEL_SECTION_ORDER]
    for key in [*_TOP_LEVEL_SECTION_ORDER, *remaining_keys]:
        if key not in config:
            continue
        value = config[key]
        if not isinstance(value, Mapping):
            raise TypeError(f"Top-level TOML key {key!r} must be a table.")
        if lines:
            lines.append("")
        _write_table(lines, (key,), value)
    return lines


def _write_table(lines: list[str], path: tuple[str, ...], table: Mapping[str, Any]) -> None:
    lines.append(f"[{'.'.join(path)}]")
    nested_entries: list[tuple[str, str, Any]] = []
    for key, value in table.items():
        if isinstance(value, Mapping):
            nested_entries.append(("table", key, value))
            continue
        if _is_array_of_tables(value):
            nested_entries.append(("array", key, value))
            continue
        lines.append(f"{key} = {_format_toml_value(value)}")
    for kind, key, value in nested_entries:
        lines.append("")
        if kind == "table":
            _write_table(lines, (*path, key), value)
        else:
            _write_array_tables(lines, (*path, key), value)


def _write_array_tables(
    lines: list[str],
    path: tuple[str, ...],
    items: Sequence[Mapping[str, Any]],
) -> None:
    for index, item in enumerate(items):
        if index > 0:
            lines.append("")
        lines.append(f"[[{'.'.join(path)}]]")
        nested_entries: list[tuple[str, str, Any]] = []
        for key, value in item.items():
            if isinstance(value, Mapping):
                nested_entries.append(("table", key, value))
                continue
            if _is_array_of_tables(value):
                nested_entries.append(("array", key, value))
                continue
            lines.append(f"{key} = {_format_toml_value(value)}")
        for kind, key, value in nested_entries:
            lines.append("")
            if kind == "table":
                _write_table(lines, (*path, key), value)
            else:
                _write_array_tables(lines, (*path, key), value)


def _is_array_of_tables(value: object) -> bool:
    return isinstance(value, list) and len(value) > 0 and all(
        isinstance(item, Mapping) for item in value
    )


def _format_toml_value(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped}"'
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return "[" + ", ".join(_format_toml_value(item) for item in value) + "]"
    raise TypeError(f"Unsupported TOML value: {value!r}")
