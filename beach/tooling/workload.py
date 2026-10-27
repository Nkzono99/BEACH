"""Particle-workload estimation for BEACH simulation configs."""

from __future__ import annotations

import copy
import csv
import math
from pathlib import Path
from typing import Any

from beach.config import ConfigError, normalize_high_level_config
from beach.config._toml import load_toml_file
from beach.config._layout import to_runtime_layout
from beach.config.schema import schema_definition_defaults, schema_definition_property_names
from beach.io.summary import load_summary_file

K_BOLTZMANN = 1.380649e-23
EV_TO_K = 1.160451812e4

# The estimator reads these [sim] and species keys. Their defaults come from the
# packaged schema, so the schema stays the single source with the Fortran types.
# The fallbacks cover keys the schema declares without a default: the box and
# field-boundary keys are resolved from [domain] and [field_boundary], one of
# the two batch-duration forms is required for dynamic sources, and the
# remaining keys are mode-specific inputs that have no neutral value.
_SIM_KEYS = (
    "batch_count",
    "batch_duration",
    "batch_duration_step",
    "particle_distribution",
    "field_solver",
    "field_bc_mode",
    "field_periodic_image_layers",
    "field_periodic_far_correction",
    "field_periodic_ewald_alpha",
    "field_periodic_ewald_layers",
    "field_periodic_cache_dir",
    "field_periodic_generation_tolerance",
    "tree_theta",
    "tree_leaf_max",
    "tree_min_nelem",
    "e0",
    "e0_abs",
    "e0_phi_xy_deg",
    "e0_phi_z_deg",
    "use_box",
    "box_min",
    "box_max",
)
_SIM_FALLBACKS: dict[str, Any] = {
    "batch_duration": 0.0,
    "batch_duration_step": 0.0,
    "field_bc_mode": "free",
    "field_periodic_cache_dir": "",
    "tree_theta": 0.5,
    "tree_leaf_max": 16,
    "e0_abs": 0.0,
    "use_box": False,
    "box_min": [-1.0, -1.0, -1.0],
    "box_max": [1.0, 1.0, 1.0],
}
_SPECIES_KEYS = (
    "species_key",
    "enabled",
    "source_mode",
    "npcls_per_step",
    "temperature_k",
    "q_particle",
    "m_particle",
    "w_particle",
    "target_macro_particles_per_batch",
    "pos_low",
    "pos_high",
    "drift_velocity",
    "emit_current_density_a_m2",
    "rays_per_batch",
    "deposit_opposite_charge_on_emit",
    "normal_drift_speed",
    "ray_direction",
    "inject_face",
    "boundary_inflow",
)
_SPECIES_FALLBACKS: dict[str, Any] = {
    "species_key": "",
    "target_macro_particles_per_batch": 0,
    "emit_current_density_a_m2": 0.0,
    "rays_per_batch": 0,
    "ray_direction": [0.0, 0.0, 0.0],
    "inject_face": "",
    "boundary_inflow": {},
}


def _defaults_from_schema(
    definition: str, keys: tuple[str, ...], fallbacks: dict[str, Any]
) -> dict[str, Any]:
    schema_defaults = schema_definition_defaults(definition)
    defaults: dict[str, Any] = {}
    for key in keys:
        if key in schema_defaults:
            defaults[key] = schema_defaults[key]
        elif key in fallbacks:
            defaults[key] = copy.deepcopy(fallbacks[key])
        else:
            raise ValueError(f"BEACH schema has no default for $defs.{definition}.{key}")
    return defaults


DEFAULT_SIM: dict[str, Any] = _defaults_from_schema("sim", _SIM_KEYS, _SIM_FALLBACKS)
DEFAULT_SPECIES: dict[str, Any] = _defaults_from_schema(
    "species", _SPECIES_KEYS, _SPECIES_FALLBACKS
)

ALLOWED_SIM_KEYS = schema_definition_property_names("sim")
_LEGACY_SIM_BOUNDARY_KEYS = frozenset(
    {
        "field_bc_mode",
        "use_box",
        "box_min",
        "box_max",
        "box_origin",
        "box_size",
        "bc_x_low",
        "bc_x_high",
        "bc_y_low",
        "bc_y_high",
        "bc_z_low",
        "bc_z_high",
    }
)
ALLOWED_SPECIES_KEYS = schema_definition_property_names("species")


class _MacroResiduals(list[float]):
    """Source residuals with optional six-face residuals from restart CSV."""

    def __init__(
        self,
        source: list[float],
        boundary: list[list[float]],
    ) -> None:
        super().__init__(source)
        self.boundary = boundary


def _rank_slot_range(total_slots: int, rank: int, n_ranks: int) -> tuple[int, int]:
    """Return the half-open global slot range a static MPI rank generates.

    The runtime splits the batch's total slot count once (``mpi_split_count``)
    and gives each rank one contiguous range (``init_particle_batch_from_config``).
    """

    base, rem = divmod(total_slots, n_ranks)
    begin = rank * base + min(rank, rem)
    return begin, begin + base + (1 if rank < rem else 0)


def _interleaved_prefix_counts(lengths: list[int], prefix: int) -> list[int]:
    """Count each segment's slots among the first ``prefix`` global slots.

    ``build_batch_injection_layout`` deals slot ``k`` of every segment that still
    holds ``k`` slots, in segment order, before slot ``k + 1`` of any segment. The
    first ``r`` complete rounds hold ``sum(min(n, r))`` slots, so the counts follow
    from the segment lengths without building the global slot list.
    """

    lo, hi = 0, max(lengths, default=0)
    while lo < hi:
        rounds = (lo + hi + 1) // 2
        if sum(min(n, rounds) for n in lengths) <= prefix:
            lo = rounds
        else:
            hi = rounds - 1
    leftover = prefix - sum(min(n, lo) for n in lengths)
    counts: list[int] = []
    for n in lengths:
        count = min(n, lo)
        if leftover > 0 and n > lo:
            count += 1
            leftover -= 1
        counts.append(count)
    return counts


def load_toml(path: Path) -> dict[str, Any]:
    """Load a TOML configuration through the shared BEACH config reader."""

    return to_runtime_layout(load_toml_file(path))


def _dot(a: list[float], b: list[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _standard_normal_pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def _standard_normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


# Box faces: normal axis, whether the face is the upper bound, and the two
# tangential axes of the face area.
_FACES: dict[str, tuple[int, bool, tuple[int, int]]] = {
    "x_low": (0, False, (1, 2)),
    "x_high": (0, True, (1, 2)),
    "y_low": (1, False, (2, 0)),
    "y_high": (1, True, (2, 0)),
    "z_low": (2, False, (0, 1)),
    "z_high": (2, True, (0, 1)),
}


def _face(inject_face: str) -> tuple[int, bool, tuple[int, int]]:
    try:
        return _FACES[inject_face]
    except KeyError:
        raise ValueError(f"unknown inject_face: {inject_face}") from None


def _inward_normal_for_face(inject_face: str) -> list[float]:
    axis, high, _ = _face(inject_face)
    normal = [0.0, 0.0, 0.0]
    normal[axis] = -1.0 if high else 1.0
    return normal


def _face_tangential_axes(inject_face: str) -> tuple[int, int]:
    return _face(inject_face)[2]


def _species_temperature_k(spec: dict[str, Any]) -> float:
    if bool(spec.get("_has_temperature_ev", False)):
        return float(spec["temperature_ev"]) * EV_TO_K
    return float(spec.get("temperature_k", DEFAULT_SPECIES["temperature_k"]))


def _species_density_m3(spec: dict[str, Any]) -> float:
    if bool(spec.get("_has_number_density_cm3", False)):
        return float(spec["number_density_cm3"]) * 1.0e6
    return float(spec["number_density_m3"])


def _check_temperature_inputs(spec: dict[str, Any]) -> None:
    has_temp_k = bool(spec.get("_has_temperature_k", False))
    has_temp_ev = bool(spec.get("_has_temperature_ev", False))
    if has_temp_k and has_temp_ev:
        raise SystemExit("Specify either temperature_ev or temperature_k, not both.")
    if has_temp_ev and float(spec["temperature_ev"]) < 0.0:
        raise SystemExit("temperature_ev must be >= 0.")
    if has_temp_k and float(spec["temperature_k"]) < 0.0:
        raise SystemExit("temperature_k must be >= 0.")


def _checked_particle_mass(spec: dict[str, Any]) -> float:
    m_particle = float(spec.get("m_particle", DEFAULT_SPECIES["m_particle"]))
    if m_particle <= 0.0:
        raise SystemExit("m_particle must be > 0.")
    return m_particle


def _compute_inflow_flux(
    number_density_m3: float,
    temperature_k: float,
    m_particle: float,
    drift_velocity: list[float],
    inward_normal: list[float],
) -> float:
    u_n = _dot(drift_velocity, inward_normal)
    sigma = math.sqrt(K_BOLTZMANN * temperature_k / m_particle)
    if sigma <= 0.0:
        return number_density_m3 * max(0.0, u_n)
    alpha = u_n / sigma
    return number_density_m3 * (
        sigma * _standard_normal_pdf(alpha) + u_n * _standard_normal_cdf(alpha)
    )


def _compute_face_area(
    inject_face: str, pos_low: list[float], pos_high: list[float]
) -> float:
    axis_t1, axis_t2 = _face_tangential_axes(inject_face)
    return (pos_high[axis_t1] - pos_low[axis_t1]) * (
        pos_high[axis_t2] - pos_low[axis_t2]
    )


def _parse_boundary_inflow_surfaces(
    sim_cfg: dict[str, Any], spec: dict[str, Any]
) -> list[dict[str, Any]]:
    """Resolve only the geometry needed for reservoir particle counts."""

    raw_inflow = spec.get("boundary_inflow", {})
    if not isinstance(raw_inflow, dict):
        raise SystemExit("particles.species.boundary_inflow must be a table.")
    unknown = set(raw_inflow) - _FACES.keys()
    if unknown:
        raise SystemExit("Unknown key in particles.species.boundary_inflow: " + sorted(unknown)[0])
    if any(value != "reservoir" for value in raw_inflow.values()):
        raise SystemExit('particles.species.boundary_inflow values must be "reservoir".')
    box_min = [float(x) for x in sim_cfg["box_min"]]
    box_max = [float(x) for x in sim_cfg["box_max"]]
    return [
        {
            "face": face,
            "normal": _inward_normal_for_face(face),
            "area": _compute_face_area(face, box_min, box_max),
        }
        for face in sorted(raw_inflow)
    ]


def _reservoir_parameters(
    sim_cfg: dict[str, Any],
    spec: dict[str, Any],
    resolved_batch_duration: float,
    species_idx: int,
    reservoir_params_by_species: list[dict[str, Any] | None],
) -> dict[str, Any]:
    """Resolve physical rates and the unique macro weight needed for division.

    Only inputs that make that arithmetic undefined or ambiguous are checked
    here. Simulation-option combinations are validated by Fortran preflight.
    """
    if not bool(sim_cfg.get("use_box", False)):
        raise SystemExit("reservoir injection requires a [domain] box.")
    surfaces = _parse_boundary_inflow_surfaces(sim_cfg, spec)

    has_w_particle = bool(spec.get("_has_w_particle", False))
    has_target_macro = bool(spec.get("_has_target_macro_particles_per_batch", False))
    if has_w_particle and has_target_macro:
        raise SystemExit(
            "reservoir injection does not allow both w_particle and target_macro_particles_per_batch."
        )
    if (not has_w_particle) and (not has_target_macro):
        raise SystemExit(
            "reservoir injection requires either w_particle or target_macro_particles_per_batch."
        )
    if has_w_particle:
        w_particle = float(spec.get("w_particle", DEFAULT_SPECIES["w_particle"]))
        if w_particle <= 0.0:
            raise SystemExit(
                "particles.species.w_particle must be > 0 for boundary_inflow"
            )
    if has_target_macro:
        target_macro = int(spec.get("target_macro_particles_per_batch", 0))
        if target_macro == 0 or target_macro < -1:
            raise SystemExit(
                "particles.species.target_macro_particles_per_batch must be > 0 or -1."
            )
        if target_macro == -1:
            if species_idx == 0:
                raise SystemExit(
                    "particles.species[1].target_macro_particles_per_batch cannot be -1."
                )
            params1 = reservoir_params_by_species[0]
            if params1 is None:
                raise SystemExit(
                    "target_macro_particles_per_batch=-1 requires species[1] reservoir injection with a resolved weight."
                )
            w_particle = float(params1["w_particle"])

    m_particle = _checked_particle_mass(spec)

    has_cm3 = bool(spec.get("_has_number_density_cm3", False))
    has_m3 = bool(spec.get("_has_number_density_m3", False))
    if has_cm3 and has_m3:
        raise SystemExit(
            "Specify either number_density_cm3 or number_density_m3, not both."
        )
    if (not has_cm3) and (not has_m3):
        raise SystemExit(
            "reservoir injection requires number_density_cm3 or number_density_m3."
        )
    if has_cm3 and float(spec["number_density_cm3"]) <= 0.0:
        raise SystemExit("number_density_cm3 must be > 0.")
    if has_m3 and float(spec["number_density_m3"]) <= 0.0:
        raise SystemExit("number_density_m3 must be > 0.")

    _check_temperature_inputs(spec)

    drift_velocity = [
        float(x) for x in spec.get("drift_velocity", DEFAULT_SPECIES["drift_velocity"])
    ]
    temperature_k = _species_temperature_k(spec)
    number_density_m3 = _species_density_m3(spec)
    surface_fluxes = [
        _compute_inflow_flux(
            number_density_m3=number_density_m3,
            temperature_k=temperature_k,
            m_particle=m_particle,
            drift_velocity=drift_velocity,
            inward_normal=list(surface["normal"]),
        )
        for surface in surfaces
    ]
    surface_rates = [
        surface_flux * float(surface["area"])
        for surface_flux, surface in zip(surface_fluxes, surfaces, strict=True)
    ]
    rate_per_s = sum(surface_rates)
    if has_target_macro and target_macro != -1:
        w_particle = rate_per_s * resolved_batch_duration / float(target_macro)
        if (not math.isfinite(w_particle)) or w_particle <= 0.0:
            raise SystemExit(
                "target_macro_particles_per_batch produced invalid w_particle."
            )

    return {
        "surfaces": surfaces,
        "w_particle": w_particle,
        "surface_rates": surface_rates,
    }


def _resolve_batch_duration(
    sim_cfg: dict[str, Any],
    sim_raw: dict[str, Any],
    has_dynamic_source_species: bool,
) -> float:
    has_batch_duration = "batch_duration" in sim_raw
    has_batch_duration_step = "batch_duration_step" in sim_raw
    if has_batch_duration and has_batch_duration_step:
        raise SystemExit(
            "sim.batch_duration and sim.batch_duration_step cannot be used together."
        )
    if has_batch_duration_step:
        batch_duration_step = float(sim_cfg.get("batch_duration_step", 0.0))
        if (not math.isfinite(batch_duration_step)) or batch_duration_step <= 0.0:
            raise SystemExit("sim.batch_duration_step must be > 0.")
        dt = float(sim_cfg.get("dt", 0.0))
        if (not math.isfinite(dt)) or dt <= 0.0:
            raise SystemExit("sim.dt must be > 0 when sim.batch_duration_step is set.")
        batch_duration = dt * batch_duration_step
        if (not math.isfinite(batch_duration)) or batch_duration <= 0.0:
            raise SystemExit(
                "sim.batch_duration_step produced invalid sim.batch_duration."
            )
        return batch_duration

    batch_duration = float(sim_cfg.get("batch_duration", 0.0))
    if not math.isfinite(batch_duration):
        raise SystemExit("sim.batch_duration must be finite.")
    if has_dynamic_source_species and batch_duration <= 0.0:
        raise SystemExit("sim.batch_duration must be > 0 for dynamic sources")
    return batch_duration


def read_macro_residuals(path: Path | None, n_species: int) -> list[float]:
    """Load resume residuals from ``macro_residuals.csv``.

    Parameters
    ----------
    path : pathlib.Path or None
        CSV path. ``None`` returns zeros.
    n_species : int
        Number of species in current configuration.

    Returns
    -------
    list[float]
        Per-species residual values.

    Raises
    ------
    SystemExit
        If ``path`` is provided but does not exist.
    """

    residuals = [0.0] * n_species
    boundary_residuals = [[0.0] * n_species for _ in range(6)]
    if path is None:
        return _MacroResiduals(residuals, boundary_residuals)
    if path.name.startswith("macro_residuals_rank") and path.suffix == ".csv":
        raise SystemExit(
            "legacy rank-local macro residual files are not supported; "
            "use the global macro_residuals.csv"
        )
    if not path.exists():
        raise SystemExit(f"macro residual file not found: {path}")

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            idx = int(row["species_idx"])
            face = int(row.get("face") or 0)
            if 1 <= idx <= n_species and face == 0:
                residuals[idx - 1] = float(row["residual"])
            elif 1 <= idx <= n_species and 1 <= face <= 6:
                boundary_residuals[face - 1][idx - 1] = float(row["residual"])
    return _MacroResiduals(residuals, boundary_residuals)


def read_summary_batches(path: Path) -> int:
    """Load completed batch count from a Fortran ``summary.txt`` file."""

    if not path.exists():
        raise SystemExit(f"summary file not found: {path}")
    try:
        values = load_summary_file(path)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    try:
        batches = int(values["batches"])
    except KeyError as exc:
        raise SystemExit(f"summary file is missing batches: {path}") from exc
    except ValueError as exc:
        raise SystemExit(f"summary batches must be an integer: {path}") from exc
    if batches < 0:
        raise SystemExit(f"summary batches must be >= 0: {path}")
    return batches


def completed_batches_from_resume_config(config: dict[str, Any]) -> int:
    """Return checkpoint batch count when config requests a restart."""

    output_raw = config.get("output", {})
    if not isinstance(output_raw, dict):
        return 0
    if not bool(output_raw.get("resume", False)):
        return 0
    checkpoint_dir_raw = output_raw.get("restart_from") or output_raw.get(
        "dir", "outputs/latest"
    )
    checkpoint_dir = Path(str(checkpoint_dir_raw))
    return read_summary_batches(checkpoint_dir / "summary.txt")


def _normalized_sections(config: dict[str, Any]) -> dict[str, Any]:
    """Check the sections the estimator reads and resolve their high-level notation."""

    sim_raw = config.get("sim", {})
    if not isinstance(sim_raw, dict):
        raise SystemExit("[sim] section must be a table.")
    unknown_sim_keys = sorted(
        set(sim_raw) - ALLOWED_SIM_KEYS - _LEGACY_SIM_BOUNDARY_KEYS
    )
    if unknown_sim_keys:
        raise SystemExit(f"Unknown key in [sim]: {unknown_sim_keys[0]}")
    domain_raw = config.get("domain")
    if domain_raw is not None and not isinstance(domain_raw, dict):
        raise SystemExit("[domain] section must be a table.")
    field_boundary_raw = config.get("field_boundary")
    if field_boundary_raw is not None and not isinstance(field_boundary_raw, dict):
        raise SystemExit("[field_boundary] section must be a table.")
    particle_boundary_raw = config.get("particle_boundary")
    if particle_boundary_raw is not None and not isinstance(
        particle_boundary_raw, dict
    ):
        raise SystemExit("[particle_boundary] section must be a table.")

    particles_raw = config.get("particles", {})
    if not isinstance(particles_raw, dict):
        raise SystemExit("[particles] section must be a table.")
    unknown_particles_keys = sorted(set(particles_raw) - {"species", "secondary_emission_max_generations"})
    if unknown_particles_keys:
        raise SystemExit(f"Unknown key in [particles]: {unknown_particles_keys[0]}")

    species_list_raw = particles_raw.get("species", [])
    if not isinstance(species_list_raw, list) or len(species_list_raw) == 0:
        raise SystemExit("At least one [[particles.species]] entry is required.")
    for raw in species_list_raw:
        if not isinstance(raw, dict):
            raise SystemExit("Each [[particles.species]] entry must be a table.")
        unknown_species_keys = sorted(set(raw) - ALLOWED_SPECIES_KEYS)
        if unknown_species_keys:
            raise SystemExit(
                f"Unknown key in [[particles.species]]: {unknown_species_keys[0]}"
            )

    workload_config = {"sim": sim_raw, "particles": particles_raw}
    if domain_raw is not None:
        workload_config["domain"] = domain_raw
    if field_boundary_raw is not None:
        workload_config["field_boundary"] = field_boundary_raw
    if particle_boundary_raw is not None:
        workload_config["particle_boundary"] = particle_boundary_raw
    try:
        return normalize_high_level_config(workload_config)
    except (ConfigError, TypeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc


def _resolved_sim_config(normalized_config: dict[str, Any]) -> dict[str, Any]:
    """Merge [sim] over schema defaults and resolve the counting geometry."""

    sim_raw = normalized_config["sim"]
    domain_raw = normalized_config.get("domain")
    sim_cfg = dict(DEFAULT_SIM)
    sim_cfg.update(sim_raw)
    if isinstance(domain_raw, dict):
        sim_cfg["use_box"] = True
        sim_cfg["box_min"] = domain_raw.get("box_min", DEFAULT_SIM["box_min"])
        sim_cfg["box_max"] = domain_raw.get("box_max", DEFAULT_SIM["box_max"])
    return sim_cfg


def _species_with_defaults(species_list_raw: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Merge each species over the defaults and record which inputs were given."""

    species_list: list[dict[str, Any]] = []
    for raw in species_list_raw:
        spec = dict(DEFAULT_SPECIES)
        spec.update(raw)
        spec["source_mode"] = (
            str(spec.get("source_mode", "volume_seed")).strip().lower()
        )
        spec["_has_number_density_cm3"] = "number_density_cm3" in raw
        spec["_has_number_density_m3"] = "number_density_m3" in raw
        spec["_has_temperature_k"] = "temperature_k" in raw
        spec["_has_temperature_ev"] = "temperature_ev" in raw
        spec["_has_w_particle"] = "w_particle" in raw
        spec["_has_target_macro_particles_per_batch"] = (
            "target_macro_particles_per_batch" in raw
        )
        species_list.append(spec)
    return species_list


def _remaining_batch_count(sim_cfg: dict[str, Any], completed_batches: int) -> tuple[int, int]:
    target_batch_count = int(sim_cfg["batch_count"])
    if target_batch_count <= 0:
        raise SystemExit("sim.batch_count must be > 0")
    if completed_batches < 0:
        raise SystemExit("completed_batches must be >= 0")
    if completed_batches > target_batch_count:
        raise SystemExit(
            "sim.batch_count must be >= completed checkpoint batches when resuming."
        )
    return target_batch_count, target_batch_count - completed_batches


def _check_parallel_layout(threads: int, mpi_ranks: int, mpi_rank: int) -> None:
    if threads <= 0:
        raise SystemExit("threads must be > 0")
    if mpi_ranks <= 0:
        raise SystemExit("mpi_ranks must be > 0")
    if mpi_rank < 0 or mpi_rank >= mpi_ranks:
        raise SystemExit("mpi_rank must satisfy 0 <= mpi_rank < mpi_ranks")


def _has_dynamic_sources(species_list: list[dict[str, Any]]) -> bool:
    """Report which sources need time to turn a physical flux into a count."""

    has_dynamic = False
    for spec in species_list:
        if not spec["enabled"]:
            continue
        mode = spec["source_mode"]
        if mode not in {"volume_seed", "photo_raycast", "secondary_emission"}:
            raise SystemExit(f"Unknown particles.species.source_mode: {mode}")
        if mode == "volume_seed" and int(spec["npcls_per_step"]) < 0:
            raise SystemExit("particles.species.npcls_per_step must be >= 0")
        has_dynamic |= mode == "photo_raycast" or bool(spec.get("boundary_inflow"))
    return has_dynamic


def _source_parameters(
    sim_cfg: dict[str, Any],
    species_list: list[dict[str, Any]],
    resolved_batch_duration: float,
) -> list[dict[str, Any] | None]:
    """Resolve the rates and weights the reservoir counter consumes.

    Photo ray counts come directly from rays_per_batch; emission phase-space,
    particle boundary actions, and source-option combinations belong to preflight.
    """

    parameters: list[dict[str, Any] | None] = [None] * len(species_list)
    for index, spec in enumerate(species_list):
        if spec["enabled"] and spec.get("boundary_inflow"):
            parameters[index] = _reservoir_parameters(
                sim_cfg, spec, resolved_batch_duration, index, parameters
            )
    return parameters


_FACE_INDEX = {"x_low": 0, "x_high": 1, "y_low": 2, "y_high": 3, "z_low": 4, "z_high": 5}


def _initial_residual_state(
    n_species: int,
    reservoir_params_by_species: list[dict[str, Any] | None],
    initial_residuals: list[float] | None,
) -> tuple[list[float], list[list[float]]]:
    """Return per-species source residuals and per-surface reservoir residuals."""

    residuals = [0.0] * n_species
    boundary_residuals = [
        [0.0] * len(params["surface_rates"]) if params is not None else []
        for params in reservoir_params_by_species
    ]
    if initial_residuals is not None:
        if len(initial_residuals) != n_species:
            raise SystemExit(
                "macro residual species count does not match config species count"
            )
        residuals = list(initial_residuals)
        restart_boundary = getattr(initial_residuals, "boundary", None)
        if restart_boundary is not None:
            for species_idx, params in enumerate(reservoir_params_by_species):
                if params is None:
                    continue
                boundary_residuals[species_idx] = [
                    float(restart_boundary[_FACE_INDEX[surface["face"]]][species_idx])
                    for surface in params["surfaces"]
                ]
    return residuals, boundary_residuals


def _reservoir_batch_counts(
    reservoir_params: dict[str, Any],
    face_residuals: list[float],
    resolved_batch_duration: float,
) -> list[int]:
    """Draw one batch of global reservoir macros per surface; update the residuals in place."""

    face_counts: list[int] = []
    for face_idx, rate_per_s in enumerate(reservoir_params["surface_rates"]):
        n_macro_expected = (
            rate_per_s * resolved_batch_duration / reservoir_params["w_particle"]
        )
        macro_budget = face_residuals[face_idx] + n_macro_expected
        if macro_budget < 0.0:
            macro_budget = 0.0
        n_face_global = math.floor(macro_budget)
        face_residuals[face_idx] = macro_budget - n_face_global
        face_counts.append(int(n_face_global))
    return face_counts


def _source_slot_count(spec: dict[str, Any]) -> int:
    """Return the global source slots of one enabled species for one batch."""

    source_mode = spec["source_mode"]
    if source_mode == "volume_seed":
        return int(spec.get("npcls_per_step", 0))
    if source_mode == "photo_raycast":
        n_rays = int(spec["rays_per_batch"])
        if n_rays < 0:
            raise SystemExit("rays_per_batch must be >= 0")
        return n_rays
    return 0


def estimate_workload(
    config: dict[str, Any],
    threads: int,
    initial_residuals: list[float] | None = None,
    mpi_ranks: int = 1,
    mpi_rank: int = 0,
    completed_batches: int = 0,
) -> dict[str, Any]:
    """Estimate particle counts for one local MPI rank without a native library.

    Use ``beachx lint`` to validate the simulation before estimating it. This
    function resolves the count arithmetic, not source/boundary/field physics
    combinations. Photo ray counts are an upper bound before ray intersections.

    Local counts follow the runtime's ``static`` particle distribution: the
    batch's global slots are interleaved round-robin over the species source and
    reservoir-face segments, and each rank takes one contiguous slot range. With
    ``particle_distribution = "shared"`` ranks claim slot chunks dynamically, so
    the same local counts are only an even-load assumption, not the exact
    per-rank numbers; ``rank_count_basis`` in the result says which applies.

    Parameters
    ----------
    config : dict[str, Any]
        Parsed TOML-like configuration dictionary.
    threads : int
        OpenMP thread count used for per-thread load bounds.
    initial_residuals : list[float] or None, default None
        Initial macro residuals for each species.
    mpi_ranks : int, default 1
        MPI world size.
    mpi_rank : int, default 0
        Local rank index whose slot range is estimated.
    completed_batches : int, default 0
        Completed checkpoint batches. ``sim.batch_count`` is treated as the
        cumulative target, so only ``sim.batch_count - completed_batches``
        batches are estimated.

    Returns
    -------
    dict[str, Any]
        Workload summary including per-batch totals, per-thread bounds, and
        resolved batch duration.

    Raises
    ------
    SystemExit
        If the particle-count arithmetic is undefined, ambiguous, or non-finite.
    """

    normalized_config = _normalized_sections(config)
    sim_cfg = _resolved_sim_config(normalized_config)
    species_list = _species_with_defaults(normalized_config["particles"]["species"])
    target_batch_count, batch_count = _remaining_batch_count(sim_cfg, completed_batches)
    _check_parallel_layout(threads, mpi_ranks, mpi_rank)
    particle_distribution = str(sim_cfg["particle_distribution"]).strip().lower()

    resolved_batch_duration = _resolve_batch_duration(
        sim_cfg=sim_cfg,
        sim_raw=normalized_config["sim"],
        has_dynamic_source_species=_has_dynamic_sources(species_list),
    )
    reservoir_params_by_species = _source_parameters(
        sim_cfg, species_list, resolved_batch_duration
    )
    residuals, boundary_residuals = _initial_residual_state(
        len(species_list), reservoir_params_by_species, initial_residuals
    )

    batch_totals: list[int] = []
    batch_thread_min: list[int] = []
    batch_thread_max: list[int] = []
    species_per_batch: list[list[int]] = []
    global_batch_totals: list[int] = []
    global_species_per_batch: list[list[int]] = []
    local_reservoir_particles = 0
    global_reservoir_particles = 0

    for _batch_idx in range(batch_count):
        global_species_counts = [0] * len(species_list)
        segment_species: list[int] = []
        segment_is_reservoir: list[bool] = []
        segment_lengths: list[int] = []
        for s_idx, spec in enumerate(species_list):
            if not bool(spec.get("enabled", True)):
                continue
            reservoir_params = reservoir_params_by_species[s_idx]
            source_count = _source_slot_count(spec)
            face_counts = (
                _reservoir_batch_counts(
                    reservoir_params, boundary_residuals[s_idx], resolved_batch_duration
                )
                if reservoir_params is not None
                else []
            )
            # Segments follow the runtime order, the species source before its
            # faces; the face order within a species does not change these counts.
            for count, is_reservoir in [(source_count, False)] + [
                (face_count, True) for face_count in face_counts
            ]:
                if count > 0:
                    segment_species.append(s_idx)
                    segment_is_reservoir.append(is_reservoir)
                    segment_lengths.append(count)
            global_species_counts[s_idx] = source_count + sum(face_counts)
            global_reservoir_particles += sum(face_counts)

        slot_begin, slot_end = _rank_slot_range(
            sum(segment_lengths), mpi_rank, mpi_ranks
        )
        species_counts = [0] * len(species_list)
        for s_idx, is_reservoir, end_count, begin_count in zip(
            segment_species,
            segment_is_reservoir,
            _interleaved_prefix_counts(segment_lengths, slot_end),
            _interleaved_prefix_counts(segment_lengths, slot_begin),
            strict=True,
        ):
            species_counts[s_idx] += end_count - begin_count
            if is_reservoir:
                local_reservoir_particles += end_count - begin_count

        batch_total = sum(species_counts)
        q, r = divmod(batch_total, threads)
        batch_totals.append(batch_total)
        global_batch_totals.append(sum(global_species_counts))
        batch_thread_min.append(q)
        batch_thread_max.append(q + 1 if r > 0 else q)
        species_per_batch.append(species_counts)
        global_species_per_batch.append(global_species_counts)

    total_particles = sum(batch_totals)
    global_total_particles = sum(global_batch_totals)

    return {
        "secondary_emission_excluded": any(
            spec.get("secondary_emission_enabled", False) and spec.get("enabled", True)
            for spec in species_list
        ),
        "batch_count": batch_count,
        "target_batch_count": target_batch_count,
        "completed_batches": completed_batches,
        "threads": threads,
        "mpi_ranks": mpi_ranks,
        "mpi_rank": mpi_rank,
        "particle_distribution": particle_distribution,
        "rank_count_basis": (
            "shared_even_load_assumption"
            if particle_distribution == "shared"
            else "static_contiguous_slots"
        ),
        "batch_totals": batch_totals,
        "batch_thread_min": batch_thread_min,
        "batch_thread_max": batch_thread_max,
        "species_per_batch": species_per_batch,
        "global_batch_totals": global_batch_totals,
        "global_species_per_batch": global_species_per_batch,
        "final_residuals": residuals,
        "total_particles": total_particles,
        "global_total_particles": global_total_particles,
        "local_reservoir_particles": local_reservoir_particles,
        "global_reservoir_particles": global_reservoir_particles,
        "resolved_batch_duration": resolved_batch_duration,
    }
