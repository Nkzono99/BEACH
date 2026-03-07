"""CLI for estimating particle workload from a Fortran TOML config."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from typing import Any, Sequence

K_BOLTZMANN = 1.380649e-23
EV_TO_K = 1.160451812e4

DEFAULT_SIM: dict[str, Any] = {
    "batch_count": 1,
    "batch_duration": 0.0,
    "batch_duration_step": 0.0,
    "use_box": False,
    "box_min": [-1.0, -1.0, -1.0],
    "box_max": [1.0, 1.0, 1.0],
}

DEFAULT_SPECIES: dict[str, Any] = {
    "enabled": True,
    "source_mode": "volume_seed",
    "npcls_per_step": 0,
    "temperature_k": 2.0e4,
    "m_particle": 9.10938356e-31,
    "w_particle": 1.0,
    "target_macro_particles_per_batch": 0,
    "pos_low": [-0.4, -0.4, 0.2],
    "pos_high": [0.4, 0.4, 0.5],
    "drift_velocity": [0.0, 0.0, -8.0e5],
    "emit_current_density_a_m2": 0.0,
    "rays_per_batch": 0,
    "deposit_opposite_charge_on_emit": False,
    "normal_drift_speed": 0.0,
    "ray_direction": [0.0, 0.0, 0.0],
    "inject_face": "",
}


def _split_count_for_rank(total_count: int, rank: int, n_ranks: int) -> int:
    if total_count < 0:
        raise SystemExit("split_count_for_rank requires total_count >= 0.")
    if n_ranks <= 0:
        raise SystemExit("split_count_for_rank requires n_ranks > 0.")
    if rank < 0 or rank >= n_ranks:
        raise SystemExit("split_count_for_rank rank out of range.")
    base, rem = divmod(total_count, n_ranks)
    return base + (1 if rank < rem else 0)


def load_toml(path: Path) -> dict[str, Any]:
    try:
        import tomllib  # py311+

        with path.open("rb") as f:
            return tomllib.load(f)
    except ModuleNotFoundError:
        try:
            import tomli  # type: ignore

            with path.open("rb") as f:
                return tomli.load(f)
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "TOML parser is missing. Use Python 3.11+ or install tomli: `python -m pip install tomli`."
            ) from exc


def _dot(a: list[float], b: list[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _standard_normal_pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def _standard_normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _inward_normal_for_face(inject_face: str) -> list[float]:
    if inject_face == "x_low":
        return [1.0, 0.0, 0.0]
    if inject_face == "x_high":
        return [-1.0, 0.0, 0.0]
    if inject_face == "y_low":
        return [0.0, 1.0, 0.0]
    if inject_face == "y_high":
        return [0.0, -1.0, 0.0]
    if inject_face == "z_low":
        return [0.0, 0.0, 1.0]
    if inject_face == "z_high":
        return [0.0, 0.0, -1.0]
    raise ValueError(f"unknown inject_face: {inject_face}")


def _face_tangential_axes(inject_face: str) -> tuple[int, int]:
    if inject_face in ("x_low", "x_high"):
        return 1, 2
    if inject_face in ("y_low", "y_high"):
        return 2, 0
    if inject_face in ("z_low", "z_high"):
        return 0, 1
    raise ValueError(f"unknown inject_face: {inject_face}")


def _axis_and_boundary_for_face(
    inject_face: str, box_min: list[float], box_max: list[float]
) -> tuple[int, float]:
    if inject_face == "x_low":
        return 0, box_min[0]
    if inject_face == "x_high":
        return 0, box_max[0]
    if inject_face == "y_low":
        return 1, box_min[1]
    if inject_face == "y_high":
        return 1, box_max[1]
    if inject_face == "z_low":
        return 2, box_min[2]
    if inject_face == "z_high":
        return 2, box_max[2]
    raise ValueError(f"unknown inject_face: {inject_face}")


def _species_temperature_k(spec: dict[str, Any]) -> float:
    if bool(spec.get("_has_temperature_ev", False)):
        return float(spec["temperature_ev"]) * EV_TO_K
    return float(spec.get("temperature_k", DEFAULT_SPECIES["temperature_k"]))


def _species_density_m3(spec: dict[str, Any]) -> float:
    if bool(spec.get("_has_number_density_cm3", False)):
        return float(spec["number_density_cm3"]) * 1.0e6
    return float(spec["number_density_m3"])


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


def _parse_reservoir_species_geometry(
    sim_cfg: dict[str, Any],
    spec: dict[str, Any],
) -> tuple[str, list[float], list[float], float]:
    inject_face = str(spec.get("inject_face", "")).strip().lower()
    _inward_normal_for_face(inject_face)
    pos_low = [float(x) for x in spec.get("pos_low", DEFAULT_SPECIES["pos_low"])]
    pos_high = [float(x) for x in spec.get("pos_high", DEFAULT_SPECIES["pos_high"])]
    box_min = [float(x) for x in sim_cfg.get("box_min", DEFAULT_SIM["box_min"])]
    box_max = [float(x) for x in sim_cfg.get("box_max", DEFAULT_SIM["box_max"])]
    axis_n, boundary_value = _axis_and_boundary_for_face(inject_face, box_min, box_max)
    if (
        abs(pos_low[axis_n] - boundary_value) > 1.0e-12
        or abs(pos_high[axis_n] - boundary_value) > 1.0e-12
    ):
        raise SystemExit(
            "reservoir_face pos_low/pos_high must lie on the selected box face."
        )
    axis_t1, axis_t2 = _face_tangential_axes(inject_face)
    if pos_high[axis_t1] < pos_low[axis_t1] or pos_high[axis_t2] < pos_low[axis_t2]:
        raise SystemExit(
            "reservoir_face tangential bounds must satisfy pos_high >= pos_low."
        )
    area = _compute_face_area(inject_face, pos_low, pos_high)
    if area <= 0.0:
        raise SystemExit("reservoir_face opening area must be positive")
    return inject_face, pos_low, pos_high, area


def _validate_reservoir_species(
    sim_cfg: dict[str, Any],
    spec: dict[str, Any],
    resolved_batch_duration: float,
    species_idx: int,
    species_list: list[dict[str, Any]],
    reservoir_params_by_species: list[dict[str, Any] | None],
) -> dict[str, Any]:
    if not bool(sim_cfg.get("use_box", False)):
        raise SystemExit("reservoir_face requires sim.use_box = true")
    if abs(float(spec.get("emit_current_density_a_m2", 0.0))) > 0.0 or int(
        spec.get("rays_per_batch", 0)
    ) != 0 or bool(spec.get("_has_ray_direction", False)) or bool(
        spec.get("_has_deposit_opposite_charge_on_emit", False)
    ):
        raise SystemExit("photo_raycast keys are not allowed for reservoir_face.")

    has_w_particle = bool(spec.get("_has_w_particle", False))
    has_target_macro = bool(spec.get("_has_target_macro_particles_per_batch", False))
    if has_w_particle and has_target_macro:
        raise SystemExit(
            "reservoir_face does not allow both w_particle and target_macro_particles_per_batch."
        )
    if (not has_w_particle) and (not has_target_macro):
        raise SystemExit(
            "reservoir_face requires either w_particle or target_macro_particles_per_batch."
        )
    if has_w_particle:
        w_particle = float(spec.get("w_particle", DEFAULT_SPECIES["w_particle"]))
        if w_particle <= 0.0:
            raise SystemExit(
                "particles.species.w_particle must be > 0 for reservoir_face"
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
            spec1 = species_list[0]
            if not bool(spec1.get("enabled", True)):
                raise SystemExit(
                    "target_macro_particles_per_batch=-1 requires particles.species[1] to be enabled."
                )
            mode1 = str(spec1.get("source_mode", "volume_seed")).strip().lower()
            if mode1 != "reservoir_face":
                raise SystemExit(
                    'target_macro_particles_per_batch=-1 requires particles.species[1].source_mode="reservoir_face".'
                )
            params1 = reservoir_params_by_species[0]
            if params1 is None or float(params1["w_particle"]) <= 0.0:
                raise SystemExit(
                    "target_macro_particles_per_batch=-1 requires species[1] to resolve a positive w_particle."
                )
            w_particle = float(params1["w_particle"])

    has_cm3 = bool(spec.get("_has_number_density_cm3", False))
    has_m3 = bool(spec.get("_has_number_density_m3", False))
    if has_cm3 and has_m3:
        raise SystemExit(
            "Specify either number_density_cm3 or number_density_m3, not both."
        )
    if (not has_cm3) and (not has_m3):
        raise SystemExit(
            "reservoir_face requires number_density_cm3 or number_density_m3."
        )
    if has_cm3 and float(spec["number_density_cm3"]) <= 0.0:
        raise SystemExit("number_density_cm3 must be > 0.")
    if has_m3 and float(spec["number_density_m3"]) <= 0.0:
        raise SystemExit("number_density_m3 must be > 0.")

    has_temp_k = bool(spec.get("_has_temperature_k", False))
    has_temp_ev = bool(spec.get("_has_temperature_ev", False))
    if has_temp_k and has_temp_ev:
        raise SystemExit("Specify either temperature_ev or temperature_k, not both.")
    if has_temp_ev and float(spec["temperature_ev"]) < 0.0:
        raise SystemExit("temperature_ev must be >= 0.")
    if has_temp_k and float(spec["temperature_k"]) < 0.0:
        raise SystemExit("temperature_k must be >= 0.")

    m_particle = float(spec.get("m_particle", DEFAULT_SPECIES["m_particle"]))
    if m_particle <= 0.0:
        raise SystemExit("m_particle must be > 0.")

    drift_velocity = [
        float(x) for x in spec.get("drift_velocity", DEFAULT_SPECIES["drift_velocity"])
    ]
    temperature_k = _species_temperature_k(spec)
    number_density_m3 = _species_density_m3(spec)
    inject_face, pos_low, pos_high, area = _parse_reservoir_species_geometry(
        sim_cfg, spec
    )
    inward_normal = _inward_normal_for_face(inject_face)
    gamma_in = _compute_inflow_flux(
        number_density_m3=number_density_m3,
        temperature_k=temperature_k,
        m_particle=m_particle,
        drift_velocity=drift_velocity,
        inward_normal=inward_normal,
    )
    if has_target_macro and target_macro != -1:
        w_particle = (
            gamma_in * area * resolved_batch_duration / float(target_macro)
        )
        if (not math.isfinite(w_particle)) or w_particle <= 0.0:
            raise SystemExit(
                "target_macro_particles_per_batch produced invalid w_particle."
            )

    return {
        "inject_face": inject_face,
        "pos_low": pos_low,
        "pos_high": pos_high,
        "area": area,
        "w_particle": w_particle,
        "drift_velocity": drift_velocity,
        "m_particle": m_particle,
        "temperature_k": temperature_k,
        "number_density_m3": number_density_m3,
        "gamma_in": gamma_in,
    }


def _validate_photo_raycast_species(
    sim_cfg: dict[str, Any],
    spec: dict[str, Any],
) -> dict[str, Any]:
    if not bool(sim_cfg.get("use_box", False)):
        raise SystemExit('particles.species.source_mode="photo_raycast" requires sim.use_box = true.')

    if abs(float(spec.get("emit_current_density_a_m2", 0.0))) <= 0.0:
        raise SystemExit("photo_raycast requires emit_current_density_a_m2 > 0.")
    rays_per_batch = int(spec.get("rays_per_batch", 0))
    if rays_per_batch <= 0:
        raise SystemExit("photo_raycast requires rays_per_batch > 0.")
    if int(spec.get("npcls_per_step", 0)) != 0:
        raise SystemExit("npcls_per_step is not allowed for photo_raycast.")

    if bool(spec.get("_has_w_particle", False)):
        raise SystemExit("w_particle is not allowed for photo_raycast.")
    if bool(spec.get("_has_target_macro_particles_per_batch", False)):
        raise SystemExit("target_macro_particles_per_batch is not allowed for photo_raycast.")
    if bool(spec.get("_has_number_density_cm3", False)) or bool(
        spec.get("_has_number_density_m3", False)
    ):
        raise SystemExit("number_density_cm3/number_density_m3 are not allowed for photo_raycast.")

    has_temp_k = bool(spec.get("_has_temperature_k", False))
    has_temp_ev = bool(spec.get("_has_temperature_ev", False))
    if has_temp_k and has_temp_ev:
        raise SystemExit("Specify either temperature_ev or temperature_k, not both.")
    if has_temp_ev and float(spec["temperature_ev"]) < 0.0:
        raise SystemExit("temperature_ev must be >= 0.")
    if has_temp_k and float(spec["temperature_k"]) < 0.0:
        raise SystemExit("temperature_k must be >= 0.")

    m_particle = float(spec.get("m_particle", DEFAULT_SPECIES["m_particle"]))
    if m_particle <= 0.0:
        raise SystemExit("m_particle must be > 0.")
    if abs(float(spec.get("q_particle", 0.0))) <= 0.0:
        raise SystemExit("q_particle must be non-zero for photo_raycast.")

    inject_face = str(spec.get("inject_face", "")).strip().lower()
    inward_normal = _inward_normal_for_face(inject_face)
    pos_low = [float(x) for x in spec.get("pos_low", DEFAULT_SPECIES["pos_low"])]
    pos_high = [float(x) for x in spec.get("pos_high", DEFAULT_SPECIES["pos_high"])]
    box_min = [float(x) for x in sim_cfg.get("box_min", DEFAULT_SIM["box_min"])]
    box_max = [float(x) for x in sim_cfg.get("box_max", DEFAULT_SIM["box_max"])]
    axis_n, boundary_value = _axis_and_boundary_for_face(inject_face, box_min, box_max)
    if (
        abs(pos_low[axis_n] - boundary_value) > 1.0e-12
        or abs(pos_high[axis_n] - boundary_value) > 1.0e-12
    ):
        raise SystemExit(
            "photo_raycast pos_low/pos_high must lie on the selected box face."
        )
    axis_t1, axis_t2 = _face_tangential_axes(inject_face)
    if pos_high[axis_t1] < pos_low[axis_t1] or pos_high[axis_t2] < pos_low[axis_t2]:
        raise SystemExit(
            "photo_raycast tangential bounds must satisfy pos_high >= pos_low."
        )
    area = _compute_face_area(inject_face, pos_low, pos_high)
    if area <= 0.0:
        raise SystemExit("photo_raycast opening area must be positive")

    has_ray_direction = bool(spec.get("_has_ray_direction", False))
    if has_ray_direction:
        ray_direction = [float(x) for x in spec["ray_direction"]]
        norm_ray = math.sqrt(_dot(ray_direction, ray_direction))
        if (not math.isfinite(norm_ray)) or norm_ray <= 0.0:
            raise SystemExit("ray_direction norm must be > 0.")
        ray_direction = [x / norm_ray for x in ray_direction]
    else:
        ray_direction = list(inward_normal)

    inward_dot = _dot(ray_direction, inward_normal)
    if inward_dot <= 0.0:
        raise SystemExit("ray_direction must point inward from inject_face.")

    return {
        "rays_per_batch": rays_per_batch,
        "emit_current_density_a_m2": float(spec.get("emit_current_density_a_m2", 0.0)),
        "inject_face": inject_face,
        "ray_direction": ray_direction,
    }


def _resolve_batch_duration(
    sim_cfg: dict[str, Any],
    sim_raw: dict[str, Any],
    has_dynamic_source_species: bool,
) -> float:
    has_batch_duration = "batch_duration" in sim_raw
    has_batch_duration_step = "batch_duration_step" in sim_raw
    if "target_npcls_species1" in sim_raw:
        raise SystemExit(
            "sim.target_npcls_species1 was removed. Use [[particles.species]].target_macro_particles_per_batch."
        )
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
    residuals = [0.0] * n_species
    if path is None:
        return residuals
    if not path.exists():
        raise SystemExit(f"macro residual file not found: {path}")

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            idx = int(row["species_idx"])
            if 1 <= idx <= n_species:
                residuals[idx - 1] = float(row["residual"])
    return residuals


def estimate_workload(
    config: dict[str, Any],
    threads: int,
    initial_residuals: list[float] | None = None,
    mpi_ranks: int = 1,
    mpi_rank: int = 0,
) -> dict[str, Any]:
    sim_raw = config.get("sim", {})
    if not isinstance(sim_raw, dict):
        raise SystemExit("[sim] section must be a table.")
    sim_cfg = dict(DEFAULT_SIM)
    sim_cfg.update(sim_raw)

    species_list_raw = config.get("particles", {}).get("species", [])
    if not isinstance(species_list_raw, list) or len(species_list_raw) == 0:
        raise SystemExit("At least one [[particles.species]] entry is required.")

    species_list: list[dict[str, Any]] = []
    for raw in species_list_raw:
        if not isinstance(raw, dict):
            raise SystemExit("Each [[particles.species]] entry must be a table.")
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
        spec["_has_ray_direction"] = "ray_direction" in raw
        spec["_has_deposit_opposite_charge_on_emit"] = (
            "deposit_opposite_charge_on_emit" in raw
        )
        species_list.append(spec)

    batch_count = int(sim_cfg["batch_count"])
    if batch_count <= 0:
        raise SystemExit("sim.batch_count must be > 0")

    if threads <= 0:
        raise SystemExit("threads must be > 0")
    if mpi_ranks <= 0:
        raise SystemExit("mpi_ranks must be > 0")
    if mpi_rank < 0 or mpi_rank >= mpi_ranks:
        raise SystemExit("mpi_rank must satisfy 0 <= mpi_rank < mpi_ranks")

    has_dynamic_source_species = False
    per_batch_volume_particles = 0
    for spec in species_list:
        if not bool(spec.get("enabled", True)):
            continue
        source_mode = str(spec.get("source_mode", "volume_seed")).strip().lower()
        if source_mode == "volume_seed":
            if bool(spec.get("_has_target_macro_particles_per_batch", False)):
                raise SystemExit(
                    "target_macro_particles_per_batch is only valid for reservoir_face."
                )
            if abs(float(spec.get("emit_current_density_a_m2", 0.0))) > 0.0 or int(
                spec.get("rays_per_batch", 0)
            ) != 0 or bool(spec.get("_has_ray_direction", False)) or bool(
                spec.get("_has_deposit_opposite_charge_on_emit", False)
            ):
                raise SystemExit(
                    'photo_raycast keys are only valid for source_mode="photo_raycast".'
                )
            n_macro = int(spec.get("npcls_per_step", 0))
            if n_macro < 0:
                raise SystemExit("particles.species.npcls_per_step must be >= 0")
            per_batch_volume_particles += n_macro
        elif source_mode == "reservoir_face":
            has_dynamic_source_species = True
        elif source_mode == "photo_raycast":
            has_dynamic_source_species = True
        else:
            raise SystemExit(f"Unknown particles.species.source_mode: {source_mode}")
    if per_batch_volume_particles <= 0 and not has_dynamic_source_species:
        raise SystemExit(
            "At least one enabled [[particles.species]] entry must have npcls_per_step > 0."
        )
    resolved_batch_duration = _resolve_batch_duration(
        sim_cfg=sim_cfg,
        sim_raw=sim_raw,
        has_dynamic_source_species=has_dynamic_source_species,
    )
    reservoir_params_by_species: list[dict[str, Any] | None] = [None] * len(
        species_list
    )
    photo_params_by_species: list[dict[str, Any] | None] = [None] * len(species_list)
    for s_idx, spec in enumerate(species_list):
        if not bool(spec.get("enabled", True)):
            continue
        source_mode = str(spec.get("source_mode", "volume_seed")).strip().lower()
        if source_mode == "reservoir_face":
            reservoir_params_by_species[s_idx] = _validate_reservoir_species(
                sim_cfg,
                spec,
                resolved_batch_duration,
                s_idx,
                species_list,
                reservoir_params_by_species,
            )
        elif source_mode == "photo_raycast":
            photo_params_by_species[s_idx] = _validate_photo_raycast_species(
                sim_cfg,
                spec,
            )

    residuals = [0.0] * len(species_list)
    if initial_residuals is not None:
        if len(initial_residuals) != len(species_list):
            raise SystemExit(
                "macro residual species count does not match config species count"
            )
        residuals = list(initial_residuals)

    batch_totals: list[int] = []
    batch_thread_min: list[int] = []
    batch_thread_max: list[int] = []
    species_per_batch: list[list[int]] = []

    for _batch_idx in range(batch_count):
        species_counts: list[int] = []
        for s_idx, spec in enumerate(species_list):
            if not bool(spec.get("enabled", True)):
                species_counts.append(0)
                continue

            source_mode = str(spec.get("source_mode", "volume_seed")).strip().lower()
            if source_mode == "volume_seed":
                n_macro_global = int(spec.get("npcls_per_step", 0))
                if n_macro_global < 0:
                    raise SystemExit("particles.species.npcls_per_step must be >= 0")
                species_counts.append(
                    _split_count_for_rank(n_macro_global, mpi_rank, mpi_ranks)
                )
                continue

            if source_mode == "reservoir_face":
                params = reservoir_params_by_species[s_idx]
                if params is None:
                    raise SystemExit(
                        "internal error: reservoir species parameters were not initialized."
                    )
                n_phys_batch = (
                    params["gamma_in"]
                    * params["area"]
                    * resolved_batch_duration
                    / float(mpi_ranks)
                )
                n_macro_expected = n_phys_batch / params["w_particle"]
                macro_budget = residuals[s_idx] + n_macro_expected
                if macro_budget < 0.0:
                    macro_budget = 0.0
                n_macro = math.floor(macro_budget)
                residuals[s_idx] = macro_budget - n_macro
                species_counts.append(int(n_macro))
                continue

            if source_mode == "photo_raycast":
                params = photo_params_by_species[s_idx]
                if params is None:
                    raise SystemExit(
                        "internal error: photo_raycast species parameters were not initialized."
                    )
                species_counts.append(
                    _split_count_for_rank(
                        int(params["rays_per_batch"]), mpi_rank, mpi_ranks
                    )
                )
                continue

            raise SystemExit(f"Unknown particles.species.source_mode: {source_mode}")

        batch_total = sum(species_counts)
        q, r = divmod(batch_total, threads)
        batch_totals.append(batch_total)
        batch_thread_min.append(q)
        batch_thread_max.append(q + 1 if r > 0 else q)
        species_per_batch.append(species_counts)

    total_particles = sum(batch_totals)

    return {
        "batch_count": batch_count,
        "threads": threads,
        "mpi_ranks": mpi_ranks,
        "mpi_rank": mpi_rank,
        "batch_totals": batch_totals,
        "batch_thread_min": batch_thread_min,
        "batch_thread_max": batch_thread_max,
        "species_per_batch": species_per_batch,
        "final_residuals": residuals,
        "total_particles": total_particles,
        "resolved_batch_duration": resolved_batch_duration,
    }


def _default_threads() -> int:
    value = os.environ.get("OMP_NUM_THREADS", "").strip()
    if value:
        try:
            parsed = int(value)
            if parsed > 0:
                return parsed
        except ValueError:
            pass
    return 1


def _default_mpi_ranks() -> int:
    for key in ("PMI_SIZE", "OMPI_COMM_WORLD_SIZE", "SLURM_NTASKS"):
        value = os.environ.get(key, "").strip()
        if not value:
            continue
        try:
            parsed = int(value)
            if parsed > 0:
                return parsed
        except ValueError:
            pass
    return 1


def _default_mpi_rank() -> int:
    for key in ("PMI_RANK", "OMPI_COMM_WORLD_RANK", "SLURM_PROCID"):
        value = os.environ.get(key, "").strip()
        if not value:
            continue
        try:
            parsed = int(value)
            if parsed >= 0:
                return parsed
        except ValueError:
            pass
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate particle workload from Fortran TOML config: "
            "per-batch local(rank) particles and per-thread particle counts."
        )
    )
    parser.add_argument("config", type=Path, help="path to beach.toml")
    parser.add_argument(
        "--threads",
        type=int,
        default=_default_threads(),
        help="OpenMP thread count for per-thread estimate (default: OMP_NUM_THREADS or 1)",
    )
    parser.add_argument(
        "--macro-residuals",
        type=Path,
        default=None,
        help="optional macro_residuals.csv to start from resume state",
    )
    parser.add_argument(
        "--mpi-ranks",
        type=int,
        default=_default_mpi_ranks(),
        help="MPI world size used for local(rank) workload estimate (default: MPI env or 1)",
    )
    parser.add_argument(
        "--mpi-rank",
        type=int,
        default=_default_mpi_rank(),
        help="rank index used for local(rank) workload estimate (default: MPI env or 0)",
    )
    parser.add_argument(
        "--show-batches",
        type=int,
        default=10,
        help="number of head batches to print in detail (default: 10)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.config.exists():
        raise SystemExit(f"config file not found: {args.config}")

    config = load_toml(args.config)
    species_raw = config.get("particles", {}).get("species", [])
    initial_residuals = read_macro_residuals(args.macro_residuals, len(species_raw))

    result = estimate_workload(
        config=config,
        threads=args.threads,
        initial_residuals=initial_residuals,
        mpi_ranks=args.mpi_ranks,
        mpi_rank=args.mpi_rank,
    )

    batch_totals = result["batch_totals"]
    batch_thread_min = result["batch_thread_min"]
    batch_thread_max = result["batch_thread_max"]
    total_particles = result["total_particles"]
    batch_count = result["batch_count"]
    threads = result["threads"]
    mpi_ranks = result["mpi_ranks"]
    mpi_rank = result["mpi_rank"]

    print(f"config={args.config}")
    print(f"threads={threads}")
    print(f"mpi_ranks={mpi_ranks}")
    print(f"mpi_rank={mpi_rank}")
    print("estimate_scope=local_rank")
    print(f"batch_count={batch_count}")
    print(f"resolved_batch_duration={result['resolved_batch_duration']}")
    print(f"total_particles={total_particles}")
    print(f"particles_per_batch_min={min(batch_totals)}")
    print(f"particles_per_batch_max={max(batch_totals)}")
    print(f"particles_per_batch_avg={total_particles / batch_count:.3f}")
    print(f"per_thread_particles_min={min(batch_thread_min)}")
    print(f"per_thread_particles_max={max(batch_thread_max)}")
    print(
        "per_thread_particles_avg=" f"{total_particles / (batch_count * threads):.3f}"
    )
    print(f"final_macro_residuals={result['final_residuals']}")

    n_detail = max(0, min(args.show_batches, batch_count))
    if n_detail > 0:
        print("batch_details=")
        for batch_idx in range(n_detail):
            species_counts = result["species_per_batch"][batch_idx]
            print(
                "  "
                f"batch={batch_idx + 1} "
                f"total={batch_totals[batch_idx]} "
                f"per_thread=[{batch_thread_min[batch_idx]},{batch_thread_max[batch_idx]}] "
                f"species={species_counts}"
            )


if __name__ == "__main__":
    main()
