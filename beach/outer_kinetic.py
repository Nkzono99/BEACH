"""Offline 1D1V Vlasov oracle for matching-plane response studies.

The implementation is intentionally independent of the BEACH runtime. It uses
positive, conservative first-order finite-volume advection and rejects results
that do not pass stationarity, far-boundary, and conservation checks.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from beach.config._toml import load_toml_file


EPS0 = 8.8541878128e-12
ELEMENTARY_CHARGE = 1.602176634e-19
ELECTRON_MASS = 9.1093837139e-31
PROTON_MASS = 1.67262192369e-27
SOLVER_SOURCE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

QUERY_HEADER = (
    "displacement_c_m2",
    "photoelectron_outward_number_flux_m2_s",
    "photoelectron_outward_mean_normal_energy_ev",
    "electron_outward_number_flux_m2_s",
    "ion_outward_number_flux_m2_s",
)
RESPONSE_HEADER = (
    *QUERY_HEADER,
    "matching_potential_v",
    "electron_inward_number_flux_m2_s",
    "ion_inward_number_flux_m2_s",
    "electron_access_potential_v",
    "ion_access_potential_v",
    "photoelectron_barrier_potential_v",
)
RAW_HEADER = (
    *RESPONSE_HEADER,
    "photoelectron_return_number_flux_m2_s",
    "photoelectron_escape_number_flux_m2_s",
    "classification",
    "failure_reason",
    "simulated_time_s",
    "averaging_window_s",
    "electron_inventory_m2",
    "ion_inventory_m2",
    "photoelectron_inventory_m2",
    "far_field_v_m",
    "far_charge_imbalance",
    "gauss_residual_c_m2",
    "charge_budget_residual",
    "mean_stationarity_metric",
    "max_drift_metric",
    "max_standard_error_metric",
    "max_velocity_boundary_loss_fraction",
    "time_steps",
)
ACCEPTED_CLASSIFICATIONS = frozenset({"steady", "stationary_average"})


class KineticConfigError(ValueError):
    """Raised when an outer-kinetic configuration is invalid."""


@dataclass(frozen=True)
class VelocityGridConfig:
    nv: int
    vmin_mps: float
    vmax_mps: float

    def validate(self, name: str) -> None:
        if self.nv < 8:
            raise KineticConfigError(f"{name}.nv must be >= 8")
        if not (
            math.isfinite(self.vmin_mps)
            and math.isfinite(self.vmax_mps)
            and self.vmin_mps < 0.0 < self.vmax_mps
        ):
            raise KineticConfigError(
                f"{name} velocity range must be finite and straddle zero"
            )


@dataclass(frozen=True)
class ReservoirConfig:
    number_density_m3: float
    temperature_ev: float
    drift_velocity_mps: float
    mass_kg: float
    charge_c: float
    grid: VelocityGridConfig

    def validate(self, name: str) -> None:
        values = (
            self.number_density_m3,
            self.temperature_ev,
            self.drift_velocity_mps,
            self.mass_kg,
            self.charge_c,
        )
        if not all(math.isfinite(value) for value in values):
            raise KineticConfigError(f"{name} values must be finite")
        if self.number_density_m3 <= 0.0 or self.temperature_ev <= 0.0:
            raise KineticConfigError(f"{name} density and temperature must be > 0")
        if self.mass_kg <= 0.0 or self.charge_c == 0.0:
            raise KineticConfigError(f"{name} mass must be > 0 and charge nonzero")
        self.grid.validate(name)


@dataclass(frozen=True)
class CertificationConfig:
    warmup_time_s: float
    averaging_window_s: float
    sample_interval_s: float
    stationarity_rtol: float
    drift_rtol: float
    sem_rtol: float
    steady_fluctuation_rtol: float
    autocorrelation_windows: float
    far_field_abs_v_m: float
    far_charge_rel: float
    gauss_abs_c_m2: float
    charge_budget_rtol: float
    velocity_loss_rtol: float

    def validate(self, max_time_s: float) -> None:
        positive = (
            self.warmup_time_s,
            self.averaging_window_s,
            self.sample_interval_s,
            self.stationarity_rtol,
            self.drift_rtol,
            self.sem_rtol,
            self.steady_fluctuation_rtol,
            self.autocorrelation_windows,
            self.far_field_abs_v_m,
            self.far_charge_rel,
            self.gauss_abs_c_m2,
            self.charge_budget_rtol,
            self.velocity_loss_rtol,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in positive):
            raise KineticConfigError("certification values must be finite and > 0")
        required_time = self.warmup_time_s + 2.0 * self.averaging_window_s
        time_tolerance = 16.0 * np.finfo(float).eps * max(max_time_s, 1.0e-300)
        if required_time > max_time_s + time_tolerance:
            raise KineticConfigError(
                "max_time_s must cover warmup plus two averaging windows"
            )
        if self.sample_interval_s > self.averaging_window_s / 4.0:
            raise KineticConfigError(
                "sample_interval_s must provide at least four samples per window"
            )


@dataclass(frozen=True)
class OuterKineticConfig:
    matching_plane_z_m: float
    z_length_m: float
    nz: int
    cfl: float
    max_time_s: float
    initial_condition: str
    electron: ReservoirConfig
    ion: ReservoirConfig
    photoelectron_grid: VelocityGridConfig
    certification: CertificationConfig

    def validate(self) -> None:
        if not math.isfinite(self.matching_plane_z_m):
            raise KineticConfigError("matching_plane_z_m must be finite")
        if not math.isfinite(self.z_length_m) or self.z_length_m <= 0.0:
            raise KineticConfigError("z_length_m must be finite and > 0")
        if self.nz < 8:
            raise KineticConfigError("nz must be >= 8")
        if not math.isfinite(self.cfl) or not 0.0 < self.cfl <= 1.0:
            raise KineticConfigError("cfl must be in (0, 1]")
        if not math.isfinite(self.max_time_s) or self.max_time_s <= 0.0:
            raise KineticConfigError("max_time_s must be finite and > 0")
        if self.initial_condition not in {"reservoir", "empty"}:
            raise KineticConfigError(
                'initial_condition must be "reservoir" or "empty"'
            )
        self.electron.validate("ambient_electron")
        self.ion.validate("ambient_ion")
        if self.electron.charge_c >= 0.0:
            raise KineticConfigError("ambient_electron charge must be negative")
        if self.ion.charge_c <= 0.0:
            raise KineticConfigError("ambient_ion charge must be positive")
        self.photoelectron_grid.validate("photoelectron")
        self.certification.validate(self.max_time_s)


@dataclass(frozen=True)
class KineticQuery:
    displacement_c_m2: float
    photoelectron_outward_number_flux_m2_s: float
    photoelectron_outward_mean_normal_energy_ev: float
    electron_outward_number_flux_m2_s: float = 0.0
    ion_outward_number_flux_m2_s: float = 0.0

    def validate(self) -> None:
        values = asdict(self)
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError("kinetic query values must be finite")
        for name, value in values.items():
            if name != "displacement_c_m2" and value < 0.0:
                raise ValueError(f"{name} must be >= 0")
        if (
            self.electron_outward_number_flux_m2_s != 0.0
            or self.ion_outward_number_flux_m2_s != 0.0
        ):
            raise ValueError(
                "kinetic oracle v1 supports only zero ambient outward fluxes"
            )
        if (
            self.photoelectron_outward_number_flux_m2_s > 0.0
            and self.photoelectron_outward_mean_normal_energy_ev <= 0.0
        ):
            raise ValueError("positive photoelectron flux requires positive energy")


@dataclass
class FieldProfile:
    electric_faces_v_m: np.ndarray
    electric_cells_v_m: np.ndarray
    potential_faces_v: np.ndarray
    potential_cells_v: np.ndarray
    rho_c_m3: np.ndarray
    gauss_residual_c_m2: float


@dataclass
class KineticResult:
    query: KineticQuery
    response: tuple[float, float, float, float, float, float]
    photoelectron_return_flux_m2_s: float
    photoelectron_escape_flux_m2_s: float
    classification: str
    failure_reason: str
    simulated_time_s: float
    averaging_window_s: float
    inventories_m2: tuple[float, float, float]
    far_field_v_m: float
    far_charge_imbalance: float
    gauss_residual_c_m2: float
    charge_budget_residual: float
    mean_stationarity_metric: float
    max_drift_metric: float
    max_standard_error_metric: float
    max_velocity_boundary_loss_fraction: float
    time_steps: int
    profile: FieldProfile
    time_history: dict[str, np.ndarray]

    def raw_row(self) -> dict[str, object]:
        values: list[object] = [
            *asdict(self.query).values(),
            *self.response,
            self.photoelectron_return_flux_m2_s,
            self.photoelectron_escape_flux_m2_s,
            self.classification,
            self.failure_reason,
            self.simulated_time_s,
            self.averaging_window_s,
            *self.inventories_m2,
            self.far_field_v_m,
            self.far_charge_imbalance,
            self.gauss_residual_c_m2,
            self.charge_budget_residual,
            self.mean_stationarity_metric,
            self.max_drift_metric,
            self.max_standard_error_metric,
            self.max_velocity_boundary_loss_fraction,
            self.time_steps,
        ]
        return dict(zip(RAW_HEADER, values, strict=True))


@dataclass
class _SpeciesState:
    name: str
    mass_kg: float
    charge_c: float
    velocity_mps: np.ndarray
    dv_mps: float
    distribution: np.ndarray
    left_inflow: np.ndarray
    right_inflow: np.ndarray
    initial_inventory_m2: float
    boundary_inventory_change_m2: float = 0.0
    velocity_inventory_change_m2: float = 0.0


def load_outer_kinetic_config(path: str | Path) -> OuterKineticConfig:
    """Load and validate one standalone outer-kinetic TOML file."""

    document = load_toml_file(path)
    allowed = {
        "outer_kinetic",
        "ambient_electron",
        "ambient_ion",
        "photoelectron",
        "certification",
    }
    unknown = set(document) - allowed
    if unknown:
        raise KineticConfigError(
            "unsupported outer-kinetic section(s): " + ", ".join(sorted(unknown))
        )
    outer = _required_table(document, "outer_kinetic")
    electron = _reservoir_from_table(
        _required_table(document, "ambient_electron"),
        default_mass=ELECTRON_MASS,
        default_charge=-ELEMENTARY_CHARGE,
    )
    ion = _reservoir_from_table(
        _required_table(document, "ambient_ion"),
        default_mass=PROTON_MASS,
        default_charge=ELEMENTARY_CHARGE,
    )
    photo = _required_table(document, "photoelectron")
    cert = _required_table(document, "certification")
    config = OuterKineticConfig(
        matching_plane_z_m=_float(outer, "matching_plane_z_m"),
        z_length_m=_float(outer, "z_length_m"),
        nz=_int(outer, "nz"),
        cfl=_float(outer, "cfl"),
        max_time_s=_float(outer, "max_time_s"),
        initial_condition=str(outer.get("initial_condition", "reservoir")),
        electron=electron,
        ion=ion,
        photoelectron_grid=VelocityGridConfig(
            nv=_int(photo, "nv"),
            vmin_mps=_float(photo, "vmin_mps"),
            vmax_mps=_float(photo, "vmax_mps"),
        ),
        certification=CertificationConfig(
            warmup_time_s=_float(cert, "warmup_time_s"),
            averaging_window_s=_float(cert, "averaging_window_s"),
            sample_interval_s=_float(cert, "sample_interval_s"),
            stationarity_rtol=_float(cert, "stationarity_rtol"),
            drift_rtol=_float(cert, "drift_rtol"),
            sem_rtol=_float(cert, "sem_rtol"),
            steady_fluctuation_rtol=_float(cert, "steady_fluctuation_rtol"),
            autocorrelation_windows=_float(cert, "autocorrelation_windows"),
            far_field_abs_v_m=_float(cert, "far_field_abs_v_m"),
            far_charge_rel=_float(cert, "far_charge_rel"),
            gauss_abs_c_m2=_float(cert, "gauss_abs_c_m2"),
            charge_budget_rtol=_float(cert, "charge_budget_rtol"),
            velocity_loss_rtol=_float(cert, "velocity_loss_rtol"),
        ),
    )
    config.validate()
    return config


def _required_table(document: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = document.get(name)
    if not isinstance(value, Mapping):
        raise KineticConfigError(f"missing [{name}] table")
    return value


def _float(table: Mapping[str, Any], key: str) -> float:
    if key not in table or isinstance(table[key], bool):
        raise KineticConfigError(f"missing numeric key {key}")
    try:
        return float(table[key])
    except (TypeError, ValueError) as exc:
        raise KineticConfigError(f"{key} must be numeric") from exc


def _int(table: Mapping[str, Any], key: str) -> int:
    if key not in table or isinstance(table[key], bool):
        raise KineticConfigError(f"missing integer key {key}")
    value = table[key]
    if not isinstance(value, int):
        raise KineticConfigError(f"{key} must be an integer")
    return value


def _reservoir_from_table(
    table: Mapping[str, Any], *, default_mass: float, default_charge: float
) -> ReservoirConfig:
    return ReservoirConfig(
        number_density_m3=_float(table, "number_density_m3"),
        temperature_ev=_float(table, "temperature_ev"),
        drift_velocity_mps=_float(table, "drift_velocity_mps"),
        mass_kg=float(table.get("mass_kg", default_mass)),
        charge_c=float(table.get("charge_c", default_charge)),
        grid=VelocityGridConfig(
            nv=_int(table, "nv"),
            vmin_mps=_float(table, "vmin_mps"),
            vmax_mps=_float(table, "vmax_mps"),
        ),
    )


def velocity_grid(config: VelocityGridConfig) -> tuple[np.ndarray, float]:
    dv = (config.vmax_mps - config.vmin_mps) / config.nv
    values = config.vmin_mps + (np.arange(config.nv) + 0.5) * dv
    return values, dv


def drifting_maxwellian(
    velocity_mps: np.ndarray,
    *,
    number_density_m3: float,
    temperature_ev: float,
    drift_velocity_mps: float,
    mass_kg: float,
) -> np.ndarray:
    sigma = math.sqrt(temperature_ev * ELEMENTARY_CHARGE / mass_kg)
    return (
        number_density_m3
        / (math.sqrt(2.0 * math.pi) * sigma)
        * np.exp(-0.5 * ((velocity_mps - drift_velocity_mps) / sigma) ** 2)
    )


def photoelectron_half_maxwellian(
    velocity_mps: np.ndarray,
    *,
    number_flux_m2_s: float,
    mean_normal_energy_ev: float,
    mass_kg: float = ELECTRON_MASS,
) -> np.ndarray:
    result = np.zeros_like(velocity_mps)
    if number_flux_m2_s == 0.0:
        return result
    sigma2 = mean_normal_energy_ev * ELEMENTARY_CHARGE / mass_kg
    positive = velocity_mps > 0.0
    result[positive] = (
        number_flux_m2_s
        / sigma2
        * np.exp(-0.5 * velocity_mps[positive] ** 2 / sigma2)
    )
    return result


def compute_field_profile(
    rho_c_m3: np.ndarray,
    *,
    displacement_c_m2: float,
    dz_m: float,
) -> FieldProfile:
    """Integrate Gauss' law from H and set the potential gauge at L."""

    rho = np.asarray(rho_c_m3, dtype=float)
    if rho.ndim != 1 or rho.size == 0:
        raise ValueError("rho_c_m3 must be a non-empty 1-D array")
    if not np.all(np.isfinite(rho)):
        raise ValueError("rho_c_m3 must be finite")
    electric_faces = np.empty(rho.size + 1)
    electric_faces[0] = displacement_c_m2 / EPS0
    electric_faces[1:] = electric_faces[0] + np.cumsum(rho) * dz_m / EPS0
    electric_cells = 0.5 * (electric_faces[:-1] + electric_faces[1:])
    potential_faces = np.zeros(rho.size + 1)
    for index in range(rho.size - 1, -1, -1):
        potential_faces[index] = (
            potential_faces[index + 1] + electric_cells[index] * dz_m
        )
    potential_cells = 0.5 * (potential_faces[:-1] + potential_faces[1:])
    gauss_residual = (
        EPS0 * (electric_faces[-1] - electric_faces[0])
        - float(np.sum(rho) * dz_m)
    )
    return FieldProfile(
        electric_faces_v_m=electric_faces,
        electric_cells_v_m=electric_cells,
        potential_faces_v=potential_faces,
        potential_cells_v=potential_cells,
        rho_c_m3=rho.copy(),
        gauss_residual_c_m2=gauss_residual,
    )


def advect_z_finite_volume(
    distribution: np.ndarray,
    velocity_mps: np.ndarray,
    *,
    dt_s: float,
    dz_m: float,
    dv_mps: float,
    left_inflow: np.ndarray,
    right_inflow: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Advance one z-advection step and return inventory change from z faces."""

    f = np.asarray(distribution, dtype=float)
    velocity = np.asarray(velocity_mps, dtype=float)
    nz, nv = f.shape
    if velocity.shape != (nv,):
        raise ValueError("velocity shape does not match distribution")
    flux = np.empty((nz + 1, nv))
    positive = velocity >= 0.0
    flux[1:nz] = velocity * np.where(positive, f[:-1], f[1:])
    flux[0] = velocity * np.where(positive, left_inflow, f[0])
    flux[nz] = velocity * np.where(positive, f[-1], right_inflow)
    advanced = f - dt_s / dz_m * (flux[1:] - flux[:-1])
    inventory_change = dt_s * float(np.sum(flux[0] - flux[-1]) * dv_mps)
    return advanced, inventory_change


def advect_v_finite_volume(
    distribution: np.ndarray,
    acceleration_m_s2: np.ndarray,
    *,
    dt_s: float,
    dv_mps: float,
    dz_m: float,
) -> tuple[np.ndarray, float]:
    """Advance one v-advection step with open, zero-inflow velocity boundaries."""

    f = np.asarray(distribution, dtype=float)
    acceleration = np.asarray(acceleration_m_s2, dtype=float)
    nz, nv = f.shape
    if acceleration.shape != (nz,):
        raise ValueError("acceleration shape does not match distribution")
    flux = np.empty((nz, nv + 1))
    positive = acceleration >= 0.0
    flux[:, 1:nv] = acceleration[:, None] * np.where(
        positive[:, None], f[:, :-1], f[:, 1:]
    )
    flux[:, 0] = np.where(positive, 0.0, acceleration * f[:, 0])
    flux[:, nv] = np.where(positive, acceleration * f[:, -1], 0.0)
    advanced = f - dt_s / dv_mps * (flux[:, 1:] - flux[:, :-1])
    inventory_change = dt_s * float(np.sum(flux[:, 0] - flux[:, -1]) * dz_m)
    return advanced, inventory_change


def number_flux_at_left(
    distribution: np.ndarray, velocity_mps: np.ndarray, dv_mps: float
) -> float:
    negative = velocity_mps < 0.0
    return float(
        np.sum((-velocity_mps[negative]) * distribution[0, negative]) * dv_mps
    )


def number_flux_at_right(
    distribution: np.ndarray, velocity_mps: np.ndarray, dv_mps: float
) -> float:
    positive = velocity_mps > 0.0
    return float(
        np.sum(velocity_mps[positive] * distribution[-1, positive]) * dv_mps
    )


def read_kinetic_queries(path: str | Path) -> list[KineticQuery]:
    queries: list[KineticQuery] = []
    with Path(path).open(newline="", encoding="utf-8") as stream:
        filtered = (
            line
            for line in stream
            if line.strip() and not line.lstrip().startswith("#")
        )
        reader = csv.DictReader(filtered)
        if tuple(reader.fieldnames or ()) != QUERY_HEADER:
            raise ValueError(
                "kinetic query CSV header must be: " + ",".join(QUERY_HEADER)
            )
        for line_number, row in enumerate(reader, start=2):
            try:
                query = KineticQuery(
                    **{key: float(row[key]) for key in QUERY_HEADER}
                )
                query.validate()
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"invalid kinetic query at CSV line {line_number}: {exc}"
                ) from exc
            queries.append(query)
    if not queries:
        raise ValueError("kinetic query CSV contains no queries")
    return queries


def config_fingerprint(config: OuterKineticConfig) -> str:
    payload = json.dumps(asdict(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


_OBSERVABLE_NAMES = (
    "matching_potential_v",
    "electron_inward_number_flux_m2_s",
    "ion_inward_number_flux_m2_s",
    "photoelectron_return_number_flux_m2_s",
    "photoelectron_escape_number_flux_m2_s",
    "total_outer_charge_c_m2",
)


def _make_reservoir_state(
    name: str,
    config: ReservoirConfig,
    *,
    nz: int,
    dz_m: float,
    initial_condition: str,
) -> _SpeciesState:
    velocity, dv = velocity_grid(config.grid)
    reservoir = drifting_maxwellian(
        velocity,
        number_density_m3=config.number_density_m3,
        temperature_ev=config.temperature_ev,
        drift_velocity_mps=config.drift_velocity_mps,
        mass_kg=config.mass_kg,
    )
    represented_density = float(np.sum(reservoir) * dv)
    if represented_density <= 0.0 or not math.isfinite(represented_density):
        raise KineticConfigError(f"{name} velocity grid resolves no reservoir density")
    reservoir *= config.number_density_m3 / represented_density
    if initial_condition == "reservoir":
        distribution = np.repeat(reservoir[None, :], nz, axis=0)
    else:
        distribution = np.zeros((nz, config.grid.nv))
    left = np.zeros_like(velocity)
    right = reservoir
    inventory = float(np.sum(distribution) * dv * dz_m)
    return _SpeciesState(
        name=name,
        mass_kg=config.mass_kg,
        charge_c=config.charge_c,
        velocity_mps=velocity,
        dv_mps=dv,
        distribution=distribution,
        left_inflow=left,
        right_inflow=right,
        initial_inventory_m2=inventory,
    )


def _make_photoelectron_state(
    config: OuterKineticConfig,
    query: KineticQuery,
    *,
    dz_m: float,
) -> _SpeciesState:
    velocity, dv = velocity_grid(config.photoelectron_grid)
    left = photoelectron_half_maxwellian(
        velocity,
        number_flux_m2_s=query.photoelectron_outward_number_flux_m2_s,
        mean_normal_energy_ev=max(
            query.photoelectron_outward_mean_normal_energy_ev, 1.0
        ),
    )
    if query.photoelectron_outward_number_flux_m2_s > 0.0:
        represented_flux = float(
            np.sum(velocity[velocity > 0.0] * left[velocity > 0.0]) * dv
        )
        if represented_flux <= 0.0 or not math.isfinite(represented_flux):
            raise KineticConfigError(
                "photoelectron velocity grid resolves no outward flux"
            )
        left *= query.photoelectron_outward_number_flux_m2_s / represented_flux
    distribution = np.zeros((config.nz, config.photoelectron_grid.nv))
    return _SpeciesState(
        name="photoelectron",
        mass_kg=ELECTRON_MASS,
        charge_c=-ELEMENTARY_CHARGE,
        velocity_mps=velocity,
        dv_mps=dv,
        distribution=distribution,
        left_inflow=left,
        right_inflow=np.zeros_like(velocity),
        initial_inventory_m2=0.0,
    )


def _inventory(state: _SpeciesState, dz_m: float) -> float:
    return float(np.sum(state.distribution) * state.dv_mps * dz_m)


def _field_from_states(
    states: tuple[_SpeciesState, ...],
    *,
    displacement_c_m2: float,
    dz_m: float,
) -> FieldProfile:
    rho = np.zeros(states[0].distribution.shape[0])
    for state in states:
        density = np.sum(state.distribution, axis=1) * state.dv_mps
        rho += state.charge_c * density
    return compute_field_profile(
        rho, displacement_c_m2=displacement_c_m2, dz_m=dz_m
    )


def _time_step(
    config: OuterKineticConfig,
    states: tuple[_SpeciesState, ...],
    field: FieldProfile,
    *,
    dz_m: float,
) -> float:
    maximum_velocity = max(
        float(np.max(np.abs(state.velocity_mps))) for state in states
    )
    dt_z = dz_m / maximum_velocity
    dt_v = math.inf
    maximum_field = float(np.max(np.abs(field.electric_cells_v_m)))
    if maximum_field > 0.0:
        for state in states:
            acceleration = abs(state.charge_c) * maximum_field / state.mass_kg
            if acceleration > 0.0:
                dt_v = min(dt_v, state.dv_mps / acceleration)
    return config.cfl * min(dt_z, dt_v)


def _advance_species_z(
    state: _SpeciesState, *, dt_s: float, dz_m: float
) -> None:
    state.distribution, change = advect_z_finite_volume(
        state.distribution,
        state.velocity_mps,
        dt_s=dt_s,
        dz_m=dz_m,
        dv_mps=state.dv_mps,
        left_inflow=state.left_inflow,
        right_inflow=state.right_inflow,
    )
    state.boundary_inventory_change_m2 += change


def _advance_species_v(
    state: _SpeciesState,
    field: FieldProfile,
    *,
    dt_s: float,
    dz_m: float,
) -> None:
    acceleration = state.charge_c * field.electric_cells_v_m / state.mass_kg
    state.distribution, change = advect_v_finite_volume(
        state.distribution,
        acceleration,
        dt_s=dt_s,
        dv_mps=state.dv_mps,
        dz_m=dz_m,
    )
    state.velocity_inventory_change_m2 += change


def _enforce_positive(states: tuple[_SpeciesState, ...]) -> None:
    for state in states:
        scale = max(float(np.max(state.distribution)), 1.0)
        minimum = float(np.min(state.distribution))
        if minimum < -1.0e-12 * scale:
            raise FloatingPointError(
                f"{state.name} distribution lost positivity: min={minimum:.6e}"
            )
        np.maximum(state.distribution, 0.0, out=state.distribution)
        if not np.all(np.isfinite(state.distribution)):
            raise FloatingPointError(
                f"{state.name} distribution became non-finite"
            )


def _sample_observables(
    states: tuple[_SpeciesState, _SpeciesState, _SpeciesState],
    field: FieldProfile,
    *,
    dz_m: float,
    reference_density_m3: float,
) -> dict[str, float]:
    electron, ion, photoelectron = states
    potential = field.potential_faces_v
    phi_h = float(potential[0])
    minimum_potential = float(np.min(potential))
    maximum_potential = float(np.max(potential))
    return {
        "matching_potential_v": phi_h,
        "electron_inward_number_flux_m2_s": number_flux_at_left(
            electron.distribution, electron.velocity_mps, electron.dv_mps
        ),
        "ion_inward_number_flux_m2_s": number_flux_at_left(
            ion.distribution, ion.velocity_mps, ion.dv_mps
        ),
        "photoelectron_return_number_flux_m2_s": number_flux_at_left(
            photoelectron.distribution,
            photoelectron.velocity_mps,
            photoelectron.dv_mps,
        ),
        "photoelectron_escape_number_flux_m2_s": number_flux_at_right(
            photoelectron.distribution,
            photoelectron.velocity_mps,
            photoelectron.dv_mps,
        ),
        "electron_access_potential_v": min(0.0, minimum_potential),
        "ion_access_potential_v": max(0.0, maximum_potential),
        "photoelectron_barrier_potential_v": min(
            phi_h, minimum_potential
        ),
        "electron_inventory_m2": _inventory(electron, dz_m),
        "ion_inventory_m2": _inventory(ion, dz_m),
        "photoelectron_inventory_m2": _inventory(photoelectron, dz_m),
        "total_outer_charge_c_m2": float(np.sum(field.rho_c_m3) * dz_m),
        "far_field_v_m": float(field.electric_faces_v_m[-1]),
        "far_charge_imbalance": float(
            abs(field.rho_c_m3[-1])
            / (ELEMENTARY_CHARGE * reference_density_m3)
        ),
        "gauss_residual_c_m2": abs(field.gauss_residual_c_m2),
    }


def _observable_scales(
    config: OuterKineticConfig, query: KineticQuery
) -> dict[str, float]:
    electron_sigma = math.sqrt(
        config.electron.temperature_ev
        * ELEMENTARY_CHARGE
        / config.electron.mass_kg
    )
    electron_flux = config.electron.number_density_m3 * electron_sigma
    ion_flux = (
        config.ion.number_density_m3
        * max(
            abs(config.ion.drift_velocity_mps),
            math.sqrt(
                config.ion.temperature_ev
                * ELEMENTARY_CHARGE
                / config.ion.mass_kg
            ),
        )
    )
    pe_flux = max(query.photoelectron_outward_number_flux_m2_s, ion_flux)
    charge_scale = max(
        abs(query.displacement_c_m2),
        ELEMENTARY_CHARGE
        * max(config.electron.number_density_m3, config.ion.number_density_m3)
        * config.z_length_m,
    )
    return {
        "matching_potential_v": max(
            config.electron.temperature_ev,
            query.photoelectron_outward_mean_normal_energy_ev,
            1.0,
        ),
        "electron_inward_number_flux_m2_s": electron_flux,
        "ion_inward_number_flux_m2_s": ion_flux,
        "photoelectron_return_number_flux_m2_s": pe_flux,
        "photoelectron_escape_number_flux_m2_s": pe_flux,
        "total_outer_charge_c_m2": charge_scale,
    }


def _window_indices(
    times: np.ndarray, *, end_s: float, width_s: float
) -> np.ndarray:
    return np.flatnonzero((times > end_s - width_s) & (times <= end_s))


def _integrated_autocorrelation_samples(values: np.ndarray) -> float | None:
    centered = values - float(np.mean(values))
    variance = float(np.dot(centered, centered) / centered.size)
    if variance <= np.finfo(float).eps:
        return 0.5
    correlation = np.correlate(centered, centered, mode="full")
    correlation = correlation[centered.size - 1 :] / (
        variance * np.arange(centered.size, 0, -1)
    )
    positive = correlation[1:]
    cutoff = np.flatnonzero(positive <= 0.0)
    count = int(cutoff[0]) if cutoff.size else positive.size
    tau = 0.5 + float(np.sum(positive[:count]))
    if not math.isfinite(tau) or tau <= 0.0:
        return None
    return tau


def _certify_time_history(
    config: OuterKineticConfig,
    query: KineticQuery,
    history: dict[str, np.ndarray],
) -> tuple[str, float, float, float, float]:
    cert = config.certification
    times = history["time_s"]
    second = _window_indices(
        times, end_s=config.max_time_s, width_s=cert.averaging_window_s
    )
    first = _window_indices(
        times,
        end_s=config.max_time_s - cert.averaging_window_s,
        width_s=cert.averaging_window_s,
    )
    if first.size < 4 or second.size < 4:
        return "unresolved_transient", math.inf, math.inf, math.inf, math.inf

    scales = _observable_scales(config, query)
    stationarity = 0.0
    drift = 0.0
    standard_error = 0.0
    fluctuation = 0.0
    autocorrelation_resolved = True
    window_duration = float(times[second[-1]] - times[second[0]])
    sample_spacing = window_duration / max(second.size - 1, 1)

    for name in _OBSERVABLE_NAMES:
        values_first = history[name][first]
        values_second = history[name][second]
        scale = max(scales[name], abs(float(np.mean(values_second))))
        stationarity = max(
            stationarity,
            abs(float(np.mean(values_second) - np.mean(values_first))) / scale,
        )
        relative_time = times[second] - float(np.mean(times[second]))
        denominator = float(np.dot(relative_time, relative_time))
        slope = (
            float(np.dot(relative_time, values_second - np.mean(values_second)))
            / denominator
            if denominator > 0.0
            else 0.0
        )
        drift = max(drift, abs(slope) * cert.averaging_window_s / scale)
        standard_deviation = float(np.std(values_second, ddof=1))
        fluctuation = max(fluctuation, standard_deviation / scale)
        tau_samples = _integrated_autocorrelation_samples(values_second)
        if tau_samples is None:
            autocorrelation_resolved = False
            continue
        effective_samples = max(
            1.0, values_second.size / max(2.0 * tau_samples, 1.0)
        )
        standard_error = max(
            standard_error,
            standard_deviation / math.sqrt(effective_samples) / scale,
        )
        if (
            tau_samples > 0.5
            and window_duration
            < cert.autocorrelation_windows * tau_samples * sample_spacing
        ):
            autocorrelation_resolved = False

    if (
        stationarity <= cert.stationarity_rtol
        and drift <= cert.drift_rtol
        and standard_error <= cert.sem_rtol
        and autocorrelation_resolved
    ):
        classification = (
            "steady"
            if fluctuation <= cert.steady_fluctuation_rtol
            else "stationary_average"
        )
    elif (
        stationarity > 3.0 * cert.stationarity_rtol
        or drift > 3.0 * cert.drift_rtol
    ):
        classification = "secular"
    else:
        classification = "unresolved_transient"
    return classification, stationarity, drift, standard_error, fluctuation


def _finalize_classification(
    time_classification: str,
    *,
    numerical_failure: str | None,
    budget_residual: float,
    velocity_loss: float,
    gauss_residual: float,
    far_field: float,
    far_charge: float,
    certification: CertificationConfig,
) -> tuple[str, str]:
    """Apply numerical and far-boundary gates without hiding time behavior."""

    if numerical_failure is not None:
        return "numerical_failure", numerical_failure
    if budget_residual > certification.charge_budget_rtol:
        return "numerical_failure", "species inventory budget tolerance exceeded"
    if velocity_loss > certification.velocity_loss_rtol:
        return "numerical_failure", "velocity-domain loss tolerance exceeded"
    if gauss_residual > certification.gauss_abs_c_m2:
        return "numerical_failure", "Gauss-law residual tolerance exceeded"
    if time_classification == "secular":
        return "secular", "time history has secular drift"
    if time_classification == "unresolved_transient":
        return (
            "unresolved_transient",
            "time history did not pass stationarity certification",
        )
    if (
        far_field > certification.far_field_abs_v_m
        or far_charge > certification.far_charge_rel
    ):
        return (
            "far_boundary_not_converged",
            "far-field or far-charge tolerance exceeded",
        )
    return time_classification, ""


def run_kinetic_query(
    config: OuterKineticConfig, query: KineticQuery
) -> KineticResult:
    """Run one fixed-input outer-sheath query and certify its mean response."""

    config.validate()
    query.validate()
    dz = config.z_length_m / config.nz
    electron = _make_reservoir_state(
        "ambient_electron",
        config.electron,
        nz=config.nz,
        dz_m=dz,
        initial_condition=config.initial_condition,
    )
    ion = _make_reservoir_state(
        "ambient_ion",
        config.ion,
        nz=config.nz,
        dz_m=dz,
        initial_condition=config.initial_condition,
    )
    photoelectron = _make_photoelectron_state(config, query, dz_m=dz)
    states = (electron, ion, photoelectron)
    history_lists: dict[str, list[float]] = {"time_s": []}
    for name in (
        *_OBSERVABLE_NAMES,
        "electron_access_potential_v",
        "ion_access_potential_v",
        "photoelectron_barrier_potential_v",
        "electron_inventory_m2",
        "ion_inventory_m2",
        "photoelectron_inventory_m2",
        "far_field_v_m",
        "far_charge_imbalance",
        "gauss_residual_c_m2",
    ):
        history_lists[name] = []

    field = _field_from_states(
        states, displacement_c_m2=query.displacement_c_m2, dz_m=dz
    )
    time_s = 0.0
    initial_sample = _sample_observables(
        states,
        field,
        dz_m=dz,
        reference_density_m3=max(
            config.electron.number_density_m3,
            config.ion.number_density_m3,
        ),
    )
    history_lists["time_s"].append(time_s)
    for name, value in initial_sample.items():
        history_lists[name].append(value)
    next_sample = config.certification.sample_interval_s
    steps = 0
    maximum_steps = 10_000_000
    numerical_failure: str | None = None

    try:
        while time_s < config.max_time_s:
            dt = _time_step(config, states, field, dz_m=dz)
            dt = min(dt, config.max_time_s - time_s)
            if next_sample > time_s:
                dt = min(dt, next_sample - time_s)
            if not math.isfinite(dt) or dt <= np.finfo(float).eps * max(
                config.max_time_s, 1.0
            ):
                raise FloatingPointError("kinetic CFL timestep underflow")
            for state in states:
                _advance_species_z(state, dt_s=0.5 * dt, dz_m=dz)
            field = _field_from_states(
                states, displacement_c_m2=query.displacement_c_m2, dz_m=dz
            )
            for state in states:
                _advance_species_v(state, field, dt_s=dt, dz_m=dz)
            for state in states:
                _advance_species_z(state, dt_s=0.5 * dt, dz_m=dz)
            _enforce_positive(states)
            time_s += dt
            steps += 1
            if steps > maximum_steps:
                raise FloatingPointError("kinetic step limit exceeded")
            field = _field_from_states(
                states, displacement_c_m2=query.displacement_c_m2, dz_m=dz
            )
            if time_s + 1.0e-15 * config.max_time_s >= next_sample:
                sample = _sample_observables(
                    states,
                    field,
                    dz_m=dz,
                    reference_density_m3=max(
                        config.electron.number_density_m3,
                        config.ion.number_density_m3,
                    ),
                )
                history_lists["time_s"].append(time_s)
                for name, value in sample.items():
                    history_lists[name].append(value)
                next_sample += config.certification.sample_interval_s
    except (FloatingPointError, OverflowError) as exc:
        numerical_failure = str(exc)

    history = {
        name: np.asarray(values, dtype=float)
        for name, values in history_lists.items()
    }
    classification, stationarity, drift, sem, _ = _certify_time_history(
        config, query, history
    )
    cert = config.certification
    final_inventories = tuple(_inventory(state, dz) for state in states)
    budget_residuals = []
    velocity_loss_fractions = []
    for state, final_inventory in zip(states, final_inventories, strict=True):
        expected = (
            state.initial_inventory_m2
            + state.boundary_inventory_change_m2
            + state.velocity_inventory_change_m2
        )
        scale = max(
            abs(state.initial_inventory_m2),
            abs(state.boundary_inventory_change_m2),
            abs(final_inventory),
            1.0,
        )
        budget_residuals.append(abs(final_inventory - expected) / scale)
        velocity_loss_fractions.append(
            abs(state.velocity_inventory_change_m2) / scale
        )
    budget_residual = max(budget_residuals)
    velocity_loss = max(velocity_loss_fractions)

    second = _window_indices(
        history["time_s"],
        end_s=config.max_time_s,
        width_s=cert.averaging_window_s,
    )
    if second.size == 0:
        second = np.arange(history["time_s"].size)
    far_field = float(np.max(np.abs(history["far_field_v_m"][second])))
    far_charge = float(np.max(history["far_charge_imbalance"][second]))
    classification, failure_reason = _finalize_classification(
        classification,
        numerical_failure=numerical_failure,
        budget_residual=budget_residual,
        velocity_loss=velocity_loss,
        gauss_residual=abs(field.gauss_residual_c_m2),
        far_field=far_field,
        far_charge=far_charge,
        certification=cert,
    )
    mean = {
        name: float(np.mean(history[name][second]))
        for name in history
        if name != "time_s"
    }
    response = (
        mean["matching_potential_v"],
        mean["electron_inward_number_flux_m2_s"],
        mean["ion_inward_number_flux_m2_s"],
        mean["electron_access_potential_v"],
        mean["ion_access_potential_v"],
        mean["photoelectron_barrier_potential_v"],
    )
    return KineticResult(
        query=query,
        response=response,
        photoelectron_return_flux_m2_s=mean[
            "photoelectron_return_number_flux_m2_s"
        ],
        photoelectron_escape_flux_m2_s=mean[
            "photoelectron_escape_number_flux_m2_s"
        ],
        classification=classification,
        failure_reason=failure_reason,
        simulated_time_s=time_s,
        averaging_window_s=cert.averaging_window_s,
        inventories_m2=tuple(float(value) for value in final_inventories),
        far_field_v_m=far_field,
        far_charge_imbalance=far_charge,
        gauss_residual_c_m2=abs(field.gauss_residual_c_m2),
        charge_budget_residual=budget_residual,
        mean_stationarity_metric=stationarity,
        max_drift_metric=drift,
        max_standard_error_metric=sem,
        max_velocity_boundary_loss_fraction=velocity_loss,
        time_steps=steps,
        profile=field,
        time_history=history,
    )


def write_kinetic_atlas(
    config: OuterKineticConfig,
    queries: list[KineticQuery],
    output_directory: str | Path,
    *,
    solver_version: str = "beach_python_fv_v1",
) -> list[KineticResult]:
    """Run queries and write raw rows, profiles, and a complete manifest."""

    output = Path(output_directory)
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"kinetic output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    profiles = output / "kinetic_response_profiles"
    profiles.mkdir(exist_ok=True)
    results: list[KineticResult] = []
    raw_path = output / "kinetic_response_raw.csv"
    temporary = raw_path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RAW_HEADER)
        writer.writeheader()
        for index, query in enumerate(queries):
            result = run_kinetic_query(config, query)
            results.append(result)
            writer.writerow(result.raw_row())
            np.savez_compressed(
                profiles / f"query_{index:06d}.npz",
                electric_faces_v_m=result.profile.electric_faces_v_m,
                potential_faces_v=result.profile.potential_faces_v,
                rho_c_m3=result.profile.rho_c_m3,
                **result.time_history,
            )
    temporary.replace(raw_path)
    raw_sha256 = _sha256_file(raw_path)
    accepted = sum(
        result.classification in ACCEPTED_CLASSIFICATIONS for result in results
    )
    manifest = {
        "schema": "beach_kinetic_response_v1",
        "solver_version": solver_version,
        "solver_source_sha256": SOLVER_SOURCE_SHA256,
        "git_commit": _git_commit(),
        "config_fingerprint_sha256": config_fingerprint(config),
        "configuration": asdict(config),
        "nz": config.nz,
        "nv_e": config.electron.grid.nv,
        "nv_i": config.ion.grid.nv,
        "nv_pe": config.photoelectron_grid.nv,
        "z_length_m": config.z_length_m,
        "z_max_m": config.matching_plane_z_m + config.z_length_m,
        "velocity_ranges_m_s": {
            "ambient_electron": [
                config.electron.grid.vmin_mps,
                config.electron.grid.vmax_mps,
            ],
            "ambient_ion": [
                config.ion.grid.vmin_mps,
                config.ion.grid.vmax_mps,
            ],
            "photoelectron": [
                config.photoelectron_grid.vmin_mps,
                config.photoelectron_grid.vmax_mps,
            ],
        },
        "cfl": config.cfl,
        "stationarity_tolerances": {
            "stationarity_rtol": config.certification.stationarity_rtol,
            "drift_rtol": config.certification.drift_rtol,
            "sem_rtol": config.certification.sem_rtol,
            "steady_fluctuation_rtol": (
                config.certification.steady_fluctuation_rtol
            ),
            "autocorrelation_windows": (
                config.certification.autocorrelation_windows
            ),
        },
        "far_boundary_tolerances": {
            "far_field_abs_v_m": config.certification.far_field_abs_v_m,
            "far_charge_rel": config.certification.far_charge_rel,
        },
        "query_count": len(results),
        "accepted_count": accepted,
        "rejected_count": len(results) - accepted,
        "classifications": {
            name: sum(result.classification == name for result in results)
            for name in sorted({result.classification for result in results})
        },
        "raw_csv": raw_path.name,
        "raw_csv_sha256": raw_sha256,
        "profiles_directory": profiles.name,
    }
    manifest_path = output / "kinetic_response_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return results


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
            cwd=Path(__file__).resolve().parents[1],
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    commit = result.stdout.strip()
    return commit if len(commit) == 40 else "unknown"


def _parse_csv_float(row: Mapping[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def convert_kinetic_table(
    raw_csv: str | Path,
    manifest_json: str | Path,
    output_csv: str | Path,
    *,
    ranges: Mapping[str, tuple[float, float]] | None = None,
    allow_stationary_average: bool = False,
) -> int:
    """Convert one fully certified Cartesian subset to response-table CSV v1."""

    manifest = json.loads(Path(manifest_json).read_text(encoding="utf-8"))
    if manifest.get("schema") != "beach_kinetic_response_v1":
        raise ValueError("kinetic manifest schema must be beach_kinetic_response_v1")
    expected_raw_hash = manifest.get("raw_csv_sha256")
    if not isinstance(expected_raw_hash, str) or len(expected_raw_hash) != 64:
        raise ValueError("kinetic manifest is missing raw_csv_sha256")
    if _sha256_file(Path(raw_csv)) != expected_raw_hash:
        raise ValueError("kinetic raw CSV does not match its manifest")
    configuration = manifest.get("configuration")
    if not isinstance(configuration, Mapping):
        raise ValueError("kinetic manifest is missing configuration")
    matching_plane_z = float(configuration["matching_plane_z_m"])
    allowed = {"steady"}
    if allow_stationary_average:
        allowed.add("stationary_average")
    rows: list[dict[str, str]] = []
    with Path(raw_csv).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != RAW_HEADER:
            raise ValueError("kinetic raw CSV header is not v1")
        for row in reader:
            for key in RESPONSE_HEADER:
                _parse_csv_float(row, key)
            include = True
            for key, bounds in (ranges or {}).items():
                if key not in QUERY_HEADER:
                    raise ValueError(f"unknown kinetic table range axis: {key}")
                value = _parse_csv_float(row, key)
                if value < bounds[0] or value > bounds[1]:
                    include = False
                    break
            if include:
                if row["classification"] not in allowed:
                    raise ValueError(
                        "kinetic table subset contains uncertified point: "
                        f"{row['classification']}"
                    )
                rows.append(row)
    if not rows:
        raise ValueError("kinetic table subset is empty")
    manifest_query_count = manifest.get("query_count")
    if not isinstance(manifest_query_count, int):
        raise ValueError("kinetic manifest is missing query_count")
    with Path(raw_csv).open(newline="", encoding="utf-8") as stream:
        raw_row_count = sum(1 for _ in csv.DictReader(stream))
    if raw_row_count != manifest_query_count:
        raise ValueError(
            "kinetic raw CSV row count does not match its manifest: "
            f"rows={raw_row_count} manifest={manifest_query_count}"
        )

    axes = {
        key: sorted({_parse_csv_float(row, key) for row in rows})
        for key in QUERY_HEADER
    }
    expected = math.prod(len(values) for values in axes.values())
    keys = {
        tuple(_parse_csv_float(row, key) for key in QUERY_HEADER) for row in rows
    }
    if len(keys) != len(rows):
        raise ValueError("kinetic table subset contains duplicate query points")
    if len(rows) != expected:
        raise ValueError(
            "kinetic table subset is not a complete Cartesian product: "
            f"rows={len(rows)} expected={expected}"
        )

    output = Path(output_csv)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        stream.write(f"# matching_plane_z_m={matching_plane_z:.17g}\n")
        writer = csv.DictWriter(stream, fieldnames=RESPONSE_HEADER)
        writer.writeheader()
        for row in sorted(
            rows,
            key=lambda item: tuple(
                _parse_csv_float(item, key) for key in QUERY_HEADER
            ),
        ):
            writer.writerow({key: row[key] for key in RESPONSE_HEADER})
    temporary.replace(output)
    return len(rows)
