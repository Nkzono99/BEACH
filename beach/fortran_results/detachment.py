"""Immutable force records and frozen-path detachment mechanics."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np


@dataclass(frozen=True)
class WrenchComponent:
    """One physical contribution to an object force and torque."""

    force_N: np.ndarray
    torque_Nm: np.ndarray
    potential_energy_J: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "force_N", _vec3(self.force_N, "force_N"))
        object.__setattr__(self, "torque_Nm", _vec3(self.torque_Nm, "torque_Nm"))
        if self.potential_energy_J is not None:
            object.__setattr__(
                self,
                "potential_energy_J",
                _finite_scalar(self.potential_energy_J, "potential_energy_J"),
            )


@dataclass(frozen=True, kw_only=True)
class ObjectWrench:
    """Immutable object-level force, torque, and component decomposition."""

    mesh_id: int
    step: int | None
    total_charge_C: float
    force_N: np.ndarray
    torque_Nm: np.ndarray
    torque_origin_m: np.ndarray
    transform: object | None = None
    transform_origin_m: np.ndarray | None = None
    components: Mapping[str, WrenchComponent] = field(default_factory=dict)
    numerical_metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "mesh_id", int(self.mesh_id))
        object.__setattr__(self, "total_charge_C", _finite_scalar(self.total_charge_C, "total_charge_C"))
        object.__setattr__(self, "force_N", _vec3(self.force_N, "force_N"))
        object.__setattr__(self, "torque_Nm", _vec3(self.torque_Nm, "torque_Nm"))
        object.__setattr__(self, "torque_origin_m", _vec3(self.torque_origin_m, "torque_origin_m"))
        if self.transform_origin_m is not None:
            object.__setattr__(
                self,
                "transform_origin_m",
                _vec3(self.transform_origin_m, "transform_origin_m"),
            )
        object.__setattr__(self, "transform", _freeze_transform(self.transform))
        component_copy: dict[str, WrenchComponent] = {}
        for name, component in self.components.items():
            if not isinstance(name, str) or not isinstance(component, WrenchComponent):
                raise TypeError("components must map strings to WrenchComponent values.")
            component_copy[name] = component
        object.__setattr__(self, "components", MappingProxyType(component_copy))
        object.__setattr__(self, "numerical_metadata", _freeze_mapping(self.numerical_metadata))


@dataclass(frozen=True)
class AdhesionProfile:
    """Finite-range object-level resisting force and its exact work."""

    kind: str
    displacement_m: np.ndarray = field(default_factory=lambda: np.empty(0))
    force_values_N: np.ndarray = field(default_factory=lambda: np.empty(0))
    constant_force_N: float = 0.0
    range_m: float = 0.0

    def __post_init__(self) -> None:
        if self.kind not in {"none", "finite_range_constant", "tabulated"}:
            raise ValueError("unknown adhesion profile kind.")
        displacement = np.asarray(self.displacement_m, dtype=np.float64)
        force = np.asarray(self.force_values_N, dtype=np.float64)
        if self.kind == "tabulated":
            if displacement.ndim != 1 or force.shape != displacement.shape:
                raise ValueError("adhesion displacement and force must be matching 1D arrays.")
            if displacement.size < 2:
                raise ValueError("tabulated adhesion requires at least two points.")
            if not np.all(np.isfinite(displacement)) or not np.all(np.isfinite(force)):
                raise ValueError("adhesion table must contain finite values.")
            if displacement[0] != 0.0:
                raise ValueError("adhesion displacement must start at zero.")
            if np.any(np.diff(displacement) <= 0.0):
                raise ValueError("adhesion displacement must be strictly increasing.")
            if np.any(force < 0.0):
                raise ValueError("adhesion force must be non-negative.")
        else:
            displacement = np.empty(0, dtype=np.float64)
            force = np.empty(0, dtype=np.float64)
        constant = _finite_scalar(self.constant_force_N, "constant_force_N")
        extent = _finite_scalar(self.range_m, "range_m")
        if constant < 0.0:
            raise ValueError("adhesion force must be non-negative.")
        if self.kind == "finite_range_constant" and extent <= 0.0:
            raise ValueError("adhesion range_m must be positive.")
        if self.kind != "finite_range_constant":
            constant = 0.0
            extent = 0.0
        object.__setattr__(self, "displacement_m", _readonly(displacement))
        object.__setattr__(self, "force_values_N", _readonly(force))
        object.__setattr__(self, "constant_force_N", constant)
        object.__setattr__(self, "range_m", extent)

    @classmethod
    def none(cls) -> "AdhesionProfile":
        return cls(kind="none")

    @classmethod
    def finite_range_constant(cls, force_N: float, range_m: float) -> "AdhesionProfile":
        return cls(
            kind="finite_range_constant",
            constant_force_N=force_N,
            range_m=range_m,
        )

    @classmethod
    def tabulated(
        cls,
        displacement_m: np.ndarray | list[float],
        force_N: np.ndarray | list[float],
    ) -> "AdhesionProfile":
        return cls(
            kind="tabulated",
            displacement_m=np.asarray(displacement_m, dtype=np.float64),
            force_values_N=np.asarray(force_N, dtype=np.float64),
        )

    @property
    def breakpoints_m(self) -> np.ndarray:
        if self.kind == "finite_range_constant":
            return _readonly(np.array([self.range_m]))
        if self.kind == "tabulated":
            return self.displacement_m
        return _readonly(np.empty(0))

    def force_N(self, displacement_m: np.ndarray | float) -> np.ndarray | float:
        h, scalar = _nonnegative_displacement(displacement_m)
        if self.kind == "none":
            result = np.zeros_like(h)
        elif self.kind == "finite_range_constant":
            result = np.where(h < self.range_m, self.constant_force_N, 0.0)
        else:
            result = np.interp(h, self.displacement_m, self.force_values_N)
            result = np.where(h > self.displacement_m[-1], 0.0, result)
        return float(result) if scalar else _readonly(result)

    def work_J(self, displacement_m: np.ndarray | float) -> np.ndarray | float:
        h, scalar = _nonnegative_displacement(displacement_m)
        if self.kind == "none":
            result = np.zeros_like(h)
        elif self.kind == "finite_range_constant":
            result = self.constant_force_N * np.minimum(h, self.range_m)
        else:
            result = _piecewise_linear_integral(
                self.displacement_m,
                self.force_values_N,
                np.minimum(h, self.displacement_m[-1]),
            )
        return float(result) if scalar else _readonly(result)


@dataclass(frozen=True, kw_only=True)
class ObjectForcePath:
    """Force and torque samples for a rigid object displacement path."""

    displacement_m: np.ndarray
    force_N: np.ndarray
    torque_Nm: np.ndarray
    electrostatic_work_J: np.ndarray
    potential_energy_J: np.ndarray | None = None
    potential_difference_work_J: np.ndarray | None = None
    component_force_N: Mapping[str, np.ndarray] = field(default_factory=dict)
    component_torque_Nm: Mapping[str, np.ndarray] = field(default_factory=dict)
    numerical_metadata: Mapping[str, object] = field(default_factory=dict)
    status: str = "converged"
    refinement_count: int = 0
    work_relative_mismatch: float | None = None

    def __post_init__(self) -> None:
        h = _path_displacement(self.displacement_m)
        npoint = h.size
        force = _path_vec(self.force_N, npoint, "force_N")
        torque = _path_vec(self.torque_Nm, npoint, "torque_Nm")
        work = _path_scalar(self.electrostatic_work_J, npoint, "electrostatic_work_J")
        if work[0] != 0.0:
            raise ValueError("electrostatic_work_J must start at zero.")
        potential = None
        potential_work = None
        if self.potential_energy_J is not None:
            potential = _path_scalar(self.potential_energy_J, npoint, "potential_energy_J")
        if self.potential_difference_work_J is not None:
            potential_work = _path_scalar(
                self.potential_difference_work_J,
                npoint,
                "potential_difference_work_J",
            )
        if self.status not in {"converged", "not_converged"}:
            raise ValueError('status must be "converged" or "not_converged".')
        if int(self.refinement_count) < 0:
            raise ValueError("refinement_count must be non-negative.")
        mismatch = self.work_relative_mismatch
        if mismatch is not None:
            mismatch = _finite_scalar(mismatch, "work_relative_mismatch")
            if mismatch < 0.0:
                raise ValueError("work_relative_mismatch must be non-negative.")

        object.__setattr__(self, "displacement_m", h)
        object.__setattr__(self, "force_N", force)
        object.__setattr__(self, "torque_Nm", torque)
        object.__setattr__(self, "electrostatic_work_J", work)
        object.__setattr__(self, "potential_energy_J", potential)
        object.__setattr__(self, "potential_difference_work_J", potential_work)
        object.__setattr__(self, "component_force_N", _freeze_array_mapping(self.component_force_N, npoint))
        object.__setattr__(self, "component_torque_Nm", _freeze_array_mapping(self.component_torque_Nm, npoint))
        object.__setattr__(self, "numerical_metadata", _freeze_mapping(self.numerical_metadata))
        object.__setattr__(self, "refinement_count", int(self.refinement_count))
        object.__setattr__(self, "work_relative_mismatch", mismatch)

    @classmethod
    def from_samples(
        cls,
        displacement_m: np.ndarray,
        force_N: np.ndarray,
        torque_Nm: np.ndarray,
        potential_energy_J: np.ndarray | None = None,
    ) -> "ObjectForcePath":
        h = _path_displacement(displacement_m)
        force = _path_vec(force_N, h.size, "force_N")
        torque = _path_vec(torque_Nm, h.size, "torque_Nm")
        work = _piecewise_linear_integral(h, force[:, 2], h)
        potential = None
        potential_work = None
        mismatch = None
        if potential_energy_J is not None:
            potential = _path_scalar(potential_energy_J, h.size, "potential_energy_J")
            potential_work = potential[0] - potential
            scale = max(float(np.max(np.abs(work))), float(np.max(np.abs(potential_work))), np.finfo(float).tiny)
            mismatch = float(np.max(np.abs(work - potential_work)) / scale)
        return cls(
            displacement_m=h,
            force_N=force,
            torque_Nm=torque,
            electrostatic_work_J=work,
            potential_energy_J=potential,
            potential_difference_work_J=potential_work,
            work_relative_mismatch=mismatch,
        )

    def evaluate_release(
        self,
        mass_kg: float,
        gravity_m_s2: float,
        adhesion: AdhesionProfile,
        eta_translation: float = 1.0,
        energy_tolerance_J: float = 1.0e-18,
        dissipation_work_J: np.ndarray | None = None,
    ) -> "DetachmentResult":
        """Evaluate from-rest accessibility for the continuous piecewise path."""

        mass = _finite_scalar(mass_kg, "mass_kg")
        gravity = _finite_scalar(gravity_m_s2, "gravity_m_s2")
        eta = _finite_scalar(eta_translation, "eta_translation")
        tolerance = _finite_scalar(energy_tolerance_J, "energy_tolerance_J")
        if mass <= 0.0:
            raise ValueError("mass_kg must be positive.")
        if gravity < 0.0:
            raise ValueError("gravity_m_s2 must be non-negative.")
        if not 0.0 <= eta <= 1.0:
            raise ValueError("eta_translation must be between zero and one.")
        if tolerance < 0.0:
            raise ValueError("energy_tolerance_J must be non-negative.")
        if not isinstance(adhesion, AdhesionProfile):
            raise TypeError("adhesion must be an AdhesionProfile.")

        h = self.displacement_m
        dissipation = _dissipation(dissipation_work_J, h.size)
        electric = self.electrostatic_work_J
        gravity_work = mass * gravity * h
        adhesion_work = np.asarray(adhesion.work_J(h), dtype=np.float64)
        gravity_corrected = electric - gravity_work
        available = gravity_corrected - adhesion_work - dissipation

        def available_at(query: np.ndarray | float) -> np.ndarray:
            x = np.asarray(query, dtype=np.float64)
            electric_x = _piecewise_linear_integral(h, self.force_N[:, 2], x)
            dissipation_x = np.interp(x, h, dissipation)
            adhesion_x = np.asarray(adhesion.work_J(x), dtype=np.float64)
            return electric_x - mass * gravity * x - adhesion_x - dissipation_x

        extra = adhesion.breakpoints_m
        extra = extra[(extra > h[0]) & (extra < h[-1])]
        mechanics_grid = np.unique(np.concatenate((h, extra)))
        critical = list(mechanics_grid)
        for left, right in zip(mechanics_grid[:-1], mechanics_grid[1:]):
            middle = 0.5 * (left + right)
            e0, em, e1 = available_at(np.array([left, middle, right]))
            a = 2.0 * (e1 + e0 - 2.0 * em)
            b = e1 - e0 - a
            scale = max(abs(e0), abs(em), abs(e1), np.finfo(float).tiny)
            if abs(a) > 32.0 * np.finfo(float).eps * scale:
                root = -b / (2.0 * a)
                if 0.0 < root < 1.0:
                    critical.append(float(left + root * (right - left)))
        critical_grid = np.unique(np.asarray(critical, dtype=np.float64))
        critical_energy = available_at(critical_grid)
        minimum = float(np.min(critical_energy))
        barrier_free = minimum >= -tolerance
        first_inaccessible = _first_inaccessible(
            critical_grid,
            critical_energy,
            available_at,
            tolerance,
        )

        electro_speed = _speed(electric, mass, eta)
        gravity_speed = _speed(gravity_corrected, mass, eta)
        speed = _speed(available, mass, eta)
        if first_inaccessible is None:
            reachable_mask = np.ones(critical_grid.size, dtype=bool)
        else:
            reachable_mask = critical_grid <= first_inaccessible
        reachable_energy = np.maximum(critical_energy[reachable_mask], 0.0)
        maximum_speed = float(_speed(reachable_energy, mass, eta).max(initial=0.0))
        initial_margin = float(
            self.force_N[0, 2]
            - mass * gravity
            - float(adhesion.force_N(0.0))
        )

        return DetachmentResult(
            mass_kg=mass,
            gravity_m_s2=gravity,
            eta_translation=eta,
            energy_tolerance_J=tolerance,
            displacement_m=h,
            electrostatic_work_J=electric,
            gravity_work_J=gravity_work,
            adhesion_work_J=adhesion_work,
            dissipation_work_J=dissipation,
            gravity_corrected_work_J=gravity_corrected,
            available_energy_J=available,
            electrostatic_only_speed_m_s=electro_speed,
            gravity_corrected_speed_m_s=gravity_speed,
            speed_m_s=speed,
            minimum_available_energy_J=minimum,
            barrier_free_from_rest=barrier_free,
            first_inaccessible_displacement_m=first_inaccessible,
            endpoint_available_energy_J=float(available[-1]),
            endpoint_positive=bool(available[-1] >= -tolerance),
            endpoint_reachable_from_rest=barrier_free,
            endpoint_speed_m_s=float(speed[-1]),
            maximum_reachable_speed_m_s=maximum_speed,
            instantaneous_force_margin_N=initial_margin,
        )


@dataclass(frozen=True, kw_only=True)
class DetachmentResult:
    """Work, speed, and continuous barrier diagnostics for one force path."""

    mass_kg: float
    gravity_m_s2: float
    eta_translation: float
    energy_tolerance_J: float
    displacement_m: np.ndarray
    electrostatic_work_J: np.ndarray
    gravity_work_J: np.ndarray
    adhesion_work_J: np.ndarray
    dissipation_work_J: np.ndarray
    gravity_corrected_work_J: np.ndarray
    available_energy_J: np.ndarray
    electrostatic_only_speed_m_s: np.ndarray
    gravity_corrected_speed_m_s: np.ndarray
    speed_m_s: np.ndarray
    minimum_available_energy_J: float
    barrier_free_from_rest: bool
    first_inaccessible_displacement_m: float | None
    endpoint_available_energy_J: float
    endpoint_positive: bool
    endpoint_reachable_from_rest: bool
    endpoint_speed_m_s: float
    maximum_reachable_speed_m_s: float
    instantaneous_force_margin_N: float

    def __post_init__(self) -> None:
        for name in (
            "displacement_m",
            "electrostatic_work_J",
            "gravity_work_J",
            "adhesion_work_J",
            "dissipation_work_J",
            "gravity_corrected_work_J",
            "available_energy_J",
            "electrostatic_only_speed_m_s",
            "gravity_corrected_speed_m_s",
            "speed_m_s",
        ):
            object.__setattr__(self, name, _readonly(np.asarray(getattr(self, name))))


def _path_displacement(value: np.ndarray) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.ndim != 1 or result.size < 2:
        raise ValueError("displacement_m must be a 1D array with at least two points.")
    if not np.all(np.isfinite(result)):
        raise ValueError("displacement_m must contain finite values.")
    if result[0] != 0.0:
        raise ValueError("displacement_m must start at zero.")
    if np.any(np.diff(result) <= 0.0):
        raise ValueError("displacement_m must be strictly increasing.")
    return _readonly(result)


def _path_vec(value: np.ndarray, npoint: int, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (npoint, 3):
        raise ValueError(f"{name} must have shape ({npoint}, 3).")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values.")
    return _readonly(result)


def _path_scalar(value: np.ndarray, npoint: int, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (npoint,):
        raise ValueError(f"{name} must have shape ({npoint},).")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values.")
    return _readonly(result)


def _piecewise_linear_integral(
    x: np.ndarray,
    y: np.ndarray,
    query: np.ndarray,
) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    query = np.asarray(query, dtype=np.float64)
    widths = np.diff(x)
    cumulative = np.concatenate(([0.0], np.cumsum(0.5 * (y[:-1] + y[1:]) * widths)))
    index = np.searchsorted(x, query, side="right") - 1
    index = np.clip(index, 0, x.size - 2)
    dx = query - x[index]
    slope = (y[index + 1] - y[index]) / widths[index]
    return cumulative[index] + y[index] * dx + 0.5 * slope * dx * dx


def _first_inaccessible(
    critical_grid: np.ndarray,
    critical_energy: np.ndarray,
    available_at,
    tolerance: float,
) -> float | None:
    shifted = critical_energy + tolerance
    if shifted[0] < 0.0:
        return float(critical_grid[0])
    for left, right, e_left, e_right in zip(
        critical_grid[:-1],
        critical_grid[1:],
        shifted[:-1],
        shifted[1:],
    ):
        if e_left < 0.0:
            return float(left)
        if e_right >= 0.0:
            continue
        if e_left == 0.0:
            return float(left)
        lo = float(left)
        hi = float(right)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if float(available_at(mid)) + tolerance >= 0.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    return None


def _dissipation(value: np.ndarray | None, npoint: int) -> np.ndarray:
    if value is None:
        return _readonly(np.zeros(npoint))
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (npoint,):
        raise ValueError(f"dissipation_work_J must have shape ({npoint},).")
    if not np.all(np.isfinite(result)):
        raise ValueError("dissipation_work_J must contain finite values.")
    if result[0] != 0.0:
        raise ValueError("dissipation_work_J must start at zero.")
    if np.any(np.diff(result) < 0.0):
        raise ValueError("dissipation_work_J must be non-decreasing.")
    return _readonly(result)


def _speed(energy: np.ndarray, mass: float, eta: float) -> np.ndarray:
    return np.sqrt(2.0 * eta * np.maximum(np.asarray(energy), 0.0) / mass)


def _nonnegative_displacement(value: np.ndarray | float) -> tuple[np.ndarray, bool]:
    result = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(result)) or np.any(result < 0.0):
        raise ValueError("displacement_m must be finite and non-negative.")
    return result, result.ndim == 0


def _vec3(value: np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (3,):
        raise ValueError(f"{name} must have shape (3,).")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values.")
    return _readonly(result)


def _readonly(value: np.ndarray) -> np.ndarray:
    result = np.array(value, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


def _finite_scalar(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _freeze_array_mapping(value: Mapping[str, np.ndarray], npoint: int) -> Mapping[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for name, array in value.items():
        if not isinstance(name, str):
            raise TypeError("component names must be strings.")
        result[name] = _path_vec(array, npoint, f"component {name}")
    return MappingProxyType(result)


def _freeze_mapping(value: Mapping[str, object]) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError("metadata must be a mapping.")
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise TypeError("metadata keys must be strings.")
        result[key] = _freeze_value(item)
    return MappingProxyType(result)


def _freeze_value(value: object) -> object:
    if isinstance(value, np.ndarray):
        return _readonly(value)
    if isinstance(value, Mapping):
        return _freeze_mapping(value)
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze_value(item) for item in value)
    if isinstance(value, (bytearray, memoryview)):
        return bytes(value)
    if isinstance(value, np.generic):
        return value.item()
    if value is None or isinstance(value, (str, bytes, bool, int, float, complex, Path)):
        return value
    raise TypeError(f"unsupported mutable metadata value: {type(value).__name__}.")


def _freeze_transform(value: object | None) -> object | None:
    if value is None:
        return None
    from .scene import RigidTransform

    if not isinstance(value, RigidTransform):
        raise TypeError("transform must be a RigidTransform or None.")
    result = RigidTransform(
        rotation=_readonly(value.rotation),
        translation_m=_readonly(value.translation_m),
    )
    result.rotation.setflags(write=False)
    result.translation_m.setflags(write=False)
    return result
