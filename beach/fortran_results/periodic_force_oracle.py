"""Finite periodic-image validation for frozen object probes."""

from __future__ import annotations

import operator
from contextlib import contextmanager
from dataclasses import dataclass, replace
from typing import Iterable, Iterator

import numpy as np

from .constants import K_COULOMB
from .detachment import ObjectForcePath, ObjectWrench, WrenchComponent
from .kernel import FieldKernel
from .object_interaction import (
    ObjectInteractionSnapshot,
    ObjectProbe,
    _aggregate_wrench,
)
from .periodic import Periodic2Config
from .scene import RigidTransform


EPS0 = 1.0 / (4.0 * np.pi * K_COULOMB)


@dataclass(frozen=True, kw_only=True)
class FiniteShellWrenchResult:
    """Raw and closure-corrected wrenches for one finite image shell."""

    image_layers: int
    closure_shift_V_m: np.ndarray
    symmetric: ObjectWrench
    e_bottom_zero: ObjectWrench
    selected: ObjectWrench
    selected_closure: str

    def __post_init__(self) -> None:
        layer = _nonnegative_integer(self.image_layers, "image_layers")
        shift = _readonly_vec3(self.closure_shift_V_m, "closure_shift_V_m")
        if self.selected_closure not in {"symmetric", "e_bottom_zero"}:
            raise ValueError(
                'selected_closure must be "symmetric" or "e_bottom_zero".'
            )
        expected = (
            self.symmetric
            if self.selected_closure == "symmetric"
            else self.e_bottom_zero
        )
        if self.selected is not expected:
            raise ValueError("selected wrench must match selected_closure.")
        object.__setattr__(self, "image_layers", layer)
        object.__setattr__(self, "closure_shift_V_m", shift)


@dataclass(frozen=True, kw_only=True)
class FiniteShellConvergenceResult:
    """Two-successive combined-gate convergence record for finite shells."""

    image_layers: np.ndarray
    symmetric_paths: tuple[ObjectForcePath, ...]
    corrected_paths: tuple[ObjectForcePath, ...]
    force_increment_error_N: np.ndarray
    work_increment_error_J: np.ndarray
    increment_converged: np.ndarray
    status: str
    selected_image_layers: int | None
    selected_path: ObjectForcePath | None
    force_tail_proxy_N: np.ndarray | None = None
    work_tail_proxy_J: np.ndarray | None = None
    reference_force_error_N: np.ndarray | None = None
    reference_work_error_J: np.ndarray | None = None
    reference_converged: np.ndarray | None = None
    reference_model: str | None = None

    def __post_init__(self) -> None:
        symmetric_paths = tuple(self.symmetric_paths)
        corrected_paths = tuple(self.corrected_paths)
        layers = np.asarray(self.image_layers)
        if layers.ndim != 1 or layers.size == 0:
            raise ValueError("image_layers must be a non-empty 1D array.")
        if not np.issubdtype(layers.dtype, np.integer) or np.any(layers < 0):
            raise ValueError("image_layers must contain non-negative integers.")
        layers = np.array(layers, dtype=np.int64, copy=True)
        if np.any(np.diff(layers) != 1):
            raise ValueError("image_layers must be consecutive.")
        if len(symmetric_paths) != layers.size or len(corrected_paths) != layers.size:
            raise ValueError("path tuples must match image_layers.")
        force_error = _readonly_nonnegative(
            self.force_increment_error_N,
            layers.size - 1,
            "force_increment_error_N",
        )
        work_error = _readonly_nonnegative(
            self.work_increment_error_J,
            layers.size - 1,
            "work_increment_error_J",
        )
        tail_factor = np.maximum(1, layers[1:]).astype(np.float64)
        force_tail = _readonly_nonnegative(
            force_error * tail_factor
            if self.force_tail_proxy_N is None
            else self.force_tail_proxy_N,
            layers.size - 1,
            "force_tail_proxy_N",
        )
        work_tail = _readonly_nonnegative(
            work_error * tail_factor
            if self.work_tail_proxy_J is None
            else self.work_tail_proxy_J,
            layers.size - 1,
            "work_tail_proxy_J",
        )
        increment = _readonly_boolean(
            self.increment_converged,
            layers.size - 1,
            "increment_converged",
        )
        if self.status not in {"converged", "not_converged"}:
            raise ValueError('status must be "converged" or "not_converged".')
        if self.status == "converged":
            if self.selected_image_layers is None or self.selected_path is None:
                raise ValueError("converged shell results require a selected path.")
            if increment.size < 2 or not np.all(increment[-2:]):
                raise ValueError(
                    "converged shell results require two successive combined gates."
                )
            if self.selected_image_layers != int(layers[-1]):
                raise ValueError("selected layer must be the final evaluated layer.")
            if self.selected_path is not corrected_paths[-1]:
                raise ValueError("selected path must be the final corrected path.")
        elif self.selected_image_layers is not None or self.selected_path is not None:
            raise ValueError("non-converged shell results cannot select a path.")
        reference_force: np.ndarray | None = None
        reference_work: np.ndarray | None = None
        reference_ok: np.ndarray | None = None
        if self.reference_model is None:
            if any(
                value is not None
                for value in (
                    self.reference_force_error_N,
                    self.reference_work_error_J,
                    self.reference_converged,
                )
            ):
                raise ValueError("reference arrays require reference_model.")
        else:
            if self.reference_model != "infinite_physical":
                raise ValueError('reference_model must be "infinite_physical" or None.')
            if any(
                value is None
                for value in (
                    self.reference_force_error_N,
                    self.reference_work_error_J,
                    self.reference_converged,
                )
            ):
                raise ValueError("reference_model requires all reference arrays.")
            reference_force = _readonly_nonnegative(
                self.reference_force_error_N,
                layers.size,
                "reference_force_error_N",
            )
            reference_work = _readonly_nonnegative(
                self.reference_work_error_J,
                layers.size,
                "reference_work_error_J",
            )
            reference_ok = _readonly_boolean(
                self.reference_converged,
                layers.size,
                "reference_converged",
            )
            if self.status == "converged" and not np.all(reference_ok[-2:]):
                raise ValueError(
                    "selected path requires two successive physical reference gates."
                )
        layers.setflags(write=False)
        object.__setattr__(self, "image_layers", layers)
        object.__setattr__(self, "symmetric_paths", symmetric_paths)
        object.__setattr__(self, "corrected_paths", corrected_paths)
        object.__setattr__(self, "force_increment_error_N", force_error)
        object.__setattr__(self, "work_increment_error_J", work_error)
        object.__setattr__(self, "force_tail_proxy_N", force_tail)
        object.__setattr__(self, "work_tail_proxy_J", work_tail)
        object.__setattr__(self, "increment_converged", increment)
        object.__setattr__(self, "reference_force_error_N", reference_force)
        object.__setattr__(self, "reference_work_error_J", reference_work)
        object.__setattr__(self, "reference_converged", reference_ok)


def finite_shell_wrench(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
    transform: RigidTransform | None,
    image_layers: int,
    closure: str,
    *,
    torque_origin: str | Iterable[float] = "geometric_area_centroid",
    components: bool = True,
) -> FiniteShellWrenchResult:
    """Evaluate one native finite shell with both explicit zero-mode closures."""

    layer = _nonnegative_integer(image_layers, "image_layers")
    selected_closure = str(closure).strip().lower()
    if selected_closure not in {"symmetric", "e_bottom_zero"}:
        raise ValueError('closure must be "symmetric" or "e_bottom_zero".')
    with _finite_probe(snapshot, probe, layer) as finite_probe:
        symmetric = finite_probe.wrench(
            transform=transform,
            torque_origin=torque_origin,
            components=components,
        )
        symmetric = _with_wrench_metadata(
            symmetric,
            image_layers=layer,
            closure="symmetric",
        )
        corrected, shift = _correct_wrench(
            snapshot,
            finite_probe,
            symmetric,
            components=components,
        )
    selected = symmetric if selected_closure == "symmetric" else corrected
    return FiniteShellWrenchResult(
        image_layers=layer,
        closure_shift_V_m=shift,
        symmetric=symmetric,
        e_bottom_zero=corrected,
        selected=selected,
        selected_closure=selected_closure,
    )


def finite_shell_convergence(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
    displacement_m: np.ndarray,
    *,
    max_layers: int = 12,
    relative_tolerance: float = 1.0e-2,
    force_floor_N: float = 1.0e-12,
    work_floor_J: float = 1.0e-18,
) -> FiniteShellConvergenceResult:
    """Increase shells until two tail/reference combined gates converge."""

    maximum = _nonnegative_integer(max_layers, "max_layers")
    relative = _nonnegative_scalar(relative_tolerance, "relative_tolerance")
    force_floor = _nonnegative_scalar(force_floor_N, "force_floor_N")
    work_floor = _nonnegative_scalar(work_floor_J, "work_floor_J")
    symmetric_paths: list[ObjectForcePath] = []
    corrected_paths: list[ObjectForcePath] = []
    force_errors: list[float] = []
    work_errors: list[float] = []
    force_tail_proxies: list[float] = []
    work_tail_proxies: list[float] = []
    increment_ok: list[bool] = []
    reference_force_errors: list[float] = []
    reference_work_errors: list[float] = []
    reference_ok: list[bool] = []
    consecutive = 0
    selected_layer: int | None = None
    selected_path: ObjectForcePath | None = None
    reference_path: ObjectForcePath | None = None
    reference_model: str | None = None
    path_grid = np.asarray(displacement_m, dtype=np.float64)
    if snapshot.periodic_model == "infinite_physical":
        reference_model = "infinite_physical"
        reference_path = probe.vertical_path(
            path_grid,
            adaptive=False,
            relative_tolerance=relative,
            force_absolute_tolerance_N=force_floor,
            work_absolute_tolerance_J=work_floor,
            components=False,
        )

    for layer in range(maximum + 1):
        with _finite_probe(snapshot, probe, layer) as finite_probe:
            symmetric = finite_probe.vertical_path(
                path_grid,
                adaptive=False,
                relative_tolerance=relative,
                force_absolute_tolerance_N=force_floor,
                work_absolute_tolerance_J=work_floor,
                components=True,
            )
            symmetric = _with_path_metadata(
                symmetric,
                image_layers=layer,
                closure="symmetric",
            )
            corrected, _ = _correct_path(snapshot, finite_probe, symmetric)
        symmetric_paths.append(symmetric)
        corrected_paths.append(corrected)
        if reference_path is not None:
            reference_force_error = float(
                np.max(
                    np.linalg.norm(
                        corrected.force_N - reference_path.force_N,
                        axis=1,
                    ),
                    initial=0.0,
                )
            )
            reference_work_error = float(
                np.max(
                    np.abs(
                        corrected.electrostatic_work_J
                        - reference_path.electrostatic_work_J
                    ),
                    initial=0.0,
                )
            )
            reference_force_scale = max(
                float(np.max(np.linalg.norm(corrected.force_N, axis=1), initial=0.0)),
                float(
                    np.max(
                        np.linalg.norm(reference_path.force_N, axis=1),
                        initial=0.0,
                    )
                ),
            )
            reference_work_scale = max(
                float(np.max(np.abs(corrected.electrostatic_work_J), initial=0.0)),
                float(
                    np.max(
                        np.abs(reference_path.electrostatic_work_J),
                        initial=0.0,
                    )
                ),
            )
            reference_force_errors.append(reference_force_error)
            reference_work_errors.append(reference_work_error)
            reference_ok.append(
                reference_force_error
                <= force_floor + relative * reference_force_scale
                and reference_work_error
                <= work_floor + relative * reference_work_scale
                and reference_path.status == "converged"
                and corrected.status == "converged"
            )
        if layer == 0:
            continue
        previous = corrected_paths[-2]
        if not np.array_equal(previous.displacement_m, corrected.displacement_m):
            raise ValueError("finite shell paths must share one displacement grid.")
        force_error = float(
            np.max(
                np.linalg.norm(corrected.force_N - previous.force_N, axis=1),
                initial=0.0,
            )
        )
        work_error = float(
            np.max(
                np.abs(
                    corrected.electrostatic_work_J
                    - previous.electrostatic_work_J
                ),
                initial=0.0,
            )
        )
        force_scale = max(
            float(np.max(np.linalg.norm(corrected.force_N, axis=1), initial=0.0)),
            float(np.max(np.linalg.norm(previous.force_N, axis=1), initial=0.0)),
        )
        work_scale = max(
            float(np.max(np.abs(corrected.electrostatic_work_J), initial=0.0)),
            float(np.max(np.abs(previous.electrostatic_work_J), initial=0.0)),
        )
        tail_factor = float(max(1, layer))
        force_tail_proxy = tail_factor * force_error
        work_tail_proxy = tail_factor * work_error
        converged = (
            force_tail_proxy <= force_floor + relative * force_scale
            and work_tail_proxy <= work_floor + relative * work_scale
            and corrected.status == "converged"
        )
        force_errors.append(force_error)
        work_errors.append(work_error)
        force_tail_proxies.append(force_tail_proxy)
        work_tail_proxies.append(work_tail_proxy)
        physical_reference_ok = reference_path is None or reference_ok[-1]
        combined_converged = converged and physical_reference_ok
        increment_ok.append(combined_converged)
        consecutive = consecutive + 1 if combined_converged else 0
        if consecutive >= 2:
            selected_layer = layer
            selected_path = corrected
            break

    status = "converged" if selected_path is not None else "not_converged"
    layers = np.arange(len(corrected_paths), dtype=np.int64)
    return FiniteShellConvergenceResult(
        image_layers=layers,
        symmetric_paths=tuple(symmetric_paths),
        corrected_paths=tuple(corrected_paths),
        force_increment_error_N=np.asarray(force_errors),
        work_increment_error_J=np.asarray(work_errors),
        force_tail_proxy_N=np.asarray(force_tail_proxies),
        work_tail_proxy_J=np.asarray(work_tail_proxies),
        increment_converged=np.asarray(increment_ok, dtype=bool),
        status=status,
        selected_image_layers=selected_layer,
        selected_path=selected_path,
        reference_force_error_N=(
            None if reference_path is None else np.asarray(reference_force_errors)
        ),
        reference_work_error_J=(
            None if reference_path is None else np.asarray(reference_work_errors)
        ),
        reference_converged=(
            None if reference_path is None else np.asarray(reference_ok, dtype=bool)
        ),
        reference_model=reference_model,
    )


@contextmanager
def _finite_probe(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
    image_layers: int,
) -> Iterator[ObjectProbe]:
    snapshot._require_usable()
    probe._require_usable()
    if probe._snapshot is not snapshot:
        raise ValueError("probe must belong to snapshot.")
    configured = snapshot._options.periodic2
    if configured is None:
        raise ValueError("finite-shell validation requires an x/y periodic snapshot.")
    periodic2 = Periodic2Config(
        axes=configured.axes,
        lengths=configured.lengths,
        origins=configured.origins,
        image_layers=image_layers,
        far_correction="none",
        ewald_alpha=configured.ewald_alpha,
        ewald_layers=configured.ewald_layers,
    )
    options = replace(snapshot._options, periodic2=periodic2)
    source_triangles = (
        snapshot._triangles_m if snapshot.source_model == "triangle_p0" else None
    )
    periodic = FieldKernel(
        snapshot._centers_m,
        snapshot._charges_C,
        options=options,
        library_path=snapshot._library_path,
        source_triangles=source_triangles,
    )
    finite_snapshot = ObjectInteractionSnapshot(
        result=snapshot.result,
        step=snapshot.step,
        triangles_m=snapshot._triangles_m,
        centers_m=snapshot._centers_m,
        charges_C=snapshot._charges_C,
        mesh_ids=snapshot._mesh_ids,
        source_model=snapshot.source_model,
        options=options,
        periodic_model=f"finite_shell_m{image_layers}",
        periodic=periodic,
        zero_mode=None,
        external_e0_V_m=snapshot._external_e0_V_m,
        library_path=snapshot._library_path,
    )
    try:
        finite_probe = finite_snapshot.object_probe(
            probe.mesh_id,
            target_integration=probe.target_integration,
            quadrature_order=probe.quadrature_order,
        )
        yield finite_probe
    finally:
        finite_snapshot.close()


def _correct_wrench(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
    symmetric: ObjectWrench,
    *,
    components: bool,
) -> tuple[ObjectWrench, np.ndarray]:
    transform = symmetric.transform
    if not isinstance(transform, RigidTransform):
        raise TypeError("finite shell wrench requires a rigid transform record.")
    assert symmetric.transform_origin_m is not None
    points = transform.apply(
        probe._target_points_m,
        origin=symmetric.transform_origin_m,
    )
    closure = _closure_wrenches(
        snapshot,
        probe,
        points,
        symmetric.torque_origin_m,
    )
    corrected_components = _correct_component_mapping(
        symmetric.components,
        closure,
        include=components,
    )
    metadata = dict(symmetric.numerical_metadata)
    metadata.update(
        {
            "finite_shell_closure": "e_bottom_zero",
            "closure_shift_V_m": closure["shift"],
            "closure_potential_gauge_m": np.array(
                [0.0, 0.0, _zero_gauge_z(snapshot)]
            ),
        }
    )
    total_closure = closure["total_external"]
    corrected = ObjectWrench(
        mesh_id=symmetric.mesh_id,
        step=symmetric.step,
        total_charge_C=symmetric.total_charge_C,
        force_N=symmetric.force_N + total_closure.force_N,
        torque_Nm=symmetric.torque_Nm + total_closure.torque_Nm,
        torque_origin_m=symmetric.torque_origin_m,
        transform=symmetric.transform,
        transform_origin_m=symmetric.transform_origin_m,
        components=corrected_components,
        numerical_metadata=metadata,
    )
    return corrected, closure["shift"]


def _with_wrench_metadata(
    wrench: ObjectWrench,
    *,
    image_layers: int,
    closure: str,
) -> ObjectWrench:
    metadata = dict(wrench.numerical_metadata)
    metadata.update(
        {
            "finite_shell_image_layers": image_layers,
            "finite_shell_closure": closure,
            "finite_shell_source_representation": "native_canonical_unwrapped",
        }
    )
    return ObjectWrench(
        mesh_id=wrench.mesh_id,
        step=wrench.step,
        total_charge_C=wrench.total_charge_C,
        force_N=wrench.force_N,
        torque_Nm=wrench.torque_Nm,
        torque_origin_m=wrench.torque_origin_m,
        transform=wrench.transform,
        transform_origin_m=wrench.transform_origin_m,
        components=wrench.components,
        numerical_metadata=metadata,
    )


def _with_path_metadata(
    path: ObjectForcePath,
    *,
    image_layers: int,
    closure: str,
) -> ObjectForcePath:
    metadata = dict(path.numerical_metadata)
    metadata.update(
        {
            "finite_shell_image_layers": image_layers,
            "finite_shell_closure": closure,
            "finite_shell_source_representation": "native_canonical_unwrapped",
        }
    )
    return ObjectForcePath(
        displacement_m=path.displacement_m,
        force_N=path.force_N,
        torque_Nm=path.torque_Nm,
        electrostatic_work_J=path.electrostatic_work_J,
        potential_energy_J=path.potential_energy_J,
        potential_difference_work_J=path.potential_difference_work_J,
        component_force_N=path.component_force_N,
        component_torque_Nm=path.component_torque_Nm,
        numerical_metadata=metadata,
        status=path.status,
        refinement_count=path.refinement_count,
        work_relative_mismatch=path.work_relative_mismatch,
        work_absolute_mismatch_J=path.work_absolute_mismatch_J,
    )


def _correct_path(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
    symmetric: ObjectForcePath,
) -> tuple[ObjectForcePath, np.ndarray]:
    h = symmetric.displacement_m
    points = np.broadcast_to(
        probe._target_points_m,
        (h.size, probe._target_points_m.shape[0], 3),
    ).copy()
    points[:, :, 2] += h[:, None]
    origins = np.broadcast_to(probe._geometric_area_centroid_m, (h.size, 3)).copy()
    origins[:, 2] += h
    closure = _closure_path_wrenches(snapshot, probe, points, origins)
    total_force = symmetric.force_N + closure["total_external"][0]
    total_torque = symmetric.torque_Nm + closure["total_external"][1]
    potential = symmetric.potential_energy_J
    if potential is None:
        raise ValueError("finite shell paths require potential energy samples.")
    corrected_potential = potential + closure["total_external"][2]
    work = _cumulative_trapezoid(h, total_force[:, 2])
    potential_work = corrected_potential[0] - corrected_potential
    absolute_mismatch = float(np.max(np.abs(work - potential_work), initial=0.0))
    scale = max(
        float(np.max(np.abs(work), initial=0.0)),
        float(np.max(np.abs(potential_work), initial=0.0)),
        np.finfo(float).tiny,
    )
    component_force = _correct_path_component_mapping(
        symmetric.component_force_N,
        closure,
        index=0,
    )
    component_torque = _correct_path_component_mapping(
        symmetric.component_torque_Nm,
        closure,
        index=1,
    )
    metadata = dict(symmetric.numerical_metadata)
    metadata.update(
        {
            "finite_shell_closure": "e_bottom_zero",
            "closure_shift_V_m": closure["shift"],
            "closure_potential_gauge_m": np.array(
                [0.0, 0.0, _zero_gauge_z(snapshot)]
            ),
        }
    )
    status = symmetric.status
    if {
        "relative_tolerance",
        "work_absolute_tolerance_J",
    }.issubset(metadata):
        relative_tolerance = float(metadata["relative_tolerance"])
        work_absolute_tolerance = float(metadata["work_absolute_tolerance_J"])
        mismatch_threshold = work_absolute_tolerance + relative_tolerance * scale
        mismatch_converged = absolute_mismatch <= mismatch_threshold
        source_reason = str(metadata.get("status_reason", ""))
        if not mismatch_converged:
            status = "not_converged"
            metadata["status_reason"] = "work_potential_mismatch"
        elif symmetric.status == "converged":
            status = "converged"
        elif source_reason == "work_potential_mismatch":
            status = "converged"
            metadata["status_reason"] = "tolerances_satisfied_after_closure"
    corrected = ObjectForcePath(
        displacement_m=h,
        force_N=total_force,
        torque_Nm=total_torque,
        electrostatic_work_J=work,
        potential_energy_J=corrected_potential,
        potential_difference_work_J=potential_work,
        component_force_N=component_force,
        component_torque_Nm=component_torque,
        numerical_metadata=metadata,
        status=status,
        refinement_count=symmetric.refinement_count,
        work_relative_mismatch=absolute_mismatch / scale,
        work_absolute_mismatch_J=absolute_mismatch,
    )
    return corrected, closure["shift"]


def _closure_wrenches(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
    points: np.ndarray,
    torque_origin_m: np.ndarray,
) -> dict[str, object]:
    shifts = _closure_source_shifts(snapshot, probe)
    gauge_z = _zero_gauge_z(snapshot)
    result: dict[str, object] = {"shift": shifts["total_external"]}
    for name in (
        "other_objects_all_images",
        "target_periodic_images",
        "total_external",
    ):
        field = np.broadcast_to(shifts[name], points.shape).copy()
        potential = -shifts[name][2] * (points[:, 2] - gauge_z)
        result[name] = _aggregate_wrench(
            points,
            probe._target_charge_weights_C,
            field,
            potential,
            torque_origin_m,
        )
    return result


def _closure_path_wrenches(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
    points: np.ndarray,
    origins: np.ndarray,
) -> dict[str, object]:
    shifts = _closure_source_shifts(snapshot, probe)
    gauge_z = _zero_gauge_z(snapshot)
    result: dict[str, object] = {"shift": shifts["total_external"]}
    for name in (
        "other_objects_all_images",
        "target_periodic_images",
        "total_external",
    ):
        field = np.broadcast_to(shifts[name], points.shape).copy()
        potential = -shifts[name][2] * (points[:, :, 2] - gauge_z)
        force_samples = probe._target_charge_weights_C[None, :, None] * field
        result[name] = (
            np.sum(force_samples, axis=1),
            np.sum(np.cross(points - origins[:, None, :], force_samples), axis=1),
            potential @ probe._target_charge_weights_C,
        )
    return result


def _closure_source_shifts(
    snapshot: ObjectInteractionSnapshot,
    probe: ObjectProbe,
) -> dict[str, np.ndarray]:
    configured = snapshot._options.periodic2
    assert configured is not None
    area = float(configured[1][0] * configured[1][1])
    target_charge = float(np.sum(snapshot._charges_C[probe._target_mask]))
    other_charge = float(np.sum(snapshot._charges_C[~probe._target_mask]))

    def shift(charge: float) -> np.ndarray:
        return np.array([0.0, 0.0, charge / (2.0 * EPS0 * area)])

    return {
        "other_objects_all_images": shift(other_charge),
        "target_periodic_images": shift(target_charge),
        "total_external": shift(other_charge + target_charge),
    }


def _correct_component_mapping(
    symmetric: dict[str, WrenchComponent] | object,
    closure: dict[str, object],
    *,
    include: bool,
) -> dict[str, WrenchComponent]:
    if not include:
        return {}
    if not hasattr(symmetric, "items"):
        raise TypeError("symmetric components must be a mapping.")
    source = dict(symmetric.items())  # type: ignore[union-attr]
    result = dict(source)
    for name in ("other_objects_all_images", "target_periodic_images"):
        result[name] = _sum_component(source[name], closure[name])
    result["total_external"] = _sum_component(
        source["total_external"],
        closure["total_external"],
    )
    return result


def _correct_path_component_mapping(
    symmetric: object,
    closure: dict[str, object],
    *,
    index: int,
) -> dict[str, np.ndarray]:
    if not hasattr(symmetric, "items"):
        return {}
    source = dict(symmetric.items())  # type: ignore[union-attr]
    if not source:
        return {}
    result = {name: np.array(values, copy=True) for name, values in source.items()}
    for name in (
        "other_objects_all_images",
        "target_periodic_images",
        "total_external",
    ):
        result[name] += closure[name][index]
    return result


def _sum_component(left: object, right: object) -> WrenchComponent:
    if not isinstance(left, WrenchComponent) or not isinstance(right, WrenchComponent):
        raise TypeError("wrench components must be WrenchComponent values.")
    potential = None
    if left.potential_energy_J is not None and right.potential_energy_J is not None:
        potential = left.potential_energy_J + right.potential_energy_J
    return WrenchComponent(
        force_N=left.force_N + right.force_N,
        torque_Nm=left.torque_Nm + right.torque_Nm,
        potential_energy_J=potential,
    )


def _zero_gauge_z(snapshot: ObjectInteractionSnapshot) -> float:
    if snapshot.source_model == "triangle_p0":
        return float(np.min(snapshot._triangles_m[:, :, 2]))
    return float(np.min(snapshot._centers_m[:, 2]))


def _cumulative_trapezoid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.concatenate(
        ([0.0], np.cumsum(0.5 * (y[:-1] + y[1:]) * np.diff(x)))
    )


def _nonnegative_integer(value: int, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a non-negative integer.") from exc
    if isinstance(value, (bool, np.bool_)) or result < 0:
        raise ValueError(f"{name} must be a non-negative integer.")
    return result


def _nonnegative_scalar(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and non-negative.")
    return result


def _readonly_vec3(value: np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain three finite values.")
    result = np.array(result, copy=True)
    result.setflags(write=False)
    return result


def _readonly_nonnegative(value: np.ndarray, size: int, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (size,) or not np.all(np.isfinite(result)) or np.any(result < 0.0):
        raise ValueError(f"{name} must contain {size} finite non-negative values.")
    result = np.array(result, copy=True)
    result.setflags(write=False)
    return result


def _readonly_boolean(value: np.ndarray, size: int, name: str) -> np.ndarray:
    result = np.asarray(value)
    if result.shape != (size,) or result.dtype != np.dtype(bool):
        raise ValueError(f"{name} must contain {size} boolean values.")
    result = np.array(result, copy=True)
    result.setflags(write=False)
    return result
