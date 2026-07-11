"""Frozen-source object force and torque under BEACH field conventions."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from threading import RLock
from typing import Iterable, Mapping

import numpy as np

from .detachment import ObjectWrench, WrenchComponent
from .kernel import (
    FieldKernel,
    FieldKernelError,
    FieldKernelOptions,
    _load_full_config,
    _options_from_result,
)
from .mesh import _triangle_centers
from .periodic_zero_mode import PeriodicZeroMode
from .scene import RigidTransform
from .selection import (
    _charges_for_step,
    _mesh_ids_or_default,
    _require_triangles,
    _resolve_result,
)
from .types import FortranRunResult


_PHYSICAL_COMPONENTS = (
    "other_objects_all_images",
    "target_periodic_images",
    "external_uniform",
    "total_external",
)


class ObjectInteractionSnapshot:
    """Immutable source geometry and charge state for object-level probes."""

    def __init__(
        self,
        *,
        result: FortranRunResult,
        step: int | None,
        triangles_m: np.ndarray,
        centers_m: np.ndarray,
        charges_C: np.ndarray,
        mesh_ids: np.ndarray,
        source_model: str,
        options: FieldKernelOptions,
        periodic_model: str,
        periodic: FieldKernel,
        zero_mode: PeriodicZeroMode | None,
        external_e0_V_m: np.ndarray,
        library_path: str | Path | None,
    ) -> None:
        self.result = result
        self.step = step
        self.periodic_model = periodic_model
        self.source_model = source_model
        self._triangles_m = _readonly(triangles_m)
        self._centers_m = _readonly(centers_m)
        self._charges_C = _readonly(charges_C)
        self._mesh_ids = _readonly_int(mesh_ids)
        self._options = options
        self._periodic = periodic
        self._zero_mode = zero_mode
        self._external_e0_V_m = _vec3(external_e0_V_m, "external_e0_V_m")
        self._library_path = library_path
        self._lock = RLock()
        self._closed = False
        self._poisoned = False
        self._probes: list[ObjectProbe] = []

    @classmethod
    def from_result(
        cls,
        result: FortranRunResult | object,
        *,
        step: int | None = -1,
        config_path: str | Path | None = None,
        periodic_model: str = "configured",
        cache_dir: str | Path | None = None,
        generation_tolerance: float | None = None,
        library_path: str | Path | None = None,
    ) -> "ObjectInteractionSnapshot":
        """Build one source snapshot from saved BEACH geometry and charges."""

        model = str(periodic_model).strip().lower()
        if model not in {"configured", "infinite_physical"}:
            raise ValueError(
                'periodic_model must be "configured" or "infinite_physical".'
            )
        resolved = _resolve_result(result)
        source_model = str(resolved.field_source_model).strip().lower()
        if source_model not in {"point", "triangle_p0"}:
            raise ValueError(
                "Object interaction requires field_source_model='point' or "
                f"'triangle_p0'; got {resolved.field_source_model!r}."
            )
        full_config = _load_full_config(resolved.directory, config_path=config_path)
        if full_config is None:
            raise ValueError(
                "Object interaction requires the run's beach.toml so boundary, "
                "uniform-field, box, and outer-plasma policies cannot silently change."
            )
        _reject_active_outer_plasma(full_config)
        _validate_full_box_config(full_config)
        triangles = np.asarray(_require_triangles(resolved), dtype=np.float64)
        centers = _triangle_centers(triangles)
        charges = np.asarray(_charges_for_step(resolved, step=step), dtype=np.float64)
        mesh_ids = np.asarray(_mesh_ids_or_default(resolved), dtype=np.int64)
        if charges.shape != (triangles.shape[0],) or not np.all(np.isfinite(charges)):
            raise ValueError("saved source charges must be finite and match the mesh.")

        options = _options_from_result(
            resolved,
            softening=None,
            periodic2=None,
            theta=None,
            leaf_max=None,
            order=4,
            config_path=config_path,
        )
        external_e0 = np.asarray(options.external_e0, dtype=np.float64)
        periodic2 = options.periodic2
        if model == "infinite_physical":
            if periodic2 is None:
                raise ValueError(
                    "periodic_model='infinite_physical' requires an x/y periodic2 run."
                )
            periodic2 = (
                periodic2[0],
                periodic2[1],
                periodic2[2],
                periodic2[3],
                "cached_kneq0",
                periodic2[5],
                periodic2[6],
            )
        far_correction = "free" if periodic2 is None else periodic2[4]
        if far_correction not in {"free", "none", "cached_kneq0"}:
            raise ValueError(
                "Object interaction supports free, configured finite images, or "
                "cached_kneq0 physical periodic fields; "
                f"got {far_correction!r}."
            )
        if periodic2 is not None and tuple(periodic2[0]) != (0, 1):
            raise ValueError("Physical periodic object interaction requires x/y axes (0, 1).")

        cache_path = options.periodic_cache_dir if cache_dir is None else str(cache_dir)
        tolerance = (
            options.periodic_generation_tolerance
            if generation_tolerance is None
            else float(generation_tolerance)
        )
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("generation_tolerance must be finite and positive.")
        options = replace(
            options,
            periodic2=periodic2,
            external_e0=(0.0, 0.0, 0.0),
            periodic_cache_dir=cache_path,
            periodic_generation_tolerance=tolerance,
        )

        source_triangles = triangles if source_model == "triangle_p0" else None
        periodic: FieldKernel | None = None
        zero_mode: PeriodicZeroMode | None = None
        try:
            periodic = FieldKernel(
                centers,
                charges,
                options=options,
                library_path=library_path,
                source_triangles=source_triangles,
            )
            if far_correction == "cached_kneq0":
                assert periodic2 is not None
                heights = (
                    triangles[:, :, 2]
                    if source_model == "triangle_p0"
                    else np.repeat(centers[:, 2, None], 3, axis=1)
                )
                zero_mode = PeriodicZeroMode(
                    heights,
                    charges,
                    float(periodic2[1][0] * periodic2[1][1]),
                    library_path=library_path,
                )
            return cls(
                result=resolved,
                step=step,
                triangles_m=triangles,
                centers_m=centers,
                charges_C=charges,
                mesh_ids=mesh_ids,
                source_model=source_model,
                options=options,
                periodic_model=model,
                periodic=periodic,
                zero_mode=zero_mode,
                external_e0_V_m=external_e0,
                library_path=library_path,
            )
        except Exception:
            if zero_mode is not None:
                zero_mode.close()
            if periodic is not None:
                periodic.close()
            raise

    @property
    def source_positions_m(self) -> np.ndarray:
        return self._centers_m

    @property
    def source_charges_C(self) -> np.ndarray:
        return self._charges_C

    def object_probe(
        self,
        mesh_id: int,
        *,
        self_policy: str = "exclude_primary_keep_images",
        target_integration: str = "auto",
        quadrature_order: int = 7,
    ) -> "ObjectProbe":
        """Create one independently movable target probe."""

        with self._lock:
            self._require_usable()
            if self_policy != "exclude_primary_keep_images":
                raise ValueError(
                    "self_policy must be 'exclude_primary_keep_images'. The legacy "
                    "exclude_target_lattice policy remains on calc_object_forces_kernel."
                )
            if target_integration != "auto":
                raise ValueError("target_integration must be 'auto' in this release phase.")
            order = int(quadrature_order)
            if order not in {3, 7}:
                raise ValueError("quadrature_order must be 3 or 7.")
            target_id = int(mesh_id)
            mask = self._mesh_ids == target_id
            if not np.any(mask):
                available = [int(value) for value in np.unique(self._mesh_ids)]
                raise ValueError(f"unknown mesh id {target_id}. available={available}")
            probe = ObjectProbe(
                snapshot=self,
                mesh_id=target_id,
                target_mask=np.asarray(mask, dtype=bool),
                target_integration=target_integration,
                quadrature_order=order,
            )
            self._probes.append(probe)
            return probe

    def close(self) -> None:
        """Close probes and native source evaluators."""

        with self._lock:
            if self._closed:
                return
            errors: list[Exception] = []
            for probe in tuple(self._probes):
                try:
                    probe.close()
                except Exception as exc:
                    errors.append(exc)
            if self._zero_mode is not None:
                try:
                    self._zero_mode.close()
                except Exception as exc:
                    errors.append(exc)
            try:
                self._periodic.close()
            except Exception as exc:
                errors.append(exc)
            if errors:
                self._poisoned = True
                raise FieldKernelError(
                    "failed to close one or more object-interaction evaluators: "
                    + "; ".join(str(error) for error in errors)
                ) from errors[0]
            self._closed = True

    def __enter__(self) -> "ObjectInteractionSnapshot":
        self._require_usable()
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _evaluate_fields(
        self,
        *,
        target_mask: np.ndarray,
        target_points_m: np.ndarray,
        primary: FieldKernel,
    ) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        with self._lock:
            self._require_usable()
            q_target = np.where(target_mask, self._charges_C, 0.0)
            q_other = np.where(target_mask, 0.0, self._charges_C)
            evaluation_error: BaseException | None = None
            try:
                self._periodic.update_charges(q_other)
                if self._zero_mode is not None:
                    self._zero_mode.update_charges(q_other)
                p_other = self._eval_periodic(target_points_m)
                z_other = self._eval_zero(target_points_m)

                self._periodic.update_charges(q_target)
                if self._zero_mode is not None:
                    self._zero_mode.update_charges(q_target)
                p_target = self._eval_periodic(target_points_m)
                z_target = self._eval_zero(target_points_m)

                direct = (
                    primary.eval_e_direct(target_points_m),
                    primary.eval_phi_direct(target_points_m),
                )
                uniform_e = np.broadcast_to(
                    self._external_e0_V_m,
                    target_points_m.shape,
                ).copy()
                gauge = np.zeros(3)
                uniform_phi = -(
                    (target_points_m - gauge[None, :]) @ self._external_e0_V_m
                )
                return {
                    "p_other": p_other,
                    "z_other": z_other,
                    "p_target": p_target,
                    "z_target": z_target,
                    "direct": direct,
                    "uniform": (uniform_e, uniform_phi),
                }
            except BaseException as exc:
                evaluation_error = exc
                raise
            finally:
                restore_errors = self._restore_full_state()
                if restore_errors:
                    self._poisoned = True
                    message = "; ".join(str(error) for error in restore_errors)
                    if evaluation_error is not None:
                        raise FieldKernelError(
                            "object field evaluation failed and source-state restore also "
                            f"failed: {message}"
                        ) from evaluation_error
                    raise FieldKernelError(
                        f"failed to restore object-interaction source state: {message}"
                    ) from restore_errors[0]

    def _eval_periodic(self, target_points_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return self._periodic.eval_e(target_points_m), self._periodic.eval_phi(
            target_points_m
        )

    def _eval_zero(self, target_points_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self._zero_mode is None:
            return np.zeros_like(target_points_m), np.zeros(target_points_m.shape[0])
        phi, ez = self._zero_mode.eval(
            target_points_m[:, 2],
            trace="principal_value",
        )
        field = np.zeros_like(target_points_m)
        field[:, 2] = ez
        return field, phi

    def _restore_full_state(self) -> list[Exception]:
        errors: list[Exception] = []
        try:
            self._periodic.update_charges(self._charges_C)
        except Exception as exc:
            errors.append(exc)
        if self._zero_mode is not None:
            try:
                self._zero_mode.update_charges(self._charges_C)
            except Exception as exc:
                errors.append(exc)
        return errors

    def _require_usable(self) -> None:
        if self._closed:
            raise FieldKernelError("object interaction snapshot is closed.")
        if self._poisoned:
            raise FieldKernelError(
                "object interaction snapshot is unavailable after state-restore failure."
            )


class ObjectProbe:
    """Selected object charge distribution movable within a frozen snapshot."""

    def __init__(
        self,
        *,
        snapshot: ObjectInteractionSnapshot,
        mesh_id: int,
        target_mask: np.ndarray,
        target_integration: str,
        quadrature_order: int,
    ) -> None:
        self._snapshot = snapshot
        self.mesh_id = mesh_id
        self.target_integration = target_integration
        self.quadrature_order = quadrature_order
        self._target_mask = np.array(target_mask, dtype=bool, copy=True)
        self._target_mask.setflags(write=False)
        self._target_centers_m = _readonly(snapshot._centers_m[self._target_mask])
        self._target_triangles_m = _readonly(snapshot._triangles_m[self._target_mask])
        self._target_charges_C = _readonly(snapshot._charges_C[self._target_mask])
        self._geometric_area_centroid_m = _geometric_area_centroid(
            self._target_triangles_m
        )
        primary_options = FieldKernelOptions(
            softening=snapshot._options.softening,
            theta=snapshot._options.theta,
            leaf_max=snapshot._options.leaf_max,
            order=snapshot._options.order,
        )
        source_triangles = (
            self._target_triangles_m
            if snapshot.source_model == "triangle_p0"
            else None
        )
        self._primary = FieldKernel(
            self._target_centers_m,
            self._target_charges_C,
            options=primary_options,
            library_path=snapshot._library_path,
            source_triangles=source_triangles,
        )
        self._closed = False

    def wrench(
        self,
        transform: RigidTransform | None = None,
        transform_origin: str | Iterable[float] = "geometric_area_centroid",
        torque_origin: str | Iterable[float] = "geometric_area_centroid",
        components: bool = True,
    ) -> ObjectWrench:
        """Evaluate force and torque with primary-only self exclusion."""

        self._require_usable()
        if self._snapshot.source_model == "triangle_p0":
            raise ValueError(
                "triangle_p0 target integration requires Task 6 quadrature; "
                "auto never falls back to centroid compatibility."
            )
        rigid = RigidTransform.identity() if transform is None else transform
        if not isinstance(rigid, RigidTransform):
            raise TypeError("transform must be a RigidTransform or None.")
        transform_origin_m = _resolve_origin(
            transform_origin,
            geometric=self._geometric_area_centroid_m,
            name="transform_origin",
        )
        target_points = rigid.apply(
            self._target_centers_m,
            origin=transform_origin_m,
        )
        _validate_target_points(target_points, self._snapshot._options)
        if isinstance(torque_origin, str) and torque_origin == "geometric_area_centroid":
            torque_origin_m = rigid.apply(
                self._geometric_area_centroid_m[None, :],
                origin=transform_origin_m,
            )[0]
        else:
            torque_origin_m = _resolve_origin(
                torque_origin,
                geometric=self._geometric_area_centroid_m,
                name="torque_origin",
            )

        fields = self._snapshot._evaluate_fields(
            target_mask=self._target_mask,
            target_points_m=target_points,
            primary=self._primary,
        )
        p_other = fields["p_other"]
        z_other = fields["z_other"]
        p_target = fields["p_target"]
        z_target = fields["z_target"]
        direct = fields["direct"]
        uniform = fields["uniform"]

        other_field = p_other[0] + z_other[0]
        other_phi = p_other[1] + z_other[1]
        images_field = p_target[0] + z_target[0] - direct[0]
        images_phi = p_target[1] + z_target[1] - direct[1]
        total_field = other_field + images_field + uniform[0]
        total_phi = other_phi + images_phi + uniform[1]
        physical = {
            "other_objects_all_images": _aggregate_wrench(
                target_points,
                self._target_charges_C,
                other_field,
                other_phi,
                torque_origin_m,
            ),
            "target_periodic_images": _aggregate_wrench(
                target_points,
                self._target_charges_C,
                images_field,
                images_phi,
                torque_origin_m,
            ),
            "external_uniform": _aggregate_wrench(
                target_points,
                self._target_charges_C,
                uniform[0],
                uniform[1],
                torque_origin_m,
            ),
            "total_external": _aggregate_wrench(
                target_points,
                self._target_charges_C,
                total_field,
                total_phi,
                torque_origin_m,
            ),
        }
        total = physical["total_external"]
        cached = (
            self._snapshot._options.periodic2 is not None
            and self._snapshot._options.periodic2[4] == "cached_kneq0"
        )
        p_full = _aggregate_wrench(
            target_points,
            self._target_charges_C,
            p_other[0] + p_target[0],
            p_other[1] + p_target[1],
            torque_origin_m,
        )
        z_full = _aggregate_wrench(
            target_points,
            self._target_charges_C,
            z_other[0] + z_target[0],
            z_other[1] + z_target[1],
            torque_origin_m,
        )
        negative_direct = _aggregate_wrench(
            target_points,
            self._target_charges_C,
            -direct[0],
            -direct[1],
            torque_origin_m,
        )
        metadata: dict[str, object] = {
            "periodic_model": self._snapshot.periodic_model,
            "effective_far_correction": (
                "free"
                if self._snapshot._options.periodic2 is None
                else self._snapshot._options.periodic2[4]
            ),
            "target_integration": "point_centroid",
            "quadrature_order": None,
            "periodic_kneq0": _wrench_metadata(p_full) if cached else None,
            "physical_k0": _wrench_metadata(z_full) if cached else None,
            "primary_free_subtraction": _wrench_metadata(negative_direct),
            "uniform_potential_gauge_m": np.zeros(3),
        }
        try:
            diagnostics = self._snapshot._periodic.diagnostics()
        except FieldKernelError:
            diagnostics = None
        if diagnostics is not None:
            metadata["periodic_cache"] = {
                "hit": diagnostics.periodic_cache_hit,
                "build_count": diagnostics.periodic_operator_build_count,
                "fingerprint": diagnostics.periodic_cache_fingerprint,
                "path": diagnostics.periodic_cache_path,
            }
        return ObjectWrench(
            mesh_id=self.mesh_id,
            step=self._snapshot.step,
            total_charge_C=float(np.sum(self._target_charges_C)),
            force_N=total.force_N,
            torque_Nm=total.torque_Nm,
            torque_origin_m=torque_origin_m,
            transform=rigid,
            transform_origin_m=transform_origin_m,
            components=physical if components else {},
            numerical_metadata=metadata,
        )

    def close(self) -> None:
        with self._snapshot._lock:
            if self._closed:
                return
            self._primary.close()
            self._closed = True

    def __enter__(self) -> "ObjectProbe":
        self._require_usable()
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _require_usable(self) -> None:
        if self._closed:
            raise FieldKernelError("object probe is closed.")
        self._snapshot._require_usable()


def _aggregate_wrench(
    points_m: np.ndarray,
    charges_C: np.ndarray,
    field_V_m: np.ndarray,
    potential_V: np.ndarray,
    origin_m: np.ndarray,
) -> WrenchComponent:
    force_samples = charges_C[:, None] * field_V_m
    force = np.sum(force_samples, axis=0)
    torque = np.sum(np.cross(points_m - origin_m[None, :], force_samples), axis=0)
    energy = float(np.dot(charges_C, potential_V))
    return WrenchComponent(force_N=force, torque_Nm=torque, potential_energy_J=energy)


def _wrench_metadata(component: WrenchComponent) -> Mapping[str, object]:
    return {
        "force_N": component.force_N,
        "torque_Nm": component.torque_Nm,
        "potential_energy_J": component.potential_energy_J,
    }


def _resolve_origin(
    value: str | Iterable[float],
    *,
    geometric: np.ndarray,
    name: str,
) -> np.ndarray:
    if isinstance(value, str):
        if value == "geometric_area_centroid":
            return geometric
        if value == "origin":
            return np.zeros(3)
        raise ValueError(
            f"{name} must be 'geometric_area_centroid', 'origin', or a 3-vector."
        )
    return _vec3(value, name)


def _validate_target_points(points_m: np.ndarray, options: FieldKernelOptions) -> None:
    points = np.asarray(points_m, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("target points must have shape (n, 3) and contain finite values.")
    if options.box_min is None or options.box_max is None:
        return
    lower = np.asarray(options.box_min, dtype=np.float64)
    upper = np.asarray(options.box_max, dtype=np.float64)
    if np.any(points < lower[None, :]) or np.any(points > upper[None, :]):
        raise ValueError("transformed target points must remain inside the configured box.")


def _geometric_area_centroid(triangles_m: np.ndarray) -> np.ndarray:
    edge1 = triangles_m[:, 1] - triangles_m[:, 0]
    edge2 = triangles_m[:, 2] - triangles_m[:, 0]
    area = 0.5 * np.linalg.norm(np.cross(edge1, edge2), axis=1)
    if np.any(~np.isfinite(area)) or np.any(area <= 0.0):
        raise ValueError("target triangles must be finite and non-degenerate.")
    centers = np.mean(triangles_m, axis=1)
    return _readonly(np.sum(area[:, None] * centers, axis=0) / np.sum(area))


def _reject_active_outer_plasma(config: Mapping[str, object] | None) -> None:
    if config is None:
        return
    outer = config.get("outer_plasma")
    if outer is None:
        return
    if not isinstance(outer, Mapping):
        raise ValueError("outer_plasma config must be a table.")
    model = str(outer.get("model", "none")).strip().lower()
    if model != "none":
        raise ValueError(
            "Object interaction does not yet support an active outer_plasma field; "
            "use outer_plasma.model='none'."
        )


def _validate_full_box_config(config: Mapping[str, object]) -> None:
    sim = config.get("sim")
    if not isinstance(sim, Mapping):
        raise ValueError("beach.toml must contain a [sim] table.")
    has_min = "box_min" in sim
    has_max = "box_max" in sim
    if not has_min or not has_max:
        raise ValueError("beach.toml must define both sim.box_min and sim.box_max.")
    lower = np.asarray(sim["box_min"], dtype=np.float64)
    upper = np.asarray(sim["box_max"], dtype=np.float64)
    if lower.shape != (3,) or upper.shape != (3,):
        raise ValueError("sim.box_min and sim.box_max must each contain three values.")
    if not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper)):
        raise ValueError("sim.box_min and sim.box_max must contain finite values.")
    if np.any(upper <= lower):
        raise ValueError("sim.box_max must be greater than sim.box_min on every axis.")


def _readonly(value: np.ndarray) -> np.ndarray:
    result = np.array(value, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


def _readonly_int(value: np.ndarray) -> np.ndarray:
    result = np.array(value, dtype=np.int64, copy=True)
    result.setflags(write=False)
    return result


def _vec3(value: Iterable[float], name: str) -> np.ndarray:
    result = np.asarray(list(value), dtype=np.float64)
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain exactly three finite values.")
    return _readonly(result)
