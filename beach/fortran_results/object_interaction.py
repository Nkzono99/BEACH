"""Frozen-source object force and torque under BEACH field conventions."""

from __future__ import annotations

from dataclasses import dataclass, replace
import operator
from pathlib import Path
from threading import RLock
from typing import Iterable, Mapping

import numpy as np

from .context import RunContext
from .detachment import ObjectForcePath, ObjectWrench, WrenchComponent
from .kernel import (
    FieldKernel,
    FieldKernelError,
    FieldKernelOptions,
    _options_from_result,
)
from .mesh import _triangle_centers, _wrap_periodic2_triangles_by_mesh_centroid
from .panel_quadrature import _quadrature_order, panel_target_quadrature
from .periodic import coerce_periodic2
from .periodic_zero_mode import PeriodicZeroMode
from .scene import RigidTransform
from .selection import (
    _charges_for_step,
    _mesh_ids_or_default,
    _require_triangle_source_model,
    _require_triangles,
)
from .types import FortranRunResult


_PHYSICAL_COMPONENTS = (
    "other_objects_all_images",
    "target_periodic_images",
    "external_uniform",
    "total_external",
)

_PATH_TARGET_BATCH_LIMIT = 131_072
_EPS0_F_M = 8.8541878128e-12


@dataclass
class _VerticalSamples:
    displacement_m: np.ndarray
    force_N: np.ndarray
    torque_Nm: np.ndarray
    torque_origin_m: np.ndarray
    potential_energy_J: np.ndarray
    component_force_N: dict[str, np.ndarray]
    component_torque_Nm: dict[str, np.ndarray]
    peak_target_batch_size: int

    def take(self, mask: np.ndarray) -> "_VerticalSamples":
        return _VerticalSamples(
            displacement_m=self.displacement_m[mask],
            force_N=self.force_N[mask],
            torque_Nm=self.torque_Nm[mask],
            torque_origin_m=self.torque_origin_m[mask],
            potential_energy_J=self.potential_energy_J[mask],
            component_force_N={
                name: values[mask] for name, values in self.component_force_N.items()
            },
            component_torque_Nm={
                name: values[mask] for name, values in self.component_torque_Nm.items()
            },
            peak_target_batch_size=self.peak_target_batch_size,
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
        options: FieldKernelOptions,
        periodic_model: str,
        periodic: FieldKernel,
        zero_mode: PeriodicZeroMode | None,
        external_e0_V_m: np.ndarray,
        library_path: str | Path | None,
        zero_mode_lower_boundary: str | None = None,
        zero_mode_area_xy_m2: float | None = None,
    ) -> None:
        self.result = result
        self.step = step
        self.periodic_model = periodic_model
        self._triangles_m = _readonly(triangles_m)
        self._centers_m = _readonly(centers_m)
        self._charges_C = _readonly(charges_C)
        self._mesh_ids = _readonly_int(mesh_ids)
        self._options = options
        self._periodic = periodic
        self._zero_mode = zero_mode
        self._zero_mode_lower_boundary = zero_mode_lower_boundary
        self._zero_mode_area_xy_m2 = zero_mode_area_xy_m2
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
        context = RunContext.from_value(result, config_path=config_path)
        resolved = context.result
        _require_triangle_source_model(resolved)
        full_config = context.config
        if full_config is None:
            raise ValueError(
                "Object interaction requires the run's beach.toml so boundary, "
                "uniform-field, box, and outer-plasma policies cannot silently change."
            )
        _reject_active_outer_plasma(full_config)
        _validate_full_box_config(full_config)
        triangles = np.asarray(_require_triangles(resolved), dtype=np.float64)
        charges = np.asarray(_charges_for_step(resolved, step=step), dtype=np.float64)
        mesh_ids = np.asarray(_mesh_ids_or_default(resolved), dtype=np.int64)
        if charges.shape != (triangles.shape[0],) or not np.all(np.isfinite(charges)):
            raise ValueError("saved source charges must be finite and match the mesh.")

        options = _options_from_result(
            resolved,
            periodic2=None,
            theta=None,
            leaf_max=None,
            order=4,
            config_path=config_path,
            context=context,
        )
        external_e0 = np.asarray(options.external_e0, dtype=np.float64)
        periodic2 = options.periodic2
        if model == "infinite_physical":
            if periodic2 is None:
                raise ValueError(
                    "periodic_model='infinite_physical' requires an x/y periodic2 run."
                )
            periodic2 = periodic2._replace(far_correction="cached_kneq0")
        far_correction = "free" if periodic2 is None else periodic2.far_correction
        if far_correction not in {"free", "none", "cached_kneq0"}:
            raise ValueError(
                "Object interaction supports free, configured finite images, or "
                "cached_kneq0 physical periodic fields; "
                f"got {far_correction!r}."
            )
        if periodic2 is not None and periodic2.axes != (0, 1):
            raise ValueError("Physical periodic object interaction requires x/y axes (0, 1).")

        centers = _triangle_centers(triangles)

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

        periodic: FieldKernel | None = None
        zero_mode: PeriodicZeroMode | None = None
        zero_mode_lower_boundary: str | None = None
        zero_mode_area_xy_m2: float | None = None
        if far_correction == "cached_kneq0":
            assert periodic2 is not None
            zero_mode_lower_boundary = _resolve_zero_mode_lower_boundary(
                full_config
            )
            zero_mode_area_xy_m2 = float(
                periodic2.lengths[0] * periodic2.lengths[1]
            )
        try:
            periodic = FieldKernel(
                triangles,
                charges,
                options=options,
                library_path=library_path,
            )
            if far_correction == "cached_kneq0":
                assert periodic2 is not None
                heights = triangles[:, :, 2]
                assert zero_mode_lower_boundary is not None
                assert zero_mode_area_xy_m2 is not None
                zero_mode = PeriodicZeroMode(
                    heights,
                    charges,
                    zero_mode_area_xy_m2,
                    e_bottom_V_m=_zero_mode_bottom_field(
                        zero_mode_lower_boundary,
                        charges,
                        zero_mode_area_xy_m2,
                    ),
                    library_path=library_path,
                )
            return cls(
                result=resolved,
                step=step,
                triangles_m=triangles,
                centers_m=centers,
                charges_C=charges,
                mesh_ids=mesh_ids,
                options=options,
                periodic_model=model,
                periodic=periodic,
                zero_mode=zero_mode,
                zero_mode_lower_boundary=zero_mode_lower_boundary,
                zero_mode_area_xy_m2=zero_mode_area_xy_m2,
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
    def source_triangles_m(self) -> np.ndarray:
        """Return immutable saved panels used by the all-source field kernel."""

        return self._triangles_m

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
            if target_integration not in {"auto", "centroid_compatibility"}:
                raise ValueError(
                    "target_integration must be 'auto' or 'centroid_compatibility'."
                )
            order = _quadrature_order(quadrature_order)
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
                    self._update_zero_mode_charges(q_other)
                p_other_native = self._eval_periodic(target_points_m)
                z_other, trace_other = self._eval_zero_mechanical(target_points_m)
                p_other = _add_field_component(p_other_native, trace_other)

                self._periodic.update_charges(q_target)
                if self._zero_mode is not None:
                    self._update_zero_mode_charges(q_target)
                p_target_native = self._eval_periodic(target_points_m)
                z_target, trace_target = self._eval_zero_mechanical(target_points_m)
                p_target = _add_field_component(p_target_native, trace_target)

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
                    "trace_other": trace_other,
                    "trace_target": trace_target,
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

    def _eval_zero_mechanical(
        self,
        target_points_m: np.ndarray,
    ) -> tuple[
        tuple[np.ndarray, np.ndarray],
        tuple[np.ndarray, np.ndarray],
    ]:
        if self._zero_mode is None:
            zeros_e = np.zeros_like(target_points_m)
            zeros_phi = np.zeros(target_points_m.shape[0])
            return (zeros_e, zeros_phi), (zeros_e.copy(), zeros_phi.copy())
        phi, ez_principal = self._zero_mode.eval(
            target_points_m[:, 2],
            trace="principal_value",
        )
        _, ez_plus = self._zero_mode.eval(
            target_points_m[:, 2],
            trace="plus",
        )
        field = np.zeros_like(target_points_m)
        field[:, 2] = ez_principal
        trace_correction = np.zeros_like(target_points_m)
        trace_correction[:, 2] = ez_plus - ez_principal
        return (field, phi), (trace_correction, np.zeros(target_points_m.shape[0]))

    def _restore_full_state(self) -> list[Exception]:
        errors: list[Exception] = []
        try:
            self._periodic.update_charges(self._charges_C)
        except Exception as exc:
            errors.append(exc)
        if self._zero_mode is not None:
            try:
                self._update_zero_mode_charges(self._charges_C)
            except Exception as exc:
                errors.append(exc)
        return errors

    def _update_zero_mode_charges(self, charges_C: np.ndarray) -> None:
        if (
            self._zero_mode is None
            or self._zero_mode_lower_boundary is None
            or self._zero_mode_area_xy_m2 is None
        ):
            raise FieldKernelError("periodic zero-mode boundary state is incomplete.")
        self._zero_mode.update_charges(
            charges_C,
            e_bottom_V_m=_zero_mode_bottom_field(
                self._zero_mode_lower_boundary,
                charges_C,
                self._zero_mode_area_xy_m2,
            ),
        )

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
        target_triangles = snapshot._triangles_m[self._target_mask]
        self._target_geometry_representation = "saved"
        if snapshot._options.periodic2 is not None:
            periodic2 = snapshot._options.periodic2
            target_triangles = _wrap_periodic2_triangles_by_mesh_centroid(
                target_triangles,
                np.full(np.count_nonzero(self._target_mask), self.mesh_id),
                axes=periodic2.axes,
                lengths=periodic2.lengths,
                origins=periodic2.origins,
            )
            self._target_geometry_representation = "periodic2_mesh_connected"
        self._target_triangles_m = _readonly(target_triangles)
        self._target_centers_m = _readonly(_triangle_centers(target_triangles))
        self._target_charges_C = _readonly(snapshot._charges_C[self._target_mask])
        self._geometric_area_centroid_m = _geometric_area_centroid(
            self._target_triangles_m
        )
        target_vertices = self._target_triangles_m.reshape(-1, 3)
        self._vertex_bounding_center_m = _readonly(
            0.5 * (np.min(target_vertices, axis=0) + np.max(target_vertices, axis=0))
        )
        self._vertex_bounding_radius_m = float(
            np.max(
                np.linalg.norm(
                    target_vertices - self._vertex_bounding_center_m[None, :],
                    axis=1,
                )
            )
        )
        if target_integration == "auto":
            (
                target_points,
                target_charge_weights,
                target_element_index,
            ) = panel_target_quadrature(
                self._target_triangles_m,
                self._target_charges_C,
                quadrature_order,
            )
            integration_label = "gauss_duffy"
            integration_order: int | None = quadrature_order
        else:
            target_points = self._target_centers_m
            target_charge_weights = self._target_charges_C
            target_element_index = np.arange(
                self._target_charges_C.size,
                dtype=np.int64,
            )
            integration_label = "centroid_compatibility"
            integration_order = None
        self._target_points_m = _readonly(target_points)
        self._target_charge_weights_C = _readonly(target_charge_weights)
        self._target_element_index = _readonly_int(target_element_index)
        self._integration_label = integration_label
        self._integration_order = integration_order
        primary_options = FieldKernelOptions(
            theta=snapshot._options.theta,
            leaf_max=snapshot._options.leaf_max,
            order=snapshot._options.order,
        )
        self._primary = FieldKernel(
            self._target_triangles_m,
            self._target_charges_C,
            options=primary_options,
            library_path=snapshot._library_path,
        )
        self._closed = False

    @property
    def geometric_area_centroid_m(self) -> np.ndarray:
        """Return the immutable area-weighted centroid used by named origins."""

        return self._geometric_area_centroid_m

    @property
    def target_geometry_representation(self) -> str:
        """Return how the selected target mesh is represented for mechanics."""

        return self._target_geometry_representation

    @property
    def target_triangles_m(self) -> np.ndarray:
        """Return the immutable target panels in their connected geometry branch."""

        return self._target_triangles_m

    @property
    def vertex_bounding_center_m(self) -> np.ndarray:
        """Return the immutable center of the target vertex axis-aligned bounds."""

        return self._vertex_bounding_center_m

    @property
    def vertex_bounding_radius_m(self) -> float:
        """Return the largest target-vertex distance from ``vertex_bounding_center_m``."""

        return self._vertex_bounding_radius_m

    def wrench(
        self,
        transform: RigidTransform | None = None,
        transform_origin: str | Iterable[float] = "geometric_area_centroid",
        torque_origin: str | Iterable[float] = "geometric_area_centroid",
        components: bool = True,
    ) -> ObjectWrench:
        """Evaluate force and torque with primary-only self exclusion."""

        self._require_usable()
        rigid = RigidTransform.identity() if transform is None else transform
        if not isinstance(rigid, RigidTransform):
            raise TypeError("transform must be a RigidTransform or None.")
        transform_origin_m = _resolve_origin(
            transform_origin,
            geometric=self._geometric_area_centroid_m,
            name="transform_origin",
        )
        target_points = rigid.apply(
            self._target_points_m,
            origin=transform_origin_m,
        )
        target_vertices = rigid.apply(
            self._target_triangles_m.reshape(-1, 3),
            origin=transform_origin_m,
        )
        _validate_target_points(target_vertices, self._snapshot._options)
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
        trace_other = fields["trace_other"]
        trace_target = fields["trace_target"]
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
                self._target_charge_weights_C,
                other_field,
                other_phi,
                torque_origin_m,
            ),
            "target_periodic_images": _aggregate_wrench(
                target_points,
                self._target_charge_weights_C,
                images_field,
                images_phi,
                torque_origin_m,
            ),
            "external_uniform": _aggregate_wrench(
                target_points,
                self._target_charge_weights_C,
                uniform[0],
                uniform[1],
                torque_origin_m,
            ),
            "total_external": _aggregate_wrench(
                target_points,
                self._target_charge_weights_C,
                total_field,
                total_phi,
                torque_origin_m,
            ),
        }
        total = physical["total_external"]
        cached = (
            self._snapshot._options.periodic2 is not None
            and self._snapshot._options.periodic2.far_correction == "cached_kneq0"
        )
        p_full = _aggregate_wrench(
            target_points,
            self._target_charge_weights_C,
            p_other[0] + p_target[0],
            p_other[1] + p_target[1],
            torque_origin_m,
        )
        z_full = _aggregate_wrench(
            target_points,
            self._target_charge_weights_C,
            z_other[0] + z_target[0],
            z_other[1] + z_target[1],
            torque_origin_m,
        )
        negative_direct = _aggregate_wrench(
            target_points,
            self._target_charge_weights_C,
            -direct[0],
            -direct[1],
            torque_origin_m,
        )
        trace_full = _aggregate_wrench(
            target_points,
            self._target_charge_weights_C,
            trace_other[0] + trace_target[0],
            trace_other[1] + trace_target[1],
            torque_origin_m,
        )
        metadata: dict[str, object] = {
            "periodic_model": self._snapshot.periodic_model,
            "effective_far_correction": (
                "free"
                if self._snapshot._options.periodic2 is None
                else self._snapshot._options.periodic2.far_correction
            ),
            "target_integration": self._integration_label,
            "quadrature_order": self._integration_order,
            "target_geometry_representation": self.target_geometry_representation,
            "vertex_bounding_center_m": self._vertex_bounding_center_m,
            "vertex_bounding_radius_m": self._vertex_bounding_radius_m,
            "torque_origin_policy": _torque_origin_policy(torque_origin),
            "torque_origin_m": torque_origin_m,
            "periodic_kneq0": _wrench_metadata(p_full) if cached else None,
            "physical_k0": _wrench_metadata(z_full) if cached else None,
            "cached_kneq0_trace_correction": (
                _wrench_metadata(trace_full) if cached else None
            ),
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
            total_charge_C=float(np.sum(self._target_charge_weights_C)),
            force_N=total.force_N,
            torque_Nm=total.torque_Nm,
            torque_origin_m=torque_origin_m,
            transform=rigid,
            transform_origin_m=transform_origin_m,
            components=physical if components else {},
            numerical_metadata=metadata,
        )

    def vertical_path(
        self,
        displacement_m: np.ndarray,
        *,
        adaptive: bool = True,
        relative_tolerance: float = 5.0e-3,
        force_absolute_tolerance_N: float = 1.0e-12,
        work_absolute_tolerance_J: float = 1.0e-18,
        max_refinement: int = 8,
        torque_origin: str | Iterable[float] = "geometric_area_centroid",
        components: bool = True,
    ) -> ObjectForcePath:
        """Evaluate a frozen-source vertical translation path."""

        self._require_usable()
        grid = _path_displacement(displacement_m)
        relative = _nonnegative_scalar(relative_tolerance, "relative_tolerance")
        force_absolute = _nonnegative_scalar(
            force_absolute_tolerance_N,
            "force_absolute_tolerance_N",
        )
        work_absolute = _nonnegative_scalar(
            work_absolute_tolerance_J,
            "work_absolute_tolerance_J",
        )
        refinement_limit = _nonnegative_integer(max_refinement, "max_refinement")
        self._validate_vertical_geometry(grid)

        samples = self._evaluate_vertical_samples(
            grid,
            torque_origin=torque_origin,
            components=components,
        )
        inserted = 0
        refinement_rounds = 0
        peak_target_batch_size = samples.peak_target_batch_size
        status = "converged"
        status_reason = "fixed_grid"
        max_force_error = 0.0
        max_work_error = 0.0
        max_work_potential_error = 0.0

        if adaptive:
            status = "not_converged"
            status_reason = "max_refinement_reached"
            for refinement_round in range(refinement_limit + 1):
                middle = 0.5 * (grid[:-1] + grid[1:])
                candidate_grid = np.empty(2 * grid.size - 1, dtype=np.float64)
                candidate_grid[0::2] = grid
                candidate_grid[1::2] = middle
                self._validate_vertical_geometry(candidate_grid)
                candidate = self._evaluate_vertical_samples(
                    candidate_grid,
                    torque_origin=torque_origin,
                    components=components,
                )
                peak_target_batch_size = max(
                    peak_target_batch_size,
                    candidate.peak_target_batch_size,
                )
                force_left = candidate.force_N[0:-2:2]
                force_middle = candidate.force_N[1::2]
                force_right = candidate.force_N[2::2]
                widths = grid[1:] - grid[:-1]
                force_error = np.linalg.norm(
                    force_middle - 0.5 * (force_left + force_right),
                    axis=1,
                )
                force_scale = np.maximum.reduce(
                    (
                        np.linalg.norm(force_left, axis=1),
                        np.linalg.norm(force_middle, axis=1),
                        np.linalg.norm(force_right, axis=1),
                    )
                )
                work_coarse = 0.5 * widths * (
                    force_left[:, 2] + force_right[:, 2]
                )
                work_fine = 0.25 * widths * (
                    force_left[:, 2]
                    + 2.0 * force_middle[:, 2]
                    + force_right[:, 2]
                )
                potential_work = (
                    candidate.potential_energy_J[0:-2:2]
                    - candidate.potential_energy_J[2::2]
                )
                work_error = np.abs(work_fine - work_coarse)
                work_potential_error = np.abs(work_fine - potential_work)
                work_scale = np.maximum(np.abs(work_coarse), np.abs(work_fine))
                work_potential_scale = np.maximum(
                    np.abs(work_fine),
                    np.abs(potential_work),
                )
                force_failing = force_error > force_absolute + relative * force_scale
                work_failing = (
                    work_error > work_absolute + relative * work_scale
                )
                potential_failing = (
                    work_potential_error
                    > work_absolute + relative * work_potential_scale
                )
                failing = force_failing | work_failing | potential_failing
                max_force_error = float(np.max(force_error, initial=0.0))
                max_work_error = float(np.max(work_error, initial=0.0))
                max_work_potential_error = float(
                    np.max(work_potential_error, initial=0.0)
                )
                if not np.any(failing):
                    status = "converged"
                    status_reason = "tolerances_satisfied"
                    break
                if refinement_round >= refinement_limit:
                    if not np.any(force_failing | work_failing) and np.any(
                        potential_failing
                    ):
                        status_reason = "work_potential_mismatch"
                    break
                keep = np.zeros(candidate_grid.size, dtype=bool)
                keep[0::2] = True
                keep[1::2] = failing
                inserted_now = int(np.count_nonzero(failing))
                inserted += inserted_now
                refinement_rounds += 1
                grid = candidate_grid[keep]
                samples = candidate.take(keep)

        work = _cumulative_trapezoid(grid, samples.force_N[:, 2])
        potential_work = samples.potential_energy_J[0] - samples.potential_energy_J
        absolute_mismatch = float(np.max(np.abs(work - potential_work), initial=0.0))
        mismatch_scale = max(
            float(np.max(np.abs(work), initial=0.0)),
            float(np.max(np.abs(potential_work), initial=0.0)),
            np.finfo(float).tiny,
        )
        relative_mismatch = absolute_mismatch / mismatch_scale
        mismatch_threshold = work_absolute + relative * mismatch_scale
        if absolute_mismatch > mismatch_threshold and status == "converged":
            status = "not_converged"
            status_reason = "work_potential_mismatch"

        metadata: dict[str, object] = {
            "periodic_model": self._snapshot.periodic_model,
            "effective_far_correction": (
                "free"
                if self._snapshot._options.periodic2 is None
                else self._snapshot._options.periodic2.far_correction
            ),
            "target_integration": self._integration_label,
            "quadrature_order": self._integration_order,
            "source_geometry_policy": "frozen",
            "target_geometry_representation": self.target_geometry_representation,
            "vertex_bounding_center_m": self._vertex_bounding_center_m,
            "vertex_bounding_radius_m": self._vertex_bounding_radius_m,
            "torque_origin_policy": _torque_origin_policy(torque_origin),
            "torque_origin_m": samples.torque_origin_m,
            "target_motion": "vertical_translation",
            "adaptive": bool(adaptive),
            "relative_tolerance": relative,
            "force_absolute_tolerance_N": force_absolute,
            "work_absolute_tolerance_J": work_absolute,
            "max_refinement": refinement_limit,
            "refinement_rounds": refinement_rounds,
            "inserted_point_count": inserted,
            "max_force_refinement_error_N": max_force_error,
            "max_work_refinement_error_J": max_work_error,
            "max_work_potential_error_J": max_work_potential_error,
            "status_reason": status_reason,
            "peak_target_batch_size": peak_target_batch_size,
        }
        return ObjectForcePath(
            displacement_m=grid,
            force_N=samples.force_N,
            torque_Nm=samples.torque_Nm,
            electrostatic_work_J=work,
            potential_energy_J=samples.potential_energy_J,
            potential_difference_work_J=potential_work,
            component_force_N=samples.component_force_N,
            component_torque_Nm=samples.component_torque_Nm,
            numerical_metadata=metadata,
            status=status,
            refinement_count=inserted,
            work_relative_mismatch=relative_mismatch,
            work_absolute_mismatch_J=absolute_mismatch,
        )

    def _validate_vertical_geometry(self, displacement_m: np.ndarray) -> None:
        translations = np.zeros((displacement_m.size, 3), dtype=np.float64)
        translations[:, 2] = displacement_m
        vertices = (
            self._target_triangles_m[None, :, :, :]
            + translations[:, None, None, :]
        )
        _validate_target_points(
            vertices.reshape(-1, 3),
            self._snapshot._options,
        )

    def _evaluate_vertical_samples(
        self,
        displacement_m: np.ndarray,
        *,
        torque_origin: str | Iterable[float],
        components: bool,
    ) -> _VerticalSamples:
        nheight = displacement_m.size
        npoint = self._target_points_m.shape[0]
        height_batch = max(1, _PATH_TARGET_BATCH_LIMIT // max(1, npoint))
        force = np.empty((nheight, 3), dtype=np.float64)
        torque = np.empty((nheight, 3), dtype=np.float64)
        potential = np.empty(nheight, dtype=np.float64)
        component_force = {
            name: np.empty((nheight, 3), dtype=np.float64)
            for name in _PHYSICAL_COMPONENTS
        }
        component_torque = {
            name: np.empty((nheight, 3), dtype=np.float64)
            for name in _PHYSICAL_COMPONENTS
        }
        if isinstance(torque_origin, str) and torque_origin == "geometric_area_centroid":
            origins = np.broadcast_to(
                self._geometric_area_centroid_m,
                (nheight, 3),
            ).copy()
            origins[:, 2] += displacement_m
        else:
            origin = _resolve_origin(
                torque_origin,
                geometric=self._geometric_area_centroid_m,
                name="torque_origin",
            )
            origins = np.broadcast_to(origin, (nheight, 3)).copy()

        peak_target_batch_size = 0
        for start in range(0, nheight, height_batch):
            stop = min(start + height_batch, nheight)
            batch_h = displacement_m[start:stop]
            points = np.broadcast_to(
                self._target_points_m,
                (batch_h.size, npoint, 3),
            ).copy()
            points[:, :, 2] += batch_h[:, None]
            flat_points = points.reshape(-1, 3)
            peak_target_batch_size = max(peak_target_batch_size, flat_points.shape[0])
            fields = self._snapshot._evaluate_fields(
                target_mask=self._target_mask,
                target_points_m=flat_points,
                primary=self._primary,
            )
            physical_fields = _physical_field_components(fields)
            for name, (field_values, potential_values) in physical_fields.items():
                aggregate = _aggregate_wrench_samples(
                    points,
                    self._target_charge_weights_C,
                    field_values.reshape(batch_h.size, npoint, 3),
                    potential_values.reshape(batch_h.size, npoint),
                    origins[start:stop],
                )
                component_force[name][start:stop] = aggregate[0]
                component_torque[name][start:stop] = aggregate[1]
                if name == "total_external":
                    force[start:stop] = aggregate[0]
                    torque[start:stop] = aggregate[1]
                    potential[start:stop] = aggregate[2]

        return _VerticalSamples(
            displacement_m=np.array(displacement_m, copy=True),
            force_N=force,
            torque_Nm=torque,
            torque_origin_m=origins,
            potential_energy_J=potential,
            component_force_N=component_force if components else {},
            component_torque_Nm=component_torque if components else {},
            peak_target_batch_size=peak_target_batch_size,
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


def _aggregate_wrench_samples(
    points_m: np.ndarray,
    charges_C: np.ndarray,
    field_V_m: np.ndarray,
    potential_V: np.ndarray,
    origin_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    force_samples = charges_C[None, :, None] * field_V_m
    force = np.sum(force_samples, axis=1)
    torque = np.sum(
        np.cross(points_m - origin_m[:, None, :], force_samples),
        axis=1,
    )
    energy = potential_V @ charges_C
    return force, torque, energy


def _physical_field_components(
    fields: Mapping[str, tuple[np.ndarray, np.ndarray]],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    p_other = fields["p_other"]
    z_other = fields["z_other"]
    p_target = fields["p_target"]
    z_target = fields["z_target"]
    direct = fields["direct"]
    uniform = fields["uniform"]
    other = p_other[0] + z_other[0], p_other[1] + z_other[1]
    images = (
        p_target[0] + z_target[0] - direct[0],
        p_target[1] + z_target[1] - direct[1],
    )
    total = other[0] + images[0] + uniform[0], other[1] + images[1] + uniform[1]
    return {
        "other_objects_all_images": other,
        "target_periodic_images": images,
        "external_uniform": uniform,
        "total_external": total,
    }


def _add_field_component(
    base: tuple[np.ndarray, np.ndarray],
    correction: tuple[np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    return base[0] + correction[0], base[1] + correction[1]


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


def _torque_origin_policy(value: str | Iterable[float]) -> str:
    if isinstance(value, str):
        if value == "geometric_area_centroid":
            return "moving_geometric_area_centroid"
        if value == "origin":
            return "fixed_origin"
    return "fixed_explicit"


def _validate_target_points(points_m: np.ndarray, options: FieldKernelOptions) -> None:
    points = np.asarray(points_m, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("target points must have shape (n, 3) and contain finite values.")
    if options.box_min is None or options.box_max is None:
        return
    lower = np.asarray(options.box_min, dtype=np.float64)
    upper = np.asarray(options.box_max, dtype=np.float64)
    bounded_axes = np.ones(3, dtype=bool)
    periodic2 = coerce_periodic2(options.periodic2, allow_cached_kneq0=True)
    if periodic2 is not None:
        bounded_axes[np.asarray(periodic2.axes, dtype=np.int64)] = False
    bounded = points[:, bounded_axes]
    if np.any(bounded < lower[bounded_axes][None, :]) or np.any(
        bounded > upper[bounded_axes][None, :]
    ):
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
    external_boundary = config.get("external_boundary")
    if external_boundary is not None:
        if not isinstance(external_boundary, Mapping):
            raise ValueError("external_boundary config must be a table.")
        field = external_boundary.get("field")
        if not isinstance(field, Mapping):
            raise ValueError("external_boundary.field config must be a table.")
        model = str(field.get("model", "")).strip().lower()
        if model != "none":
            raise ValueError(
                "Object interaction does not yet support an active "
                "external_boundary.field; use "
                "external_boundary.field.model='none'."
            )
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


def _resolve_zero_mode_lower_boundary(config: Mapping[str, object]) -> str:
    periodic2 = config.get("periodic2")
    if not isinstance(periodic2, Mapping):
        raise ValueError(
            "cached_kneq0 object interaction requires an explicit [periodic2] "
            "table with zero_mode_policy and lower_boundary_model."
        )
    nonzero_backend = str(
        periodic2.get("nonzero_mode_backend", "")
    ).strip().lower()
    if nonzero_backend != "cached_kneq0":
        raise ValueError(
            "cached_kneq0 object interaction requires "
            "periodic2.nonzero_mode_backend='cached_kneq0'."
        )
    zero_mode_policy = str(periodic2.get("zero_mode_policy", "")).strip().lower()
    if zero_mode_policy != "exclude_k0":
        raise ValueError(
            "cached_kneq0 object interaction requires "
            "periodic2.zero_mode_policy='exclude_k0'."
        )
    lower_boundary = str(
        periodic2.get("lower_boundary_model", "")
    ).strip().lower()
    if lower_boundary not in {"e_bottom_zero", "symmetric_vacuum"}:
        raise ValueError(
            "cached_kneq0 object interaction requires "
            "periodic2.lower_boundary_model='e_bottom_zero' or "
            "'symmetric_vacuum'."
        )
    return lower_boundary


def _zero_mode_bottom_field(
    lower_boundary: str,
    charges_C: np.ndarray,
    area_xy_m2: float,
) -> float:
    if lower_boundary == "e_bottom_zero":
        return 0.0
    if lower_boundary == "symmetric_vacuum":
        return -float(np.sum(charges_C)) / (2.0 * _EPS0_F_M * area_xy_m2)
    raise ValueError(f"unsupported periodic2 lower boundary: {lower_boundary!r}.")


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
    return np.array(result, copy=True)


def _nonnegative_scalar(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and non-negative.")
    return result


def _nonnegative_integer(value: int, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a non-negative integer.") from exc
    if isinstance(value, (bool, np.bool_)) or result < 0:
        raise ValueError(f"{name} must be a non-negative integer.")
    return result


def _cumulative_trapezoid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    increments = 0.5 * (y[:-1] + y[1:]) * np.diff(x)
    return np.concatenate(([0.0], np.cumsum(increments)))


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
