"""ctypes bindings for the BEACH Fortran field-kernel shared library."""

from __future__ import annotations

import ctypes
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .context import RunContext, load_config_for_output
from .mesh import _triangle_centers
from .panel_quadrature import panel_target_quadrature
from .periodic import (
    Periodic2Input,
    coerce_periodic2 as _coerce_periodic2,
    periodic2_from_sim as _periodic2_from_sim,
)
from .selection import (
    _charges_for_step,
    _mesh_ids_or_default,
    _require_triangle_source_model,
    _require_triangles,
)
from .types import FortranRunResult


_STATUS_MESSAGES = {
    0: "ok",
    1: "invalid kernel handle",
    2: "invalid kernel argument",
    3: "field kernel is not ready",
}

_FAR_CORRECTION_CODES = {
    "auto": 0,
    "none": 1,
    "cached_kneq0": 3,
}

FIELD_KERNEL_ABI_MAJOR = 2
FIELD_KERNEL_ABI_MINOR = 0


class FieldKernelError(RuntimeError):
    """Raised when the shared field kernel cannot be used."""


@dataclass(frozen=True)
class FieldKernelOptions:
    """Options passed to the Fortran Coulomb FMM kernel."""

    theta: float = 0.5
    leaf_max: int = 16
    order: int = 4
    periodic2: Periodic2Input | None = None
    box_min: tuple[float, float, float] | None = None
    box_max: tuple[float, float, float] | None = None
    external_e0: tuple[float, float, float] = (0.0, 0.0, 0.0)
    periodic_cache_dir: str = ".beach_cache/periodic2"
    periodic_generation_tolerance: float = 1.0e-8


@dataclass(frozen=True)
class FieldKernelDiagnostics:
    """Periodic operator-cache diagnostics for a built field kernel."""

    periodic_cache_hit: bool | None
    periodic_operator_build_count: int
    periodic_cache_fingerprint: str | None
    periodic_cache_path: Path | None


@dataclass(frozen=True)
class FieldKernelBuildInfo:
    """Source identity embedded in a BEACH field-kernel library."""

    schema_version: int
    version: str
    version_mode: str
    source_commit: str
    build_id: str


@dataclass(frozen=True)
class KernelObjectForceRecord:
    """Object-level net force/torque computed by the Fortran field kernel."""

    mesh_id: int
    step: int | None
    total_charge_C: float
    center_m: np.ndarray
    force_N: np.ndarray
    torque_Nm: np.ndarray


def field_kernel_build_info(
    library_path: str | Path | None = None,
) -> FieldKernelBuildInfo:
    """Read and validate build-origin metadata from a field-kernel library."""

    lib = _load_kernel_library(library_path)
    _configure_library(lib)
    getter = getattr(lib, "beach_kernel_get_build_info", None)
    if getter is None:
        raise FieldKernelError(
            "field-kernel build-info ABI requires shared-library symbol "
            "beach_kernel_get_build_info."
        )

    text_length = ctypes.c_int()

    def read_info(buffer: ctypes.Array[ctypes.c_char]) -> int:
        return int(
            getter(
                ctypes.cast(buffer, ctypes.c_void_p),
                ctypes.c_int(len(buffer)),
                ctypes.byref(text_length),
            )
        )

    buffer = ctypes.create_string_buffer(1)
    probe_status = read_info(buffer)
    if probe_status not in (0, 2):
        _check_status(probe_status, "beach_kernel_get_build_info")
    if not 0 < text_length.value <= 65535:
        raise FieldKernelError("field-kernel build-info returned an invalid text length.")

    buffer = ctypes.create_string_buffer(text_length.value + 1)
    status = read_info(buffer)
    _check_status(status, "beach_kernel_get_build_info")
    if not 0 < text_length.value < len(buffer):
        raise FieldKernelError("field-kernel build-info returned an invalid text length.")
    if buffer.raw[text_length.value] != 0:
        raise FieldKernelError("field-kernel build-info is not NUL terminated.")
    try:
        payload = buffer.raw[: text_length.value].decode("utf-8")
    except UnicodeDecodeError as exc:
        raise FieldKernelError("field-kernel build-info is not valid UTF-8.") from exc

    fields: dict[str, str] = {}
    for line in payload.split("\n"):
        key, separator, value = line.partition("=")
        if not separator or not key or not value or key in fields:
            raise FieldKernelError("field-kernel build-info payload is malformed.")
        fields[key] = value

    expected_keys = (
        "build_info_schema_version",
        "build_version",
        "build_version_mode",
        "build_source_commit",
        "build_id",
    )
    if tuple(fields) != expected_keys:
        raise FieldKernelError("field-kernel build-info payload has unexpected fields.")
    if fields["build_info_schema_version"] != "1":
        raise FieldKernelError("field-kernel build-info schema is not supported.")
    return FieldKernelBuildInfo(
        schema_version=1,
        version=fields["build_version"],
        version_mode=fields["build_version_mode"],
        source_commit=fields["build_source_commit"],
        build_id=fields["build_id"],
    )


class FieldKernel:
    """Thin Python wrapper around ``libbeach_field_kernel``.

    The kernel uses the same Fortran FMM core as the simulator. Source geometry
    is fixed after ``build``; charges can be refreshed cheaply with
    :meth:`update_charges`. A ``cached_kneq0`` periodic plan evaluates only the
    nonzero in-plane Fourier component; it is not the simulator's total field
    unless the physical zero mode and any outer-boundary contribution are
    composed separately.
    """

    def __init__(
        self,
        source_triangles: np.ndarray,
        source_charges: np.ndarray,
        *,
        options: FieldKernelOptions | None = None,
        library_path: str | Path | None = None,
    ) -> None:
        self._lib = _load_kernel_library(library_path)
        _configure_library(self._lib)
        self._handle = ctypes.c_void_p()
        self._closed = False
        self._source_count = 0
        self._options = FieldKernelOptions()
        self._cached_target_free_axis: int | None = None
        self._cached_target_free_bounds: tuple[float, float] | None = None
        status = self._lib.beach_kernel_create(ctypes.byref(self._handle))
        _check_status(status, "beach_kernel_create")
        try:
            self.build(source_triangles, options=options)
            self.update_charges(source_charges)
        except Exception:
            self.close()
            raise

    @classmethod
    def from_result(
        cls,
        result: FortranRunResult | object,
        *,
        step: int | None = -1,
        periodic2: Mapping[str, object] | None = None,
        theta: float | None = None,
        leaf_max: int | None = None,
        order: int = 4,
        config_path: str | Path | None = None,
        library_path: str | Path | None = None,
    ) -> "FieldKernel":
        """Build a kernel from one BEACH output directory."""

        context = RunContext.from_value(result, config_path=config_path)
        resolved = context.result
        _require_triangle_source_model(resolved)
        triangles = _require_triangles(resolved)
        charges = _charges_for_step(resolved, step=step)
        options = _options_from_result(
            resolved,
            periodic2=periodic2,
            theta=theta,
            leaf_max=leaf_max,
            order=order,
            config_path=config_path,
            context=context,
        )
        return cls(
            triangles,
            charges,
            options=options,
            library_path=library_path,
        )

    @staticmethod
    def is_available(library_path: str | Path | None = None) -> bool:
        """Return whether the shared kernel library can be loaded."""

        try:
            lib = _load_kernel_library(library_path)
            _configure_library(lib)
        except FieldKernelError:
            return False
        return True

    def build(
        self,
        source_triangles: np.ndarray,
        *,
        options: FieldKernelOptions | None = None,
    ) -> None:
        """Build or rebuild a constant-density triangle-panel plan."""

        self._require_open()
        opts = options or FieldKernelOptions()
        vertices = _triangles_to_fortran_vertices(source_triangles)
        self._build_geometry(vertices, opts)

    def _build_geometry(
        self,
        panel_vertices: tuple[np.ndarray, np.ndarray, np.ndarray],
        opts: FieldKernelOptions,
    ) -> None:
        nsrc = panel_vertices[0].shape[1]
        if nsrc <= 0:
            raise ValueError("source_triangles must contain at least one triangle.")
        if opts.order < 1:
            raise ValueError("FMM expansion order must be >= 1.")

        periodic_cfg = _coerce_periodic2(opts.periodic2, allow_cached_kneq0=True)
        periodic_axes = _null_ptr()
        periodic_len = _null_ptr()
        box_min = _null_ptr()
        box_max = _null_ptr()
        use_periodic2 = 0
        image_layers = 1
        far_correction = 0
        ewald_alpha = 0.0
        ewald_layers = 4
        cached_target_free_axis: int | None = None
        cached_target_free_bounds: tuple[float, float] | None = None
        keepalive: list[np.ndarray] = [*panel_vertices]

        far_key = "none"
        if periodic_cfg is not None:
            axes = periodic_cfg.axes
            lengths = periodic_cfg.lengths
            image_layers = periodic_cfg.image_layers
            far_key = periodic_cfg.far_correction
            ewald_alpha = periodic_cfg.ewald_alpha
            ewald_layers = periodic_cfg.ewald_layers
            box_min_vec, box_max_vec = _periodic_box_vectors(
                axes=axes,
                lengths=lengths,
                origins=periodic_cfg.origins,
                box_min=opts.box_min,
                box_max=opts.box_max,
                source_vertices_3xn=panel_vertices,
            )
            axes_1based = np.ascontiguousarray(np.asarray(axes, dtype=np.int32) + 1)
            lengths_vec = np.ascontiguousarray(np.asarray(lengths, dtype=np.float64))
            box_min_arr = np.ascontiguousarray(np.asarray(box_min_vec, dtype=np.float64))
            box_max_arr = np.ascontiguousarray(np.asarray(box_max_vec, dtype=np.float64))
            keepalive.extend([axes_1based, lengths_vec, box_min_arr, box_max_arr])
            periodic_axes = axes_1based.ctypes.data_as(ctypes.c_void_p)
            periodic_len = lengths_vec.ctypes.data_as(ctypes.c_void_p)
            box_min = box_min_arr.ctypes.data_as(ctypes.c_void_p)
            box_max = box_max_arr.ctypes.data_as(ctypes.c_void_p)
            use_periodic2 = 1
            far_correction = _far_correction_code(far_key)
            if far_key == "cached_kneq0":
                cached_target_free_axis = next(
                    axis for axis in range(3) if axis not in axes
                )
                cached_target_free_bounds = (
                    float(box_min_vec[cached_target_free_axis]),
                    float(box_max_vec[cached_target_free_axis]),
                )

        self._set_periodic_cache_options(opts, far_key=far_key)

        status = self._lib.beach_kernel_build(
            self._handle,
            ctypes.c_int(nsrc),
            *[
                vertex.ctypes.data_as(ctypes.c_void_p)
                for vertex in panel_vertices
            ],
            ctypes.c_double(opts.theta),
            ctypes.c_int(opts.leaf_max),
            ctypes.c_int(opts.order),
            ctypes.c_int(use_periodic2),
            periodic_axes,
            periodic_len,
            ctypes.c_int(image_layers),
            ctypes.c_int(far_correction),
            ctypes.c_double(ewald_alpha),
            ctypes.c_int(ewald_layers),
            box_min,
            box_max,
        )
        _check_status(status, "beach_kernel_build")
        self._source_count = nsrc
        self._options = opts
        self._cached_target_free_axis = cached_target_free_axis
        self._cached_target_free_bounds = cached_target_free_bounds
        self._keepalive = keepalive

    def _set_periodic_cache_options(
        self, opts: FieldKernelOptions, *, far_key: str
    ) -> None:
        setter = getattr(self._lib, "beach_kernel_set_periodic_cache", None)
        requires_cached_model = far_key == "cached_kneq0"
        requires_nondefault_options = (
            opts.periodic_cache_dir != FieldKernelOptions.periodic_cache_dir
            or opts.periodic_generation_tolerance
            != FieldKernelOptions.periodic_generation_tolerance
        )
        if setter is None:
            if requires_cached_model:
                raise FieldKernelError(
                    "cached_kneq0 field kernels require shared-library symbol "
                    "beach_kernel_set_periodic_cache."
                )
            if requires_nondefault_options:
                raise FieldKernelError(
                    "non-default periodic cache configuration requires shared-library symbol "
                    "beach_kernel_set_periodic_cache."
                )
            return

        try:
            path_bytes = str(opts.periodic_cache_dir).encode("utf-8")
        except UnicodeEncodeError as exc:
            raise ValueError("periodic_cache_dir must be valid UTF-8.") from exc
        path_buffer = ctypes.create_string_buffer(path_bytes)
        status = setter(
            self._handle,
            ctypes.cast(path_buffer, ctypes.c_void_p),
            ctypes.c_int(len(path_bytes)),
            ctypes.c_double(opts.periodic_generation_tolerance),
        )
        _check_status(status, "beach_kernel_set_periodic_cache")

    def update_charges(self, source_charges: np.ndarray) -> None:
        """Refresh source charges without rebuilding source geometry."""

        self._require_open()
        q = _charges_1d(source_charges, expected=self._source_count, name="source_charges")
        status = self._lib.beach_kernel_update_charges(
            self._handle,
            ctypes.c_int(self._source_count),
            q.ctypes.data_as(ctypes.c_void_p),
        )
        _check_status(status, "beach_kernel_update_charges")
        self._charges_keepalive = q

    def eval_e(self, points: np.ndarray) -> np.ndarray:
        """Evaluate electric field vectors at points with shape ``(n, 3)``.

        For ``cached_kneq0``, this returns the ``k!=0`` component only.
        """

        self._require_open()
        target_pos = _points_to_fortran_3xn(points, name="points")
        self._validate_cached_target_points(target_pos)
        ntarget = target_pos.shape[1]
        e = np.zeros((3, ntarget), dtype=np.float64, order="F")
        if ntarget == 0:
            return np.empty((0, 3), dtype=np.float64)
        status = self._lib.beach_kernel_eval_e(
            self._handle,
            ctypes.c_int(ntarget),
            target_pos.ctypes.data_as(ctypes.c_void_p),
            e.ctypes.data_as(ctypes.c_void_p),
        )
        _check_status(status, "beach_kernel_eval_e")
        field = np.ascontiguousarray(e.T)
        external_e0 = np.asarray(self._options.external_e0, dtype=np.float64)
        if np.any(external_e0 != 0.0):
            field = field + external_e0
        return field

    def eval_phi(self, points: np.ndarray) -> np.ndarray:
        """Evaluate electric potential at points with shape ``(n, 3)``.

        For ``cached_kneq0``, this returns the ``k!=0`` component only.
        """

        self._require_open()
        target_pos = _points_to_fortran_3xn(points, name="points")
        self._validate_cached_target_points(target_pos)
        ntarget = target_pos.shape[1]
        phi = np.zeros(ntarget, dtype=np.float64)
        if ntarget == 0:
            return phi
        status = self._lib.beach_kernel_eval_phi(
            self._handle,
            ctypes.c_int(ntarget),
            target_pos.ctypes.data_as(ctypes.c_void_p),
            phi.ctypes.data_as(ctypes.c_void_p),
        )
        _check_status(status, "beach_kernel_eval_phi")
        external_e0 = np.asarray(self._options.external_e0, dtype=np.float64)
        if np.any(external_e0 != 0.0):
            phi = phi - target_pos.T @ external_e0
        return phi

    def eval_e_direct(self, points: np.ndarray) -> np.ndarray:
        """Evaluate the exact non-periodic all-source electric field."""

        self._require_open()
        evaluator = self._direct_evaluator("beach_kernel_eval_e_direct")
        target_pos = _points_to_fortran_3xn(points, name="points")
        ntarget = target_pos.shape[1]
        e = np.zeros((3, ntarget), dtype=np.float64, order="F")
        if ntarget == 0:
            return np.empty((0, 3), dtype=np.float64)
        status = evaluator(
            self._handle,
            ctypes.c_int(ntarget),
            target_pos.ctypes.data_as(ctypes.c_void_p),
            e.ctypes.data_as(ctypes.c_void_p),
        )
        _check_status(status, "beach_kernel_eval_e_direct")
        return np.ascontiguousarray(e.T)

    def eval_phi_direct(self, points: np.ndarray) -> np.ndarray:
        """Evaluate the exact non-periodic all-source electric potential."""

        self._require_open()
        evaluator = self._direct_evaluator("beach_kernel_eval_phi_direct")
        target_pos = _points_to_fortran_3xn(points, name="points")
        ntarget = target_pos.shape[1]
        phi = np.zeros(ntarget, dtype=np.float64)
        if ntarget == 0:
            return phi
        status = evaluator(
            self._handle,
            ctypes.c_int(ntarget),
            target_pos.ctypes.data_as(ctypes.c_void_p),
            phi.ctypes.data_as(ctypes.c_void_p),
        )
        _check_status(status, "beach_kernel_eval_phi_direct")
        return phi

    def force_on_charges(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        *,
        origin: Iterable[float] = (0.0, 0.0, 0.0),
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return net force and torque on target charges in the current field.

        For ``cached_kneq0``, the current field is the ``k!=0`` component only.
        """

        self._require_open()
        target_pos = _points_to_fortran_3xn(positions, name="positions")
        self._validate_cached_target_points(target_pos)
        ntarget = target_pos.shape[1]
        target_q = _charges_1d(charges, expected=ntarget, name="charges")
        origin_arr = _vec3(origin, name="origin")
        force = np.zeros(3, dtype=np.float64)
        torque = np.zeros(3, dtype=np.float64)
        status = self._lib.beach_kernel_force_on_charges(
            self._handle,
            ctypes.c_int(ntarget),
            target_pos.ctypes.data_as(ctypes.c_void_p),
            target_q.ctypes.data_as(ctypes.c_void_p),
            origin_arr.ctypes.data_as(ctypes.c_void_p),
            force.ctypes.data_as(ctypes.c_void_p),
            torque.ctypes.data_as(ctypes.c_void_p),
        )
        _check_status(status, "beach_kernel_force_on_charges")
        external_e0 = np.asarray(self._options.external_e0, dtype=np.float64)
        if np.any(external_e0 != 0.0) and ntarget > 0:
            force_i = target_q[:, None] * external_e0[None, :]
            force += np.sum(force_i, axis=0)
            rel_pos = target_pos.T - origin_arr[None, :]
            torque += np.sum(np.cross(rel_pos, force_i), axis=0)
        return force, torque

    def diagnostics(self) -> FieldKernelDiagnostics:
        """Return periodic operator-cache diagnostics for this built plan."""

        self._require_open()
        getter = getattr(self._lib, "beach_kernel_get_periodic_cache_info", None)
        if getter is None:
            raise FieldKernelError(
                "field-kernel diagnostics require shared-library symbol "
                "beach_kernel_get_periodic_cache_info."
            )

        hit = ctypes.c_int()
        build_count = ctypes.c_int()
        fingerprint_length = ctypes.c_int()
        path_length = ctypes.c_int()
        def read_info(
            fingerprint_buffer: ctypes.Array[ctypes.c_char],
            path_buffer: ctypes.Array[ctypes.c_char],
        ) -> int:
            return int(
                getter(
                    self._handle,
                    ctypes.byref(hit),
                    ctypes.byref(build_count),
                    ctypes.cast(fingerprint_buffer, ctypes.c_void_p),
                    ctypes.c_int(len(fingerprint_buffer)),
                    ctypes.byref(fingerprint_length),
                    ctypes.cast(path_buffer, ctypes.c_void_p),
                    ctypes.c_int(len(path_buffer)),
                    ctypes.byref(path_length),
                )
            )

        fingerprint_buffer = ctypes.create_string_buffer(1)
        path_buffer = ctypes.create_string_buffer(1)
        probe_status = read_info(fingerprint_buffer, path_buffer)
        if probe_status not in (0, 2):
            _check_status(probe_status, "beach_kernel_get_periodic_cache_info")
        if fingerprint_length.value < 0 or path_length.value < 0:
            raise FieldKernelError("field-kernel diagnostics returned a negative text length.")

        fingerprint_buffer = ctypes.create_string_buffer(fingerprint_length.value + 1)
        path_buffer = ctypes.create_string_buffer(path_length.value + 1)
        status = read_info(fingerprint_buffer, path_buffer)
        _check_status(status, "beach_kernel_get_periodic_cache_info")
        if not 0 <= fingerprint_length.value < len(fingerprint_buffer):
            raise FieldKernelError(
                "field-kernel diagnostics returned an invalid fingerprint length."
            )
        if not 0 <= path_length.value < len(path_buffer):
            raise FieldKernelError(
                "field-kernel diagnostics returned an invalid path length."
            )
        if hit.value not in (0, 1) or build_count.value < 0:
            raise FieldKernelError(
                "field-kernel diagnostics returned invalid numeric metadata."
            )

        try:
            fingerprint_text = fingerprint_buffer.raw[
                : fingerprint_length.value
            ].decode("utf-8")
            path_text = path_buffer.raw[: path_length.value].decode("utf-8")
        except UnicodeDecodeError as exc:
            raise FieldKernelError(
                "field-kernel diagnostics returned invalid UTF-8 text."
            ) from exc
        has_cache_identity = bool(fingerprint_text or path_text)
        return FieldKernelDiagnostics(
            periodic_cache_hit=bool(hit.value) if has_cache_identity else None,
            periodic_operator_build_count=build_count.value,
            periodic_cache_fingerprint=fingerprint_text or None,
            periodic_cache_path=Path(path_text) if path_text else None,
        )

    def close(self) -> None:
        """Release the Fortran kernel handle."""

        if self._closed:
            return
        if getattr(self, "_handle", None) is not None and self._handle.value:
            status = self._lib.beach_kernel_destroy(self._handle)
            _check_status(status, "beach_kernel_destroy")
            self._handle = ctypes.c_void_p()
        self._closed = True

    def __enter__(self) -> "FieldKernel":
        self._require_open()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # type: ignore[no-untyped-def]
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _require_open(self) -> None:
        if self._closed or not self._handle.value:
            raise FieldKernelError("field kernel is closed.")

    def _validate_cached_target_points(self, target_pos: np.ndarray) -> None:
        axis = self._cached_target_free_axis
        bounds = self._cached_target_free_bounds
        if axis is None or bounds is None or target_pos.shape[1] == 0:
            return
        lower, upper = bounds
        coordinates = target_pos[axis]
        outside = (coordinates < lower) | (coordinates > upper)
        if not np.any(outside):
            return
        target_index = int(np.flatnonzero(outside)[0])
        axis_name = ("x", "y", "z")[axis]
        raise ValueError(
            "cached_kneq0 target points must lie inside the configured target "
            f"box on the non-periodic {axis_name} axis; target[{target_index}] "
            f"has {axis_name}={coordinates[target_index]:.17g}, expected "
            f"{lower:.17g} <= {axis_name} <= {upper:.17g}. Periodic-axis "
            "coordinates are wrapped automatically."
        )

    def _direct_evaluator(self, symbol: str):  # type: ignore[no-untyped-def]
        if self._options.periodic2 is not None:
            raise FieldKernelError(
                "exact-direct evaluation requires a non-periodic field-kernel plan."
            )
        evaluator = getattr(self._lib, symbol, None)
        if evaluator is None:
            raise FieldKernelError(
                f"exact-direct evaluation requires shared-library symbol {symbol}."
            )
        return evaluator


def calc_object_forces_kernel(
    result: FortranRunResult | object,
    *,
    step: int | None = -1,
    target_mesh_ids: int | Iterable[int] | None = None,
    periodic2: Mapping[str, object] | None = None,
    theta: float | None = None,
    leaf_max: int | None = None,
    order: int = 4,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
) -> tuple[KernelObjectForceRecord, ...]:
    """Compute object-wise net force using the Fortran FMM field kernel.

    For each target object, its own source charges are zeroed before evaluating
    the target panel lattice is integrated with seventh-order Gauss-Duffy
    quadrature. A configured uniform ``sim.e0`` is added to the target field,
    matching the simulator pusher semantics while avoiding central-cell
    primary self-force contamination.
    """

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    _require_triangle_source_model(resolved)
    triangles = _require_triangles(resolved)
    centers = _triangle_centers(triangles)
    charges = _charges_for_step(resolved, step=step)
    mesh_ids = _mesh_ids_or_default(resolved)
    available_ids = tuple(int(v) for v in np.unique(mesh_ids))
    if target_mesh_ids is None:
        target_ids = available_ids
    elif isinstance(target_mesh_ids, (int, np.integer)):
        target_ids = (int(target_mesh_ids),)
    else:
        target_ids = tuple(dict.fromkeys(int(v) for v in target_mesh_ids))
    missing = [mid for mid in target_ids if mid not in available_ids]
    if missing:
        raise ValueError(f"unknown mesh id(s): {missing}. available={list(available_ids)}")

    _require_total_field_config(
        context,
        operation="calc_object_forces_kernel",
    )
    options = _options_from_result(
        resolved,
        periodic2=periodic2,
        theta=theta,
        leaf_max=leaf_max,
        order=order,
        config_path=config_path,
        context=context,
    )
    _require_total_field_reconstruction(
        context,
        options,
        operation="calc_object_forces_kernel",
    )
    records: list[KernelObjectForceRecord] = []
    with FieldKernel(triangles, charges, options=options, library_path=library_path) as kernel:
        for mesh_id in target_ids:
            mask = mesh_ids == mesh_id
            if not np.any(mask):
                continue
            source_q = np.asarray(charges, dtype=np.float64).copy()
            source_q[mask] = 0.0
            kernel.update_charges(source_q)
            target_q = np.asarray(charges[mask], dtype=np.float64)
            target_triangles = triangles[mask]
            target_centers = centers[mask]
            target_points, target_weights, _ = panel_target_quadrature(
                target_triangles,
                target_q,
                7,
            )
            center = target_centers.mean(axis=0)
            force, torque = kernel.force_on_charges(
                target_points,
                target_weights,
                origin=center,
            )
            records.append(
                KernelObjectForceRecord(
                    mesh_id=mesh_id,
                    step=step,
                    total_charge_C=float(np.sum(target_q)),
                    center_m=center,
                    force_N=force,
                    torque_Nm=torque,
                )
            )
    return tuple(records)


def field_kernel_options_from_result(
    result: FortranRunResult | object,
    *,
    periodic2: Mapping[str, object] | None = None,
    theta: float | None = None,
    leaf_max: int | None = None,
    order: int = 4,
    config_path: str | Path | None = None,
) -> FieldKernelOptions:
    """Resolve field-kernel options from a BEACH result and optional config."""

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    _require_triangle_source_model(resolved)
    return _options_from_result(
        resolved,
        periodic2=periodic2,
        theta=theta,
        leaf_max=leaf_max,
        order=order,
        config_path=config_path,
        context=context,
    )


def _options_from_result(
    resolved: FortranRunResult,
    *,
    periodic2: Mapping[str, object] | None,
    theta: float | None,
    leaf_max: int | None,
    order: int,
    config_path: str | Path | None,
    context: RunContext | None = None,
) -> FieldKernelOptions:
    run_context = context or RunContext.from_value(resolved, config_path=config_path)
    sim = run_context.sim
    periodic_cfg = _coerce_periodic2(periodic2, allow_cached_kneq0=True)
    if periodic_cfg is None and sim is not None:
        allow_historical_root_oracle = run_context.requested_config_path is None
        periodic_cfg = _coerce_periodic2(
            _periodic2_from_sim(
                sim,
                allow_cached_kneq0=True,
                allow_historical_root_oracle=allow_historical_root_oracle,
            ),
            allow_cached_kneq0=True,
            allow_historical_root_oracle=allow_historical_root_oracle,
        )
    resolved_theta = float(theta if theta is not None else (sim or {}).get("tree_theta", 0.5))
    resolved_leaf_max = int(leaf_max if leaf_max is not None else (sim or {}).get("tree_leaf_max", 16))
    box_min: tuple[float, float, float] | None = None
    box_max: tuple[float, float, float] | None = None
    if sim is not None and "box_min" in sim and "box_max" in sim:
        box_min = tuple(float(v) for v in sim["box_min"])  # type: ignore[index]
        box_max = tuple(float(v) for v in sim["box_max"])  # type: ignore[index]
    return FieldKernelOptions(
        external_e0=_external_e0_from_sim(sim),
        theta=resolved_theta,
        leaf_max=resolved_leaf_max,
        order=int(order),
        periodic2=periodic_cfg,
        box_min=box_min,
        box_max=box_max,
        periodic_cache_dir=str(
            (sim or {}).get("field_periodic_cache_dir", ".beach_cache/periodic2")
        ),
        periodic_generation_tolerance=float(
            (sim or {}).get("field_periodic_generation_tolerance", 1.0e-8)
        ),
    )


def _require_total_field_reconstruction(
    context: RunContext,
    options: FieldKernelOptions,
    *,
    operation: str,
) -> None:
    """Reject component-only kernels where a public API promises total fields."""

    _require_total_field_config(context, operation=operation)
    _require_complete_total_kernel(options, operation=operation)


def _require_total_field_config(
    context: RunContext,
    *,
    operation: str,
) -> None:
    """Validate configuration before resolving a total-field kernel."""

    config = context.config
    if config is None:
        raise ValueError(
            f"{operation} requires the run's beach.toml to reconstruct the "
            "simulator total field; pass config_path or keep beach.toml near "
            "the output directory."
        )
    sim = context.sim
    if (
        context.requested_config_path is not None
        and sim is not None
        and str(sim.get("field_periodic_far_correction", "")).strip().lower()
        == "m2l_root_oracle"
    ):
        raise ValueError(
            'periodic2.far_correction "m2l_root_oracle" was removed; use '
            '"none" for total-field post-processing.'
        )
    if _config_has_active_outer_field(config):
        raise ValueError(
            f"{operation} cannot reconstruct the simulator total field while "
            "an outer field model is active. Use saved simulator diagnostics "
            "or a post-processing API that explicitly composes the outer state."
        )


def _require_complete_total_kernel(
    options: FieldKernelOptions,
    *,
    operation: str,
) -> None:
    """Reject a native kernel option set that represents only one component."""

    periodic_cfg = _coerce_periodic2(
        options.periodic2,
        allow_cached_kneq0=True,
    )
    if periodic_cfg is not None and periodic_cfg.far_correction == "cached_kneq0":
        raise ValueError(
            f"{operation} cannot treat cached_kneq0 as the simulator total "
            "field: the native kernel returns only the k!=0 component and does "
            "not compose the physical periodic zero mode. Use saved "
            "mesh_potential.csv where applicable or ObjectInteractionSnapshot "
            "for supported force workflows."
        )


def _config_has_active_outer_field(config: Mapping[str, object]) -> bool:
    legacy = config.get("outer_plasma")
    if isinstance(legacy, Mapping) and _field_model_is_active(legacy.get("model")):
        return True

    external = config.get("external_boundary")
    if not isinstance(external, Mapping):
        return False
    field = external.get("field")
    return isinstance(field, Mapping) and _field_model_is_active(field.get("model"))


def _field_model_is_active(value: object) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() not in {"", "none", "not_applicable"}


def _load_sim_config(
    output_dir: Path,
    *,
    config_path: str | Path | None,
) -> Mapping[str, object] | None:
    config = _load_full_config(output_dir, config_path=config_path)
    if config is None:
        return None
    sim = config.get("sim")
    return sim if isinstance(sim, Mapping) else None


def _load_full_config(
    output_dir: Path,
    *,
    config_path: str | Path | None,
) -> Mapping[str, object] | None:
    return load_config_for_output(output_dir, config_path=config_path)


def _load_sim_config_near_output(output_dir: Path) -> Mapping[str, object] | None:
    return _load_sim_config(output_dir, config_path=None)


def _external_e0_from_sim(sim: Mapping[str, object] | None) -> tuple[float, float, float]:
    if sim is None:
        return (0.0, 0.0, 0.0)
    has_vector = "e0" in sim
    has_abs = "e0_abs" in sim
    has_phi_xy = "e0_phi_xy_deg" in sim
    has_phi_z = "e0_phi_z_deg" in sim
    if has_vector and (has_abs or has_phi_xy or has_phi_z):
        raise ValueError(
            "sim.e0 cannot be combined with sim.e0_abs/e0_phi_xy_deg/e0_phi_z_deg."
        )
    if has_vector:
        e0 = _vec3(sim["e0"], name="e0")
        return tuple(float(v) for v in e0)
    if (has_phi_xy or has_phi_z) and not has_abs:
        raise ValueError("sim.e0_phi_xy_deg/e0_phi_z_deg require sim.e0_abs.")
    if not has_abs:
        return (0.0, 0.0, 0.0)

    e0_abs = float(sim.get("e0_abs", 0.0))
    phi_xy = np.deg2rad(float(sim.get("e0_phi_xy_deg", 0.0)))
    phi_z = np.deg2rad(float(sim.get("e0_phi_z_deg", 0.0)))
    if not np.isfinite(e0_abs) or e0_abs < 0.0:
        raise ValueError("sim.e0_abs must be finite and >= 0.")
    if not np.isfinite(phi_xy) or not np.isfinite(phi_z):
        raise ValueError("sim.e0_phi_xy_deg/e0_phi_z_deg must be finite.")
    return (
        float(e0_abs * np.cos(phi_z) * np.cos(phi_xy)),
        float(e0_abs * np.cos(phi_z) * np.sin(phi_xy)),
        float(e0_abs * np.sin(phi_z)),
    )


def _periodic_box_vectors(
    *,
    axes: tuple[int, int],
    lengths: tuple[float, float],
    origins: tuple[float, float],
    box_min: tuple[float, float, float] | None,
    box_max: tuple[float, float, float] | None,
    source_vertices_3xn: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    if box_min is None:
        all_vertices = np.concatenate(source_vertices_3xn, axis=1)
        mins = np.min(all_vertices, axis=1)
        maxs = np.max(all_vertices, axis=1)
        span = np.maximum(maxs - mins, 1.0)
        box_min_arr = mins - 0.5 * span
        box_max_arr = maxs + 0.5 * span
    else:
        box_min_arr = np.asarray(box_min, dtype=np.float64)
        if box_max is None:
            raise ValueError("box_max is required when box_min is set.")
        box_max_arr = np.asarray(box_max, dtype=np.float64)

    box_min_arr = box_min_arr.copy()
    box_max_arr = box_max_arr.copy()
    for axis, origin, length in zip(axes, origins, lengths):
        box_min_arr[axis] = origin
        box_max_arr[axis] = origin + length
    if np.any(box_max_arr <= box_min_arr):
        raise ValueError("periodic2 target box must satisfy box_max > box_min.")
    return box_min_arr, box_max_arr


def _load_kernel_library(library_path: str | Path | None) -> ctypes.CDLL:
    errors: list[str] = []
    for path in _candidate_library_paths(library_path):
        try:
            return ctypes.CDLL(str(path))
        except OSError as exc:
            errors.append(f"{path}: {exc}")
    hints = [
        "Set BEACH_FIELD_KERNEL_LIB to the shared library path,",
        "or build it with `make build-kernel`.",
    ]
    detail = "\n".join(errors[-4:])
    raise FieldKernelError("BEACH field kernel library is not available. " + " ".join(hints) + ("\n" + detail if detail else ""))


def _candidate_library_paths(library_path: str | Path | None) -> list[Path]:
    if library_path is not None:
        return [Path(library_path)]
    names = _library_names()
    paths: list[Path] = []
    env_path = os.environ.get("BEACH_FIELD_KERNEL_LIB")
    if env_path:
        paths.append(Path(env_path))
    here = Path(__file__).resolve()
    package_root = here.parents[1]
    repo_root = here.parents[2]
    for name in names:
        paths.extend(
            [
                package_root / "lib" / name,
                repo_root / "build" / name,
                repo_root / ".local" / "lib" / name,
            ]
        )
    return paths


def _library_names() -> tuple[str, ...]:
    if sys.platform == "darwin":
        return ("libbeach_field_kernel.dylib", "libbeach_field_kernel.so")
    if sys.platform.startswith("win"):
        return ("beach_field_kernel.dll", "libbeach_field_kernel.dll")
    return ("libbeach_field_kernel.so",)


def _configure_library(lib: ctypes.CDLL) -> None:
    if getattr(lib, "_beach_kernel_ctypes_configured", False):
        return
    c_void_p = ctypes.c_void_p
    c_int = ctypes.c_int
    c_double = ctypes.c_double

    abi_getter = getattr(lib, "beach_kernel_get_abi_version", None)
    if abi_getter is not None:
        abi_getter.argtypes = [c_void_p, c_void_p]
        abi_getter.restype = c_int

    lib.beach_kernel_create.argtypes = [ctypes.POINTER(c_void_p)]
    lib.beach_kernel_create.restype = c_int
    lib.beach_kernel_destroy.argtypes = [c_void_p]
    lib.beach_kernel_destroy.restype = c_int
    lib.beach_kernel_build.argtypes = [
        c_void_p,
        c_int,
        c_void_p,
        c_void_p,
        c_void_p,
        c_double,
        c_int,
        c_int,
        c_int,
        c_void_p,
        c_void_p,
        c_int,
        c_int,
        c_double,
        c_int,
        c_void_p,
        c_void_p,
    ]
    lib.beach_kernel_build.restype = c_int
    lib.beach_kernel_update_charges.argtypes = [c_void_p, c_int, c_void_p]
    lib.beach_kernel_update_charges.restype = c_int
    lib.beach_kernel_eval_e.argtypes = [c_void_p, c_int, c_void_p, c_void_p]
    lib.beach_kernel_eval_e.restype = c_int
    lib.beach_kernel_eval_phi.argtypes = [c_void_p, c_int, c_void_p, c_void_p]
    lib.beach_kernel_eval_phi.restype = c_int
    direct_e = getattr(lib, "beach_kernel_eval_e_direct", None)
    if direct_e is not None:
        direct_e.argtypes = [c_void_p, c_int, c_void_p, c_void_p]
        direct_e.restype = c_int
    direct_phi = getattr(lib, "beach_kernel_eval_phi_direct", None)
    if direct_phi is not None:
        direct_phi.argtypes = [c_void_p, c_int, c_void_p, c_void_p]
        direct_phi.restype = c_int
    lib.beach_kernel_force_on_charges.argtypes = [
        c_void_p,
        c_int,
        c_void_p,
        c_void_p,
        c_void_p,
        c_void_p,
        c_void_p,
    ]
    lib.beach_kernel_force_on_charges.restype = c_int
    setter = getattr(lib, "beach_kernel_set_periodic_cache", None)
    if setter is not None:
        setter.argtypes = [c_void_p, c_void_p, c_int, c_double]
        setter.restype = c_int
    getter = getattr(lib, "beach_kernel_get_periodic_cache_info", None)
    if getter is not None:
        getter.argtypes = [
            c_void_p,
            c_void_p,
            c_void_p,
            c_void_p,
            c_int,
            c_void_p,
            c_void_p,
            c_int,
            c_void_p,
        ]
        getter.restype = c_int
    build_info_getter = getattr(lib, "beach_kernel_get_build_info", None)
    if build_info_getter is not None:
        build_info_getter.argtypes = [c_void_p, c_int, c_void_p]
        build_info_getter.restype = c_int
    _validate_library_abi(lib)
    lib._beach_kernel_ctypes_configured = True


def _validate_library_abi(lib: ctypes.CDLL) -> None:
    getter = getattr(lib, "beach_kernel_get_abi_version", None)
    if getter is None:
        raise FieldKernelError(
            "field-kernel ABI version attestation is required; "
            "the loaded library predates the triangle-only ABI v2 contract."
        )

    major = ctypes.c_int()
    minor = ctypes.c_int()
    status = getter(ctypes.byref(major), ctypes.byref(minor))
    _check_status(status, "beach_kernel_get_abi_version")
    if major.value != FIELD_KERNEL_ABI_MAJOR or minor.value < FIELD_KERNEL_ABI_MINOR:
        raise FieldKernelError(
            "field-kernel ABI is incompatible: "
            f"library={major.value}.{minor.value}, "
            f"required={FIELD_KERNEL_ABI_MAJOR}.{FIELD_KERNEL_ABI_MINOR}."
        )


def _check_status(status: int, operation: str) -> None:
    if int(status) == 0:
        return
    message = _STATUS_MESSAGES.get(int(status), f"unknown status {status}")
    raise FieldKernelError(f"{operation} failed: {message}.")


def _points_to_fortran_3xn(points: np.ndarray, *, name: str) -> np.ndarray:
    arr = np.asarray(points, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n_points, 3).")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return np.asfortranarray(arr.T)


def _triangles_to_fortran_vertices(
    triangles: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = np.asarray(triangles, dtype=np.float64)
    if arr.ndim != 3 or arr.shape[1:] != (3, 3):
        raise ValueError("source_triangles must have shape (n_triangles, 3, 3).")
    if arr.shape[0] <= 0:
        raise ValueError("source_triangles must contain at least one triangle.")
    if not np.all(np.isfinite(arr)):
        raise ValueError("source_triangles must contain finite values.")
    edge1 = arr[:, 1] - arr[:, 0]
    edge2 = arr[:, 2] - arr[:, 0]
    double_area = np.linalg.norm(np.cross(edge1, edge2), axis=1)
    scale = np.maximum(1.0, np.max(np.abs(arr), axis=(1, 2)))
    if np.any(double_area <= 64.0 * np.finfo(np.float64).eps * scale * scale):
        raise ValueError("source_triangles must contain non-degenerate triangles.")
    return tuple(np.asfortranarray(arr[:, vertex, :].T) for vertex in range(3))  # type: ignore[return-value]


def _charges_1d(charges: np.ndarray, *, expected: int, name: str) -> np.ndarray:
    arr = np.asarray(charges, dtype=np.float64)
    if arr.ndim != 1 or arr.shape[0] != expected:
        raise ValueError(f"{name} must have shape ({expected},).")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return np.ascontiguousarray(arr)


def _vec3(value: Iterable[float], *, name: str) -> np.ndarray:
    arr = np.asarray(list(value), dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"{name} must contain exactly 3 values.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return np.ascontiguousarray(arr)


def _far_correction_code(value: str) -> int:
    key = str(value).strip().lower()
    if key == "m2l_root_oracle":
        raise ValueError(
            'periodic far correction "m2l_root_oracle" was removed; '
            'use "cached_kneq0".'
        )
    if key not in _FAR_CORRECTION_CODES:
        raise ValueError(
            'periodic far correction must be "auto", "none", or "cached_kneq0".'
        )
    return _FAR_CORRECTION_CODES[key]


def _null_ptr() -> ctypes.c_void_p:
    return ctypes.c_void_p()
