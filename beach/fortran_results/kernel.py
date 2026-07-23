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
from .periodic import Periodic2Input
from .potential import (
    _coerce_periodic2,
    _periodic2_from_sim,
    _resolve_softening,
)
from .selection import (
    _charges_for_step,
    _mesh_ids_or_default,
    _require_point_source_model,
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

FIELD_KERNEL_ABI_MAJOR = 1
FIELD_KERNEL_ABI_MINOR = 0


class FieldKernelError(RuntimeError):
    """Raised when the shared field kernel cannot be used."""


@dataclass(frozen=True)
class FieldKernelOptions:
    """Options passed to the Fortran Coulomb FMM kernel."""

    softening: float = 0.0
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
    :meth:`update_charges`.
    """

    def __init__(
        self,
        source_positions: np.ndarray,
        source_charges: np.ndarray,
        *,
        options: FieldKernelOptions | None = None,
        library_path: str | Path | None = None,
        source_triangles: np.ndarray | None = None,
    ) -> None:
        self._lib = _load_kernel_library(library_path)
        _configure_library(self._lib)
        self._handle = ctypes.c_void_p()
        self._closed = False
        self._source_count = 0
        self._options = FieldKernelOptions()
        status = self._lib.beach_kernel_create(ctypes.byref(self._handle))
        _check_status(status, "beach_kernel_create")
        try:
            if source_triangles is None:
                self.build(source_positions, options=options)
            else:
                self.build_panel(source_triangles, options=options)
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
        softening: float | None = None,
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
        triangles = _require_triangles(resolved)
        centers = _triangle_centers(triangles)
        charges = _charges_for_step(resolved, step=step)
        options = _options_from_result(
            resolved,
            softening=softening,
            periodic2=periodic2,
            theta=theta,
            leaf_max=leaf_max,
            order=order,
            config_path=config_path,
            context=context,
        )
        source_triangles = triangles if resolved.field_source_model == "triangle_p0" else None
        return cls(
            centers,
            charges,
            options=options,
            library_path=library_path,
            source_triangles=source_triangles,
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
        source_positions: np.ndarray,
        *,
        options: FieldKernelOptions | None = None,
    ) -> None:
        """Build or rebuild the source-geometry plan."""

        self._require_open()
        opts = options or FieldKernelOptions()
        src_pos = _points_to_fortran_3xn(source_positions, name="source_positions")
        self._build_geometry(src_pos, None, opts)

    def build_panel(
        self,
        source_triangles: np.ndarray,
        *,
        options: FieldKernelOptions | None = None,
    ) -> None:
        """Build or rebuild a constant-density triangle-panel plan."""

        self._require_open()
        opts = options or FieldKernelOptions()
        if opts.softening != 0.0:
            raise ValueError("triangle-panel field kernels require softening=0.")
        vertices = _triangles_to_fortran_vertices(source_triangles)
        src_pos = (vertices[0] + vertices[1] + vertices[2]) / 3.0
        self._build_geometry(src_pos, vertices, opts)

    def _build_geometry(
        self,
        src_pos: np.ndarray,
        panel_vertices: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
        opts: FieldKernelOptions,
    ) -> None:
        nsrc = src_pos.shape[1]
        if nsrc <= 0:
            raise ValueError("source_positions must contain at least one point.")

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
        keepalive: list[np.ndarray] = [src_pos]
        if panel_vertices is not None:
            keepalive.extend(panel_vertices)

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
                source_positions_3xn=src_pos,
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

        self._set_periodic_cache_options(opts, far_key=far_key)

        geometry_args: list[object]
        if panel_vertices is None:
            build_function = self._lib.beach_kernel_build
            geometry_args = [src_pos.ctypes.data_as(ctypes.c_void_p)]
            operation = "beach_kernel_build"
        else:
            build_function = self._lib.beach_kernel_build_panel
            geometry_args = [vertex.ctypes.data_as(ctypes.c_void_p) for vertex in panel_vertices]
            operation = "beach_kernel_build_panel"

        status = build_function(
            self._handle,
            ctypes.c_int(nsrc),
            *geometry_args,
            ctypes.c_double(opts.theta),
            ctypes.c_int(opts.leaf_max),
            ctypes.c_int(opts.order),
            ctypes.c_double(opts.softening),
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
        _check_status(status, operation)
        self._source_count = nsrc
        self._options = opts
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
        """Evaluate electric field vectors at points with shape ``(n, 3)``."""

        self._require_open()
        target_pos = _points_to_fortran_3xn(points, name="points")
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
        """Evaluate electric potential at points with shape ``(n, 3)``."""

        self._require_open()
        target_pos = _points_to_fortran_3xn(points, name="points")
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
        """Return net force and torque on target charges in the current field."""

        self._require_open()
        target_pos = _points_to_fortran_3xn(positions, name="positions")
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
    softening: float | None = None,
    periodic2: Mapping[str, object] | None = None,
    theta: float | None = None,
    leaf_max: int | None = None,
    order: int = 4,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
) -> tuple[KernelObjectForceRecord, ...]:
    """Compute object-wise net force using the Fortran FMM field kernel.

    For each target object, its own source charges are zeroed before evaluating
    ``sum(q_i E_not_self(r_i))``. A configured uniform ``sim.e0`` is added to
    the target field, matching the simulator pusher semantics while avoiding
    Coulomb self-force contamination.
    """

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    _require_point_source_model(resolved)
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

    options = _options_from_result(
        resolved,
        softening=softening,
        periodic2=periodic2,
        theta=theta,
        leaf_max=leaf_max,
        order=order,
        config_path=config_path,
        context=context,
    )
    records: list[KernelObjectForceRecord] = []
    with FieldKernel(centers, charges, options=options, library_path=library_path) as kernel:
        for mesh_id in target_ids:
            mask = mesh_ids == mesh_id
            if not np.any(mask):
                continue
            source_q = np.asarray(charges, dtype=np.float64).copy()
            source_q[mask] = 0.0
            kernel.update_charges(source_q)
            target_centers = centers[mask]
            target_q = np.asarray(charges[mask], dtype=np.float64)
            center = target_centers.mean(axis=0)
            force, torque = kernel.force_on_charges(target_centers, target_q, origin=center)
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
    softening: float | None = None,
    periodic2: Mapping[str, object] | None = None,
    theta: float | None = None,
    leaf_max: int | None = None,
    order: int = 4,
    config_path: str | Path | None = None,
) -> FieldKernelOptions:
    """Resolve field-kernel options from a BEACH result and optional config."""

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    _require_point_source_model(resolved)
    return _options_from_result(
        resolved,
        softening=softening,
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
    softening: float | None,
    periodic2: Mapping[str, object] | None,
    theta: float | None,
    leaf_max: int | None,
    order: int,
    config_path: str | Path | None,
    context: RunContext | None = None,
) -> FieldKernelOptions:
    run_context = context or RunContext.from_value(resolved, config_path=config_path)
    sim = run_context.sim
    resolved_softening = _resolve_kernel_softening(
        resolved,
        sim=sim,
        softening=softening,
        context=run_context,
    )
    periodic_cfg = _coerce_periodic2(periodic2, allow_cached_kneq0=True)
    if periodic_cfg is None and sim is not None:
        periodic_cfg = _coerce_periodic2(
            _periodic2_from_sim(sim, allow_cached_kneq0=True),
            allow_cached_kneq0=True,
        )
    resolved_theta = float(theta if theta is not None else (sim or {}).get("tree_theta", 0.5))
    resolved_leaf_max = int(leaf_max if leaf_max is not None else (sim or {}).get("tree_leaf_max", 16))
    box_min: tuple[float, float, float] | None = None
    box_max: tuple[float, float, float] | None = None
    if sim is not None and "box_min" in sim and "box_max" in sim:
        box_min = tuple(float(v) for v in sim["box_min"])  # type: ignore[index]
        box_max = tuple(float(v) for v in sim["box_max"])  # type: ignore[index]
    return FieldKernelOptions(
        softening=resolved_softening,
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


def _resolve_kernel_softening(
    resolved: FortranRunResult,
    *,
    sim: Mapping[str, object] | None,
    softening: float | None,
    context: RunContext | None = None,
) -> float:
    if softening is not None or sim is None:
        return _resolve_softening(resolved, softening, context=context)
    value = float(sim.get("softening", 0.0))
    if not np.isfinite(value) or value < 0.0:
        raise ValueError("softening must be finite and >= 0.")
    return value


def _periodic_box_vectors(
    *,
    axes: tuple[int, int],
    lengths: tuple[float, float],
    origins: tuple[float, float],
    box_min: tuple[float, float, float] | None,
    box_max: tuple[float, float, float] | None,
    source_positions_3xn: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if box_min is None:
        mins = np.min(source_positions_3xn, axis=1)
        maxs = np.max(source_positions_3xn, axis=1)
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
        c_double,
        c_int,
        c_int,
        c_double,
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
    lib.beach_kernel_build_panel.argtypes = [
        c_void_p,
        c_int,
        c_void_p,
        c_void_p,
        c_void_p,
        c_double,
        c_int,
        c_int,
        c_double,
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
    lib.beach_kernel_build_panel.restype = c_int
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
        return

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
