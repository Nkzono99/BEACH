"""ctypes bindings for the BEACH Fortran field-kernel shared library."""

from __future__ import annotations

import ctypes
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .mesh import _triangle_centers
from .potential import (
    _auto_periodic2_from_result,
    _coerce_periodic2,
    _find_config_path_near_output,
    _load_toml,
    _periodic2_from_sim,
    _resolve_softening,
)
from .selection import _charges_for_step, _mesh_ids_or_default, _require_triangles, _resolve_result
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
    "m2l_root_oracle": 2,
}


class FieldKernelError(RuntimeError):
    """Raised when the shared field kernel cannot be used."""


@dataclass(frozen=True)
class FieldKernelOptions:
    """Options passed to the Fortran Coulomb FMM kernel."""

    softening: float = 0.0
    theta: float = 0.5
    leaf_max: int = 16
    order: int = 4
    periodic2: tuple[
        tuple[int, int],
        tuple[float, float],
        tuple[float, float],
        int,
        str,
        float,
        int,
    ] | None = None
    box_min: tuple[float, float, float] | None = None
    box_max: tuple[float, float, float] | None = None


@dataclass(frozen=True)
class KernelObjectForceRecord:
    """Object-level net force/torque computed by the Fortran field kernel."""

    mesh_id: int
    step: int | None
    total_charge_C: float
    center_m: np.ndarray
    force_N: np.ndarray
    torque_Nm: np.ndarray


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
    ) -> None:
        self._lib = _load_kernel_library(library_path)
        _configure_library(self._lib)
        self._handle = ctypes.c_void_p()
        self._closed = False
        self._source_count = 0
        status = self._lib.beach_kernel_create(ctypes.byref(self._handle))
        _check_status(status, "beach_kernel_create")
        try:
            self.build(source_positions, options=options)
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

        resolved = _resolve_result(result)
        centers = _triangle_centers(_require_triangles(resolved))
        charges = _charges_for_step(resolved, step=step)
        options = _options_from_result(
            resolved,
            softening=softening,
            periodic2=periodic2,
            theta=theta,
            leaf_max=leaf_max,
            order=order,
            config_path=config_path,
        )
        return cls(centers, charges, options=options, library_path=library_path)

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
        nsrc = src_pos.shape[1]
        if nsrc <= 0:
            raise ValueError("source_positions must contain at least one point.")

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

        if opts.periodic2 is not None:
            axes, lengths, origins, image_layers, far_key, ewald_alpha, ewald_layers = opts.periodic2
            box_min_vec, box_max_vec = _periodic_box_vectors(
                axes=axes,
                lengths=lengths,
                origins=origins,
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

        status = self._lib.beach_kernel_build(
            self._handle,
            ctypes.c_int(nsrc),
            src_pos.ctypes.data_as(ctypes.c_void_p),
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
        _check_status(status, "beach_kernel_build")
        self._source_count = nsrc
        self._keepalive = keepalive

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
        return np.ascontiguousarray(e.T)

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
        return force, torque

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
    ``sum(q_i E_not_self(r_i))``. This avoids self-force contamination while
    preserving the simulator's FMM/periodic2 far-correction semantics.
    """

    resolved = _resolve_result(result)
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

    resolved = _resolve_result(result)
    return _options_from_result(
        resolved,
        softening=softening,
        periodic2=periodic2,
        theta=theta,
        leaf_max=leaf_max,
        order=order,
        config_path=config_path,
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
) -> FieldKernelOptions:
    sim = _load_sim_config(resolved.directory, config_path=config_path)
    resolved_softening = _resolve_kernel_softening(resolved, sim=sim, softening=softening)
    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        if sim is None and config_path is None:
            periodic_cfg = _auto_periodic2_from_result(resolved)
        elif sim is not None:
            periodic_cfg = _coerce_periodic2(_periodic2_from_sim(sim))
    resolved_theta = float(theta if theta is not None else (sim or {}).get("tree_theta", 0.5))
    resolved_leaf_max = int(leaf_max if leaf_max is not None else (sim or {}).get("tree_leaf_max", 16))
    box_min: tuple[float, float, float] | None = None
    box_max: tuple[float, float, float] | None = None
    if sim is not None and "box_min" in sim and "box_max" in sim:
        box_min = tuple(float(v) for v in sim["box_min"])  # type: ignore[index]
        box_max = tuple(float(v) for v in sim["box_max"])  # type: ignore[index]
    return FieldKernelOptions(
        softening=resolved_softening,
        theta=resolved_theta,
        leaf_max=resolved_leaf_max,
        order=int(order),
        periodic2=periodic_cfg,
        box_min=box_min,
        box_max=box_max,
    )


def _load_sim_config(
    output_dir: Path,
    *,
    config_path: str | Path | None,
) -> Mapping[str, object] | None:
    if config_path is None:
        path = _find_config_path_near_output(output_dir)
    else:
        path = Path(config_path)
        if not path.exists():
            raise ValueError(f'config file is not found: "{path}".')
    if path is None:
        return None
    config = _load_toml(path)
    sim = config.get("sim")
    return sim if isinstance(sim, Mapping) else None


def _load_sim_config_near_output(output_dir: Path) -> Mapping[str, object] | None:
    return _load_sim_config(output_dir, config_path=None)


def _resolve_kernel_softening(
    resolved: FortranRunResult,
    *,
    sim: Mapping[str, object] | None,
    softening: float | None,
) -> float:
    if softening is not None or sim is None:
        return _resolve_softening(resolved, softening)
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
    c_void_p = ctypes.c_void_p
    c_int = ctypes.c_int
    c_double = ctypes.c_double

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
    lib.beach_kernel_update_charges.argtypes = [c_void_p, c_int, c_void_p]
    lib.beach_kernel_update_charges.restype = c_int
    lib.beach_kernel_eval_e.argtypes = [c_void_p, c_int, c_void_p, c_void_p]
    lib.beach_kernel_eval_e.restype = c_int
    lib.beach_kernel_eval_phi.argtypes = [c_void_p, c_int, c_void_p, c_void_p]
    lib.beach_kernel_eval_phi.restype = c_int
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
    if key not in _FAR_CORRECTION_CODES:
        raise ValueError('periodic far correction must be "auto", "none", or "m2l_root_oracle".')
    return _FAR_CORRECTION_CODES[key]


def _null_ptr() -> ctypes.c_void_p:
    return ctypes.c_void_p()
