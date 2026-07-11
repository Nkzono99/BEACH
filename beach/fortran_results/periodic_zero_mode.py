"""Thin binding for BEACH's native physical periodic zero mode."""

from __future__ import annotations

import ctypes
from pathlib import Path

import numpy as np

from .kernel import FieldKernelError, _load_kernel_library


_STATUS_MESSAGES = {
    0: "ok",
    1: "invalid zero-mode handle",
    2: "invalid zero-mode argument",
    3: "periodic zero mode is not ready",
}

_TRACE_CODES = {
    "minus": -1,
    "principal_value": 0,
    "plus": 1,
}


class PeriodicZeroMode:
    """Evaluate the simulator's physical x/y-periodic zero mode.

    ``source_heights_m`` stores the three vertex heights for each source. Point
    sources use the same height in all three columns.
    """

    def __init__(
        self,
        source_heights_m: np.ndarray,
        source_charges_C: np.ndarray,
        area_xy_m2: float,
        *,
        e_bottom_V_m: float = 0.0,
        z_gauge_m: float | None = None,
        phi_gauge_V: float = 0.0,
        library_path: str | Path | None = None,
    ) -> None:
        heights = np.asarray(source_heights_m, dtype=np.float64)
        if heights.ndim != 2 or heights.shape[1] != 3 or heights.shape[0] == 0:
            raise ValueError("source_heights_m must have shape (n_sources, 3).")
        if not np.all(np.isfinite(heights)):
            raise ValueError("source_heights_m must contain finite values.")
        charges = _charges(source_charges_C, heights.shape[0])
        area = float(area_xy_m2)
        if not np.isfinite(area) or area <= 0.0:
            raise ValueError("area_xy_m2 must be finite and positive.")
        e_bottom = _finite_scalar(e_bottom_V_m, "e_bottom_V_m")
        phi_gauge = _finite_scalar(phi_gauge_V, "phi_gauge_V")
        if z_gauge_m is None:
            z_gauge = float(np.min(heights))
        else:
            z_gauge = _finite_scalar(z_gauge_m, "z_gauge_m")

        self._lib = _load_kernel_library(library_path)
        _configure_zero_mode_library(self._lib)
        self._handle = ctypes.c_void_p()
        self._closed = False
        self._source_count = heights.shape[0]
        self._e_bottom_V_m = e_bottom
        self._z_gauge_m = z_gauge
        self._phi_gauge_V = phi_gauge

        status = self._lib.beach_zero_mode_create(ctypes.byref(self._handle))
        _check_status(status, "beach_zero_mode_create")
        try:
            native_heights = np.asfortranarray(heights.T)
            status = self._lib.beach_zero_mode_build(
                self._handle,
                ctypes.c_int(self._source_count),
                ctypes.c_void_p(native_heights.ctypes.data),
                ctypes.c_double(area),
            )
            _check_status(status, "beach_zero_mode_build")
            self.update_charges(charges)
        except Exception:
            self.close()
            raise

    def update_charges(self, source_charges_C: np.ndarray) -> None:
        """Refresh source charges without rebuilding the height plan."""

        self._require_open()
        charges = _charges(source_charges_C, self._source_count)
        status = self._lib.beach_zero_mode_update(
            self._handle,
            ctypes.c_int(self._source_count),
            ctypes.c_void_p(charges.ctypes.data),
            ctypes.c_double(self._e_bottom_V_m),
            ctypes.c_double(self._z_gauge_m),
            ctypes.c_double(self._phi_gauge_V),
        )
        _check_status(status, "beach_zero_mode_update")

    def eval(
        self,
        z_m: np.ndarray | float,
        trace: str = "principal_value",
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return native potential and vertical field at the requested heights."""

        self._require_open()
        if trace not in _TRACE_CODES:
            allowed = ", ".join(sorted(_TRACE_CODES))
            raise ValueError(f"trace must be one of: {allowed}.")
        z = np.asarray(z_m, dtype=np.float64)
        if not np.all(np.isfinite(z)):
            raise ValueError("z_m must contain finite values.")
        shape = z.shape
        native_z = np.ascontiguousarray(z.reshape(-1))
        phi = np.empty(native_z.size, dtype=np.float64)
        ez = np.empty(native_z.size, dtype=np.float64)
        if native_z.size:
            status = self._lib.beach_zero_mode_eval(
                self._handle,
                ctypes.c_int(native_z.size),
                ctypes.c_void_p(native_z.ctypes.data),
                ctypes.c_int(_TRACE_CODES[trace]),
                ctypes.c_void_p(phi.ctypes.data),
                ctypes.c_void_p(ez.ctypes.data),
            )
            _check_status(status, "beach_zero_mode_eval")
        phi = _readonly(phi.reshape(shape))
        ez = _readonly(ez.reshape(shape))
        return phi, ez

    def close(self) -> None:
        """Release the native handle. Repeated calls are harmless."""

        if self._closed:
            return
        self._closed = True
        if self._handle.value is not None:
            status = self._lib.beach_zero_mode_destroy(self._handle)
            self._handle = ctypes.c_void_p()
            _check_status(status, "beach_zero_mode_destroy")

    def __enter__(self) -> "PeriodicZeroMode":
        self._require_open()
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _require_open(self) -> None:
        if self._closed:
            raise FieldKernelError("periodic zero-mode handle is closed.")


def _configure_zero_mode_library(lib: ctypes.CDLL) -> None:
    if getattr(lib, "_beach_zero_mode_ctypes_configured", False):
        return
    names = (
        "beach_zero_mode_create",
        "beach_zero_mode_destroy",
        "beach_zero_mode_build",
        "beach_zero_mode_update",
        "beach_zero_mode_eval",
    )
    try:
        create, destroy, build, update, evaluate = (getattr(lib, name) for name in names)
    except AttributeError as exc:
        raise FieldKernelError(
            "periodic zero-mode symbols are unavailable; rebuild with `make build-kernel`."
        ) from exc

    c_void_p = ctypes.c_void_p
    c_int = ctypes.c_int
    c_double = ctypes.c_double
    create.argtypes = [ctypes.POINTER(c_void_p)]
    create.restype = c_int
    destroy.argtypes = [c_void_p]
    destroy.restype = c_int
    build.argtypes = [c_void_p, c_int, c_void_p, c_double]
    build.restype = c_int
    update.argtypes = [c_void_p, c_int, c_void_p, c_double, c_double, c_double]
    update.restype = c_int
    evaluate.argtypes = [c_void_p, c_int, c_void_p, c_int, c_void_p, c_void_p]
    evaluate.restype = c_int
    lib._beach_zero_mode_ctypes_configured = True


def _check_status(status: int, operation: str) -> None:
    value = int(status)
    if value == 0:
        return
    message = _STATUS_MESSAGES.get(value, f"unknown status {value}")
    raise FieldKernelError(f"{operation} failed: {message}.")


def _charges(value: np.ndarray, expected: int) -> np.ndarray:
    charges = np.asarray(value, dtype=np.float64)
    if charges.shape != (expected,):
        raise ValueError(f"source_charges_C must have shape ({expected},).")
    if not np.all(np.isfinite(charges)):
        raise ValueError("source_charges_C must contain finite values.")
    return np.ascontiguousarray(charges)


def _finite_scalar(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _readonly(value: np.ndarray) -> np.ndarray:
    result = np.array(value, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result
