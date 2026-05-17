"""What-if scene editing for BEACH output meshes."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Literal, Mapping

import numpy as np

from .kernel import (
    FieldKernel,
    field_kernel_options_from_result,
    KernelObjectForceRecord,
)
from .mesh import _triangle_centers
from .selection import (
    _charges_for_step,
    _mesh_ids_or_default,
    _require_triangles,
    _resolve_result,
)
from .types import FortranRunResult


TransformBackend = Literal["numpy", "numba", "auto"]
RotationOrigin = Literal["object_center", "origin"]

_NUMBA_TRANSFORM_KERNEL = None
_NUMBA_TRANSFORM_ERROR: str | None = None


@dataclass(frozen=True)
class RigidTransform:
    """Rigid transform for row-major point arrays.

    The transform is applied as ``origin + R @ (point - origin) + translation``.
    """

    rotation: np.ndarray
    translation_m: np.ndarray

    def __post_init__(self) -> None:
        rotation = np.asarray(self.rotation, dtype=np.float64)
        translation = np.asarray(self.translation_m, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise ValueError("rotation must have shape (3, 3).")
        if translation.shape != (3,):
            raise ValueError("translation_m must have shape (3,).")
        if not np.all(np.isfinite(rotation)):
            raise ValueError("rotation must contain finite values.")
        if not np.all(np.isfinite(translation)):
            raise ValueError("translation_m must contain finite values.")
        object.__setattr__(self, "rotation", np.ascontiguousarray(rotation))
        object.__setattr__(self, "translation_m", np.ascontiguousarray(translation))

    @classmethod
    def identity(cls) -> "RigidTransform":
        """Return an identity transform."""

        return cls(np.eye(3, dtype=np.float64), np.zeros(3, dtype=np.float64))

    @classmethod
    def translation(cls, by: Iterable[float]) -> "RigidTransform":
        """Return a pure translation transform."""

        return cls(np.eye(3, dtype=np.float64), _vec3(by, name="by"))

    @classmethod
    def from_axis_angle(
        cls,
        axis: Iterable[float],
        *,
        angle_rad: float | None = None,
        angle_deg: float | None = None,
        translation_m: Iterable[float] = (0.0, 0.0, 0.0),
    ) -> "RigidTransform":
        """Build a transform from an axis-angle rotation."""

        if (angle_rad is None) == (angle_deg is None):
            raise ValueError("provide exactly one of angle_rad or angle_deg.")
        angle = float(np.deg2rad(angle_deg) if angle_deg is not None else angle_rad)
        axis_arr = _vec3(axis, name="axis")
        norm = float(np.linalg.norm(axis_arr))
        if norm == 0.0:
            raise ValueError("axis must be non-zero.")
        x, y, z = axis_arr / norm
        c = float(np.cos(angle))
        s = float(np.sin(angle))
        one_c = 1.0 - c
        rotation = np.array(
            [
                [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
                [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
                [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
            ],
            dtype=np.float64,
        )
        return cls(rotation, _vec3(translation_m, name="translation_m"))

    def apply(
        self,
        points: np.ndarray,
        *,
        origin: Iterable[float] = (0.0, 0.0, 0.0),
        backend: TransformBackend = "numpy",
    ) -> np.ndarray:
        """Apply this transform to ``(..., 3)`` points."""

        return _apply_transform_points(
            points,
            rotation=self.rotation,
            translation=self.translation_m,
            origin=_vec3(origin, name="origin"),
            backend=backend,
        )


@dataclass(frozen=True)
class BeachScene:
    """Editable what-if view of one BEACH charge snapshot.

    Geometry edits are immutable: ``move``/``rotate`` return a new scene, keeping
    charges attached to their original mesh elements.
    """

    result: FortranRunResult
    step: int | None
    centers_m: np.ndarray
    triangles_m: np.ndarray
    charges_C: np.ndarray
    mesh_ids: np.ndarray
    config_path: Path | None = None
    transform_backend: TransformBackend = "numpy"

    def __post_init__(self) -> None:
        centers = np.asarray(self.centers_m, dtype=np.float64)
        triangles = np.asarray(self.triangles_m, dtype=np.float64)
        charges = np.asarray(self.charges_C, dtype=np.float64)
        mesh_ids = np.asarray(self.mesh_ids, dtype=np.int64)
        if centers.ndim != 2 or centers.shape[1] != 3:
            raise ValueError("centers_m must have shape (n_elem, 3).")
        if triangles.ndim != 3 or triangles.shape[1:] != (3, 3):
            raise ValueError("triangles_m must have shape (n_elem, 3, 3).")
        if charges.ndim != 1 or charges.shape[0] != centers.shape[0]:
            raise ValueError("charges_C must have shape (n_elem,).")
        if mesh_ids.ndim != 1 or mesh_ids.shape[0] != centers.shape[0]:
            raise ValueError("mesh_ids must have shape (n_elem,).")
        if triangles.shape[0] != centers.shape[0]:
            raise ValueError("triangles_m and centers_m must have the same length.")
        if not np.all(np.isfinite(centers)):
            raise ValueError("centers_m must contain finite values.")
        if not np.all(np.isfinite(triangles)):
            raise ValueError("triangles_m must contain finite values.")
        if not np.all(np.isfinite(charges)):
            raise ValueError("charges_C must contain finite values.")
        backend = _normalize_transform_backend(self.transform_backend)
        object.__setattr__(self, "centers_m", _readonly_array(centers, np.float64))
        object.__setattr__(self, "triangles_m", _readonly_array(triangles, np.float64))
        object.__setattr__(self, "charges_C", _readonly_array(charges, np.float64))
        object.__setattr__(self, "mesh_ids", _readonly_array(mesh_ids, np.int64))
        object.__setattr__(self, "transform_backend", backend)
        if self.config_path is not None:
            object.__setattr__(self, "config_path", Path(self.config_path))

    @classmethod
    def from_result(
        cls,
        result: FortranRunResult | object,
        *,
        step: int | None = -1,
        config_path: str | Path | None = None,
        transform_backend: TransformBackend = "numpy",
    ) -> "BeachScene":
        """Create a scene from a BEACH output snapshot."""

        resolved = _resolve_result(result)
        triangles = np.asarray(_require_triangles(resolved), dtype=np.float64)
        return cls(
            result=resolved,
            step=step,
            centers_m=_triangle_centers(triangles),
            triangles_m=triangles,
            charges_C=_charges_for_step(resolved, step=step),
            mesh_ids=_mesh_ids_or_default(resolved),
            config_path=None if config_path is None else Path(config_path),
            transform_backend=transform_backend,
        )

    @property
    def centers(self) -> np.ndarray:
        """Alias for element centers in meters."""

        return self.centers_m

    @property
    def triangles(self) -> np.ndarray:
        """Alias for triangle vertices in meters."""

        return self.triangles_m

    @property
    def charges(self) -> np.ndarray:
        """Alias for per-element charges in coulombs."""

        return self.charges_C

    @property
    def object_ids(self) -> tuple[int, ...]:
        """Return sorted mesh ids present in this scene."""

        return tuple(int(v) for v in np.unique(self.mesh_ids))

    def object_center(self, mesh_ids: int | Iterable[int]) -> np.ndarray:
        """Return the centroid of one or more objects."""

        ids = _mesh_id_tuple(mesh_ids)
        mask = self._mask_for_mesh_ids(ids)
        return np.asarray(self.centers_m[mask].mean(axis=0), dtype=np.float64)

    def with_transform_backend(self, backend: TransformBackend) -> "BeachScene":
        """Return a scene using another point-transform backend."""

        return replace(self, transform_backend=backend)

    def move(
        self,
        mesh_ids: int | Iterable[int],
        *,
        by: Iterable[float],
    ) -> "BeachScene":
        """Translate one or more objects by ``by`` meters."""

        return self.transform(
            mesh_ids,
            RigidTransform.translation(by),
            origin="origin",
        )

    def rotate(
        self,
        mesh_ids: int | Iterable[int],
        *,
        axis: Iterable[float],
        angle_rad: float | None = None,
        angle_deg: float | None = None,
        origin: RotationOrigin | Iterable[float] = "object_center",
    ) -> "BeachScene":
        """Rotate one or more objects around an origin."""

        transform = RigidTransform.from_axis_angle(
            axis,
            angle_rad=angle_rad,
            angle_deg=angle_deg,
        )
        return self.transform(mesh_ids, transform, origin=origin)

    def transform(
        self,
        mesh_ids: int | Iterable[int],
        transform: RigidTransform,
        *,
        origin: RotationOrigin | Iterable[float] = "object_center",
    ) -> "BeachScene":
        """Apply a rigid transform to one or more objects."""

        ids = _mesh_id_tuple(mesh_ids)
        mask = self._mask_for_mesh_ids(ids)
        origin_arr = self._resolve_origin(ids, origin)
        centers = np.asarray(self.centers_m).copy()
        triangles = np.asarray(self.triangles_m).copy()
        centers[mask] = transform.apply(
            centers[mask],
            origin=origin_arr,
            backend=self.transform_backend,
        )
        triangles[mask] = transform.apply(
            triangles[mask],
            origin=origin_arr,
            backend=self.transform_backend,
        )
        return replace(self, centers_m=centers, triangles_m=triangles)

    def field_kernel(
        self,
        *,
        softening: float | None = None,
        periodic2: Mapping[str, object] | None = None,
        theta: float | None = None,
        leaf_max: int | None = None,
        order: int = 4,
        config_path: str | Path | None = None,
        library_path: str | Path | None = None,
    ) -> FieldKernel:
        """Build a Fortran field kernel for the edited scene geometry."""

        options = field_kernel_options_from_result(
            self.result,
            softening=softening,
            periodic2=periodic2,
            theta=theta,
            leaf_max=leaf_max,
            order=order,
            config_path=self._resolve_config_path(config_path),
        )
        return FieldKernel(
            self.centers_m,
            self.charges_C,
            options=options,
            library_path=library_path,
        )

    def calc_object_forces_kernel(
        self,
        *,
        target_mesh_ids: int | Iterable[int] | None = None,
        softening: float | None = None,
        periodic2: Mapping[str, object] | None = None,
        theta: float | None = None,
        leaf_max: int | None = None,
        order: int = 4,
        config_path: str | Path | None = None,
        library_path: str | Path | None = None,
    ) -> tuple[KernelObjectForceRecord, ...]:
        """Compute object-wise force/torque for the edited scene geometry."""

        target_ids = self._target_ids(target_mesh_ids)
        options = field_kernel_options_from_result(
            self.result,
            softening=softening,
            periodic2=periodic2,
            theta=theta,
            leaf_max=leaf_max,
            order=order,
            config_path=self._resolve_config_path(config_path),
        )
        records: list[KernelObjectForceRecord] = []
        with FieldKernel(
            self.centers_m,
            self.charges_C,
            options=options,
            library_path=library_path,
        ) as kernel:
            for mesh_id in target_ids:
                mask = self.mesh_ids == mesh_id
                source_q = np.asarray(self.charges_C, dtype=np.float64).copy()
                source_q[mask] = 0.0
                kernel.update_charges(source_q)
                target_centers = self.centers_m[mask]
                target_q = self.charges_C[mask]
                center = np.asarray(target_centers.mean(axis=0), dtype=np.float64)
                force, torque = kernel.force_on_charges(
                    target_centers,
                    target_q,
                    origin=center,
                )
                records.append(
                    KernelObjectForceRecord(
                        mesh_id=mesh_id,
                        step=self.step,
                        total_charge_C=float(np.sum(target_q)),
                        center_m=center,
                        force_N=force,
                        torque_Nm=torque,
                    )
                )
        return tuple(records)

    def _mask_for_mesh_ids(self, mesh_ids: tuple[int, ...]) -> np.ndarray:
        available = set(self.object_ids)
        missing = [mesh_id for mesh_id in mesh_ids if mesh_id not in available]
        if missing:
            raise ValueError(f"unknown mesh id(s): {missing}. available={sorted(available)}")
        return np.isin(self.mesh_ids, np.asarray(mesh_ids, dtype=np.int64))

    def _target_ids(self, target_mesh_ids: int | Iterable[int] | None) -> tuple[int, ...]:
        if target_mesh_ids is None:
            return self.object_ids
        target_ids = _mesh_id_tuple(target_mesh_ids)
        self._mask_for_mesh_ids(target_ids)
        return target_ids

    def _resolve_origin(
        self,
        mesh_ids: tuple[int, ...],
        origin: RotationOrigin | Iterable[float],
    ) -> np.ndarray:
        if isinstance(origin, str):
            if origin == "object_center":
                return self.object_center(mesh_ids)
            if origin == "origin":
                return np.zeros(3, dtype=np.float64)
            raise ValueError('origin must be "object_center", "origin", or a 3-vector.')
        return _vec3(origin, name="origin")

    def _resolve_config_path(self, config_path: str | Path | None) -> Path | None:
        if config_path is not None:
            return Path(config_path)
        return self.config_path


def _mesh_id_tuple(mesh_ids: int | Iterable[int]) -> tuple[int, ...]:
    if isinstance(mesh_ids, (int, np.integer)):
        return (int(mesh_ids),)
    ids = tuple(dict.fromkeys(int(v) for v in mesh_ids))
    if len(ids) == 0:
        raise ValueError("at least one mesh id must be provided.")
    return ids


def _vec3(value: Iterable[float], *, name: str) -> np.ndarray:
    arr = np.asarray(list(value), dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"{name} must contain exactly 3 values.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return np.ascontiguousarray(arr)


def _readonly_array(value: np.ndarray, dtype: object) -> np.ndarray:
    arr = np.ascontiguousarray(np.asarray(value, dtype=dtype).copy())
    arr.setflags(write=False)
    return arr


def _normalize_transform_backend(backend: str) -> TransformBackend:
    if backend not in {"numpy", "numba", "auto"}:
        raise ValueError('transform_backend must be "numpy", "numba", or "auto".')
    return backend  # type: ignore[return-value]


def _apply_transform_points(
    points: np.ndarray,
    *,
    rotation: np.ndarray,
    translation: np.ndarray,
    origin: np.ndarray,
    backend: TransformBackend,
) -> np.ndarray:
    arr = np.asarray(points, dtype=np.float64)
    if arr.shape[-1:] != (3,):
        raise ValueError("points must have shape (..., 3).")
    original_shape = arr.shape
    flat = np.ascontiguousarray(arr.reshape(-1, 3))
    backend = _normalize_transform_backend(backend)
    if backend == "numba" or (backend == "auto" and flat.shape[0] >= 8192):
        kernel = _get_numba_transform_kernel()
        if kernel is not None:
            return np.asarray(
                kernel(flat, rotation, translation, origin),
                dtype=np.float64,
            ).reshape(original_shape)
        if backend == "numba":
            raise RuntimeError(_numba_backend_error())
    return _apply_transform_points_numpy(
        flat,
        rotation=rotation,
        translation=translation,
        origin=origin,
    ).reshape(original_shape)


def _apply_transform_points_numpy(
    points: np.ndarray,
    *,
    rotation: np.ndarray,
    translation: np.ndarray,
    origin: np.ndarray,
) -> np.ndarray:
    return (points - origin) @ rotation.T + origin + translation


def _apply_transform_points_loop(
    points: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
    origin: np.ndarray,
) -> np.ndarray:
    out = np.empty_like(points)
    for i in range(points.shape[0]):
        for j in range(3):
            value = origin[j] + translation[j]
            for k in range(3):
                value += rotation[j, k] * (points[i, k] - origin[k])
            out[i, j] = value
    return out


def _get_numba_transform_kernel():  # type: ignore[no-untyped-def]
    global _NUMBA_TRANSFORM_ERROR, _NUMBA_TRANSFORM_KERNEL
    if _NUMBA_TRANSFORM_KERNEL is not None:
        return _NUMBA_TRANSFORM_KERNEL
    if _NUMBA_TRANSFORM_ERROR is not None:
        return None
    try:
        import numba as nb  # type: ignore[import-not-found]

        _NUMBA_TRANSFORM_KERNEL = nb.njit(cache=True)(_apply_transform_points_loop)
    except Exception as exc:  # pragma: no cover - depends on optional dependency
        _NUMBA_TRANSFORM_ERROR = str(exc)
        return None
    return _NUMBA_TRANSFORM_KERNEL


def _numba_backend_error() -> str:
    if _NUMBA_TRANSFORM_ERROR:
        return f"numba transform backend is unavailable: {_NUMBA_TRANSFORM_ERROR}"
    return "numba transform backend is unavailable. Install BEACH with the accel extra."
