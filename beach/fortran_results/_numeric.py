"""Shared numeric input conversion and immutable arrays for result analysis."""

from __future__ import annotations

import operator
from typing import Iterable

import numpy as np


def _vec3(value: Iterable[float], *, name: str) -> np.ndarray:
    arr = np.asarray(list(value), dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"{name} must contain exactly 3 values.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return np.ascontiguousarray(arr)


def _points(value: np.ndarray) -> np.ndarray:
    points = np.asarray(value, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n_points, 3).")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must contain finite values.")
    return points


def _finite_scalar(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _readonly(value: np.ndarray) -> np.ndarray:
    result = np.array(value, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


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
