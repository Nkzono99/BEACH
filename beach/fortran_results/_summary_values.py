"""Parse shared scalar summary values without changing persisted formats."""

from __future__ import annotations

import numpy as np


def _parse_nonnegative_int(value: str, *, key: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise ValueError(f"summary.txt {key} must be >= 0.")
    return parsed


def _parse_nonnegative_finite_float(value: str, *, key: str) -> float:
    parsed = float(value)
    if not np.isfinite(parsed) or parsed < 0.0:
        raise ValueError(f"summary.txt {key} must be finite and >= 0.")
    return parsed


def _parse_optional_nonnegative_int(data: dict[str, str], key: str) -> int | None:
    if key not in data:
        return None
    return _parse_nonnegative_int(data[key], key=key)


def _parse_optional_finite_float(data: dict[str, str], key: str) -> float | None:
    if key not in data:
        return None
    parsed = float(data[key])
    if not np.isfinite(parsed):
        raise ValueError(f"summary.txt {key} must be finite.")
    return parsed


def _parse_optional_nonnegative_finite_float(
    data: dict[str, str], key: str
) -> float | None:
    if key not in data:
        return None
    return _parse_nonnegative_finite_float(data[key], key=key)


def _parse_optional_bool(data: dict[str, str], key: str) -> bool | None:
    if key not in data:
        return None
    value = data[key].strip().lower()
    if value in {"t", "true"}:
        return True
    if value in {"f", "false"}:
        return False
    raise ValueError(f"summary.txt {key} must be true or false.")
