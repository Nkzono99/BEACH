"""Readers for auxiliary diagnostic CSV products (light: saved values only)."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from .types import (
    FixedCurrentHistoryEntry,
    MatchingPlaneSpectrumHistoryEntry,
    TopReferenceHistoryEntry,
)


def _rows(
    path: str | Path,
    fields: dict[str, tuple[str, type]],
    *,
    defaults: dict[str, float] | None = None,
    allow_nan: tuple[str, ...] = (),
) -> list[dict[str, Any]] | None:
    path = Path(path)
    if not path.exists():
        return None
    defaults = defaults or {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = set(fields) - set(defaults)
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path.name} is missing required columns.")
        result = []
        for line_number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"{path.name} line {line_number} has an invalid column count."
                )
            values = {}
            for key, (attribute, kind) in fields.items():
                raw = row.get(key)
                value = defaults[key] if raw is None else kind(raw.strip())
                if isinstance(value, float) and not math.isfinite(value):
                    if key not in allow_nan or not math.isnan(value):
                        raise ValueError(f"{path.name} {key} must be finite.")
                values[attribute] = value
            result.append(values)
    return result


def load_top_reference_history(
    path: str | Path,
) -> tuple[TopReferenceHistoryEntry, ...] | None:
    """Read z-high plane samples; absent files return ``None``, header-only files ``()``.

    The released 3.0 eight-column format defaults the two switch mismatch values
    to zero, matching the simulator's persisted-history migration.
    """

    rows = _rows(
        path,
        {
            "batch": ("batch", int),
            "simulated_time_s": ("simulated_time_s", float),
            "z_high_m": ("z_high_m", float),
            "sample_n": ("sample_n", int),
            "potential_mean_V": ("potential_mean_v", float),
            "potential_std_V": ("potential_std_v", float),
            "potential_min_V": ("potential_min_v", float),
            "potential_max_V": ("potential_max_v", float),
            "switch_phi_rms_mismatch_V": ("switch_phi_rms_mismatch_v", float),
            "switch_Ez_rms_mismatch_V_m": ("switch_ez_rms_mismatch_v_m", float),
        },
        defaults={"switch_phi_rms_mismatch_V": 0.0, "switch_Ez_rms_mismatch_V_m": 0.0},
    )
    return (
        None if rows is None else tuple(TopReferenceHistoryEntry(**row) for row in rows)
    )


def load_matching_plane_spectrum_history(
    path: str | Path,
) -> tuple[MatchingPlaneSpectrumHistoryEntry, ...] | None:
    """Read nonempty photoelectron energy bins; absent files return ``None``."""

    rows = _rows(
        path,
        {
            "batch": ("batch", int),
            "simulated_time_s": ("simulated_time_s", float),
            "bin": ("bin", int),
            "energy_low_eV": ("energy_low_ev", float),
            "energy_high_eV": ("energy_high_ev", float),
            "photoelectron_outward_flux_m2_s": (
                "photoelectron_outward_flux_m2_s",
                float,
            ),
        },
    )
    return (
        None
        if rows is None
        else tuple(MatchingPlaneSpectrumHistoryEntry(**row) for row in rows)
    )


def load_fixed_current_history(
    path: str | Path,
) -> tuple[FixedCurrentHistoryEntry, ...] | None:
    """Read fixed-current window diagnostics, including the defined NaN ratio for zero tracking."""

    rows = _rows(
        path,
        {
            "batch": ("batch", int),
            "simulated_time_s": ("simulated_time_s", float),
            "window_batches": ("window_batches", int),
            "species_idx": ("species_idx", int),
            "channel": ("channel", str),
            "tracked_charge_C": ("tracked_charge_c", float),
            "target_charge_C": ("target_charge_c", float),
            "target_over_tracked": ("target_over_tracked", float),
        },
        allow_nan=("target_over_tracked",),
    )
    return (
        None if rows is None else tuple(FixedCurrentHistoryEntry(**row) for row in rows)
    )
