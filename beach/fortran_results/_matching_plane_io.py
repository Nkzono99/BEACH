"""Read and validate committed matching-plane state and coupling history."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ._summary_values import (
    _parse_nonnegative_finite_float,
    _parse_nonnegative_int,
    _parse_optional_bool,
    _parse_optional_finite_float,
)
from .types import MatchingPlaneHistoryEntry, MatchingPlaneState

_MATCHING_PLANE_HISTORY_HEADER = (
    "batch",
    "simulated_time_s",
    "D_H_C_m2",
    "phi_H_V",
    "electron_inward_flux_m2_s",
    "ion_inward_flux_m2_s",
    "electron_access_potential_V",
    "ion_access_potential_V",
    "photoelectron_barrier_potential_V",
    "photoelectron_outward_flux_m2_s",
    "photoelectron_mean_normal_energy_eV",
    "electron_outward_flux_m2_s",
    "ion_outward_flux_m2_s",
    "photoelectron_return_flux_m2_s",
    "photoelectron_escape_flux_m2_s",
    "iterations",
    "residual",
)


def _parse_matching_plane_state(data: dict[str, str]) -> MatchingPlaneState | None:
    valid = _parse_optional_bool(data, "matching_plane_state_valid")
    if valid is None or not valid:
        return None

    required = (
        "matching_plane_displacement_C_m2",
        "matching_plane_phi_V",
        "matching_plane_electron_inward_flux_m2_s",
        "matching_plane_ion_inward_flux_m2_s",
        "matching_plane_electron_access_potential_V",
        "matching_plane_ion_access_potential_V",
        "matching_plane_photoelectron_barrier_potential_V",
        "matching_plane_photoelectron_outward_flux_m2_s",
        "matching_plane_photoelectron_mean_normal_energy_eV",
        "matching_plane_electron_outward_flux_m2_s",
        "matching_plane_ion_outward_flux_m2_s",
        "matching_plane_photoelectron_return_flux_m2_s",
        "matching_plane_photoelectron_escape_flux_m2_s",
        "matching_plane_iterations",
        "matching_plane_residual",
    )
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(
            "summary.txt valid matching-plane state requires "
            + ", ".join(missing)
            + "."
        )

    state = MatchingPlaneState(
        displacement_c_m2=_parse_optional_finite_float(
            data, "matching_plane_displacement_C_m2"
        ),
        phi_v=_parse_optional_finite_float(data, "matching_plane_phi_V"),
        electron_inward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_electron_inward_flux_m2_s"],
            key="matching_plane_electron_inward_flux_m2_s",
        ),
        ion_inward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_ion_inward_flux_m2_s"],
            key="matching_plane_ion_inward_flux_m2_s",
        ),
        electron_access_potential_v=_parse_optional_finite_float(
            data, "matching_plane_electron_access_potential_V"
        ),
        ion_access_potential_v=_parse_optional_finite_float(
            data, "matching_plane_ion_access_potential_V"
        ),
        photoelectron_barrier_potential_v=_parse_optional_finite_float(
            data, "matching_plane_photoelectron_barrier_potential_V"
        ),
        photoelectron_outward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_outward_flux_m2_s"],
            key="matching_plane_photoelectron_outward_flux_m2_s",
        ),
        photoelectron_mean_normal_energy_ev=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_mean_normal_energy_eV"],
            key="matching_plane_photoelectron_mean_normal_energy_eV",
        ),
        electron_outward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_electron_outward_flux_m2_s"],
            key="matching_plane_electron_outward_flux_m2_s",
        ),
        ion_outward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_ion_outward_flux_m2_s"],
            key="matching_plane_ion_outward_flux_m2_s",
        ),
        photoelectron_return_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_return_flux_m2_s"],
            key="matching_plane_photoelectron_return_flux_m2_s",
        ),
        photoelectron_escape_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_escape_flux_m2_s"],
            key="matching_plane_photoelectron_escape_flux_m2_s",
        ),
        iterations=_parse_nonnegative_int(
            data["matching_plane_iterations"], key="matching_plane_iterations"
        ),
        residual=_parse_nonnegative_finite_float(
            data["matching_plane_residual"], key="matching_plane_residual"
        ),
    )
    _validate_matching_plane_state(state, context="summary.txt")
    return state


def _load_matching_plane_history_if_exists(
    path: Path,
) -> tuple[MatchingPlaneHistoryEntry, ...] | None:
    if not path.exists():
        return None

    entries: list[MatchingPlaneHistoryEntry] = []
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != _MATCHING_PLANE_HISTORY_HEADER:
            raise ValueError(
                "matching_plane_history.csv header does not match the 17-column contract."
            )
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(row.get(key) is None for key in _MATCHING_PLANE_HISTORY_HEADER):
                raise ValueError(
                    f"matching_plane_history.csv line {line_number} must have exactly 17 columns."
                )
            state = MatchingPlaneState(
                displacement_c_m2=_parse_history_finite(row, "D_H_C_m2", line_number),
                phi_v=_parse_history_finite(row, "phi_H_V", line_number),
                electron_inward_flux_m2_s=_parse_history_nonnegative(
                    row, "electron_inward_flux_m2_s", line_number
                ),
                ion_inward_flux_m2_s=_parse_history_nonnegative(
                    row, "ion_inward_flux_m2_s", line_number
                ),
                electron_access_potential_v=_parse_history_finite(
                    row, "electron_access_potential_V", line_number
                ),
                ion_access_potential_v=_parse_history_finite(
                    row, "ion_access_potential_V", line_number
                ),
                photoelectron_barrier_potential_v=_parse_history_finite(
                    row, "photoelectron_barrier_potential_V", line_number
                ),
                photoelectron_outward_flux_m2_s=_parse_history_nonnegative(
                    row, "photoelectron_outward_flux_m2_s", line_number
                ),
                photoelectron_mean_normal_energy_ev=_parse_history_nonnegative(
                    row, "photoelectron_mean_normal_energy_eV", line_number
                ),
                electron_outward_flux_m2_s=_parse_history_nonnegative(
                    row, "electron_outward_flux_m2_s", line_number
                ),
                ion_outward_flux_m2_s=_parse_history_nonnegative(
                    row, "ion_outward_flux_m2_s", line_number
                ),
                photoelectron_return_flux_m2_s=_parse_history_nonnegative(
                    row, "photoelectron_return_flux_m2_s", line_number
                ),
                photoelectron_escape_flux_m2_s=_parse_history_nonnegative(
                    row, "photoelectron_escape_flux_m2_s", line_number
                ),
                iterations=_parse_history_nonnegative_int(row, "iterations", line_number),
                residual=_parse_history_nonnegative(row, "residual", line_number),
            )
            _validate_matching_plane_state(
                state, context=f"matching_plane_history.csv line {line_number}"
            )
            batch = _parse_history_nonnegative_int(row, "batch", line_number)
            if batch <= 0 or (entries and batch <= entries[-1].batch):
                raise ValueError(
                    "matching_plane_history.csv batch values must be positive and strictly increasing."
                )
            entries.append(
                MatchingPlaneHistoryEntry(
                    batch=batch,
                    simulated_time_s=_parse_history_nonnegative(
                        row, "simulated_time_s", line_number
                    ),
                    state=state,
                )
            )
    return tuple(entries)


def _parse_history_finite(row: dict[str, str | None], key: str, line: int) -> float:
    raw = row.get(key)
    assert raw is not None
    parsed = float(raw)
    if not np.isfinite(parsed):
        raise ValueError(f"matching_plane_history.csv line {line} {key} must be finite.")
    return parsed


def _parse_history_nonnegative(
    row: dict[str, str | None], key: str, line: int
) -> float:
    parsed = _parse_history_finite(row, key, line)
    if parsed < 0.0:
        raise ValueError(
            f"matching_plane_history.csv line {line} {key} must be >= 0."
        )
    return parsed


def _parse_history_nonnegative_int(
    row: dict[str, str | None], key: str, line: int
) -> int:
    raw = row.get(key)
    assert raw is not None
    try:
        parsed = int(raw)
    except ValueError as exc:
        raise ValueError(
            f"matching_plane_history.csv line {line} {key} must be an integer."
        ) from exc
    if parsed < 0:
        raise ValueError(
            f"matching_plane_history.csv line {line} {key} must be >= 0."
        )
    return parsed


def _validate_matching_plane_state(
    state: MatchingPlaneState, *, context: str
) -> None:
    if state.iterations <= 0:
        raise ValueError(f"{context} matching-plane iterations must be > 0.")
    budget_scale = max(
        1.0,
        state.photoelectron_outward_flux_m2_s,
        state.photoelectron_return_flux_m2_s,
        state.photoelectron_escape_flux_m2_s,
    )
    budget_error = abs(
        state.photoelectron_outward_flux_m2_s
        - state.photoelectron_return_flux_m2_s
        - state.photoelectron_escape_flux_m2_s
    )
    if budget_error > np.sqrt(np.finfo(float).eps) * budget_scale:
        raise ValueError(f"{context} matching-plane photoelectron budget is inconsistent.")
