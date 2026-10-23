"""Shared parser and core contract for Fortran ``summary.txt`` files."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from dataclasses import dataclass
from math import isfinite
from typing import Literal


def parse_summary_text(text: str, *, source: str = "summary.txt") -> dict[str, str]:
    """Parse key/value summary text and reject duplicate normalized keys."""

    values: dict[str, str] = {}
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        normalized_key = key.strip()
        if not normalized_key:
            raise ValueError(f"{source} contains an empty key on line {line_number}.")
        if normalized_key in values:
            raise ValueError(
                f"{source} contains duplicate key {normalized_key!r} on line {line_number}."
            )
        values[normalized_key] = value.strip()
    return values


def load_summary_file(
    path: Path,
    *,
    required_keys: Iterable[str] = (),
) -> dict[str, str]:
    """Load one summary file and enforce its required-key subset."""

    values = parse_summary_text(path.read_text(encoding="utf-8"), source=str(path))
    missing = [key for key in required_keys if key not in values]
    if missing:
        raise ValueError(f"Missing keys in summary: {missing}")
    return values


@dataclass(frozen=True)
class SummaryField:
    """One persisted summary scalar and its destination in :class:`RunResult`.

    ``kind`` states the scalar contract; ``nonnegative_float`` and
    ``nonnegative_int`` reject negative values. Missing optional values use
    ``default``. Complex field/matching-plane receipts have their own schema.
    """

    key: str
    attribute: str
    group: Literal["summary", "receipts", "diagnostics", "field"]
    kind: Literal[
        "nonnegative_int", "nonnegative_float", "float", "bool", "str", "model"
    ]
    default: str | int | float | bool | None = None
    required: bool = False


SUMMARY_FIELDS = (
    SummaryField(
        "mesh_nelem", "mesh_nelem", "summary", "nonnegative_int", required=True
    ),
    SummaryField(
        "processed_particles",
        "processed_particles",
        "summary",
        "nonnegative_int",
        required=True,
    ),
    SummaryField("absorbed", "absorbed", "summary", "nonnegative_int", required=True),
    SummaryField("escaped", "escaped", "summary", "nonnegative_int", required=True),
    SummaryField("batches", "batches", "summary", "nonnegative_int", required=True),
    SummaryField(
        "last_rel_change",
        "last_rel_change",
        "summary",
        "nonnegative_float",
        required=True,
    ),
    SummaryField(
        "escaped_boundary", "escaped_boundary", "summary", "nonnegative_int", 0
    ),
    SummaryField(
        "survived_max_step", "survived_max_step", "summary", "nonnegative_int", 0
    ),
    SummaryField(
        "simulated_time_s", "simulated_time_s", "summary", "nonnegative_float"
    ),
    SummaryField(
        "checkpoint_schema_version",
        "checkpoint_schema_version",
        "receipts",
        "nonnegative_int",
    ),
    SummaryField("model_fingerprint", "model_fingerprint", "receipts", "str"),
    SummaryField("mesh_fingerprint", "mesh_fingerprint", "receipts", "str"),
    SummaryField("species_fingerprint", "species_fingerprint", "receipts", "str"),
    SummaryField("field_source_model", "source_model", "field", "model", "unknown"),
    SummaryField("field_kernel_id", "kernel_id", "field", "str"),
    SummaryField(
        "charge_ledger_residual_C", "charge_ledger_residual_c", "diagnostics", "float"
    ),
    SummaryField("periodic2_cache_hit", "periodic2_cache_hit", "diagnostics", "bool"),
    SummaryField(
        "periodic2_operator_build_count",
        "periodic2_operator_build_count",
        "diagnostics",
        "nonnegative_int",
    ),
    SummaryField(
        "periodic2_cache_fingerprint",
        "periodic2_cache_fingerprint",
        "diagnostics",
        "str",
    ),
    SummaryField("periodic2_cache_path", "periodic2_cache_path", "diagnostics", "str"),
    SummaryField(
        "periodic2_generation_tolerance",
        "periodic2_generation_tolerance",
        "diagnostics",
        "float",
    ),
    SummaryField(
        "multiple_box_events_retry_attempted",
        "multiple_box_events_retry_attempted",
        "diagnostics",
        "nonnegative_int",
        0,
    ),
    SummaryField(
        "multiple_box_events_retry_resolved",
        "multiple_box_events_retry_resolved",
        "diagnostics",
        "nonnegative_int",
        0,
    ),
    SummaryField(
        "multiple_box_events_soft_discarded",
        "multiple_box_events_soft_discarded",
        "diagnostics",
        "nonnegative_int",
        0,
    ),
    SummaryField(
        "multiple_box_events_soft_discarded_abs_charge_C",
        "multiple_box_events_soft_discarded_abs_charge_c",
        "diagnostics",
        "nonnegative_float",
        0.0,
    ),
    SummaryField(
        "multiple_box_events_soft_discard_fraction",
        "multiple_box_events_soft_discard_fraction",
        "diagnostics",
        "nonnegative_float",
        0.0,
    ),
    SummaryField(
        "adaptive_nonzero_mode_rejected_trials",
        "adaptive_nonzero_mode_rejected_trials",
        "diagnostics",
        "nonnegative_int",
        0,
    ),
    SummaryField(
        "adaptive_nonzero_mode_last_batch_duration_s",
        "adaptive_nonzero_mode_last_batch_duration_s",
        "diagnostics",
        "nonnegative_float",
    ),
    SummaryField(
        "adaptive_nonzero_mode_last_potential_step_V",
        "adaptive_nonzero_mode_last_potential_step_v",
        "diagnostics",
        "nonnegative_float",
    ),
    SummaryField(
        "adaptive_nonzero_mode_omp_threads",
        "adaptive_nonzero_mode_omp_threads",
        "diagnostics",
        "nonnegative_int",
    ),
    SummaryField(
        "periodic2_max_nonzero_mode_potential_step_V",
        "periodic2_max_nonzero_mode_potential_step_v",
        "diagnostics",
        "nonnegative_float",
    ),
)

# The core required subset is derived from the same table used by the reader.
CORE_SUMMARY_REQUIRED_KEYS = tuple(
    field.key for field in SUMMARY_FIELDS if field.required
)


def _parse_summary_fields(data: dict[str, str]) -> dict[str, dict[str, object]]:
    groups: dict[str, dict[str, object]] = {
        name: {} for name in ("summary", "receipts", "diagnostics", "field")
    }
    for field in SUMMARY_FIELDS:
        value: object = field.default
        raw = data.get(field.key)
        if raw is None:
            if field.required:
                raise ValueError(f"Missing keys in summary: {[field.key]}")
        elif field.kind == "nonnegative_int":
            value = int(raw)
            if value < 0:
                raise ValueError(f"summary.txt {field.key} must be >= 0.")
        elif field.kind in {"float", "nonnegative_float"}:
            value = float(raw)
            if not isfinite(value) or (
                field.kind == "nonnegative_float" and value < 0.0
            ):
                qualifier = (
                    "finite and >= 0" if field.kind == "nonnegative_float" else "finite"
                )
                raise ValueError(f"summary.txt {field.key} must be {qualifier}.")
        elif field.kind == "bool":
            if raw.strip().lower() not in {"t", "true", "f", "false"}:
                raise ValueError(f"summary.txt {field.key} must be true or false.")
            value = raw.strip().lower() in {"t", "true"}
        else:
            value = raw.strip().lower() if field.kind == "model" else raw
        groups[field.group][field.attribute] = value
    return groups


def summary_fields_markdown() -> str:
    """Generate the summary scalar reference table used in OutputReference.md."""

    lines = [
        "| summary.txt key | RunResult attribute | Type | Required / default |",
        "|---|---|---|---|",
    ]
    for field in SUMMARY_FIELDS:
        default = "required" if field.required else f"`{field.default!r}`"
        lines.append(
            f"| `{field.key}` | `{field.group}.{field.attribute}` | `{field.kind}` | {default} |"
        )
    return "\n".join(lines)
