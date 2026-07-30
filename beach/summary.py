"""Shared parser and core contract for Fortran ``summary.txt`` files."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path


CURRENT_CHECKPOINT_SCHEMA_VERSION = 5
CORE_SUMMARY_REQUIRED_KEYS = (
    "mesh_nelem",
    "processed_particles",
    "absorbed",
    "escaped",
    "batches",
    "last_rel_change",
)


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
