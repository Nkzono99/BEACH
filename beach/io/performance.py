"""Read the saved performance profile independently of plotting or CLI imports."""

from __future__ import annotations

import csv
from pathlib import Path


def load_performance_profile(
    path: str | Path,
) -> tuple[dict[str, str], list[dict[str, str | float]]]:
    """Load one ``performance_profile.csv`` file."""

    path = Path(path)
    metadata: dict[str, str] = {}
    csv_lines: list[str] = []

    if not path.exists():
        raise FileNotFoundError(path)

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            comment = stripped[1:].strip()
            if "=" in comment:
                key, value = comment.split("=", 1)
                metadata[key.strip()] = value.strip()
            continue
        csv_lines.append(line)

    if not csv_lines:
        raise ValueError("CSV body is empty.")

    reader = csv.DictReader(csv_lines)
    rows: list[dict[str, str | float]] = []
    for row in reader:
        if row is None:
            continue
        rows.append(
            {
                "region": _required_field(row, "region"),
                "calls_sum": float(_required_field(row, "calls_sum")),
                "calls_mean": float(_required_field(row, "calls_mean")),
                "rank_min_s": float(_required_field(row, "rank_min_s")),
                "rank_mean_s": float(_required_field(row, "rank_mean_s")),
                "rank_max_s": float(_required_field(row, "rank_max_s")),
                "imbalance_ratio": float(_required_field(row, "imbalance_ratio")),
            }
        )

    if not rows:
        raise ValueError("No profile rows were found.")
    return metadata, rows


def _required_field(row: dict[str, str], key: str) -> str:
    value = row.get(key)
    if value is None or value.strip() == "":
        raise ValueError(f"Missing required column: {key}")
    return value.strip()
