"""Standard-output report format shared by ``beachx`` subcommands.

A command builds one ordered mapping and passes it to :func:`emit`. Without
``--json`` every scalar prints as one ``key=value`` line and every table (a list of
mappings) prints one line per row of space-separated ``field=value`` pairs. With
``--json`` the same mapping prints as one JSON object. An empty table prints as
``key=[]``.
"""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Mapping
from pathlib import Path

import numpy as np


def emit(args: argparse.Namespace, report: Mapping[str, object]) -> None:
    """Print ``report`` as ``key=value`` lines, or as JSON with ``--json``."""

    if getattr(args, "json", False):
        print(json.dumps(json_value(report), allow_nan=False))
        return
    for key, value in report.items():
        if _is_table(value):
            for row in value:
                print(" ".join(f"{field}={_text(item)}" for field, item in row.items()))
        else:
            print(f"{key}={_text(value)}")


def json_value(value: object) -> object:
    """Convert numpy values, paths, and non-finite floats to JSON values.

    Non-finite floats become the strings ``"inf"``, ``"-inf"``, and ``"nan"``.
    """

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return [json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    if isinstance(value, Mapping):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    raise TypeError(f"unsupported report value: {type(value).__name__}")


def _is_table(value: object) -> bool:
    return (
        isinstance(value, list)
        and bool(value)
        and all(isinstance(row, Mapping) for row in value)
    )


def _text(value: object) -> str:
    value = json_value(value)
    if isinstance(value, float):
        return f"{value:.6e}"
    if isinstance(value, str):
        return value
    return json.dumps(value, separators=(",", ":"))
