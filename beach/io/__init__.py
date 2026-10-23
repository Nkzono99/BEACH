"""Readers and data types for BEACH output directories (numpy only when accessed)."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS = {
    "history": ("BatchHistory",),
    "results": ("load_run", "list_runs"),
    "summary": (
        "SUMMARY_FIELDS",
        "SummaryField",
        "summary_fields_markdown",
        "load_summary_file",
        "parse_summary_text",
    ),
    "csv_history": (
        "load_top_reference_history",
        "load_matching_plane_spectrum_history",
        "load_fixed_current_history",
    ),
    "performance": ("load_performance_profile",),
    "types": (
        "RunResult",
        "RunSummary",
        "RunMesh",
        "RunHistory",
        "RunReceipts",
        "RunField",
        "RunDiagnostics",
        "MeshSource",
        "MeshSelection",
        "ChargeLedgerEntry",
        "FieldReconstructionReceipt",
        "MatchingPlaneState",
        "MatchingPlaneHistoryEntry",
        "MatchingPlaneSpectrumHistoryEntry",
        "TopReferenceHistoryEntry",
        "FixedCurrentHistoryEntry",
        "SecondaryEmissionHistoryEntry",
    ),
}
_MODULE_OF = {name: module for module, names in _EXPORTS.items() for name in names}
__all__ = list(_MODULE_OF)


def __getattr__(name: str) -> Any:
    module = _MODULE_OF.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f"{__name__}.{module}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted({*globals(), *__all__})
