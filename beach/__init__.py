"""BEACH (BEM + Accumulated CHarge) Python package.

The public names below are defined in the layer subpackages :mod:`beach.io`,
:mod:`beach.native`, :mod:`beach.analysis`, and :mod:`beach.plot`, plus the
:class:`~beach.facade.Beach` facade. They are loaded on first access (PEP 562),
so lightweight subpackages such as :mod:`beach.config` and :mod:`beach.tooling`
import without numpy or matplotlib.
"""

from __future__ import annotations

import importlib
from typing import Any

_EXPORTS = {
    "beach.io.csv_history": (
        "load_top_reference_history", "load_fixed_current_history", "load_matching_plane_spectrum_history",
    ),
    "beach.io.performance": ("load_performance_profile",),
    "beach.io.summary": ("SUMMARY_FIELDS", "SummaryField", "summary_fields_markdown"),
    "beach.analysis.context": ("RunContext",),
    "beach.io.constants": ("K_COULOMB",),
    "beach.io.history": ("BatchHistory",),
    "beach.io.results": ("list_runs", "load_run"),
    "beach.io.types": (
        "ChargeLedgerEntry",
        "CoulombInteraction",
        "CoulombMobilityAnalysis",
        "CoulombMobilityRecord",
        "FieldReconstructionReceipt",
        "RunResult",
        "RunSummary", "RunMesh", "RunHistory", "RunReceipts", "RunField", "RunDiagnostics",
        "TopReferenceHistoryEntry", "FixedCurrentHistoryEntry", "MatchingPlaneSpectrumHistoryEntry",
        "MatchingPlaneHistoryEntry",
        "MatchingPlaneState",
        "MeshSelection",
        "MeshSource",
        "PotentialSlice2D",
        "SecondaryEmissionHistoryEntry",
    ),
    "beach.native.kernel": (
        "FieldKernel",
        "FieldKernelBuildInfo",
        "FieldKernelDiagnostics",
        "FieldKernelError",
        "FieldKernelOptions",
        "field_kernel_build_info",
    ),
    "beach.native.periodic_zero_mode": ("PeriodicZeroMode",),
    "beach.analysis.coulomb": ("calc_coulomb",),
    "beach.analysis.detachment": (
        "AdhesionProfile",
        "DetachmentResult",
        "ObjectForcePath",
        "ObjectWrench",
        "WrenchComponent",
    ),
    "beach.analysis.field_kernel": (
        "KernelObjectForceRecord",
        "calc_object_forces_kernel",
        "field_kernel_options_from_result",
    ),
    "beach.analysis.field_lines": ("compute_electric_field_points", "trace_field_lines"),
    "beach.analysis.mobility": ("analyze_coulomb_mobility",),
    "beach.analysis.object_interaction": ("ObjectInteractionSnapshot", "ObjectProbe"),
    "beach.analysis.periodic_force_oracle": (
        "FiniteShellConvergenceResult",
        "FiniteShellWrenchResult",
        "finite_shell_convergence",
        "finite_shell_wrench",
    ),
    "beach.analysis.potential": (
        "compute_potential_mesh",
        "compute_potential_points",
        "compute_potential_slices",
    ),
    "beach.analysis.scene": ("BeachScene", "RigidTransform"),
    "beach.plot.animation": ("animate_history_mesh",),
    "beach.plot.boxplot": ("plot_mesh_source_boxplot",),
    "beach.plot.field_lines": ("plot_field_lines_3d",),
    "beach.plot.matrix": ("plot_coulomb_force_matrix",),
    "beach.plot.mesh": ("plot_charge_mesh", "plot_charges", "plot_potential_mesh"),
    "beach.plot.slices": ("plot_potential_slices",),
    "beach.facade": ("Beach",),
}
_MODULE_OF = {name: module for module, names in _EXPORTS.items() for name in names}

__all__ = list(_MODULE_OF)


def __getattr__(name: str) -> Any:
    module = _MODULE_OF.get(name)
    if module is not None:
        value = getattr(importlib.import_module(module), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted({*globals(), *__all__})
