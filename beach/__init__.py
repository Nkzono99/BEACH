"""BEACH (BEM + Accumulated CHarge) Python package.

The public names below come from :mod:`beach.fortran_results`. They are loaded
on first access (PEP 562), so lightweight submodules such as :mod:`beach.config`
and :mod:`beach.runops_adapter` import without numpy or matplotlib.
"""

from __future__ import annotations

import importlib
from typing import Any

__all__ = [
    "animate_history_mesh",
    "analyze_coulomb_mobility",
    "AdhesionProfile",
    "Beach",
    "BeachScene",
    "calc_coulomb",
    "calc_object_forces_kernel",
    "compute_electric_field_points",
    "ChargeLedgerEntry",
    "K_COULOMB",
    "CoulombInteraction",
    "CoulombMobilityAnalysis",
    "CoulombMobilityRecord",
    "DetachmentResult",
    "compute_potential_mesh",
    "compute_potential_points",
    "compute_potential_slices",
    "FieldReconstructionReceipt",
    "FortranChargeHistory",
    "FortranRunResult",
    "field_kernel_build_info",
    "field_kernel_options_from_result",
    "FieldKernel",
    "FieldKernelBuildInfo",
    "FieldKernelDiagnostics",
    "FieldKernelError",
    "FieldKernelOptions",
    "FiniteShellConvergenceResult",
    "FiniteShellWrenchResult",
    "KernelObjectForceRecord",
    "MatchingPlaneHistoryEntry",
    "MatchingPlaneState",
    "ObjectForcePath",
    "ObjectInteractionSnapshot",
    "ObjectProbe",
    "ObjectWrench",
    "PeriodicZeroMode",
    "MeshSelection",
    "MeshSource",
    "PotentialSlice2D",
    "RigidTransform",
    "load_fortran_result",
    "list_fortran_runs",
    "plot_charges",
    "plot_charge_mesh",
    "plot_coulomb_force_matrix",
    "plot_field_lines_3d",
    "plot_mesh_source_boxplot",
    "plot_potential_mesh",
    "plot_potential_slices",
    "trace_field_lines",
    "SecondaryEmissionHistoryEntry",
    "WrenchComponent",
    "finite_shell_convergence",
    "finite_shell_wrench",
]


def __getattr__(name: str) -> Any:
    if name in __all__:
        value = getattr(importlib.import_module(".fortran_results", __name__), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted({*globals(), *__all__})
