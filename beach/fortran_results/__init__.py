"""Utilities for reading and visualizing Fortran simulation outputs."""

from .animation import _select_frame_columns, animate_history_mesh
from .constants import K_COULOMB
from .coulomb import calc_coulomb
from .facade import Beach
from .field_lines import (
    compute_electric_field_points,
    plot_field_lines_3d,
    trace_field_lines,
)
from .history import FortranChargeHistory
from .io import list_fortran_runs, load_fortran_result
from .mesh import _surface_charge_density
from .mobility import analyze_coulomb_mobility
from .plotting import (
    plot_charge_mesh,
    plot_charges,
    plot_coulomb_force_matrix,
    plot_mesh_source_boxplot,
    plot_potential_mesh,
    plot_potential_slices,
)
from .potential import (
    compute_potential_mesh,
    compute_potential_points,
    compute_potential_slices,
)
from .types import (
    CoulombInteraction,
    CoulombMobilityAnalysis,
    CoulombMobilityRecord,
    FortranRunResult,
    MeshSelection,
    MeshSource,
    PotentialSlice2D,
)

__all__ = [
    "K_COULOMB",
    "FortranChargeHistory",
    "MeshSource",
    "MeshSelection",
    "CoulombInteraction",
    "CoulombMobilityRecord",
    "CoulombMobilityAnalysis",
    "FortranRunResult",
    "PotentialSlice2D",
    "load_fortran_result",
    "list_fortran_runs",
    "Beach",
    "calc_coulomb",
    "analyze_coulomb_mobility",
    "compute_electric_field_points",
    "trace_field_lines",
    "plot_field_lines_3d",
    "plot_charges",
    "plot_charge_mesh",
    "plot_coulomb_force_matrix",
    "plot_mesh_source_boxplot",
    "compute_potential_mesh",
    "compute_potential_points",
    "compute_potential_slices",
    "plot_potential_slices",
    "plot_potential_mesh",
    "animate_history_mesh",
    "_select_frame_columns",
    "_surface_charge_density",
]
