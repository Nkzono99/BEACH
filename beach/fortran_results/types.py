"""Public data types for Fortran output handling."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .history import FortranChargeHistory


@dataclass(frozen=True)
class MeshSource:
    """Metadata for one source mesh record.

    Attributes
    ----------
    mesh_id : int
        Mesh identifier written in ``mesh_sources.csv``.
    source_kind : str
        Mesh source category (for example ``"template"`` or ``"obj"``).
    template_kind : str
        Template subtype when ``source_kind`` is ``"template"``.
    elem_count : int
        Number of triangle elements belonging to the mesh.
    """

    mesh_id: int
    source_kind: str
    template_kind: str
    elem_count: int


@dataclass(frozen=True)
class MeshSelection:
    """Selection of one or more meshes extracted from one run.

    Attributes
    ----------
    directory : pathlib.Path
        Output directory this selection belongs to.
    mesh_ids : tuple of int
        Mesh identifiers included in this selection.
    elem_indices : numpy.ndarray
        Original element indices in the full mesh arrays.
    triangles : numpy.ndarray
        Selected triangle coordinates with shape ``(n_elem, 3, 3)``.
    charges : numpy.ndarray
        Selected per-element charges with shape ``(n_elem,)``.
    step : int or None, default None
        History batch step used for the charge snapshot.
    """

    directory: Path
    mesh_ids: tuple[int, ...]
    elem_indices: np.ndarray
    triangles: np.ndarray
    charges: np.ndarray
    step: int | None = None


@dataclass(frozen=True)
class CoulombInteraction:
    """Coulomb force/torque summary between two mesh groups.

    Attributes
    ----------
    group_a_mesh_ids : tuple of int
        Mesh ids used as the target group.
    group_b_mesh_ids : tuple of int
        Mesh ids used as the source group.
    step : int or None
        History step used for the interaction snapshot.
    softening : float
        Softening length used in the force computation.
    torque_origin_m : numpy.ndarray
        Torque origin in meters with shape ``(3,)``.
    force_on_a_N : numpy.ndarray
        Net force on group A in newtons with shape ``(3,)``.
    force_on_b_N : numpy.ndarray
        Net force on group B in newtons with shape ``(3,)``.
    torque_on_a_Nm : numpy.ndarray
        Net torque on group A in newton-meters with shape ``(3,)``.
    torque_on_b_Nm : numpy.ndarray
        Net torque on group B in newton-meters with shape ``(3,)``.
    mean_force_on_a_per_element_N : numpy.ndarray
        Mean force per target element in newtons with shape ``(3,)``.
    mean_torque_on_a_per_element_Nm : numpy.ndarray
        Mean torque per target element in newton-meters with shape ``(3,)``.
    """

    group_a_mesh_ids: tuple[int, ...]
    group_b_mesh_ids: tuple[int, ...]
    step: int | None
    softening: float
    torque_origin_m: np.ndarray
    force_on_a_N: np.ndarray
    force_on_b_N: np.ndarray
    torque_on_a_Nm: np.ndarray
    torque_on_b_Nm: np.ndarray
    mean_force_on_a_per_element_N: np.ndarray
    mean_torque_on_a_per_element_Nm: np.ndarray


@dataclass(frozen=True)
class FortranRunResult:
    """Container for one Fortran simulation output directory.

    Attributes
    ----------
    directory : pathlib.Path
        Output directory path.
    mesh_nelem : int
        Number of mesh elements.
    processed_particles : int
        Number of processed particles.
    absorbed : int
        Number of absorbed particles.
    escaped : int
        Number of escaped particles.
    batches : int
        Number of processed batches.
    escaped_boundary : int
        Number of particles escaped by boundary condition.
    survived_max_step : int
        Number of particles that reached max-step limit.
    last_rel_change : float
        Last relative charge-change metric.
    charges : numpy.ndarray
        Final per-element charge array with shape ``(mesh_nelem,)``.
    triangles : numpy.ndarray or None, default None
        Triangle vertices with shape ``(mesh_nelem, 3, 3)``.
    mesh_ids : numpy.ndarray or None, default None
        Per-element mesh id array with shape ``(mesh_nelem,)``.
    mesh_sources : dict[int, MeshSource] or None, default None
        Mesh-source metadata indexed by mesh id.
    history : FortranChargeHistory or None, default None
        Lazy charge-history accessor.
    """

    directory: Path
    mesh_nelem: int
    processed_particles: int
    absorbed: int
    escaped: int
    batches: int
    escaped_boundary: int
    survived_max_step: int
    last_rel_change: float
    charges: np.ndarray
    triangles: np.ndarray | None = None
    mesh_ids: np.ndarray | None = None
    mesh_sources: dict[int, MeshSource] | None = None
    history: FortranChargeHistory | None = None

    def history_at(self, step: int = -1) -> np.ndarray:
        """Return per-element charges at one history batch step.

        Parameters
        ----------
        step : int, default -1
            Batch step to read. ``-1`` selects the latest history step.

        Returns
        -------
        numpy.ndarray
            Per-element charge array with shape ``(mesh_nelem,)``.

        Raises
        ------
        ValueError
            If history data is missing or empty.
        """

        if self.history is None or not self.history.has_data:
            raise ValueError(
                "charge_history.csv is not found or empty. Enable history output and rerun."
            )
        return self.history.get_step(step)


@dataclass(frozen=True)
class PotentialSlice2D:
    """Electric potential sampled on one axis-aligned 2D slice.

    Attributes
    ----------
    plane : str
        Plane name (``"xy"``, ``"yz"``, or ``"xz"``).
    axis_u : str
        First in-plane axis label.
    axis_v : str
        Second in-plane axis label.
    fixed_axis : str
        Axis held constant for the slice.
    fixed_value_m : float
        Fixed-axis coordinate in meters.
    u_values_m : numpy.ndarray
        Sample coordinates along ``axis_u``.
    v_values_m : numpy.ndarray
        Sample coordinates along ``axis_v``.
    potential_V : numpy.ndarray
        Potential grid in volts with shape ``(n_v, n_u)``.
    """

    plane: str
    axis_u: str
    axis_v: str
    fixed_axis: str
    fixed_value_m: float
    u_values_m: np.ndarray
    v_values_m: np.ndarray
    potential_V: np.ndarray
