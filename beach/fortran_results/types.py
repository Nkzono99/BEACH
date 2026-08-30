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
    surface_model : str
        Surface response model for this mesh group.
    epsilon_r : float
        Relative permittivity associated with this mesh group.
    """

    mesh_id: int
    source_kind: str
    template_kind: str
    elem_count: int
    surface_model: str = "insulator"
    epsilon_r: float = 1.0


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
    torque_origin_m: np.ndarray
    force_on_a_N: np.ndarray
    force_on_b_N: np.ndarray
    torque_on_a_Nm: np.ndarray
    torque_on_b_Nm: np.ndarray
    mean_force_on_a_per_element_N: np.ndarray
    mean_torque_on_a_per_element_Nm: np.ndarray


@dataclass(frozen=True)
class CoulombMobilityRecord:
    """Per-object Coulomb mobility summary.

    Attributes
    ----------
    mesh_id : int
        Mesh identifier of the analyzed object.
    label : str
        Human-readable object label.
    kind : str
        Object kind resolved from ``beach.toml`` or ``mesh_sources.csv``.
    step : int or None
        History step used for the analysis snapshot.
    center_m : numpy.ndarray
        Object centroid in meters with shape ``(3,)``.
    force_N : numpy.ndarray
        Net Coulomb force in newtons with shape ``(3,)``.
    torque_Nm : numpy.ndarray
        Net Coulomb torque about object center in newton-meters with shape ``(3,)``.
    force_normal_N : float
        Force component along support normal. Positive means lift direction.
    force_tangent_N : float
        Tangential-force magnitude relative to the support normal.
    torque_normal_Nm : float
        Torque component along support normal.
    torque_tangent_Nm : float
        Tangential-torque magnitude relative to the support normal.
    mass_kg : float or None
        Estimated object mass. ``None`` when geometry or density is unavailable.
    characteristic_radius_m : float or None
        Characteristic radius used for rolling metrics when available.
    weight_support_N : float or None
        Weight component resisted by the support normal.
    resisting_normal_N : float or None
        Baseline normal resistance used for lift checks
        (typically weight + adhesion).
    effective_normal_load_N : float or None
        Remaining compressive normal load after Coulomb lift is applied.
    lift_ratio : float or None
        ``max(force_normal, 0) / resisting_normal`` when available.
    slide_ratio : float or None
        ``force_tangent / (mu_static * effective_normal_load)`` when available.
    roll_ratio : float or None
        ``torque_tangent / (mu_roll * effective_normal_load * radius)`` when available.
    notes : tuple of str
        Diagnostic notes about unavailable assumptions or approximations.
    """

    mesh_id: int
    label: str
    kind: str
    step: int | None
    center_m: np.ndarray
    force_N: np.ndarray
    torque_Nm: np.ndarray
    force_normal_N: float
    force_tangent_N: float
    torque_normal_Nm: float
    torque_tangent_Nm: float
    mass_kg: float | None
    characteristic_radius_m: float | None
    weight_support_N: float | None
    resisting_normal_N: float | None
    effective_normal_load_N: float | None
    lift_ratio: float | None
    slide_ratio: float | None
    roll_ratio: float | None
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class CoulombMobilityAnalysis:
    """Collection of per-object Coulomb mobility records.

    Attributes
    ----------
    step : int or None
        History step used for the analysis snapshot.
    gravity_m_s2 : numpy.ndarray
        Gravity vector in meters per second squared with shape ``(3,)``.
    support_normal_m : numpy.ndarray
        Unit support normal with shape ``(3,)``.
    support_kinds : tuple of str
        Object kinds treated as supports and excluded by default from targets.
    density_kg_m3 : float or None
        Bulk density used for mass estimation.
    mu_static : float or None
        Static-friction coefficient used for slide checks.
    mu_roll : float or None
        Rolling-resistance coefficient used for roll checks.
    adhesion_force_N : float
        Additional normal adhesion resisting lift.
    records : tuple of CoulombMobilityRecord
        Per-object analysis records.
    """

    step: int | None
    gravity_m_s2: np.ndarray
    support_normal_m: np.ndarray
    support_kinds: tuple[str, ...]
    density_kg_m3: float | None
    mu_static: float | None
    mu_roll: float | None
    adhesion_force_N: float
    records: tuple[CoulombMobilityRecord, ...]


@dataclass(frozen=True)
class ChargeLedgerEntry:
    """One species row from ``charge_ledger.csv``."""

    batch: int
    species_idx: int
    injected_from_remote_c: float
    emitted_from_surface_c: float
    absorbed_on_surface_c: float
    escaped_to_infinity_c: float
    discarded_unresolved_c: float
    injected_count: int
    emitted_count: int
    absorbed_count: int
    escaped_count: int
    discarded_unresolved_count: int
    neutral_return_correction_c: float = 0.0
    neutral_return_weight_scale: float = 1.0
    neutral_return_unresolved_fraction: float = 0.0
    fixed_absorbed_target_charge_c: float = 0.0
    fixed_absorbed_applied_charge_c: float = 0.0
    fixed_absorbed_weight_scale: float = 1.0
    fixed_emission_target_charge_c: float = 0.0
    fixed_emission_applied_charge_c: float = 0.0
    fixed_emission_weight_scale: float = 1.0
    fixed_escape_target_charge_c: float = 0.0
    fixed_escape_applied_charge_c: float = 0.0
    fixed_escape_correction_c: float = 0.0
    fixed_current_correction_c: float = 0.0


@dataclass(frozen=True)
class FieldReconstructionReceipt:
    """Resolved simulator settings required to reconstruct its electric field."""

    schema_version: int
    resolved_field_solver: str
    fmm_expansion_order: int
    field_bc_mode: str
    tree_theta: float
    tree_leaf_max: int
    external_e0_v_m: tuple[float, float, float]
    use_box: bool
    box_min_m: tuple[float, float, float]
    box_max_m: tuple[float, float, float]
    boundary_low: tuple[int, int, int]
    boundary_high: tuple[int, int, int]
    periodic_image_layers: int
    periodic_far_correction: str
    periodic_nonzero_mode_backend: str
    periodic_zero_mode_policy: str
    periodic_lower_boundary_model: str
    periodic_reference_mode_layers: int
    periodic_panel_quadrature_order: int
    periodic_ewald_alpha: float
    periodic_ewald_layers: int
    periodic_cache_dir: str
    periodic_generation_tolerance: float

    def as_legacy_sim_mapping(self) -> dict[str, object]:
        """Return the combined simulator mapping consumed by field utilities."""

        sim: dict[str, object] = {
            "field_solver": self.resolved_field_solver,
            "field_bc_mode": self.field_bc_mode,
            "tree_theta": self.tree_theta,
            "tree_leaf_max": self.tree_leaf_max,
            "e0": self.external_e0_v_m,
            "field_periodic_image_layers": self.periodic_image_layers,
            "field_periodic_far_correction": self.periodic_far_correction,
            "field_periodic_ewald_alpha": self.periodic_ewald_alpha,
            "field_periodic_ewald_layers": self.periodic_ewald_layers,
            "field_periodic_cache_dir": self.periodic_cache_dir,
            "field_periodic_generation_tolerance": (
                self.periodic_generation_tolerance
            ),
        }
        if self.use_box:
            sim["box_min"] = self.box_min_m
            sim["box_max"] = self.box_max_m
        for axis, name in enumerate(("x", "y", "z")):
            sim[f"bc_{name}_low"] = (
                "periodic" if self.boundary_low[axis] == 2 else "open"
            )
            sim[f"bc_{name}_high"] = (
                "periodic" if self.boundary_high[axis] == 2 else "open"
            )
        return sim


@dataclass(frozen=True)
class MatchingPlaneState:
    """One accepted quasistatic matching-plane coupling receipt."""

    displacement_c_m2: float
    phi_v: float
    electron_inward_flux_m2_s: float
    ion_inward_flux_m2_s: float
    electron_access_potential_v: float
    ion_access_potential_v: float
    photoelectron_barrier_potential_v: float
    photoelectron_outward_flux_m2_s: float
    photoelectron_mean_normal_energy_ev: float
    electron_outward_flux_m2_s: float
    ion_outward_flux_m2_s: float
    photoelectron_return_flux_m2_s: float
    photoelectron_escape_flux_m2_s: float
    iterations: int
    residual: float


@dataclass(frozen=True)
class MatchingPlaneHistoryEntry:
    """Matching-plane receipt selected by ``output.history_stride``."""

    batch: int
    simulated_time_s: float
    state: MatchingPlaneState


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
    multiple_box_events_retry_attempted : int, default 0
        Number of failed boundary steps replayed by the configured retry backend.
    multiple_box_events_retry_resolved : int, default 0
        Number of replayed boundary steps completed successfully.
    multiple_box_events_soft_discard_fraction : float, default 0.0
        Derived fraction of processed particles removed by the soft-discard fallback.
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
    mesh_potential_v : numpy.ndarray or None, default None
        Optional per-element centroid potential loaded from ``mesh_potential.csv``.
    history : FortranChargeHistory or None, default None
        Lazy charge-history accessor.
    simulated_time_s : float or None, default None
        Sum of accepted batch durations reported by the runtime.
    adaptive_nonzero_mode_rejected_trials : int, default 0
        Cumulative rejected adaptive ``k != 0`` trial count.
    adaptive_nonzero_mode_last_batch_duration_s : float or None, default None
        Duration of the last accepted adaptive trial.
    adaptive_nonzero_mode_last_potential_step_v : float or None, default None
        Last accepted maximum ``k != 0`` potential change.
    adaptive_nonzero_mode_omp_threads : int or None, default None
        OpenMP team size recorded and enforced for an adaptive restart.
    periodic2_max_nonzero_mode_potential_step_v : float or None, default None
        Configured adaptive ``k != 0`` potential-change limit.
    field_reconstruction : FieldReconstructionReceipt or None, default None
        Resolved field-model receipt embedded in new ``summary.txt`` outputs.
    matching_plane_state : MatchingPlaneState or None, default None
        Last accepted matching-plane state when the summary marks it valid.
    matching_plane_history : tuple[MatchingPlaneHistoryEntry, ...] or None, default None
        Matching-plane history. A present header-only CSV produces an empty tuple.
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
    mesh_potential_v: np.ndarray | None = None
    history: FortranChargeHistory | None = None
    checkpoint_schema_version: int | None = None
    model_fingerprint: str | None = None
    mesh_fingerprint: str | None = None
    species_fingerprint: str | None = None
    charge_ledger_residual_c: float | None = None
    charge_ledger: tuple[ChargeLedgerEntry, ...] | None = None
    field_source_model: str = "triangle_p0"
    field_kernel_id: str | None = None
    periodic2_cache_hit: bool | None = None
    periodic2_operator_build_count: int | None = None
    periodic2_cache_fingerprint: str | None = None
    periodic2_cache_path: str | None = None
    periodic2_generation_tolerance: float | None = None
    multiple_box_events_retry_attempted: int = 0
    multiple_box_events_retry_resolved: int = 0
    multiple_box_events_soft_discarded: int = 0
    multiple_box_events_soft_discarded_abs_charge_c: float = 0.0
    simulated_time_s: float | None = None
    adaptive_nonzero_mode_rejected_trials: int = 0
    adaptive_nonzero_mode_last_batch_duration_s: float | None = None
    adaptive_nonzero_mode_last_potential_step_v: float | None = None
    adaptive_nonzero_mode_omp_threads: int | None = None
    periodic2_max_nonzero_mode_potential_step_v: float | None = None
    field_reconstruction: FieldReconstructionReceipt | None = None
    matching_plane_state: MatchingPlaneState | None = None
    matching_plane_history: tuple[MatchingPlaneHistoryEntry, ...] | None = None
    multiple_box_events_soft_discard_fraction: float = 0.0

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
