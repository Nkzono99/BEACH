"""Output fixtures and the deterministic panel kernel used by result tests."""

from pathlib import Path
import numpy as np
import pytest
import beach.analysis.coulomb as coulomb_module
import beach.analysis.field_kernel as field_kernel_module
import beach.analysis.potential as potential_module
from beach import FortranChargeHistory, K_COULOMB
from beach.analysis.panel_quadrature import panel_target_quadrature


class _UnitPanelKernel:
    """Small deterministic stand-in; native panel accuracy is tested separately."""

    def __init__(self, triangles, charges, *, options, **_kwargs) -> None:
        panels = np.asarray(triangles, dtype=float)
        self.centers = panels.mean(axis=1)
        (
            self.source_quadrature,
            self.source_quadrature_weights,
            self.source_element_index,
        ) = panel_target_quadrature(
            panels,
            np.ones(panels.shape[0], dtype=float),
            order=7,
        )
        self.charges = np.asarray(charges, dtype=float).copy()
        self.options = options
        self.periodic2 = options.periodic2
        self.closed = False

    def __enter__(self):
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def close(self) -> None:
        self.closed = True

    @property
    def source_charges(self) -> np.ndarray:
        return self.charges.copy()

    def update_charges(self, charges: np.ndarray) -> None:
        self.charges = np.asarray(charges, dtype=float).copy()

    def _wrapped_targets(self, points: np.ndarray) -> np.ndarray:
        targets = np.asarray(points, dtype=float)
        if self.periodic2 is None:
            return targets
        axes, lengths, origins, *_rest = self.periodic2
        targets = targets.copy()
        for axis, length, origin in zip(axes, lengths, origins):
            targets[:, axis] = origin + np.mod(targets[:, axis] - origin, length)
        return targets

    def eval_phi(self, points: np.ndarray) -> np.ndarray:
        targets = self._wrapped_targets(points)
        values = np.zeros(targets.shape[0], dtype=float)
        shifts = [np.zeros(3)]
        if self.periodic2 is not None:
            axes, lengths, _origins, layers, *_rest = self.periodic2
            shifts = []
            for first in range(-layers, layers + 1):
                for second in range(-layers, layers + 1):
                    shift = np.zeros(3)
                    shift[axes[0]] = first * lengths[0]
                    shift[axes[1]] = second * lengths[1]
                    shifts.append(shift)
        for shift in shifts:
            delta = targets[:, None, :] - (self.source_quadrature + shift)[None, :, :]
            distance = np.linalg.norm(delta, axis=2)
            inverse_distance = np.zeros_like(distance)
            np.divide(1.0, distance, out=inverse_distance, where=distance > 0.0)
            sample_charges = (
                self.charges[self.source_element_index] * self.source_quadrature_weights
            )
            values += K_COULOMB * np.sum(
                inverse_distance * sample_charges[None, :],
                axis=1,
            )
        return values

    def eval_e(self, points: np.ndarray) -> np.ndarray:
        targets = self._wrapped_targets(points)
        field = np.zeros_like(targets)
        shifts = [np.zeros(3)]
        if self.periodic2 is not None:
            axes, lengths, _origins, layers, *_rest = self.periodic2
            shifts = []
            for first in range(-layers, layers + 1):
                for second in range(-layers, layers + 1):
                    shift = np.zeros(3)
                    shift[axes[0]] = first * lengths[0]
                    shift[axes[1]] = second * lengths[1]
                    shifts.append(shift)
        for shift in shifts:
            delta = targets[:, None, :] - (self.source_quadrature + shift)[None, :, :]
            distance = np.linalg.norm(delta, axis=2)
            inverse_cube = np.zeros_like(distance)
            np.divide(1.0, distance**3, out=inverse_cube, where=distance > 0.0)
            sample_charges = (
                self.charges[self.source_element_index] * self.source_quadrature_weights
            )
            field += K_COULOMB * np.sum(
                sample_charges[None, :, None] * delta * inverse_cube[:, :, None],
                axis=1,
            )
        return field

    def force_on_charges(
        self,
        points: np.ndarray,
        charges: np.ndarray,
        *,
        origin: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        target_points = np.asarray(points, dtype=float)
        force_samples = np.asarray(charges, dtype=float)[:, None] * self.eval_e(
            target_points
        )
        return (
            np.sum(force_samples, axis=0),
            np.sum(
                np.cross(target_points - np.asarray(origin)[None, :], force_samples),
                axis=0,
            ),
        )


def _write_complete_free_field_config(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "beach.toml").write_text(
        '[sim]\nfield_bc_mode = "free"\n\n[external_boundary.field]\nmodel = "none"\n',
        encoding="utf-8",
    )


@pytest.fixture(autouse=True)
def _unit_panel_kernel(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(potential_module, "FieldKernel", _UnitPanelKernel)
    monkeypatch.setattr(coulomb_module, "FieldKernel", _UnitPanelKernel)
    monkeypatch.setattr(field_kernel_module, "FieldKernel", _UnitPanelKernel)


def _make_history(
    *,
    mesh_nelem: int,
    history: np.ndarray,
    batch_indices: np.ndarray | None = None,
    processed_particles_by_batch: np.ndarray | None = None,
    rel_change_by_batch: np.ndarray | None = None,
) -> FortranChargeHistory:
    n_snapshots = int(history.shape[1])
    batches = (
        np.arange(1, n_snapshots + 1, dtype=np.int64)
        if batch_indices is None
        else batch_indices.astype(np.int64, copy=False)
    )
    processed = (
        np.zeros(n_snapshots, dtype=np.int64)
        if processed_particles_by_batch is None
        else processed_particles_by_batch.astype(np.int64, copy=False)
    )
    rel = (
        np.zeros(n_snapshots, dtype=float)
        if rel_change_by_batch is None
        else rel_change_by_batch.astype(float, copy=False)
    )
    return FortranChargeHistory.from_arrays(
        Path("dummy_charge_history.csv"),
        mesh_nelem=mesh_nelem,
        history=history,
        processed_particles_by_batch=processed,
        rel_change_by_batch=rel,
        batch_indices=batches,
    )


def _write_charge_history(
    path: Path,
    rows: list[str],
    *,
    mesh_nelem: int = 3,
) -> FortranChargeHistory:
    path.write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )
    return FortranChargeHistory(path, mesh_nelem=mesh_nelem)


def _write_three_mesh_fixture(out: Path) -> None:
    _write_complete_free_field_config(out)
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=3",
                "processed_particles=12",
                "absorbed=9",
                "escaped=3",
                "batches=10",
                "last_rel_change=1.0e-6",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-9\n2,2.0e-9\n3,-3.0e-9\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,1.0e-9,1\n"
        "2,1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,2.0e-9,2\n"
        "3,2.0,0.0,0.0,2.0,1.0,0.0,2.0,0.0,1.0,-3.0e-9,3\n",
        encoding="utf-8",
    )
    (out / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,elem_count\n"
        "1,template,plane,1\n"
        "2,template,box,1\n"
        "3,template,sphere,1\n",
        encoding="utf-8",
    )
    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        "1,2,1.0e-1,1,5.0e-10\n"
        "1,2,1.0e-1,2,1.0e-9\n"
        "1,2,1.0e-1,3,-1.5e-9\n"
        "10,12,1.0e-6,1,1.0e-9\n"
        "10,12,1.0e-6,2,2.0e-9\n"
        "10,12,1.0e-6,3,-3.0e-9\n",
        encoding="utf-8",
    )


def _write_coulomb_matrix_fixture(out: Path) -> None:
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=3",
                "processed_particles=12",
                "absorbed=9",
                "escaped=3",
                "batches=10",
                "last_rel_change=1.0e-6",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-9\n2,2.0e-9\n3,-3.0e-9\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,1.0e-9,1\n"
        "2,1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,2.0e-9,2\n"
        "3,2.0,0.0,0.0,2.0,1.0,0.0,2.0,0.0,1.0,-3.0e-9,3\n",
        encoding="utf-8",
    )
    (out / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,elem_count\n"
        "1,template,plane,1\n"
        "2,template,sphere,1\n"
        "3,template,sphere,1\n",
        encoding="utf-8",
    )
    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        "1,2,1.0e-1,1,5.0e-10\n"
        "1,2,1.0e-1,2,1.0e-9\n"
        "1,2,1.0e-1,3,-1.5e-9\n"
        "10,12,1.0e-6,1,1.0e-9\n"
        "10,12,1.0e-6,2,2.0e-9\n"
        "10,12,1.0e-6,3,-3.0e-9\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[mesh]",
                'mode = "template"',
                "",
                "[[mesh.templates]]",
                'kind = "plane"',
                "enabled = true",
                "",
                "[[mesh.templates]]",
                'kind = "sphere"',
                "enabled = true",
                "",
                "[[mesh.templates]]",
                'kind = "sphere"',
                "enabled = true",
            ]
        ),
        encoding="utf-8",
    )


def _write_mobility_fixture(out: Path) -> None:
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=2",
                "absorbed=2",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-4\n2,2.0e-4\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,1.0e-4,1\n"
        "2,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,2.0e-4,2\n",
        encoding="utf-8",
    )
    (out / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,plane,insulator,1.0,1\n"
        "2,template,sphere,conductor,2.5,1\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[mesh]",
                'mode = "template"',
                "",
                "[[mesh.templates]]",
                'kind = "plane"',
                "enabled = true",
                "",
                "[[mesh.templates]]",
                'kind = "sphere"',
                "enabled = true",
                "radius = 0.1",
                "center = [0.3333333333, 0.3333333333, 1.0]",
                "n_lon = 8",
                "n_lat = 4",
            ]
        ),
        encoding="utf-8",
    )


def _write_minimal_result_fixture(
    out: Path,
    *,
    processed_particles: int = 10,
    absorbed: int = 7,
    escaped: int = 3,
    charges_text: str | None = None,
    mesh_triangles_text: str | None = None,
    mesh_sources_text: str | None = None,
    mesh_potential_text: str | None = None,
    summary_extra: list[str] | None = None,
    field_source_model: str | None = "triangle_p0",
) -> None:
    summary_lines = [
        "mesh_nelem=2",
        f"processed_particles={processed_particles}",
        f"absorbed={absorbed}",
        f"escaped={escaped}",
        "batches=1",
        "last_rel_change=1.0e-8",
    ]
    if field_source_model is not None:
        summary_lines.append(f"field_source_model={field_source_model}")
    if summary_extra is not None:
        summary_lines.extend(summary_extra)
    (out / "summary.txt").write_text("\n".join(summary_lines), encoding="utf-8")
    (out / "charges.csv").write_text(
        charges_text or "elem_idx,charge_C\n1,1.0e-10\n2,-2.0e-10\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        mesh_triangles_text
        or (
            "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
            "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
            "2,0,0,1,1,0,1,0,1,1,-2.0e-10,2\n"
        ),
        encoding="utf-8",
    )
    if mesh_sources_text is not None:
        (out / "mesh_sources.csv").write_text(mesh_sources_text, encoding="utf-8")
    if mesh_potential_text is not None:
        (out / "mesh_potential.csv").write_text(mesh_potential_text, encoding="utf-8")


def _matching_plane_summary_extra(
    *,
    photoelectron_escape_flux_m2_s: float = 1.5,
) -> list[str]:
    return [
        "matching_plane_state_valid=T",
        "matching_plane_displacement_C_m2=-2.0e-8",
        "matching_plane_phi_V=1.5",
        "matching_plane_electron_inward_flux_m2_s=4.0",
        "matching_plane_ion_inward_flux_m2_s=5.0",
        "matching_plane_electron_access_potential_V=-1.0",
        "matching_plane_ion_access_potential_V=2.0",
        "matching_plane_photoelectron_barrier_potential_V=3.0",
        "matching_plane_photoelectron_outward_flux_m2_s=2.0",
        "matching_plane_photoelectron_mean_normal_energy_eV=6.0",
        "matching_plane_electron_outward_flux_m2_s=7.0",
        "matching_plane_ion_outward_flux_m2_s=8.0",
        "matching_plane_photoelectron_return_flux_m2_s=0.5",
        f"matching_plane_photoelectron_escape_flux_m2_s={photoelectron_escape_flux_m2_s}",
        "matching_plane_iterations=3",
        "matching_plane_residual=1.0e-3",
    ]


def _matching_history_row(
    *,
    batch: int = 1,
    electron_inward_flux: str = "4.0",
    photoelectron_outward_flux: str = "2.0",
    photoelectron_return_flux: str = "0.5",
    photoelectron_escape_flux: str = "1.5",
) -> str:
    return ",".join(
        (
            str(batch),
            "1.0e-6",
            "-2.0e-8",
            "1.5",
            electron_inward_flux,
            "5.0",
            "-1.0",
            "2.0",
            "3.0",
            photoelectron_outward_flux,
            "6.0",
            "7.0",
            "8.0",
            photoelectron_return_flux,
            photoelectron_escape_flux,
            "3",
            "1.0e-3",
        )
    )


def _charge_ledger_row(
    *,
    species_idx: str = "1",
    injected_charge: str = "-3.0",
) -> str:
    return ",".join(
        (
            "1",
            species_idx,
            injected_charge,
            "0.0",
            "-2.0",
            "-1.0",
            "-0.1",
            "3",
            "0",
            "2",
            "1",
            "1",
        )
    )
