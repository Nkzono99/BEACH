from dataclasses import fields
from pathlib import Path

import numpy as np
import pytest

import beach.fortran_results.coulomb as coulomb_module
import beach.fortran_results.potential as potential_module
from beach.fortran_results import (
    analyze_coulomb_mobility,
    Beach,
    calc_coulomb,
    CoulombMobilityAnalysis,
    FortranChargeHistory,
    K_COULOMB,
    FortranRunResult,
    _select_frame_columns,
    animate_history_mesh,
    compute_potential_points,
    compute_potential_slices,
    compute_potential_mesh,
    compute_electric_field_points,
    _surface_charge_density,
    list_fortran_runs,
    load_fortran_result,
    plot_mesh_source_boxplot,
    plot_charge_mesh,
    plot_coulomb_force_matrix,
    plot_potential_slices,
    plot_potential_mesh,
)
from beach.fortran_results.mesh import (
    _wrap_periodic2_triangles_by_centroid,
    _wrap_periodic2_triangles_by_mesh_centroid,
)
from beach.fortran_results.objects import resolve_object_specs
from beach.fortran_results.panel_quadrature import panel_target_quadrature
from beach.fortran_results.plotting import _periodic2_for_coulomb_matrix
from beach.fortran_results.potential import (
    _auto_periodic2_from_result,
    _coerce_periodic2,
    _potential_history,
)


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
        self.periodic2 = options.periodic2

    def __enter__(self):
        return self

    def __exit__(self, *_args: object) -> None:
        pass

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
            delta = (
                targets[:, None, :]
                - (self.source_quadrature + shift)[None, :, :]
            )
            distance = np.linalg.norm(delta, axis=2)
            inverse_distance = np.zeros_like(distance)
            np.divide(1.0, distance, out=inverse_distance, where=distance > 0.0)
            sample_charges = (
                self.charges[self.source_element_index]
                * self.source_quadrature_weights
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
            delta = (
                targets[:, None, :]
                - (self.source_quadrature + shift)[None, :, :]
            )
            distance = np.linalg.norm(delta, axis=2)
            inverse_cube = np.zeros_like(distance)
            np.divide(1.0, distance**3, out=inverse_cube, where=distance > 0.0)
            sample_charges = (
                self.charges[self.source_element_index]
                * self.source_quadrature_weights
            )
            field += K_COULOMB * np.sum(
                sample_charges[None, :, None]
                * delta
                * inverse_cube[:, :, None],
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
        '[sim]\nfield_bc_mode = "free"\n\n'
        '[external_boundary.field]\nmodel = "none"\n',
        encoding="utf-8",
    )


@pytest.fixture(autouse=True)
def _unit_panel_kernel(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(potential_module, "FieldKernel", _UnitPanelKernel)
    monkeypatch.setattr(coulomb_module, "FieldKernel", _UnitPanelKernel)


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
    charges_text: str | None = None,
    mesh_triangles_text: str | None = None,
    mesh_sources_text: str | None = None,
    mesh_potential_text: str | None = None,
    summary_extra: list[str] | None = None,
    field_source_model: str | None = "triangle_p0",
) -> None:
    summary_lines = [
        "mesh_nelem=2",
        "processed_particles=10",
        "absorbed=7",
        "escaped=3",
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


def test_load_fortran_result(tmp_path: Path) -> None:
    out = tmp_path / "run1"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=10",
                "absorbed=7",
                "escaped=3",
                "batches=1",
                "escaped_boundary=2",
                "survived_max_step=1",
                "multiple_box_events_soft_discarded=4",
                "multiple_box_events_soft_discarded_abs_charge_C=2.5e-15",
                "last_rel_change=1.0e-8",
                "simulated_time_s=7.5e-6",
                "adaptive_nonzero_mode_rejected_trials=3",
                "adaptive_nonzero_mode_last_batch_duration_s=2.5e-6",
                "adaptive_nonzero_mode_last_potential_step_V=8.0e-3",
                "adaptive_nonzero_mode_omp_threads=16",
                "periodic2_max_nonzero_mode_potential_step_V=1.0e-2",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-10\n2,-2.0e-10\n", encoding="utf-8"
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
        "2,0,0,1,1,0,1,0,1,1,-2.0e-10,2\n",
        encoding="utf-8",
    )
    (out / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,plane,insulator,1.0,1\n"
        "2,template,sphere,conductor,2.5,1\n",
        encoding="utf-8",
    )
    (out / "mesh_potential.csv").write_text(
        "elem_idx,potential_V\n1,1.5\n2,-2.5\n",
        encoding="utf-8",
    )

    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        "1,5,3.0e-1,1,2.0e-11\n"
        "1,5,3.0e-1,2,-1.0e-11\n"
        "3,10,1.0e-8,1,1.0e-10\n"
        "3,10,1.0e-8,2,-2.0e-10\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.mesh_nelem == 2
    assert result.absorbed == 7
    assert result.escaped_boundary == 2
    assert result.survived_max_step == 1
    assert result.multiple_box_events_soft_discarded == 4
    assert result.multiple_box_events_soft_discarded_abs_charge_c == 2.5e-15
    assert result.simulated_time_s == 7.5e-6
    assert result.adaptive_nonzero_mode_rejected_trials == 3
    assert result.adaptive_nonzero_mode_last_batch_duration_s == 2.5e-6
    assert result.adaptive_nonzero_mode_last_potential_step_v == 8.0e-3
    assert result.adaptive_nonzero_mode_omp_threads == 16
    assert result.periodic2_max_nonzero_mode_potential_step_v == 1.0e-2
    assert result.triangles is not None
    assert result.triangles.shape == (2, 3, 3)
    assert result.history is not None
    assert result.mesh_ids is not None
    np.testing.assert_array_equal(result.mesh_ids, np.array([1, 2]))
    assert result.mesh_sources is not None
    assert result.mesh_sources[1].template_kind == "plane"
    assert result.mesh_sources[1].surface_model == "insulator"
    assert result.mesh_sources[1].epsilon_r == 1.0
    assert result.mesh_sources[2].surface_model == "conductor"
    assert result.mesh_sources[2].epsilon_r == 2.5
    np.testing.assert_allclose(result.charges, np.array([1.0e-10, -2.0e-10]))
    assert result.mesh_potential_v is not None
    np.testing.assert_allclose(result.mesh_potential_v, np.array([1.5, -2.5]))
    np.testing.assert_allclose(
        result.history.as_array(),
        np.array([[2.0e-11, 1.0e-10], [-1.0e-11, -2.0e-10]]),
    )
    np.testing.assert_array_equal(
        result.history.processed_particles_by_batch, np.array([5, 10])
    )
    np.testing.assert_allclose(
        result.history.rel_change_by_batch, np.array([3.0e-1, 1.0e-8])
    )
    np.testing.assert_array_equal(result.history.batch_indices, np.array([1, 3]))


def test_panel_estimators_reject_removed_point_output(tmp_path: Path) -> None:
    out = tmp_path / "run_removed_point"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "field_kernel_id=softened_point",
        ],
        field_source_model="point",
    )
    result = load_fortran_result(out)

    assert result.field_source_model == "point"
    with pytest.raises(ValueError, match="triangle_p0"):
        compute_potential_mesh(result)
    with pytest.raises(ValueError, match="triangle_p0"):
        compute_electric_field_points(result, np.zeros((1, 3)))
    with pytest.raises(ValueError, match="triangle_p0"):
        calc_coulomb(result, 1, 2, step=None)


def test_panel_estimators_reject_output_without_source_model_receipt(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_missing_source_model"
    out.mkdir()
    _write_minimal_result_fixture(out, field_source_model=None)

    result = load_fortran_result(out)

    assert result.field_source_model == "unknown"
    with pytest.raises(ValueError, match="triangle_p0"):
        compute_potential_mesh(result)
    with pytest.raises(ValueError, match="triangle_p0"):
        compute_electric_field_points(result, np.zeros((1, 3)))
    with pytest.raises(ValueError, match="triangle_p0"):
        calc_coulomb(result, 1, 2, step=None)


def test_load_fortran_result_rejects_charges_row_count_mismatch(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_charges_count"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        charges_text="elem_idx,charge_C\n1,1.0e-10\n",
    )

    with pytest.raises(ValueError, match="charges.csv row count"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_charges_invalid_elem_idx(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_charges_elem"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        charges_text="elem_idx,charge_C\n1,1.0e-10\n3,-2.0e-10\n",
    )

    with pytest.raises(ValueError, match="charges.csv elem_idx"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_nonfinite_charge(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_charge_nan"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        charges_text="elem_idx,charge_C\n1,1.0e-10\n2,nan\n",
    )

    with pytest.raises(ValueError, match="charge_C"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_mesh_triangle_row_count_mismatch(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_bad_triangles_count"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_triangles_text=(
            "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
            "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
        ),
    )

    with pytest.raises(ValueError, match="mesh_triangles.csv row count"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_nonfinite_triangle_vertex(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_triangle_nan"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_triangles_text=(
            "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
            "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
            "2,0,0,1,1,nan,1,0,1,1,-2.0e-10,2\n"
        ),
    )

    with pytest.raises(ValueError, match="vertex values"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_bad_mesh_source_material(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_mesh_source_material"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_sources_text=(
            "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
            "1,template,plane,insulator,1.0,1\n"
            "2,template,sphere,dielectric,nan,1\n"
        ),
    )

    with pytest.raises(ValueError, match="epsilon_r"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_nonfinite_mesh_potential(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_mesh_potential"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_potential_text="elem_idx,potential_V\n1,1.5\n2,inf\n",
    )

    with pytest.raises(ValueError, match="potential_V"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_negative_summary_stat(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_summary"
    out.mkdir()
    _write_minimal_result_fixture(out, summary_extra=["escaped_boundary=-1"])

    with pytest.raises(ValueError, match="escaped_boundary"):
        load_fortran_result(out)


def test_resolve_object_specs_prefers_mesh_source_materials(tmp_path: Path) -> None:
    out = tmp_path / "run_materials"
    out.mkdir()
    _write_mobility_fixture(out)

    result = load_fortran_result(out)
    specs = resolve_object_specs(result, config_path=out / "beach.toml")

    assert [spec.kind for spec in specs] == ["plane", "sphere"]
    assert [spec.surface_model for spec in specs] == ["insulator", "conductor"]
    assert [spec.epsilon_r for spec in specs] == [1.0, 2.5]


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*missing elem_idx=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,2,2.0e-9",
                "1,2,1.0e-1,2,9.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*duplicate elem_idx=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,2,2.0e-9",
                "1,2,1.0e-1,4,4.0e-9",
            ],
            r"batch=1.*out-of-range elem_idx=4",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,2,nan",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*elem_idx=2.*non-finite charge_C",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,3,1.0e-1,2,2.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*inconsistent processed_particles.*elem_idx=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,2.0e-1,2,2.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*inconsistent rel_change.*elem_idx=2",
        ),
        (
            [
                "2,2,1.0e-1,1,1.0e-9",
                "2,2,1.0e-1,2,2.0e-9",
                "2,2,1.0e-1,3,3.0e-9",
                "1,3,5.0e-2,1,4.0e-9",
                "1,3,5.0e-2,2,5.0e-9",
                "1,3,5.0e-2,3,6.0e-9",
            ],
            r"batch=1.*previous batch=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,nan,2,2.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*elem_idx=2.*non-finite rel_change",
        ),
    ],
    ids=[
        "missing-element",
        "duplicate-element",
        "out-of-range-element",
        "non-finite-charge",
        "inconsistent-processed-particles",
        "inconsistent-rel-change",
        "decreasing-batch",
        "non-finite-rel-change",
    ],
)
def test_charge_history_rejects_corrupt_batch_when_indexed(
    tmp_path: Path,
    rows: list[str],
    message: str,
) -> None:
    history = _write_charge_history(tmp_path / "charge_history.csv", rows)

    for _ in range(2):
        with pytest.raises(ValueError, match=message):
            _ = history.batch_indices


def test_load_fortran_result_defers_corrupt_history_validation(tmp_path: Path) -> None:
    out = tmp_path / "run_corrupt_lazy_history"
    out.mkdir()
    _write_minimal_result_fixture(out)
    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n1,2,1.0e-1,1,1.0e-9\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.history is not None
    with pytest.raises(ValueError, match=r"batch=1.*missing elem_idx=2"):
        _ = result.history.batch_indices


def test_charge_history_accepts_permuted_dense_rows(tmp_path: Path) -> None:
    history = _write_charge_history(
        tmp_path / "charge_history.csv",
        [
            "1,2,1.0e-1,3,3.0e-9",
            "1,2,1.0e-1,1,1.0e-9",
            "1,2,1.0e-1,2,2.0e-9",
        ],
    )

    np.testing.assert_array_equal(history.batch_indices, np.array([1]))
    np.testing.assert_allclose(history.get_step(1), np.array([1.0e-9, 2.0e-9, 3.0e-9]))


def test_load_fortran_result_history_supports_step_access(tmp_path: Path) -> None:
    out = tmp_path / "run_lazy"
    out.mkdir()
    _write_three_mesh_fixture(out)

    result = load_fortran_result(out)

    assert result.history is not None
    np.testing.assert_allclose(result.history[1], np.array([5.0e-10, 1.0e-9, -1.5e-9]))
    np.testing.assert_allclose(
        result.history_at(10), np.array([1.0e-9, 2.0e-9, -3.0e-9])
    )
    history = result.history.as_array()
    np.testing.assert_allclose(
        history,
        np.array(
            [
                [5.0e-10, 1.0e-9],
                [1.0e-9, 2.0e-9],
                [-1.5e-9, -3.0e-9],
            ]
        ),
    )
    np.testing.assert_array_equal(result.history.batch_indices, np.array([1, 10]))


def test_list_fortran_runs(tmp_path: Path) -> None:
    valid = tmp_path / "valid"
    valid.mkdir()
    (valid / "summary.txt").write_text(
        "mesh_nelem=1\nprocessed_particles=1\nabsorbed=1\nescaped=0\nbatches=1\nlast_rel_change=0.0\n",
        encoding="utf-8",
    )
    (valid / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    invalid = tmp_path / "invalid"
    invalid.mkdir()
    (invalid / "summary.txt").write_text("mesh_nelem=1\n", encoding="utf-8")

    runs = list_fortran_runs(tmp_path)
    assert runs == [valid]


def test_beach_uses_outputs_latest_by_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "outputs" / "latest"
    out.mkdir(parents=True)
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=1",
                "absorbed=1",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    monkeypatch.chdir(tmp_path)

    beach = Beach()
    assert beach.output_dir == Path("outputs/latest")
    assert beach.result.mesh_nelem == 1


def test_beach_default_lazy_history_keeps_step_access(tmp_path: Path) -> None:
    out = tmp_path / "run_beach_lazy"
    out.mkdir()
    _write_three_mesh_fixture(out)

    beach = Beach(out)
    result = beach.result

    assert result.history is not None
    np.testing.assert_allclose(beach.get_mesh_charge(2, step=10), np.array([2.0e-9]))
    np.testing.assert_allclose(
        result.history.get_step(10), np.array([1.0e-9, 2.0e-9, -3.0e-9])
    )


def test_load_fortran_result_without_new_summary_keys(tmp_path: Path) -> None:
    out = tmp_path / "legacy"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=4",
                "absorbed=2",
                "escaped=2",
                "batches=1",
                "last_rel_change=0.5",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    result = load_fortran_result(out)

    assert result.escaped_boundary == 0
    assert result.survived_max_step == 0
    assert result.checkpoint_schema_version is None
    assert result.charge_ledger is None
    assert result.simulated_time_s is None
    assert result.adaptive_nonzero_mode_rejected_trials == 0
    assert result.adaptive_nonzero_mode_last_batch_duration_s is None
    assert result.adaptive_nonzero_mode_last_potential_step_v is None
    assert result.adaptive_nonzero_mode_omp_threads is None
    assert result.periodic2_max_nonzero_mode_potential_step_v is None


@pytest.mark.parametrize(
    ("extra_line", "match"),
    [
        (
            None,
            "adaptive nonzero-mode output requires adaptive_nonzero_mode_omp_threads",
        ),
        (
            "adaptive_nonzero_mode_omp_threads=0",
            "adaptive_nonzero_mode_omp_threads must be > 0",
        ),
    ],
)
def test_load_fortran_result_rejects_incomplete_adaptive_summary(
    tmp_path: Path, extra_line: str | None, match: str
) -> None:
    out = tmp_path / "incomplete_adaptive"
    out.mkdir()
    lines = [
        "mesh_nelem=1",
        "processed_particles=4",
        "absorbed=2",
        "escaped=2",
        "batches=1",
        "last_rel_change=0.5",
        "field_source_model=triangle_p0",
        "periodic2_max_nonzero_mode_potential_step_V=1.0e-2",
        "simulated_time_s=1.0e-3",
        "adaptive_nonzero_mode_rejected_trials=1",
        "adaptive_nonzero_mode_last_batch_duration_s=5.0e-4",
        "adaptive_nonzero_mode_last_potential_step_V=8.0e-3",
    ]
    if extra_line is not None:
        lines.append(extra_line)
    (out / "summary.txt").write_text("\n".join(lines), encoding="utf-8")
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,0.0\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


def test_load_fortran_result_rejects_duplicate_summary_keys(tmp_path: Path) -> None:
    out = tmp_path / "duplicate_summary"
    out.mkdir()
    _write_minimal_result_fixture(out)
    summary_path = out / "summary.txt"
    summary_path.write_text(
        summary_path.read_text(encoding="utf-8") + "\nbatches=2\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate key 'batches'"):
        load_fortran_result(out)


def test_load_fortran_result_model_contract_and_charge_ledger(tmp_path: Path) -> None:
    out = tmp_path / "contract"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "checkpoint_schema_version=4",
            "model_fingerprint=0123456789ABCDEF",
            "mesh_fingerprint=1111111122222222",
            "species_fingerprint=3333333344444444",
            "charge_ledger_residual_C=1.0e-18",
            "periodic2_cache_hit=T",
            "periodic2_operator_build_count=0",
            "periodic2_cache_fingerprint=ABCDEF0123456789",
            "periodic2_cache_path=.beach_cache/periodic2/operator.bin",
            "periodic2_generation_tolerance=1.0e-8",
            "max_outer_flight_time_s=3.0e-6",
            "max_outer_frozen_field_ratio=3.0e-2",
            "max_outer_energy_relative_error=4.0e-5",
            "implicit_mean_last_returned_outer_flight_time_mean_s=2.5e-6",
            "implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2=7.5e-12",
            "coupling_outer_queue_enabled=T",
            "outer_photoelectron_population_fraction=1.25",
            "outer_photoelectron_column_per_area_m2=2.0e12",
            "outer_photoelectron_column_target_per_area_m2=2.1e12",
            "outer_photoelectron_column_residual_per_area_m2=-1.0e11",
            "outer_queue_event_count=7",
            "outer_queue_signed_charge_C=-3.5e-16",
            "outer_queue_fingerprint=FEDCBA9876543210",
        ],
    )
    (out / "charge_ledger.csv").write_text(
        "batch,species_idx,injected_from_remote_C,emitted_from_surface_C,"
        "absorbed_on_surface_C,escaped_to_infinity_C,discarded_unresolved_C,"
        "interface_outward_gross_C,interface_returned_gross_C,"
        "neutral_return_correction_C,neutral_return_weight_scale,"
        "neutral_return_unresolved_fraction,injected_count,"
        "emitted_count,absorbed_count,escaped_count,discarded_unresolved_count\n"
        "1,1,-3,0,-2,-1,-0.1,0,0,-0.1,1.05,0.05,3,0,2,1,1\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.checkpoint_schema_version == 4
    assert result.model_fingerprint == "0123456789ABCDEF"
    assert result.mesh_fingerprint == "1111111122222222"
    assert result.species_fingerprint == "3333333344444444"
    assert result.charge_ledger_residual_c == pytest.approx(1.0e-18)
    assert result.charge_ledger is not None
    assert result.charge_ledger[0].species_idx == 1
    assert result.charge_ledger[0].injected_from_remote_c == pytest.approx(-3.0)
    assert result.charge_ledger[0].escaped_count == 1
    assert result.charge_ledger[0].neutral_return_correction_c == pytest.approx(-0.1)
    assert result.charge_ledger[0].neutral_return_weight_scale == pytest.approx(1.05)
    assert result.charge_ledger[0].neutral_return_unresolved_fraction == pytest.approx(
        0.05
    )
    assert result.periodic2_cache_hit is True
    assert result.periodic2_operator_build_count == 0
    assert result.periodic2_cache_fingerprint == "ABCDEF0123456789"
    assert result.periodic2_cache_path == ".beach_cache/periodic2/operator.bin"
    assert result.periodic2_generation_tolerance == pytest.approx(1.0e-8)
    assert result.max_outer_flight_time_s == pytest.approx(3.0e-6)
    assert result.max_outer_frozen_field_ratio == pytest.approx(3.0e-2)
    assert result.max_outer_energy_relative_error == pytest.approx(4.0e-5)
    assert result.implicit_mean_last_returned_outer_flight_time_mean_s == pytest.approx(
        2.5e-6
    )
    assert (
        result.implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_c_m2
        == pytest.approx(7.5e-12)
    )
    assert result.coupling_outer_queue_enabled is True
    assert result.outer_photoelectron_population_fraction == pytest.approx(1.25)
    assert result.outer_photoelectron_column_per_area_m2 == pytest.approx(2.0e12)
    assert result.outer_photoelectron_column_target_per_area_m2 == pytest.approx(2.1e12)
    assert result.outer_photoelectron_column_residual_per_area_m2 == pytest.approx(
        -1.0e11
    )
    assert result.outer_queue_event_count == 7
    assert result.outer_queue_signed_charge_c == pytest.approx(-3.5e-16)
    assert result.outer_queue_fingerprint == "FEDCBA9876543210"


def test_fortran_run_result_appends_implicit_mean_diagnostics() -> None:
    field_names = tuple(field.name for field in fields(FortranRunResult))
    legacy_tail = field_names.index(
        "multiple_box_events_soft_discarded_abs_charge_c"
    )

    assert (
        field_names.index("implicit_mean_last_returned_outer_flight_time_mean_s")
        > legacy_tail
    )
    assert (
        field_names.index(
            "implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_c_m2"
        )
        > legacy_tail
    )


@pytest.mark.parametrize(
    "present_key",
    [
        "implicit_mean_last_returned_outer_flight_time_mean_s",
        "implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2",
    ],
)
def test_load_fortran_result_requires_complete_implicit_mean_shadow_pair(
    tmp_path: Path, present_key: str
) -> None:
    out = tmp_path / present_key
    out.mkdir()
    values = {
        "implicit_mean_last_returned_outer_flight_time_mean_s": "2.5e-6",
        "implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2": "7.5e-12",
    }
    _write_minimal_result_fixture(
        out,
        summary_extra=[f"{present_key}={values[present_key]}"],
    )

    with pytest.raises(
        ValueError,
        match="implicit_mean shadow diagnostics must appear together",
    ):
        load_fortran_result(out)


def test_load_fortran_result_queue_disabled_has_no_queue_state(tmp_path: Path) -> None:
    out = tmp_path / "queue_disabled"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=["coupling_outer_queue_enabled=F"],
    )

    result = load_fortran_result(out)

    assert result.coupling_outer_queue_enabled is False
    assert result.outer_photoelectron_population_fraction is None
    assert result.outer_photoelectron_column_per_area_m2 is None
    assert result.outer_photoelectron_column_target_per_area_m2 is None
    assert result.outer_photoelectron_column_residual_per_area_m2 is None
    assert result.outer_queue_event_count is None
    assert result.outer_queue_signed_charge_c is None
    assert result.outer_queue_fingerprint is None


def test_load_fortran_result_queue_disabled_rejects_queue_state(
    tmp_path: Path,
) -> None:
    out = tmp_path / "queue_disabled_with_state"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "coupling_outer_queue_enabled=F",
            "outer_queue_event_count=0",
        ],
    )

    with pytest.raises(ValueError, match="queue_enabled=false forbids"):
        load_fortran_result(out)


@pytest.mark.parametrize(
    "missing_key",
    [
        "outer_photoelectron_population_fraction",
        "outer_photoelectron_column_per_area_m2",
        "outer_photoelectron_column_target_per_area_m2",
        "outer_photoelectron_column_residual_per_area_m2",
        "outer_queue_event_count",
        "outer_queue_signed_charge_C",
        "outer_queue_fingerprint",
    ],
)
def test_load_fortran_result_queue_enabled_requires_complete_state(
    tmp_path: Path, missing_key: str
) -> None:
    out = tmp_path / missing_key
    out.mkdir()
    queue_values = {
        "outer_photoelectron_population_fraction": "1.25",
        "outer_photoelectron_column_per_area_m2": "2.0e12",
        "outer_photoelectron_column_target_per_area_m2": "2.1e12",
        "outer_photoelectron_column_residual_per_area_m2": "-1.0e11",
        "outer_queue_event_count": "7",
        "outer_queue_signed_charge_C": "-3.5e-16",
        "outer_queue_fingerprint": "FEDCBA9876543210",
    }
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "coupling_outer_queue_enabled=T",
            *[
                f"{key}={value}"
                for key, value in queue_values.items()
                if key != missing_key
            ],
        ],
    )

    with pytest.raises(ValueError, match="queue_enabled=true requires"):
        load_fortran_result(out)


@pytest.mark.parametrize("fingerprint", ["fedcba9876543210", "ABCDEF"])
def test_load_fortran_result_rejects_invalid_queue_fingerprint(
    tmp_path: Path, fingerprint: str
) -> None:
    out = tmp_path / fingerprint
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "coupling_outer_queue_enabled=T",
            "outer_photoelectron_population_fraction=1.25",
            "outer_photoelectron_column_per_area_m2=2.0e12",
            "outer_photoelectron_column_target_per_area_m2=2.1e12",
            "outer_photoelectron_column_residual_per_area_m2=-1.0e11",
            "outer_queue_event_count=7",
            "outer_queue_signed_charge_C=-3.5e-16",
            f"outer_queue_fingerprint={fingerprint}",
        ],
    )

    with pytest.raises(ValueError, match="16 uppercase hexadecimal"):
        load_fortran_result(out)


def test_beach_get_mesh_supports_step_selection(tmp_path: Path) -> None:
    out = tmp_path / "run_mesh_step"
    out.mkdir()
    _write_three_mesh_fixture(out)

    beach = Beach(out)
    mesh1 = beach.get_mesh(1, step=10)

    assert mesh1.mesh_ids == (1,)
    assert mesh1.step == 10
    np.testing.assert_allclose(mesh1.charges, np.array([1.0e-9]))
    np.testing.assert_allclose(beach.get_mesh_charge(1, step=1), np.array([5.0e-10]))


def test_get_mesh_defaults_to_latest_history_step(tmp_path: Path) -> None:
    out = tmp_path / "run_mesh_default_latest"
    out.mkdir()
    _write_three_mesh_fixture(out)
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,9.0e-9\n2,9.0e-9\n3,9.0e-9\n",
        encoding="utf-8",
    )

    beach = Beach(out)
    mesh1 = beach.get_mesh(1)

    np.testing.assert_allclose(mesh1.charges, np.array([1.0e-9]))
    np.testing.assert_allclose(beach.get_mesh_charge(3), np.array([-3.0e-9]))


def test_beach_get_mesh_returns_tuple_for_multiple_ids(tmp_path: Path) -> None:
    out = tmp_path / "run_mesh_tuple"
    out.mkdir()
    _write_three_mesh_fixture(out)

    beach = Beach(out)
    mesh2, mesh3 = beach.get_mesh(2, 3)

    assert mesh2.mesh_ids == (2,)
    assert mesh3.mesh_ids == (3,)
    np.testing.assert_allclose(mesh2.charges, np.array([2.0e-9]))
    np.testing.assert_allclose(mesh3.charges, np.array([-3.0e-9]))


def test_calc_coulomb_defaults_to_latest_history_step(tmp_path: Path) -> None:
    out = tmp_path / "run_calc_coulomb_default_latest"
    out.mkdir()
    _write_complete_free_field_config(out)
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=2",
                "absorbed=2",
                "escaped=0",
                "batches=2",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,0.0\n2,0.0\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,1\n"
        "2,1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,0.0,2\n",
        encoding="utf-8",
    )
    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        "1,1,1.0e-1,1,5.0e-10\n"
        "1,1,1.0e-1,2,-1.0e-9\n"
        "2,2,1.0e-2,1,1.0e-9\n"
        "2,2,1.0e-2,2,-2.0e-9\n",
        encoding="utf-8",
    )

    interaction = calc_coulomb(Beach(out), 1, 2)
    explicit_latest = calc_coulomb(Beach(out), 1, 2, step=-1)
    np.testing.assert_allclose(interaction.force_on_a_N, explicit_latest.force_on_a_N)
    assert interaction.force_on_a_N[0] > 0.0
    np.testing.assert_allclose(interaction.force_on_a_N[1:], np.zeros(2), atol=1.0e-18)


def test_calc_coulomb_accepts_target_source_keywords(tmp_path: Path) -> None:
    out = tmp_path / "run_calc_coulomb_keywords"
    out.mkdir()
    _write_three_mesh_fixture(out)
    beach = Beach(out)

    interaction = calc_coulomb(
        beach,
        target=1,
        source=2,
        step=10,
        torque_origin="target_center",
    )
    interaction_legacy_origin = calc_coulomb(
        beach,
        target=1,
        source=2,
        step=10,
        torque_origin="group_a_center",
    )

    assert interaction.group_a_mesh_ids == (1,)
    assert interaction.group_b_mesh_ids == (2,)
    np.testing.assert_allclose(
        interaction.torque_origin_m, interaction_legacy_origin.torque_origin_m
    )
    np.testing.assert_allclose(
        interaction.force_on_a_N, interaction_legacy_origin.force_on_a_N
    )


def test_calc_coulomb_accepts_composite_groups(tmp_path: Path) -> None:
    out = tmp_path / "run_calc_coulomb"
    out.mkdir()
    _write_three_mesh_fixture(out)
    beach = Beach(out)

    mesh1 = beach.get_mesh(1)
    mesh2, mesh3 = beach.get_mesh(2, 3)
    interaction = beach.calc_coulomb([mesh1, mesh2], [mesh3], step=10)

    assert interaction.group_a_mesh_ids == (1, 2)
    assert interaction.group_b_mesh_ids == (3,)
    assert interaction.step == 10
    assert interaction.force_on_a_N.shape == (3,)
    assert interaction.torque_on_a_Nm.shape == (3,)
    np.testing.assert_allclose(
        interaction.force_on_a_N + interaction.force_on_b_N, np.zeros(3), atol=1.0e-18
    )


def test_calc_coulomb_panel_force_obeys_action_reaction(tmp_path: Path) -> None:
    out = tmp_path / "run_calc_coulomb_expected"
    out.mkdir()
    _write_complete_free_field_config(out)
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=1",
                "absorbed=1",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-9\n2,-2.0e-9\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,1.0e-9,1\n"
        "2,1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,-2.0e-9,2\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    interaction = calc_coulomb(result, 1, 2)

    assert interaction.force_on_a_N[0] > 0.0
    np.testing.assert_allclose(
        interaction.force_on_a_N + interaction.force_on_b_N,
        np.zeros(3),
        atol=1.0e-18,
    )
    np.testing.assert_allclose(interaction.force_on_a_N[1:], np.zeros(2), atol=1.0e-18)
    np.testing.assert_allclose(
        interaction.torque_on_a_Nm + interaction.torque_on_b_Nm,
        np.zeros(3),
        atol=1.0e-18,
    )
    assert np.linalg.norm(interaction.torque_on_a_Nm) > 0.0


def test_surface_charge_density_uses_triangle_area() -> None:
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    charges = np.array([6.0, -1.5])

    density = _surface_charge_density(charges, triangles)

    np.testing.assert_allclose(density, np.array([6.0, -3.0]))


def test_compute_potential_mesh_uses_precomputed_output_when_available() -> None:
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9, -1.0e-9]),
        mesh_potential_v=np.array([3.0, -4.0]),
    )

    potential = compute_potential_mesh(result)

    np.testing.assert_allclose(potential, np.array([3.0, -4.0]))


def test_compute_potential_mesh_requires_triangles_when_precomputed_output_is_incompatible() -> (
    None
):
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        mesh_potential_v=np.array([2.0]),
    )

    with pytest.raises(ValueError, match="mesh_triangles.csv"):
        compute_potential_mesh(result, reference_point=[0.0, 0.0, 1.0])


def test_compute_potential_mesh_is_linear_in_panel_charges(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    charges = np.array([2.0e-9, -1.0e-9])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
    )

    scaled_result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=-0.25 * charges,
        triangles=triangles,
    )

    potential = compute_potential_mesh(result)
    scaled = compute_potential_mesh(scaled_result)

    np.testing.assert_allclose(scaled, -0.25 * potential)


def test_compute_potential_mesh_rejects_removed_softening_keyword() -> None:
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9, 1.0e-9]),
        triangles=triangles,
    )

    with pytest.raises(TypeError, match="softening"):
        compute_potential_mesh(result, softening=1.0)  # type: ignore[call-arg]


def test_compute_potential_mesh_includes_finite_panel_self_term(
    tmp_path: Path,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    potential = compute_potential_mesh(result)

    assert potential.shape == (1,)
    assert np.isfinite(potential[0])
    assert potential[0] > 0.0


def test_compute_potential_mesh_uses_panel_centroids_as_targets(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    charges = np.array([2.0e-9, -1.0e-9])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
    )

    captured: dict[str, np.ndarray] = {}

    class _CapturePanelKernel(_UnitPanelKernel):
        def __init__(self, source_triangles, source_charges, **kwargs) -> None:
            captured["source_triangles"] = np.asarray(source_triangles).copy()
            super().__init__(source_triangles, source_charges, **kwargs)

        def eval_phi(self, points: np.ndarray) -> np.ndarray:
            captured.setdefault("targets", np.asarray(points).copy())
            return super().eval_phi(points)

    monkeypatch.setattr(potential_module, "FieldKernel", _CapturePanelKernel)

    potential = compute_potential_mesh(result)

    np.testing.assert_array_equal(captured["source_triangles"], triangles)
    np.testing.assert_allclose(captured["targets"], triangles.mean(axis=1))
    assert np.all(np.isfinite(potential))


def test_compute_potential_mesh_supports_reference_point_difference(
    tmp_path: Path,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    charges = np.array([2.0e-9, -1.0e-9])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
    )

    reference = np.array([0.0, 0.0, 6.0])
    phi = compute_potential_mesh(
        result,
        reference_point=reference,
    )
    phi_abs = compute_potential_mesh(result)
    phi_ref = compute_potential_points(
        result,
        reference.reshape(1, 3),
    )[0]

    np.testing.assert_allclose(phi, phi_abs - phi_ref)


def test_compute_potential_mesh_reads_panel_solver_options_from_config(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    out = tmp_path / "run_panel_options"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,2.0e-9\n2,-1.0e-9\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,3.0,0.0,0.0,0.0,3.0,0.0,2.0e-9,1\n"
        "2,3.0,0.0,0.0,6.0,0.0,0.0,3.0,3.0,0.0,-1.0e-9,2\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                "tree_theta = 0.5",
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    captured: dict[str, object] = {}

    class _CaptureOptionsKernel(_UnitPanelKernel):
        def __init__(self, source_triangles, source_charges, *, options, **kwargs) -> None:
            captured["options"] = options
            super().__init__(
                source_triangles,
                source_charges,
                options=options,
                **kwargs,
            )

    monkeypatch.setattr(potential_module, "FieldKernel", _CaptureOptionsKernel)
    potential = compute_potential_mesh(result)

    assert np.all(np.isfinite(potential))
    assert captured["options"].theta == pytest.approx(0.5)


def test_compute_potential_mesh_supports_species1_injection_reference_from_config(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_species1_ref"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,2.0e-9\n", encoding="utf-8")
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,2.0e-9,1\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                "tree_theta = 0.5",
                "",
                "[particles]",
                "",
                "[[particles.species]]",
                'inject_face = "z_high"',
                "pos_low = [0.0, 0.0, 10.0]",
                "pos_high = [2.0, 2.0, 10.0]",
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    absolute = compute_potential_mesh(result)
    phi = compute_potential_mesh(
        result,
        reference_point="species1_injection_center",
    )

    reference = compute_potential_points(
        result,
        np.array([[1.0, 1.0, 10.0]]),
    )[0]
    np.testing.assert_allclose(phi, absolute - reference)


def test_compute_potential_mesh_returns_finite_panel_values(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9, -1.0e-9]),
        triangles=triangles,
    )

    potential = compute_potential_mesh(result)

    assert np.all(np.isfinite(potential))


def test_compute_potential_mesh_requires_triangles() -> None:
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0]),
        triangles=None,
    )

    with pytest.raises(ValueError, match="mesh_triangles.csv"):
        compute_potential_mesh(result)


def test_compute_potential_mesh_rejects_removed_softening_argument() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0]),
        triangles=triangles,
    )

    with pytest.raises(TypeError, match="softening"):
        compute_potential_mesh(result, softening=-1.0)  # type: ignore[call-arg]


def test_compute_potential_mesh_rejects_removed_self_term_argument() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0]),
        triangles=triangles,
    )

    with pytest.raises(TypeError, match="self_term"):
        compute_potential_mesh(result, self_term="invalid")  # type: ignore[call-arg]


def test_compute_potential_mesh_has_no_softened_point_compatibility_mode() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0]),
        triangles=triangles,
    )

    with pytest.raises(TypeError, match="self_term"):
        compute_potential_mesh(  # type: ignore[call-arg]
            result,
            self_term="softened_point",
        )


def test_compute_potential_mesh_rejects_degenerate_panel(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    with pytest.raises(ValueError, match="non-degenerate"):
        compute_potential_mesh(result)


def test_compute_potential_points_returns_finite_panel_potential(
    tmp_path: Path,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charge,
        triangles=triangles,
    )
    points = np.array([[0.0, 0.0, 2.0], [0.0, 0.0, 4.0]])

    potential = compute_potential_points(result, points)

    assert np.all(np.isfinite(potential))
    assert potential[0] > potential[1] > 0.0


def test_compute_potential_points_supports_periodic2_image_sum(
    tmp_path: Path,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charge,
        triangles=triangles,
    )
    points = np.array([[0.0, 0.0, 2.0]])
    periodic2 = {
        "axes": (0, 1),
        "lengths": (1.0, 1.0),
        "image_layers": 1,
        "far_correction": "none",
        "ewald_layers": 4,
    }

    potential = compute_potential_points(result, points, periodic2=periodic2)
    free = compute_potential_points(result, points)

    assert np.all(np.isfinite(potential))
    assert potential[0] > free[0] > 0.0


def test_compute_potential_mesh_supports_periodic2_image_sum(
    tmp_path: Path,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charge,
        triangles=triangles,
    )
    periodic2 = {
        "axes": (0, 1),
        "lengths": (1.0, 1.0),
        "image_layers": 1,
        "far_correction": "none",
        "ewald_layers": 4,
    }

    potential = compute_potential_mesh(
        result,
        periodic2=periodic2,
    )
    free = compute_potential_mesh(result)

    assert np.all(np.isfinite(potential))
    assert potential[0] > free[0] > 0.0


@pytest.mark.parametrize("far_correction", ["auto", "m2l_root_oracle"])
def test_compute_potential_points_auto_detects_periodic2_from_config(
    tmp_path: Path, far_correction: str
) -> None:
    out = tmp_path / "run_periodic"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,2.0e-9\n", encoding="utf-8")
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,-1.0,-1.0,0.0,1.0,-1.0,0.0,0.0,2.0,0.0,2.0e-9,1\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [1.0, 1.0, 1.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                "field_periodic_image_layers = 1",
                f'field_periodic_far_correction = "{far_correction}"',
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    points = np.array([[0.0, 0.0, 0.5]])
    potential = compute_potential_points(result, points)
    explicit = compute_potential_points(
        result,
        points,
        periodic2={
            "axes": (0, 1),
            "lengths": (1.0, 1.0),
            "image_layers": 1,
            "far_correction": "none",
        },
    )

    np.testing.assert_allclose(potential, explicit)


def test_coulomb_matrix_auto_reader_accepts_historical_root_oracle(
    tmp_path: Path,
) -> None:
    (tmp_path / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [1.0, 1.0, 1.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                'field_periodic_far_correction = "m2l_root_oracle"',
                "field_periodic_ewald_layers = 4",
            ]
        ),
        encoding="utf-8",
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0]),
        triangles=np.array(
            [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
        ),
    )

    periodic2 = _periodic2_for_coulomb_matrix(result, config_path=None)

    assert periodic2 is not None
    assert periodic2.far_correction == "none"
    with pytest.raises(ValueError, match='was removed; use "none"'):
        _periodic2_for_coulomb_matrix(
            result,
            config_path=tmp_path / "beach.toml",
        )


def test_composed_analyses_normalize_auto_loaded_historical_root_oracle(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    out = tmp_path / "run_historical_periodic"
    out.mkdir()
    _write_coulomb_matrix_fixture(out)
    config_path = out / "beach.toml"
    config_path.write_text(
        config_path.read_text(encoding="utf-8")
        + "\n".join(
            [
                "",
                "",
                "[sim]",
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [3.0, 2.0, 2.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                "field_periodic_image_layers = 0",
                'field_periodic_far_correction = "m2l_root_oracle"',
                "",
            ]
        ),
        encoding="utf-8",
    )

    beach = Beach(out)
    slices = beach.compute_potential_slices(
        box_min=[0.0, 0.0, -1.0],
        box_max=[3.0, 2.0, 2.0],
        grid_n=2,
    )
    assert all(np.all(np.isfinite(item.potential_V)) for item in slices.values())

    analysis = beach.analyze_coulomb_mobility()
    assert analysis.records

    fig, ax = beach.plot_coulomb_force_matrix(component="x")
    assert np.all(np.isfinite(getattr(ax, "_beach_coulomb_matrix")["matrix"]))
    fig.clf()

    explicit_beach = Beach(out, config_path=config_path)
    with pytest.raises(ValueError, match='was removed; use "none"'):
        explicit_beach.compute_potential_points(np.array([[0.5, 0.5, 0.5]]))
    with pytest.raises(ValueError, match='was removed; use "none"'):
        explicit_beach.analyze_coulomb_mobility()


def test_compute_potential_points_wraps_periodic2_points_to_fundamental_cell(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_periodic_wrap"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,2.0e-9\n", encoding="utf-8")
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,0.75,0.0,0.0,0.0,0.75,0.0,2.0e-9,1\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [1.0, 1.0, 1.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                "field_periodic_image_layers = 1",
                'field_periodic_far_correction = "none"',
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    points = np.array(
        [
            [0.25, 0.25, 1.0],
            [1.25, -0.75, 1.0],
        ]
    )

    potential = compute_potential_points(result, points)

    assert np.all(np.isfinite(potential))
    np.testing.assert_allclose(potential[0], potential[1])


def test_auto_periodic2_from_result_defaults_far_correction_to_none(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_periodic_default_far_correction"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,2.0e-9\n", encoding="utf-8")
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,-1.0,-1.0,0.0,1.0,-1.0,0.0,0.0,2.0,0.0,2.0e-9,1\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_solver = "fmm"',
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [1.0, 1.0, 1.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                "field_periodic_image_layers = 1",
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    periodic2 = _auto_periodic2_from_result(result)

    assert periodic2 is not None
    assert periodic2[4] == "none"
    assert periodic2[6] == 4


def test_auto_periodic2_from_result_preserves_none_far_correction(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_periodic_none_far_correction"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,2.0e-9\n", encoding="utf-8")
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,-1.0,-1.0,0.0,1.0,-1.0,0.0,0.0,2.0,0.0,2.0e-9,1\n",
        encoding="utf-8",
    )
    (out / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_solver = "fmm"',
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [1.0, 1.0, 1.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                "field_periodic_image_layers = 1",
                'field_periodic_far_correction = "none"',
                "field_periodic_ewald_layers = 4",
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    periodic2 = _auto_periodic2_from_result(result)

    assert periodic2 is not None
    assert periodic2[4] == "none"
    assert periodic2[6] == 4


def test_potential_history_reuses_panel_kernel_with_periodic2() -> None:
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charges_history = np.array([[2.0e-9, -1.0e-9]])
    periodic2 = _coerce_periodic2(
        {"axes": (0, 1), "lengths": (1.0, 1.0), "image_layers": 1}
    )
    assert periodic2 is not None

    potential = _potential_history(
        charges_history,
        triangles,
        periodic2=periodic2,
    )

    free = _potential_history(charges_history, triangles)

    assert np.all(np.isfinite(potential))
    np.testing.assert_allclose(potential[:, 1], -0.5 * potential[:, 0])
    assert abs(potential[0, 0]) > abs(free[0, 0])


def test_coerce_periodic2_rejects_legacy_ewald_modes() -> None:
    with pytest.raises(
        ValueError,
        match='periodic2.far_correction must be "auto" or "none"',
    ):
        _coerce_periodic2(
            {
                "axes": (0, 1),
                "lengths": (1.0, 1.0),
                "image_layers": 1,
                "far_correction": "ewald",
                "ewald_alpha": 1.2,
                "ewald_layers": 4,
            }
        )


def test_coerce_periodic2_accepts_auto_default() -> None:
    periodic2 = _coerce_periodic2(
        {
            "axes": (0, 1),
            "lengths": (1.0, 1.0),
            "image_layers": 1,
            "far_correction": "auto",
            "ewald_layers": 4,
        }
    )

    assert periodic2 is not None
    assert periodic2[4] == "none"
    assert periodic2[6] == 4


def test_coerce_periodic2_preserves_none() -> None:
    periodic2 = _coerce_periodic2(
        {
            "axes": (0, 1),
            "lengths": (1.0, 1.0),
            "image_layers": 1,
            "far_correction": "none",
            "ewald_layers": 4,
        }
    )

    assert periodic2 is not None
    assert periodic2[4] == "none"
    assert periodic2[6] == 4


def test_coerce_periodic2_rejects_removed_root_oracle() -> None:
    with pytest.raises(
        ValueError,
        match='periodic2.far_correction "m2l_root_oracle" was removed',
    ):
        _coerce_periodic2(
            {
                "axes": (0, 1),
                "lengths": (1.0, 1.0),
                "image_layers": 1,
                "far_correction": "m2l_root_oracle",
                "ewald_layers": 4,
            }
        )


@pytest.mark.parametrize("far_correction", ["m2l_root", "m2l_root_trunc"])
def test_coerce_periodic2_rejects_removed_far_correction_aliases(
    far_correction: str,
) -> None:
    with pytest.raises(
        ValueError,
        match='periodic2.far_correction must be "auto" or "none"',
    ):
        _coerce_periodic2(
            {
                "axes": (0, 1),
                "lengths": (1.0, 1.0),
                "image_layers": 1,
                "far_correction": far_correction,
                "ewald_layers": 4,
            }
        )


def test_potential_history_supports_reference_point_difference() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    charges_history = np.array([[2.0e-9, 4.0e-9]])
    reference_point = np.array([0.0, 0.0, 2.0])

    potential = _potential_history(
        charges_history,
        triangles,
        reference_point=reference_point,
    )
    absolute = _potential_history(charges_history, triangles)

    assert np.all(np.isfinite(potential))
    np.testing.assert_allclose(potential[:, 1], 2.0 * potential[:, 0])
    reference_offset = absolute - potential
    assert np.all(reference_offset > 0.0)
    np.testing.assert_allclose(reference_offset[:, 1], 2.0 * reference_offset[:, 0])


def test_potential_animation_rejects_cached_kneq0_component(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    (output_dir / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -2.0]",
                "box_max = [4.0, 5.0, 12.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                'field_periodic_far_correction = "cached_kneq0"',
            ]
        ),
        encoding="utf-8",
    )
    triangles = np.array(
        [[[0.2, 0.2, 0.0], [0.8, 0.2, 0.0], [0.2, 0.8, 0.0]]]
    )
    result = FortranRunResult(
        directory=output_dir,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=1,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([2.0e-9]),
        triangles=triangles,
        history=_make_history(
            mesh_nelem=1,
            history=np.array([[1.0e-9, 2.0e-9]]),
        ),
    )
    with pytest.raises(ValueError, match="only the k!=0 component"):
        animate_history_mesh(
            result,
            quantity="potential",
            reference_point=[2.0, 2.5, 10.0],
        )


def test_compute_potential_points_validates_shape_and_chunk_size() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    with pytest.raises(ValueError, match="shape"):
        compute_potential_points(result, np.array([0.0, 0.0, 1.0]))
    with pytest.raises(ValueError, match="chunk_size"):
        compute_potential_points(result, np.array([[0.0, 0.0, 1.0]]), chunk_size=0)
    with pytest.raises(ValueError, match="axes"):
        compute_potential_points(
            result,
            np.array([[0.0, 0.0, 1.0]]),
            periodic2={"axes": (0, 0), "lengths": (1.0, 1.0)},
        )


def test_compute_potential_slices_returns_xy_yz_xz_grids(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charge,
        triangles=triangles,
    )

    slices = compute_potential_slices(
        result,
        box_min=[-1.0, -1.0, -1.0],
        box_max=[1.0, 1.0, 1.0],
        grid_n=5,
        xy_z=0.5,
        yz_x=0.5,
        xz_y=0.5,
    )

    assert set(slices.keys()) == {"xy", "yz", "xz"}
    for plane in ("xy", "yz", "xz"):
        slc = slices[plane]
        assert slc.potential_V.shape == (5, 5)
        assert slc.u_values_m.shape == (5,)
        assert slc.v_values_m.shape == (5,)

    center_values = np.array(
        [
            slices["xy"].potential_V[2, 2],
            slices["yz"].potential_V[2, 2],
            slices["xz"].potential_V[2, 2],
        ]
    )
    assert np.all(np.isfinite(center_values))
    assert np.all(center_values > 0.0)
    assert np.ptp(center_values) > 0.0


def test_compute_potential_slices_rejects_invalid_grid_or_coordinate() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    with pytest.raises(ValueError, match="grid_n"):
        compute_potential_slices(
            result,
            box_min=[0.0, 0.0, 0.0],
            box_max=[1.0, 1.0, 1.0],
            grid_n=1,
        )
    with pytest.raises(ValueError, match="xy_z"):
        compute_potential_slices(
            result,
            box_min=[0.0, 0.0, 0.0],
            box_max=[1.0, 1.0, 1.0],
            grid_n=8,
            xy_z=2.0,
        )


def test_plot_potential_slices_returns_figure_and_axes(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    fig, axes = plot_potential_slices(
        result,
        box_min=[-1.0, -1.0, -1.0],
        box_max=[1.0, 1.0, 1.0],
        grid_n=10,
        xy_z=0.3,
        yz_x=0.3,
        xz_y=0.3,
    )

    assert fig is not None
    assert len(axes) == 3
    fig.clf()


def test_plot_potential_slices_applies_vmin_vmax(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    fig, axes = plot_potential_slices(
        result,
        box_min=[-1.0, -1.0, -1.0],
        box_max=[1.0, 1.0, 1.0],
        grid_n=8,
        vmin=-5.0,
        vmax=12.0,
    )

    for ax in axes:
        assert ax.collections
        assert ax.collections[0].get_clim() == (-5.0, 12.0)
    fig.clf()


def test_plot_potential_slices_rejects_invalid_vmin_vmax(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    with pytest.raises(ValueError, match="vmin"):
        plot_potential_slices(
            result,
            box_min=[-1.0, -1.0, -1.0],
            box_max=[1.0, 1.0, 1.0],
            grid_n=8,
            vmin=1.0,
            vmax=1.0,
        )


def test_wrap_periodic2_triangles_by_centroid_keeps_face_shape() -> None:
    triangles = np.array(
        [
            [[1.10, 0.10, 0.0], [1.20, 0.10, 0.0], [1.10, 0.20, 0.0]],
            [[-0.20, 1.05, 1.0], [-0.10, 1.05, 1.0], [-0.20, 1.15, 1.0]],
        ]
    )

    wrapped = _wrap_periodic2_triangles_by_centroid(
        triangles,
        axes=(0, 1),
        lengths=(1.0, 1.0),
        origins=(0.0, 0.0),
    )

    centers = wrapped.mean(axis=1)
    assert np.all(centers[:, 0] >= 0.0)
    assert np.all(centers[:, 0] < 1.0)
    assert np.all(centers[:, 1] >= 0.0)
    assert np.all(centers[:, 1] < 1.0)
    np.testing.assert_allclose(
        wrapped[:, 1, :] - wrapped[:, 0, :],
        triangles[:, 1, :] - triangles[:, 0, :],
    )
    np.testing.assert_allclose(
        wrapped[:, 2, :] - wrapped[:, 0, :],
        triangles[:, 2, :] - triangles[:, 0, :],
    )


def test_wrap_periodic2_triangles_by_mesh_centroid_keeps_object_intact() -> None:
    triangles = np.array(
        [
            [[0.90, 0.10, 0.0], [0.95, 0.10, 0.0], [0.90, 0.15, 0.0]],
            [[1.10, 0.10, 0.0], [1.15, 0.10, 0.0], [1.10, 0.15, 0.0]],
        ]
    )

    face_wrapped = _wrap_periodic2_triangles_by_centroid(
        triangles,
        axes=(0, 1),
        lengths=(1.0, 1.0),
        origins=(0.0, 0.0),
    )
    mesh_wrapped = _wrap_periodic2_triangles_by_mesh_centroid(
        triangles,
        np.array([2, 2]),
        axes=(0, 1),
        lengths=(1.0, 1.0),
        origins=(0.0, 0.0),
    )

    original_dx = triangles[1].mean(axis=0)[0] - triangles[0].mean(axis=0)[0]
    face_dx = face_wrapped[1].mean(axis=0)[0] - face_wrapped[0].mean(axis=0)[0]
    mesh_dx = mesh_wrapped[1].mean(axis=0)[0] - mesh_wrapped[0].mean(axis=0)[0]
    assert abs(face_dx) > 0.5
    assert mesh_dx == pytest.approx(original_dx)


def test_plot_charge_mesh_can_apply_periodic2_mesh_wrapping() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    triangles = np.array([[[1.10, 0.10, 0.0], [1.20, 0.10, 0.0], [1.10, 0.20, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    fig, ax = plot_charge_mesh(
        result,
        periodic2={"axes": (0, 1), "lengths": (1.0, 1.0), "origins": (0.0, 0.0)},
        apply_periodic2_mesh=True,
    )

    x_center = 0.5 * sum(ax.get_xlim())
    y_center = 0.5 * sum(ax.get_ylim())
    assert x_center == pytest.approx(0.15)
    assert y_center == pytest.approx(0.15)
    fig.clf()


def test_plot_potential_mesh_returns_figure_and_axes(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9, -1.0e-9]),
        triangles=triangles,
    )

    fig, ax = plot_potential_mesh(result, reference_point=None)

    assert fig is not None
    assert ax is not None
    fig.clf()

    fig, ax = plot_potential_mesh(
        result,
        reference_point=None,
        axis_unit="um",
    )
    assert ax.get_xlabel() == "x [um]"
    fig.clf()


def test_plot_potential_mesh_can_apply_periodic2_mesh_wrapping(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array([[[1.10, 0.10, 0.0], [1.20, 0.10, 0.0], [1.10, 0.20, 0.0]]])
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    fig, ax = plot_potential_mesh(
        result,
        reference_point=None,
        periodic2={"axes": (0, 1), "lengths": (1.0, 1.0), "origins": (0.0, 0.0)},
        apply_periodic2_mesh=True,
    )

    x_center = 0.5 * sum(ax.get_xlim())
    y_center = 0.5 * sum(ax.get_ylim())
    assert x_center == pytest.approx(0.15)
    assert y_center == pytest.approx(0.15)
    fig.clf()


def test_plot_potential_mesh_accepts_custom_view_angles(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9, -1.0e-9]),
        triangles=triangles,
    )

    fig, ax = plot_potential_mesh(
        result,
        reference_point=None,
        view_elev=11.0,
        view_azim=37.0,
    )

    assert ax.elev == pytest.approx(11.0)
    assert ax.azim == pytest.approx(37.0)
    ax.view_init(elev=5.0, azim=15.0)
    assert ax.elev == pytest.approx(5.0)
    assert ax.azim == pytest.approx(15.0)
    fig.clf()


def test_plot_mesh_source_boxplot_charge_uses_area_weighted_statistics() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [2.0, 0.0, 1.0], [0.0, 2.0, 1.0]],
        ]
    )
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0, 10.0]),
        triangles=triangles,
        mesh_ids=np.array([11, 11], dtype=np.int64),
    )

    fig, ax = plot_mesh_source_boxplot(result, quantity="charge")

    stats = getattr(ax, "_beach_box_stats")
    assert len(stats) == 1
    assert stats[0]["med"] == pytest.approx(10.0)
    fig.clf()


def test_plot_mesh_source_boxplot_supports_potential_quantity(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 0.8, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9, -2.0e-9]),
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    fig, ax = plot_mesh_source_boxplot(
        result,
        quantity="potential",
    )

    stats = getattr(ax, "_beach_box_stats")
    assert len(stats) == 2
    assert all(np.isfinite(float(item["med"])) for item in stats)
    fig.clf()


def test_plot_mesh_source_boxplot_rejects_invalid_quantity() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
    )

    with pytest.raises(ValueError, match="quantity"):
        plot_mesh_source_boxplot(result, quantity="invalid")


def test_beach_plot_mesh_source_boxplot_uses_mesh_source_labels(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    out = tmp_path / "run_source_boxplot"
    out.mkdir()
    _write_three_mesh_fixture(out)

    beach = Beach(out)
    fig, ax = beach.plot_mesh_source_boxplot(quantity="charge")

    stats = getattr(ax, "_beach_box_stats")
    labels = [str(item["label"]) for item in stats]
    assert "id=1 (template/plane/insulator/eps=1)" in labels
    assert "id=2 (template/box/insulator/eps=1)" in labels
    assert "id=3 (template/sphere/insulator/eps=1)" in labels
    fig.clf()


def test_plot_coulomb_force_matrix_auto_labels_targets_from_config(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    out = tmp_path / "run_coulomb_matrix"
    out.mkdir()
    _write_coulomb_matrix_fixture(out)

    beach = Beach(out)
    fig, ax = beach.plot_coulomb_force_matrix(component="x")

    matrix_info = getattr(ax, "_beach_coulomb_matrix")
    assert matrix_info["target_labels"] == ("plane", "sphere1", "sphere2")
    assert matrix_info["source_labels"] == ("plane", "sphere1", "sphere2", "net")
    matrix = matrix_info["matrix"]
    assert np.all(np.isfinite(matrix))
    np.testing.assert_allclose(np.diag(matrix[:3]), np.zeros(3))
    np.testing.assert_allclose(
        matrix[-1],
        np.sum(matrix[:-1], axis=0),
    )
    fig.clf()

    fig, ax = beach.plot_coulomb_force_matrix(component="x", workers=2)
    worker_info = getattr(ax, "_beach_coulomb_matrix")
    assert worker_info["workers"] == 2
    np.testing.assert_allclose(worker_info["matrix"], matrix_info["matrix"])
    fig.clf()


def test_plot_coulomb_force_matrix_accepts_explicit_target_kinds(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    out = tmp_path / "run_coulomb_matrix_plane"
    out.mkdir()
    _write_coulomb_matrix_fixture(out)

    fig, ax = plot_coulomb_force_matrix(
        Beach(out),
        component="x",
        target_kinds=("plane",),
        annotate=False,
    )

    matrix_info = getattr(ax, "_beach_coulomb_matrix")
    assert matrix_info["target_labels"] == ("plane",)
    assert matrix_info["source_labels"] == ("plane", "sphere1", "sphere2", "net")
    matrix = matrix_info["matrix"]
    assert np.all(np.isfinite(matrix))
    assert matrix[0, 0] == pytest.approx(0.0)
    np.testing.assert_allclose(
        matrix[-1, 0],
        np.sum(matrix[:-1, 0]),
    )
    fig.clf()


def test_analyze_coulomb_mobility_defaults_to_non_support_objects(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_mobility_default"
    out.mkdir()
    _write_coulomb_matrix_fixture(out)

    analysis = analyze_coulomb_mobility(Beach(out))

    assert isinstance(analysis, CoulombMobilityAnalysis)
    assert analysis.support_kinds == ("plane",)
    assert [record.label for record in analysis.records] == ["sphere1", "sphere2"]
    assert all(record.lift_ratio is None for record in analysis.records)


def test_analyze_coulomb_mobility_computes_lift_and_slide_ratios(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_mobility_ratios"
    out.mkdir()
    _write_mobility_fixture(out)

    analysis = analyze_coulomb_mobility(
        Beach(out),
        density_kg_m3=1000.0,
        mu_static=0.5,
    )

    assert len(analysis.records) == 1
    record = analysis.records[0]
    expected_interaction = calc_coulomb(Beach(out), target=2, source=1)
    expected_force = expected_interaction.force_on_a_N[2]
    expected_mass = 1000.0 * (4.0 * np.pi * 0.1**3 / 3.0)
    expected_weight = expected_mass * 9.81

    assert record.label == "sphere"
    np.testing.assert_allclose(
        record.force_N,
        np.array([0.0, 0.0, expected_force]),
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        record.torque_Nm,
        expected_interaction.torque_on_a_Nm,
        atol=1.0e-12,
    )
    assert record.force_normal_N == pytest.approx(expected_force)
    assert record.force_tangent_N == pytest.approx(0.0, abs=1.0e-12)
    assert record.mass_kg == pytest.approx(expected_mass)
    assert record.weight_support_N == pytest.approx(expected_weight)
    assert record.resisting_normal_N == pytest.approx(expected_weight)
    assert record.effective_normal_load_N == pytest.approx(0.0)
    assert record.lift_ratio == pytest.approx(expected_force / expected_weight)
    assert record.slide_ratio == pytest.approx(0.0, abs=1.0e-12)
    assert record.notes == ()


def test_select_frame_columns_with_frame_stride() -> None:
    cols = _select_frame_columns(10, frame_stride=3, total_frames=None)
    np.testing.assert_array_equal(cols, np.array([0, 3, 6, 9], dtype=np.int64))


def test_select_frame_columns_with_total_frames() -> None:
    cols = _select_frame_columns(10, frame_stride=1, total_frames=4)
    np.testing.assert_array_equal(cols, np.array([0, 3, 6, 9], dtype=np.int64))


def test_select_frame_columns_with_total_frames_larger_than_snapshots() -> None:
    cols = _select_frame_columns(3, frame_stride=1, total_frames=10)
    np.testing.assert_array_equal(cols, np.array([0, 1, 2], dtype=np.int64))


def test_select_frame_columns_with_one_total_frame() -> None:
    cols = _select_frame_columns(7, frame_stride=1, total_frames=1)
    np.testing.assert_array_equal(cols, np.array([0], dtype=np.int64))


def test_animate_history_mesh_requires_charge_history(tmp_path: Path) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
        history=None,
    )

    with pytest.raises(ValueError, match="charge_history.csv"):
        animate_history_mesh(result, tmp_path / "history.gif")


def test_animate_history_mesh_rejects_invalid_quantity(tmp_path: Path) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
        history=_make_history(
            mesh_nelem=1,
            history=np.array([[1.0e-9, 1.2e-9]]),
        ),
    )

    with pytest.raises(ValueError, match="quantity"):
        animate_history_mesh(result, tmp_path / "history.gif", quantity="invalid")


@pytest.mark.parametrize("source_model", ["point", "unknown"])
def test_potential_animation_rejects_legacy_or_missing_source_receipt(
    tmp_path: Path,
    source_model: str,
) -> None:
    pytest.importorskip("matplotlib")
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
        dtype=float,
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
        history=_make_history(
            mesh_nelem=1,
            history=np.array([[1.0e-9, 1.2e-9]]),
        ),
        field_source_model=source_model,
    )

    with pytest.raises(ValueError, match="field_source_model='triangle_p0'"):
        animate_history_mesh(result, quantity="potential")


def test_animate_history_mesh_rejects_invalid_total_frames(tmp_path: Path) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
        history=_make_history(
            mesh_nelem=1,
            history=np.array([[1.0e-9, 1.2e-9]]),
        ),
    )

    with pytest.raises(ValueError, match="total_frames"):
        animate_history_mesh(result, tmp_path / "history.gif", total_frames=0)


def test_animate_history_mesh_rejects_stride_and_total_frames_combination(
    tmp_path: Path,
) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = FortranRunResult(
        directory=Path("dummy"),
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([1.0e-9]),
        triangles=triangles,
        history=_make_history(
            mesh_nelem=1,
            history=np.array([[1.0e-9, 1.2e-9]]),
        ),
    )

    with pytest.raises(ValueError, match="cannot be used together"):
        animate_history_mesh(
            result,
            tmp_path / "history.gif",
            frame_stride=2,
            total_frames=2,
        )


@pytest.mark.parametrize("quantity", ["charge", "potential"])
def test_animate_history_mesh_writes_gif(tmp_path: Path, quantity: str) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    pytest.importorskip("PIL")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [2.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    charge_history = np.array(
        [
            [1.0e-10, 2.0e-10, 3.0e-10],
            [-2.0e-10, -1.0e-10, -0.5e-10],
        ]
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charge_history[:, -1],
        triangles=triangles,
        history=_make_history(
            mesh_nelem=2,
            history=charge_history,
            batch_indices=np.array([1, 2, 3]),
            processed_particles_by_batch=np.array([10, 20, 30]),
            rel_change_by_batch=np.array([1.0e-1, 1.0e-2, 1.0e-3]),
        ),
    )
    out = tmp_path / f"{quantity}.gif"

    written = animate_history_mesh(
        result,
        out,
        quantity=quantity,
        fps=4,
        frame_stride=1,
    )

    assert written == out
    assert out.exists()
    assert out.stat().st_size > 0
