from pathlib import Path

import numpy as np
import pytest

from beach.fortran_results import (
    Beach,
    calc_coulomb,
    FortranChargeHistory,
    K_COULOMB,
    FortranRunResult,
    _select_frame_columns,
    animate_history_mesh,
    compute_potential_points,
    compute_potential_slices,
    compute_potential_mesh,
    _surface_charge_density,
    list_fortran_runs,
    load_fortran_result,
    plot_mesh_source_boxplot,
    plot_coulomb_force_matrix,
    plot_potential_slices,
    plot_potential_mesh,
)
from beach.fortran_results.potential import _coerce_periodic2, _potential_history


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


def _write_three_mesh_fixture(out: Path) -> None:
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=3",
                "processed_particles=12",
                "absorbed=9",
                "escaped=3",
                "batches=10",
                "last_rel_change=1.0e-6",
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
                "last_rel_change=1.0e-8",
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
        "mesh_id,source_kind,template_kind,elem_count\n"
        "1,template,plane,1\n"
        "2,template,sphere,1\n",
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
    assert result.triangles is not None
    assert result.triangles.shape == (2, 3, 3)
    assert result.history is not None
    assert result.mesh_ids is not None
    np.testing.assert_array_equal(result.mesh_ids, np.array([1, 2]))
    assert result.mesh_sources is not None
    assert result.mesh_sources[1].template_kind == "plane"
    np.testing.assert_allclose(result.charges, np.array([1.0e-10, -2.0e-10]))
    np.testing.assert_allclose(
        result.history.as_array(),
        np.array([[2.0e-11, 1.0e-10], [-1.0e-11, -2.0e-10]]),
    )
    np.testing.assert_array_equal(result.history.processed_particles_by_batch, np.array([5, 10]))
    np.testing.assert_allclose(result.history.rel_change_by_batch, np.array([3.0e-1, 1.0e-8]))
    np.testing.assert_array_equal(result.history.batch_indices, np.array([1, 3]))


def test_load_fortran_result_history_supports_step_access(tmp_path: Path) -> None:
    out = tmp_path / "run_lazy"
    out.mkdir()
    _write_three_mesh_fixture(out)

    result = load_fortran_result(out)

    assert result.history is not None
    np.testing.assert_allclose(result.history[1], np.array([5.0e-10, 1.0e-9, -1.5e-9]))
    np.testing.assert_allclose(result.history_at(10), np.array([1.0e-9, 2.0e-9, -3.0e-9]))
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
    np.testing.assert_allclose(result.history.get_step(10), np.array([1.0e-9, 2.0e-9, -3.0e-9]))


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
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    result = load_fortran_result(out)

    assert result.escaped_boundary == 0
    assert result.survived_max_step == 0


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
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=2",
                "absorbed=2",
                "escaped=0",
                "batches=2",
                "last_rel_change=0.0",
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
    expected_fx = K_COULOMB * 2.0e-18
    np.testing.assert_allclose(interaction.force_on_a_N, np.array([expected_fx, 0.0, 0.0]))


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
    np.testing.assert_allclose(interaction.force_on_a_N, interaction_legacy_origin.force_on_a_N)


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


def test_calc_coulomb_matches_two_charge_expected_force(tmp_path: Path) -> None:
    out = tmp_path / "run_calc_coulomb_expected"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=1",
                "absorbed=1",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
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

    expected_fx = K_COULOMB * 2.0e-18
    np.testing.assert_allclose(interaction.force_on_a_N, np.array([expected_fx, 0.0, 0.0]))
    np.testing.assert_allclose(interaction.torque_on_a_Nm, np.zeros(3))


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


def test_compute_potential_mesh_matches_area_equivalent_expected_values() -> None:
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    charges = np.array([2.0e-9, -1.0e-9])
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
        charges=charges,
        triangles=triangles,
    )

    distance = 3.0
    self_coeff = 2.0 * np.sqrt(np.pi) / np.sqrt(4.5)
    expected = K_COULOMB * np.array(
        [
            charges[0] * self_coeff + charges[1] / distance,
            charges[0] / distance + charges[1] * self_coeff,
        ]
    )

    potential = compute_potential_mesh(result, self_term="area_equivalent")

    np.testing.assert_allclose(potential, expected)


def test_compute_potential_mesh_changes_with_softening_for_offdiag_terms() -> None:
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

    potential_small = compute_potential_mesh(
        result, softening=1.0e-6, self_term="exclude"
    )
    potential_large = compute_potential_mesh(result, softening=1.0, self_term="exclude")

    assert np.all(np.isfinite(potential_small))
    assert np.all(np.isfinite(potential_large))
    assert np.all(potential_small > potential_large)


def test_compute_potential_mesh_excludes_self_term_for_single_triangle() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]])
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

    potential = compute_potential_mesh(result, self_term="exclude")

    np.testing.assert_allclose(potential, np.array([0.0]))


def test_compute_potential_mesh_softened_point_matches_legacy_behavior() -> None:
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    charges = np.array([2.0e-9, -1.0e-9])
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
        charges=charges,
        triangles=triangles,
    )

    softening = 0.5
    distance = 3.0
    expected = K_COULOMB * np.array(
        [
            charges[0] / softening + charges[1] / np.sqrt(distance**2 + softening**2),
            charges[0] / np.sqrt(distance**2 + softening**2) + charges[1] / softening,
        ]
    )

    potential = compute_potential_mesh(
        result, softening=softening, self_term="softened_point"
    )

    np.testing.assert_allclose(potential, expected)


def test_compute_potential_mesh_supports_reference_point_difference() -> None:
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [3.0, 3.0, 0.0]],
        ]
    )
    charges = np.array([2.0e-9, -1.0e-9])
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
        charges=charges,
        triangles=triangles,
    )

    reference = np.array([0.0, 0.0, 6.0])
    phi = compute_potential_mesh(
        result,
        softening=0.5,
        self_term="softened_point",
        reference_point=reference,
    )
    phi_abs = compute_potential_mesh(
        result,
        softening=0.5,
        self_term="softened_point",
    )
    phi_ref = compute_potential_points(
        result,
        reference.reshape(1, 3),
        softening=0.5,
    )[0]

    np.testing.assert_allclose(phi, phi_abs - phi_ref)


def test_compute_potential_mesh_auto_uses_softening_from_config(tmp_path: Path) -> None:
    out = tmp_path / "run_softening_auto"
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
                "softening = 0.5",
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    distance = 3.0
    softening = 0.5
    charges = np.array([2.0e-9, -1.0e-9])
    expected = K_COULOMB * np.array(
        [
            charges[0] / softening + charges[1] / np.sqrt(distance**2 + softening**2),
            charges[0] / np.sqrt(distance**2 + softening**2) + charges[1] / softening,
        ]
    )

    potential = compute_potential_mesh(result)

    np.testing.assert_allclose(potential, expected)


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
                "softening = 0.0",
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
    phi = compute_potential_mesh(
        result,
        self_term="exclude",
        reference_point="species1_injection_center",
    )

    distance = np.sqrt((2.0 / 3.0) ** 2 + (2.0 / 3.0) ** 2 + 10.0**2)
    expected = np.array([-K_COULOMB * 2.0e-9 / distance])
    np.testing.assert_allclose(phi, expected)


def test_compute_potential_mesh_allows_zero_softening_for_non_softened_self_terms() -> (
    None
):
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
        charges=np.array([1.0e-9, -1.0e-9]),
        triangles=triangles,
    )

    potential_area = compute_potential_mesh(
        result, softening=0.0, self_term="area_equivalent"
    )
    potential_exclude = compute_potential_mesh(
        result, softening=0.0, self_term="exclude"
    )

    assert np.all(np.isfinite(potential_area))
    assert np.all(np.isfinite(potential_exclude))


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


def test_compute_potential_mesh_rejects_negative_softening() -> None:
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

    with pytest.raises(ValueError, match="softening"):
        compute_potential_mesh(result, softening=-1.0)


def test_compute_potential_mesh_rejects_invalid_self_term() -> None:
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

    with pytest.raises(ValueError, match="self_term"):
        compute_potential_mesh(result, self_term="invalid")


def test_compute_potential_mesh_requires_positive_softening_for_softened_point() -> (
    None
):
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

    with pytest.raises(ValueError, match="softening"):
        compute_potential_mesh(result, softening=0.0, self_term="softened_point")


def test_compute_potential_mesh_handles_degenerate_triangle_in_area_equivalent() -> (
    None
):
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]])
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

    potential = compute_potential_mesh(result, self_term="area_equivalent")

    assert np.isfinite(potential[0])


def test_compute_potential_points_matches_centroid_point_charge_model() -> None:
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
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
        charges=charge,
        triangles=triangles,
    )
    points = np.array([[0.0, 0.0, 2.0], [0.0, 0.0, 4.0]])

    potential = compute_potential_points(result, points)

    expected = K_COULOMB * charge[0] / np.array([2.0, 4.0])
    np.testing.assert_allclose(potential, expected)


def test_compute_potential_points_supports_periodic2_image_sum() -> None:
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
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
        charges=charge,
        triangles=triangles,
    )
    points = np.array([[0.0, 0.0, 2.0]])
    periodic2 = {"axes": (0, 1), "lengths": (1.0, 1.0), "image_layers": 1}

    potential = compute_potential_points(result, points, periodic2=periodic2)

    expected_sum = 0.0
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            radius = np.sqrt(float(ix * ix + iy * iy) + 4.0)
            expected_sum += charge[0] / radius
    np.testing.assert_allclose(potential, np.array([K_COULOMB * expected_sum]))


def test_compute_potential_mesh_supports_periodic2_image_sum() -> None:
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
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
        charges=charge,
        triangles=triangles,
    )
    periodic2 = {"axes": (0, 1), "lengths": (1.0, 1.0), "image_layers": 1}

    potential = compute_potential_mesh(
        result,
        self_term="exclude",
        periodic2=periodic2,
    )

    expected_sum = 0.0
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            if ix == 0 and iy == 0:
                continue
            radius = np.sqrt(float(ix * ix + iy * iy))
            expected_sum += charge[0] / radius
    np.testing.assert_allclose(potential, np.array([K_COULOMB * expected_sum]))


def test_compute_potential_points_auto_detects_periodic2_from_config(tmp_path: Path) -> None:
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
                'field_periodic_far_correction = "none"',
            ]
        ),
        encoding="utf-8",
    )

    result = load_fortran_result(out)
    points = np.array([[0.0, 0.0, 2.0]])
    potential = compute_potential_points(result, points)

    expected_sum = 0.0
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            radius = np.sqrt(float(ix * ix + iy * iy) + 4.0)
            expected_sum += 2.0e-9 / radius
    np.testing.assert_allclose(potential, np.array([K_COULOMB * expected_sum]))


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

    expected_sum = 0.0
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            expected_sum += 2.0e-9 / np.sqrt(float(ix * ix + iy * iy) + 1.0)
    expected = np.array([K_COULOMB * expected_sum, K_COULOMB * expected_sum])
    np.testing.assert_allclose(potential, expected)


def test_potential_history_supports_periodic2_image_sum() -> None:
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charges_history = np.array([[2.0e-9, -1.0e-9]])
    periodic2 = _coerce_periodic2(
        {"axes": (0, 1), "lengths": (1.0, 1.0), "image_layers": 1}
    )
    assert periodic2 is not None

    potential = _potential_history(
        charges_history,
        triangles,
        softening=0.0,
        self_term="exclude",
        periodic2=periodic2,
    )

    image_sum = 0.0
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            if ix == 0 and iy == 0:
                continue
            image_sum += 1.0 / np.sqrt(float(ix * ix + iy * iy))
    expected = K_COULOMB * charges_history * image_sum
    np.testing.assert_allclose(potential, expected)


def test_potential_history_supports_reference_point_difference() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    charges_history = np.array([[2.0e-9, 4.0e-9]])
    reference_point = np.array([0.0, 0.0, 2.0])

    potential = _potential_history(
        charges_history,
        triangles,
        softening=0.0,
        self_term="exclude",
        reference_point=reference_point,
    )

    distance = np.sqrt((1.0 / 3.0) ** 2 + (1.0 / 3.0) ** 2 + 2.0**2)
    expected = -K_COULOMB * charges_history / distance
    np.testing.assert_allclose(potential, expected)


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


def test_compute_potential_slices_returns_xy_yz_xz_grids() -> None:
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
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

    center_expected = K_COULOMB * charge[0] / 0.5
    np.testing.assert_allclose(slices["xy"].potential_V[2, 2], center_expected)
    np.testing.assert_allclose(slices["yz"].potential_V[2, 2], center_expected)
    np.testing.assert_allclose(slices["xz"].potential_V[2, 2], center_expected)


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


def test_plot_potential_slices_returns_figure_and_axes() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

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


def test_plot_potential_slices_applies_vmin_vmax() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

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


def test_plot_potential_slices_rejects_invalid_vmin_vmax() -> None:
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

    with pytest.raises(ValueError, match="vmin"):
        plot_potential_slices(
            result,
            box_min=[-1.0, -1.0, -1.0],
            box_max=[1.0, 1.0, 1.0],
            grid_n=8,
            vmin=1.0,
            vmax=1.0,
        )


def test_plot_potential_mesh_returns_figure_and_axes() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
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
        charges=np.array([1.0e-9, -1.0e-9]),
        triangles=triangles,
    )

    fig, ax = plot_potential_mesh(result, softening=0.5, self_term="softened_point")

    assert fig is not None
    assert ax is not None
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


def test_plot_mesh_source_boxplot_supports_potential_quantity() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 0.8, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
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
        charges=np.array([1.0e-9, -2.0e-9]),
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    fig, ax = plot_mesh_source_boxplot(
        result,
        quantity="potential",
        softening=0.1,
        self_term="softened_point",
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
    assert "id=1 (template/plane)" in labels
    assert "id=2 (template/box)" in labels
    assert "id=3 (template/sphere)" in labels
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
    np.testing.assert_allclose(
        matrix_info["matrix"],
        K_COULOMB
        * np.array(
            [
                [0.0, 2.0e-18, -0.75e-18],
                [-2.0e-18, 0.0, -6.0e-18],
                [0.75e-18, 6.0e-18, 0.0],
                [-1.25e-18, 8.0e-18, -6.75e-18],
            ]
        ),
    )
    fig.clf()


def test_plot_coulomb_force_matrix_accepts_explicit_target_kinds(tmp_path: Path) -> None:
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
    np.testing.assert_allclose(
        matrix_info["matrix"][:, 0],
        K_COULOMB * np.array([0.0, -2.0e-18, 0.75e-18, -1.25e-18]),
    )
    fig.clf()


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
        directory=Path("dummy"),
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
        softening=0.5,
        self_term="softened_point",
    )

    assert written == out
    assert out.exists()
    assert out.stat().st_size > 0
