from pathlib import Path

import numpy as np
import pytest

from beach.fortran_results import (
    K_COULOMB,
    FortranRunResult,
    compute_potential_mesh,
    _surface_charge_density,
    list_fortran_runs,
    load_fortran_result,
    plot_potential_mesh,
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
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C\n"
        "1,0,0,0,1,0,0,0,1,0,1.0e-10\n"
        "2,0,0,1,1,0,1,0,1,1,-2.0e-10\n",
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
    assert result.charge_history is not None
    np.testing.assert_allclose(result.charges, np.array([1.0e-10, -2.0e-10]))
    np.testing.assert_allclose(
        result.charge_history,
        np.array([[2.0e-11, 1.0e-10], [-1.0e-11, -2.0e-10]]),
    )
    np.testing.assert_array_equal(
        result.processed_particles_by_batch, np.array([5, 10])
    )
    np.testing.assert_allclose(result.rel_change_by_batch, np.array([3.0e-1, 1.0e-8]))
    np.testing.assert_array_equal(result.batch_indices, np.array([1, 3]))


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

    potential = compute_potential_mesh(result)

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
