"""Result analysis contracts."""

from beach.io.types import RunMesh, RunSummary

from pathlib import Path
import numpy as np
import pytest
import beach.analysis.potential as potential_module
from beach import (
    analyze_coulomb_mobility,
    Beach,
    calc_coulomb,
    CoulombMobilityAnalysis,
    FieldKernelOptions,
    RunResult,
    compute_potential_points,
    compute_potential_slices,
    compute_potential_mesh,
    compute_electric_field_points,
    load_run,
)
from beach.analysis.mesh import (
    _surface_charge_density,
    _wrap_periodic2_triangles_by_centroid,
    _wrap_periodic2_triangles_by_mesh_centroid,
)
from beach.analysis.objects import resolve_object_specs
from beach.analysis.potential import _potential_history
from beach.native.periodic2 import coerce_periodic2
from result_fixtures import (
    _UnitPanelKernel,
    _write_complete_free_field_config,
    _write_coulomb_matrix_fixture,
    _write_minimal_result_fixture,
    _write_mobility_fixture,
    _write_three_mesh_fixture,
    _unit_panel_kernel as _unit_panel_kernel,
)



def _record_potential_kernels(monkeypatch: pytest.MonkeyPatch) -> list[_UnitPanelKernel]:
    """Install a unit panel kernel that remembers every kernel the potential API builds."""

    built: list[_UnitPanelKernel] = []

    class _RecordingKernel(_UnitPanelKernel):
        def __init__(self, *args, **kwargs) -> None:
            super().__init__(*args, **kwargs)
            built.append(self)

    monkeypatch.setattr(potential_module, "FieldKernel", _RecordingKernel)
    return built

@pytest.mark.parametrize(
    ("field_source_model", "summary_extra", "expected_source_model"),
    [
        pytest.param("point", ["field_kernel_id=softened_point"], "point", id="point"),
        pytest.param(None, None, "unknown", id="missing-receipt"),
    ],
)
def test_panel_estimators_reject_incompatible_source_model_receipt(
    tmp_path: Path,
    field_source_model: str | None,
    summary_extra: list[str] | None,
    expected_source_model: str,
) -> None:
    out = tmp_path / "run_incompatible_source_model"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=summary_extra,
        field_source_model=field_source_model,
    )
    result = load_run(out)

    assert result.field.source_model == expected_source_model
    with pytest.raises(ValueError, match="triangle_p0"):
        compute_potential_mesh(result)
    with pytest.raises(ValueError, match="triangle_p0"):
        compute_electric_field_points(result, np.zeros((1, 3)))
    with pytest.raises(ValueError, match="triangle_p0"):
        calc_coulomb(result, 1, 2, step=None)


def test_resolve_object_specs_prefers_mesh_source_materials(tmp_path: Path) -> None:
    out = tmp_path / "run_materials"
    out.mkdir()
    _write_mobility_fixture(out)

    result = load_run(out)
    specs = resolve_object_specs(result, config_path=out / "beach.toml")

    assert [spec.kind for spec in specs] == ["plane", "sphere"]
    assert [spec.surface_model for spec in specs] == ["insulator", "conductor"]
    assert [spec.epsilon_r for spec in specs] == [1.0, 2.5]


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

    result = load_run(out)
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
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9, -1.0e-9]),
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(mesh_potential_v=np.array([3.0, -4.0])),
    )

    potential = compute_potential_mesh(result)

    np.testing.assert_allclose(potential, np.array([3.0, -4.0]))


def test_compute_potential_mesh_source_and_kernel_choose_the_evaluation(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_potential_source"
    out.mkdir()
    _write_complete_free_field_config(out)
    _write_minimal_result_fixture(
        out, mesh_potential_text="elem_idx,potential_V\n1,3.0\n2,-4.0\n"
    )
    result = load_run(out)
    kernel = _UnitPanelKernel(
        result.mesh.triangles, result.charges, options=FieldKernelOptions()
    )

    saved = compute_potential_mesh(result)
    recomputed = compute_potential_mesh(result, source="kernel")
    borrowed = compute_potential_mesh(result, kernel=kernel)

    np.testing.assert_allclose(saved, np.array([3.0, -4.0]))
    np.testing.assert_allclose(compute_potential_mesh(result, source="saved"), saved)
    assert not np.allclose(recomputed, saved)
    np.testing.assert_allclose(borrowed, recomputed)
    assert not kernel.closed
    with pytest.raises(ValueError, match='source="saved" cannot be combined'):
        compute_potential_mesh(result, source="saved", reference_point=[0.0, 0.0, 2.0])
    (out / "mesh_potential.csv").unlink()
    with pytest.raises(ValueError, match="requires mesh_potential.csv"):
        compute_potential_mesh(load_run(out), source="saved")


def test_calc_coulomb_with_run_kernel_matches_source_only_plan_and_restores_it(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_coulomb_kernel"
    out.mkdir()
    _write_three_mesh_fixture(out)
    result = load_run(out)
    kernel = _UnitPanelKernel(
        result.mesh.triangles, result.charges, options=FieldKernelOptions()
    )

    expected = calc_coulomb(result, 1, [2, 3], step=1)
    borrowed = calc_coulomb(result, 1, [2, 3], step=1, kernel=kernel)

    np.testing.assert_allclose(borrowed.force_on_a_N, expected.force_on_a_N)
    np.testing.assert_allclose(borrowed.torque_on_a_Nm, expected.torque_on_a_Nm)
    np.testing.assert_allclose(kernel.source_charges, result.charges)
    assert not kernel.closed


def test_compute_potential_mesh_requires_triangles_when_precomputed_output_is_incompatible() -> (
    None
):
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(mesh_potential_v=np.array([2.0])),
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
    result = RunResult(
        directory=tmp_path,
        charges=charges,
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    scaled_result = RunResult(
        directory=tmp_path,
        charges=-0.25 * charges,
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    potential = compute_potential_mesh(result)
    scaled = compute_potential_mesh(scaled_result)

    np.testing.assert_allclose(scaled, -0.25 * potential)


@pytest.mark.parametrize(
    ("removed_keyword", "value"),
    [
        pytest.param("softening", 1.0, id="softening"),
        pytest.param("self_term", "softened_point", id="self-term"),
    ],
)
def test_compute_potential_mesh_rejects_removed_keywords(
    removed_keyword: str,
    value: object,
) -> None:
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
    )
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([0.0]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    with pytest.raises(TypeError, match=removed_keyword):
        compute_potential_mesh(  # type: ignore[call-arg]
            result,
            **{removed_keyword: value},
        )


def test_compute_potential_mesh_includes_finite_panel_self_term(
    tmp_path: Path,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]])
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    potential = compute_potential_mesh(result)

    assert potential.shape == (1,)
    assert np.isfinite(potential[0])
    assert potential[0] > 0.0


def test_compute_potential_mesh_uses_panel_centroids_as_targets(
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
    result = RunResult(
        directory=tmp_path,
        charges=charges,
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    mesh_potential = compute_potential_mesh(result)
    centroid_potential = compute_potential_points(result, triangles.mean(axis=1))

    np.testing.assert_allclose(mesh_potential, centroid_potential)


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
    result = RunResult(
        directory=tmp_path,
        charges=charges,
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
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
                "tree_theta = 0.37",
            ]
        ),
        encoding="utf-8",
    )

    result = load_run(out)

    class _ThetaProbeKernel(_UnitPanelKernel):
        def __init__(
            self, source_triangles, source_charges, *, options, **kwargs
        ) -> None:
            self.theta = float(options.theta)
            super().__init__(
                source_triangles,
                source_charges,
                options=options,
                **kwargs,
            )

        def eval_phi(self, points: np.ndarray) -> np.ndarray:
            return np.full(np.asarray(points).shape[0], self.theta)

    monkeypatch.setattr(potential_module, "FieldKernel", _ThetaProbeKernel)
    potential = compute_potential_mesh(result)

    np.testing.assert_allclose(potential, np.array([0.37, 0.37]))


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

    result = load_run(out)
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


def test_compute_potential_mesh_requires_triangles() -> None:
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([0.0]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=None),
    )

    with pytest.raises(ValueError, match="mesh_triangles.csv"):
        compute_potential_mesh(result)


def test_compute_potential_mesh_rejects_degenerate_panel(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]])
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    with pytest.raises(ValueError, match="non-degenerate"):
        compute_potential_mesh(result)


def test_compute_potential_points_returns_finite_panel_potential(
    tmp_path: Path,
) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charge = np.array([2.0e-9])
    result = RunResult(
        directory=tmp_path,
        charges=charge,
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
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
    result = RunResult(
        directory=tmp_path,
        charges=charge,
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
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
    result = RunResult(
        directory=tmp_path,
        charges=charge,
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
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


@pytest.mark.parametrize(
    "far_correction",
    [None, "none", "auto", "m2l_root_oracle"],
)
def test_compute_potential_points_normalizes_periodic2_config_far_correction(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    far_correction: str | None,
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
                *(
                    [f'field_periodic_far_correction = "{far_correction}"']
                    if far_correction is not None
                    else []
                ),
            ]
        ),
        encoding="utf-8",
    )

    result = load_run(out)
    points = np.array([[0.0, 0.0, 0.5]])
    built = _record_potential_kernels(monkeypatch)
    compute_potential_points(result, points)
    compute_potential_points(
        result,
        points,
        periodic2={
            "axes": (0, 1),
            "lengths": (1.0, 1.0),
            "image_layers": 1,
            "far_correction": "none",
        },
    )

    configured, explicit = (kernel.options.periodic2 for kernel in built)
    assert configured.far_correction == "none"
    assert configured == explicit


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
    with pytest.raises(ValueError, match='was removed; use "none"'):
        explicit_beach.plot_coulomb_force_matrix(component="x")


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

    result = load_run(out)
    points = np.array(
        [
            [0.25, 0.25, 1.0],
            [1.25, -0.75, 1.0],
        ]
    )

    potential = compute_potential_points(result, points)

    assert np.all(np.isfinite(potential))
    np.testing.assert_allclose(potential[0], potential[1])


def test_potential_history_reuses_panel_kernel_with_periodic2(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    charges_history = np.array([[2.0e-9, -1.0e-9]])
    periodic2 = coerce_periodic2(
        {"axes": (0, 1), "lengths": (1.0, 1.0), "image_layers": 1}
    )
    assert periodic2 is not None
    built = _record_potential_kernels(monkeypatch)

    potential = _potential_history(
        charges_history,
        triangles,
        periodic2=periodic2,
    )

    # One plan serves every frame and is closed when the history is done.
    assert len(built) == 1
    assert built[0].periodic2 == periodic2
    assert built[0].closed

    free = _potential_history(charges_history, triangles)

    assert np.all(np.isfinite(potential))
    np.testing.assert_allclose(potential[:, 1], -0.5 * potential[:, 0])
    assert abs(potential[0, 0]) > abs(free[0, 0])


def test_compute_potential_points_normalizes_auto_far_correction_to_none(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _write_complete_free_field_config(tmp_path)
    result = RunResult(
        directory=tmp_path,
        charges=np.array([2.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(
            triangles=np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
        ),
    )
    points = np.array([[0.0, 0.0, 2.0]])
    periodic2 = {
        "axes": (0, 1),
        "lengths": (1.0, 1.0),
        "image_layers": 1,
        "far_correction": "auto",
    }

    built = _record_potential_kernels(monkeypatch)
    compute_potential_points(result, points, periodic2=periodic2)
    compute_potential_points(
        result,
        points,
        periodic2={**periodic2, "far_correction": "none"},
    )

    auto, explicit = (kernel.options.periodic2 for kernel in built)
    assert auto.far_correction == "none"
    assert auto == explicit


@pytest.mark.parametrize(
    ("far_correction", "match"),
    [
        pytest.param("ewald", "periodic2.far_correction must be", id="ewald"),
        pytest.param("m2l_root", "periodic2.far_correction must be", id="m2l-root"),
        pytest.param(
            "m2l_root_trunc",
            "periodic2.far_correction must be",
            id="m2l-root-trunc",
        ),
        pytest.param(
            "m2l_root_oracle",
            "periodic2.far_correction .* was removed",
            id="m2l-root-oracle",
        ),
    ],
)
def test_compute_potential_points_rejects_removed_periodic2_far_correction(
    tmp_path: Path,
    far_correction: str,
    match: str,
) -> None:
    _write_complete_free_field_config(tmp_path)
    result = RunResult(
        directory=tmp_path,
        charges=np.array([2.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(
            triangles=np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
        ),
    )

    with pytest.raises(ValueError, match=match):
        compute_potential_points(
            result,
            np.array([[0.0, 0.0, 2.0]]),
            periodic2={
                "axes": (0, 1),
                "lengths": (1.0, 1.0),
                "image_layers": 1,
                "far_correction": far_correction,
            },
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


def test_compute_potential_points_validates_shape_and_chunk_size() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
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


class _AffinePotentialKernel(_UnitPanelKernel):
    """phi = x + 10 y + 100 z, so every sample identifies its own coordinates."""

    def eval_phi(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(points, dtype=float) @ np.array([1.0, 10.0, 100.0])


def test_compute_potential_slices_returns_xy_yz_xz_grids(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(potential_module, "FieldKernel", _AffinePotentialKernel)
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 2.0, 0.0]]])
    result = RunResult(
        directory=tmp_path,
        charges=np.array([2.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )
    box_min = np.array([-1.0, -2.0, 0.0])
    box_max = np.array([3.0, 2.0, 8.0])
    grid_n = 5

    slices = compute_potential_slices(
        result,
        box_min=box_min,
        box_max=box_max,
        grid_n=grid_n,
        xy_z=6.0,
        yz_x=-0.5,
        xz_y=1.5,
    )

    gradient = {"x": 1.0, "y": 10.0, "z": 100.0}
    axes = {
        name: np.linspace(box_min[index], box_max[index], grid_n)
        for index, name in enumerate("xyz")
    }
    expected_layout = {
        "xy": ("x", "y", "z", 6.0),
        "yz": ("y", "z", "x", -0.5),
        "xz": ("x", "z", "y", 1.5),
    }
    assert set(slices) == set(expected_layout)
    for plane, (axis_u, axis_v, fixed_axis, fixed_value) in expected_layout.items():
        slc = slices[plane]
        assert (slc.plane, slc.axis_u, slc.axis_v, slc.fixed_axis) == (
            plane,
            axis_u,
            axis_v,
            fixed_axis,
        )
        assert slc.fixed_value_m == fixed_value
        np.testing.assert_allclose(slc.u_values_m, axes[axis_u])
        np.testing.assert_allclose(slc.v_values_m, axes[axis_v])
        # potential_V has shape (n_v, n_u): rows follow v and columns follow u.
        expected = (
            gradient[axis_u] * axes[axis_u][None, :]
            + gradient[axis_v] * axes[axis_v][:, None]
            + gradient[fixed_axis] * fixed_value
        )
        np.testing.assert_allclose(slc.potential_V, expected)


def test_compute_potential_slices_rejects_invalid_grid_or_coordinate() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
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


def test_compute_potential_mesh_accepts_a_single_use_reference(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    result = RunResult(
        directory=tmp_path,
        summary=RunSummary(
            mesh_nelem=1, processed_particles=0, absorbed=0, escaped=0, batches=0
        ),
        charges=np.array([2.0e-9]),
        mesh=RunMesh(
            triangles=np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]]], dtype=float)
        ),
    )
    reference = [0.0, 0.0, 2.0]
    expected = compute_potential_mesh(result, reference_point=reference)
    actual = compute_potential_mesh(result, reference_point=(value for value in reference))
    np.testing.assert_array_equal(actual, expected)


def test_saved_potential_plot_rejects_an_explicit_reference_without_a_config(tmp_path: Path) -> None:
    from beach import plot_potential_mesh

    _write_minimal_result_fixture(tmp_path, mesh_potential_text="elem_idx,potential_V\n1,1.0\n2,2.0\n")
    with pytest.raises(ValueError, match='source="saved" cannot be combined'):
        plot_potential_mesh(load_run(tmp_path), source="saved", reference_point="species1_injection_center")
