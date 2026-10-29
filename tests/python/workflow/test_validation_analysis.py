"""Analysis of staged archives: estimator audits, artifacts, and release inputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from validation_support import (
    LEGACY_NATIVE_KEYS,
    _load_tool,
    _wrench_component,
    _write_release_mechanics_fixture,
)

from beach import RunDiagnostics, RunField, RunHistory, RunMesh, RunReceipts, RunSummary

pytestmark = pytest.mark.workflow


def _legacy_current_rows() -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    wrenches: list[dict[str, object]] = []
    curves: list[dict[str, object]] = []
    paths: list[dict[str, object]] = []
    radius = 3.5e-5
    for batch, mesh_id in LEGACY_NATIVE_KEYS:
        charge = (mesh_id - 5) * 1.0e-15
        wrenches.append(
            {
                "case": "archived_v1_3",
                "periodic_model": "configured",
                "resolved_batch": batch,
                "mesh_id": mesh_id,
                "component": "total_external",
                "total_charge_C": charge,
                "status": "available",
            }
        )
        for displacement in (0.0, 2.0 * radius):
            for component, scale in (
                ("other_objects_all_images", 0.5),
                ("target_periodic_images", 0.25),
                ("total_external", 0.75),
            ):
                curves.append(
                    {
                        "case": "archived_v1_3",
                        "periodic_model": "configured",
                        "resolved_batch": batch,
                        "mesh_id": mesh_id,
                        "component": component,
                        "displacement_m": displacement,
                        "force_z_N": scale * (mesh_id * 1.0e-12 + batch * 1.0e-18),
                        "status": "converged",
                    }
                )
        paths.append(
            {
                "case": "archived_v1_3",
                "periodic_model": "configured",
                "resolved_batch": batch,
                "mesh_id": mesh_id,
                "radius_m": radius,
                "endpoint_work_J": 4.0e-16,
                "status": "available",
            }
        )
    return wrenches, curves, paths


def test_legacy_estimator_audit_compares_exact_native_archive_keys(
    archive_run: Path,
) -> None:
    tool = _load_tool()
    wrenches, curves, paths = _legacy_current_rows()

    rows, summary = tool._legacy_estimator_audit(
        archive_run,
        wrenches,
        curves,
        paths,
        strict=True,
    )

    assert summary["status"] == "complete"
    assert summary["covered_native_keys"] == [
        [149001, 7],
        [180001, 6],
        [279001, 6],
        [279001, 7],
    ]
    assert len(rows) == 36
    assert {int(row["batch"]) for row in rows} == {149001, 180001, 279001}
    assert 280000 not in {int(row["batch"]) for row in rows}
    assert {row["estimator"] for row in rows} == {
        "direct_object_z",
        "local_pair_proxy",
        "moving_top_mesh_pairwise_coulomb",
    }
    assert {
        row["component"]
        for row in rows
        if row["comparison_kind"] == "legacy_moving_sphere_force"
    } == {
        "other_objects_all_images",
        "target_periodic_images",
        "total_external",
    }
    assert all(row["status"] == "computed" for row in rows)


@pytest.mark.parametrize(
    "mutation",
    [
        "missing",
        "schema",
        "duplicate",
        "nonfinite",
        "charge",
        "radius",
        "curve_coverage",
        "endpoint",
        "work",
    ],
)
def test_legacy_estimator_audit_strictly_validates_inputs_and_mapping(
    archive_run: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    analysis = archive_run / "analysis/local_release"
    wrenches, curves, paths = _legacy_current_rows()
    if mutation == "missing":
        (analysis / "force_timeseries.csv").unlink()
    elif mutation == "schema":
        path = analysis / "moving_sphere_force_curves.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text(
            lines[0].replace(",net_force_z_n", "") + "\n" + "\n".join(lines[1:]) + "\n",
            encoding="utf-8",
        )
    elif mutation == "duplicate":
        path = analysis / "moving_sphere_force_curves.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("\n".join([*lines, lines[1]]) + "\n", encoding="utf-8")
    elif mutation == "nonfinite":
        path = analysis / "force_timeseries.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        fields = lines[1].split(",")
        fields[10] = "nan"
        lines[1] = ",".join(fields)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    elif mutation == "charge":
        wrenches[0]["total_charge_C"] = 9.0e-15
    elif mutation == "radius":
        paths[0]["radius_m"] = 9.0e-5
    elif mutation == "curve_coverage":
        curves.pop()
    else:
        path = analysis / "moving_sphere_release_summary.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        fields = lines[1].split(",")
        fields[12 if mutation == "endpoint" else 19] = "9.0e-6"
        lines[1] = ",".join(fields)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="legacy estimator"):
        tool._legacy_estimator_audit(
            archive_run,
            wrenches,
            curves,
            paths,
            strict=True,
        )


def test_legacy_estimator_audit_non_strict_reports_unavailable(
    archive_run: Path,
) -> None:
    tool = _load_tool()
    (archive_run / "analysis/local_release/force_timeseries.csv").unlink()

    rows, summary = tool._legacy_estimator_audit(
        archive_run,
        [],
        [],
        [],
        strict=False,
    )

    assert summary["status"] == "unavailable"
    assert rows[0]["status"] == "unavailable"


def test_analyze_archive_only_writes_stable_explicit_missing_artifacts(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    (archive_run / "work/latest/summary.txt").write_text(
        "mesh_nelem=3392\nprocessed_particles=100\nabsorbed=60\n"
        "escaped=40\nbatches=280000\nlast_rel_change=1.0e-3\n",
        encoding="utf-8",
    )
    (archive_run / "work/latest/charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-15\n",
        encoding="utf-8",
    )

    class FakeBeach:
        def __init__(self, output_dir, *, config_path=None):
            assert Path(output_dir) == archive_run / "work/latest"
            assert Path(config_path) == archive_run / "input/beach.toml"
            self.result = SimpleNamespace(
                charges=[1.0e-15, -2.0e-15],
                ledger=None,
                summary=RunSummary(
                    mesh_nelem=3392,
                    processed_particles=100,
                    absorbed=60,
                    escaped=40,
                    batches=280000,
                    escaped_boundary=39,
                    survived_max_step=1,
                    last_rel_change=1.0e-3,
                ),
                mesh=RunMesh(mesh_ids=None, mesh_potential_v=None),
                diagnostics=RunDiagnostics(charge_ledger_residual_c=None),
                receipts=RunReceipts(
                    model_fingerprint=None,
                    mesh_fingerprint=None,
                    species_fingerprint=None,
                ),
                history=RunHistory(),
                field=RunField(),
            )

    monkeypatch.setattr(tool, "Beach", FakeBeach)
    report = tool.analyze_validation(
        archive_run,
        validation_root,
        library=binary,
    )

    analysis = validation_root / "analysis"
    expected = {
        "run_summary.csv",
        "charge_history_pair.csv",
        "particle_ledger_pair.csv",
        "mesh_potential_pair.csv",
        "snapshot_manifest.csv",
        "object_wrench.csv",
        "object_path_curves.csv",
        "object_path_summary.csv",
        "finite_shell_convergence.csv",
        "comparison_matrix.csv",
        "legacy_estimator_comparison.csv",
        "analysis_summary.json",
        "review_ja.md",
        "artifacts.json",
    }
    assert expected.issubset({path.name for path in analysis.iterdir()})
    assert report["cases"]["archived_v1_3"]["status"] == "available"
    assert report["cases"]["new_finite_configured"]["status"] == "missing"
    assert report["cases"]["new_infinite_physical"]["status"] == "missing"
    assert report["schema_version"] == 2
    qualification = report["numerical_qualification_for_local_frozen_model"]
    assert qualification["status"] == "not_qualified"
    assert qualification["status_semantics"] == (
        "path_work_shell_on_fixed_saved_discretization"
    )
    assert qualification["saved_sphere_mesh_refinement_status"] == "not_evaluated"
    assert qualification["source_discretization_refinement_status"] == "not_evaluated"
    assert qualification["plane_oracle_used_as_sphere_error_bar"] is False
    review = (analysis / "review_ja.md").read_text(encoding="utf-8")
    assert "旧結果は二つの実行バイナリ" in review
    assert "新 finite: 未実行" in review
    assert "fixed-discretization path/work/shell gate" in review
    assert "保存済み sphere mesh/source refinement: not_evaluated" in review
    assert "平面 oracle の誤差を sphere error bar に流用しません" in review
    with (analysis / "comparison_matrix.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    assert rows[0]["status"] == "missing_comparison"


def test_analyze_available_archive_uses_public_wrench_path_and_shell_apis(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    cache_path = validation_root / "cache/periodic2/operator.bin"
    cache_path.write_bytes(b"analysis-cache\n")
    (archive_run / "work/latest/summary.txt").write_text(
        "mesh_nelem=2\nprocessed_particles=100\nabsorbed=60\n"
        "escaped=40\nbatches=280000\nlast_rel_change=1.0e-3\n",
        encoding="utf-8",
    )
    (archive_run / "work/latest/charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-15\n2,2.0e-15\n",
        encoding="utf-8",
    )
    calls: dict[str, list[object]] = {
        "models": [],
        "steps": [],
        "shell_meshes": [],
        "path_grids": [],
    }

    class FakeRelease:
        def __init__(self, grid: np.ndarray):
            self.gravity_work_J = grid * 0.1
            self.adhesion_work_J = grid * 0.2
            self.available_energy_J = grid * 0.3
            self.electrostatic_only_speed_m_s = grid + 1.0
            self.gravity_corrected_speed_m_s = grid + 0.5
            self.speed_m_s = grid + 0.25
            self.barrier_free_from_rest = True
            self.endpoint_positive = True
            self.endpoint_reachable_from_rest = True
            self.endpoint_speed_m_s = 2.0
            self.maximum_reachable_speed_m_s = 3.0
            self.instantaneous_force_margin_N = 4.0
            self.minimum_available_energy_J = -5.0e-19
            self.first_inaccessible_displacement_m = None
            self.endpoint_available_energy_J = 6.0e-18

    class FakePath:
        def __init__(self, grid: np.ndarray, *, qualified: bool):
            self.displacement_m = grid
            self.force_N = np.column_stack((grid * 0.0, grid * 0.0, grid + 3.0))
            self.torque_Nm = np.column_stack((grid + 1.0, grid * 0.0, grid * 0.0))
            self.electrostatic_work_J = grid * 3.0
            self.potential_difference_work_J = grid * 3.0 if qualified else None
            self.component_force_N = {
                "other_objects_all_images": self.force_N * 0.25,
                "target_periodic_images": self.force_N * 0.5,
                "external_uniform": self.force_N * 0.25,
                "total_external": self.force_N,
            }
            self.component_torque_Nm = {
                "other_objects_all_images": self.torque_Nm * 0.25,
                "target_periodic_images": self.torque_Nm * 0.5,
                "external_uniform": self.torque_Nm * 0.25,
                "total_external": self.torque_Nm,
            }
            self.status = "converged" if qualified else "not_converged"
            self.work_relative_mismatch = 0.0
            torque_origins = np.broadcast_to(
                np.array([1.0, 2.0, 3.0]), (grid.size, 3)
            ).copy()
            torque_origins[:, 2] += grid
            self.numerical_metadata = {
                "target_geometry_representation": "periodic2_mesh_connected",
                "vertex_bounding_radius_m": 3.5e-5,
                "torque_origin_policy": "moving_geometric_area_centroid",
                "torque_origin_m": torque_origins,
            }

        def evaluate_release(self, **_kwargs):
            return FakeRelease(self.displacement_m)

    class FakeProbe:
        def __init__(self, mesh_id: int, periodic_model: str):
            self.mesh_id = mesh_id
            self.periodic_model = periodic_model
            self.vertex_bounding_radius_m = 3.5e-5
            self.target_geometry_representation = "periodic2_mesh_connected"
            self.geometric_area_centroid_m = np.array([1.0, 2.0, 3.0])

        def wrench(self, **_kwargs):
            total = SimpleNamespace(
                force_N=np.array([1.0, 2.0, 3.0]),
                torque_Nm=np.array([4.0, 5.0, 6.0]),
                potential_energy_J=7.0,
            )
            other = _wrench_component(
                total.force_N * 0.25,
                total.torque_Nm * 0.25,
                total.potential_energy_J * 0.25,
            )
            images = _wrench_component(
                total.force_N * 0.5,
                total.torque_Nm * 0.5,
                total.potential_energy_J * 0.5,
            )
            uniform = _wrench_component(
                total.force_N * 0.25,
                total.torque_Nm * 0.25,
                total.potential_energy_J * 0.25,
            )
            return SimpleNamespace(
                total_charge_C=8.0e-15,
                force_N=total.force_N,
                torque_Nm=total.torque_Nm,
                torque_origin_m=np.array([1.0, 2.0, 3.0]),
                components={
                    "other_objects_all_images": other,
                    "target_periodic_images": images,
                    "external_uniform": uniform,
                    "total_external": total,
                },
                numerical_metadata={
                    "effective_far_correction": (
                        "none"
                        if self.periodic_model == "configured"
                        else "cached_kneq0"
                    ),
                    "target_geometry_representation": ("periodic2_mesh_connected"),
                    "vertex_bounding_radius_m": self.vertex_bounding_radius_m,
                    "torque_origin_policy": "moving_geometric_area_centroid",
                    "torque_origin_m": np.array([1.0, 2.0, 3.0]),
                    "periodic_cache": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "hit": True,
                            "build_count": 0,
                            "fingerprint": "analysis-cache-fingerprint",
                            "path": cache_path,
                        }
                    ),
                    "periodic_kneq0": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "force_N": np.array([0.1, 0.2, 0.3]),
                            "torque_Nm": np.array([0.4, 0.5, 0.6]),
                            "potential_energy_J": 0.7,
                        }
                    ),
                    "physical_k0": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "force_N": np.array([0.0, 0.0, 0.8]),
                            "torque_Nm": np.zeros(3),
                            "potential_energy_J": 0.9,
                        }
                    ),
                    "cached_kneq0_trace_correction": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "force_N": np.array([9.0, 8.0, 7.0]),
                            "torque_Nm": np.array([6.0, 5.0, 4.0]),
                            "potential_energy_J": 3.0,
                        }
                    ),
                    "primary_free_subtraction": {
                        "force_N": np.array([-0.1, -0.2, -0.3]),
                        "torque_Nm": np.array([-0.4, -0.5, -0.6]),
                        "potential_energy_J": -0.7,
                    },
                },
            )

        def vertical_path(self, displacement_m, **_kwargs):
            calls["path_grids"].append(np.asarray(displacement_m, dtype=float))
            return FakePath(
                np.asarray(displacement_m, dtype=float),
                qualified=self.periodic_model == "configured",
            )

    class FakeSnapshot:
        def __init__(self, model: str):
            self.model = model

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def object_probe(self, mesh_id: int):
            return FakeProbe(mesh_id, self.model)

    triangles = np.array(
        [
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]],
        ]
    )

    class FakeBeach:
        def __init__(self, output_dir, *, config_path=None):
            assert Path(output_dir) == archive_run / "work/latest"
            assert Path(config_path) == archive_run / "input/beach.toml"
            self.result = SimpleNamespace(
                charges=np.array([1.0e-15, 2.0e-15]),
                ledger=None,
                summary=RunSummary(
                    mesh_nelem=2,
                    processed_particles=100,
                    absorbed=60,
                    escaped=40,
                    batches=280000,
                    escaped_boundary=39,
                    survived_max_step=1,
                    last_rel_change=1.0e-3,
                ),
                mesh=RunMesh(
                    mesh_ids=np.array([6, 7]),
                    triangles=triangles,
                    mesh_potential_v=None,
                ),
                diagnostics=RunDiagnostics(charge_ledger_residual_c=None),
                receipts=RunReceipts(
                    model_fingerprint=None,
                    mesh_fingerprint=None,
                    species_fingerprint=None,
                ),
                history=RunHistory(
                    charge=SimpleNamespace(
                        has_data=True,
                        batch_indices=np.array([149001, 180001, 279001]),
                    )
                ),
                field=RunField(),
            )

        def object_interaction_snapshot(self, *, periodic_model, step, **_kwargs):
            calls["models"].append(periodic_model)
            calls["steps"].append(step)
            return FakeSnapshot(periodic_model)

    def fake_shell(_snapshot, probe, displacement_m, **_kwargs):
        calls["shell_meshes"].append(probe.mesh_id)
        assert np.asarray(displacement_m)[0] == 0.0
        return SimpleNamespace(
            image_layers=np.array([0, 1, 2, 3]),
            force_increment_error_N=np.array([1.0, 0.1, 0.01]),
            work_increment_error_J=np.array([2.0, 0.2, 0.02]),
            force_tail_proxy_N=np.array([10.0, 1.0, 0.1]),
            work_tail_proxy_J=np.array([20.0, 2.0, 0.2]),
            increment_converged=np.array([False, True, True]),
            reference_model="infinite_physical",
            reference_force_error_N=np.array([3.0, 2.0, 1.0, 0.1]),
            reference_work_error_J=np.array([6.0, 4.0, 2.0, 0.2]),
            reference_converged=np.array([False, False, True, True]),
            status="converged",
            selected_image_layers=3,
        )

    monkeypatch.setattr(tool, "Beach", FakeBeach)
    monkeypatch.setattr(tool, "finite_shell_convergence", fake_shell)
    report = tool.analyze_validation(
        archive_run,
        validation_root,
        library=binary,
    )

    model_steps = list(zip(calls["models"], calls["steps"]))
    expected_model_steps = {
        (model, step)
        for model in ("configured", "infinite_physical")
        for step in (149001, 180001, 279001, None)
    }
    assert len(model_steps) == len(expected_model_steps)
    assert set(model_steps) == expected_model_steps
    assert len(calls["shell_meshes"]) == 2
    assert set(calls["shell_meshes"]) == {6, 7}
    assert calls["path_grids"]
    assert all(grid[0] == 0.0 for grid in calls["path_grids"])
    assert all(grid[-1] == pytest.approx(7.0e-5) for grid in calls["path_grids"])
    with (validation_root / "analysis/object_wrench.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        wrench_rows = list(csv.DictReader(stream))
    assert any(
        row["mesh_id"] == "6"
        and row["resolved_batch"] == "280000"
        and row["effective_far_correction"] == "cached_kneq0"
        and row["zero_mode_policy"] == "exclude_k0"
        and row["lower_boundary_model"] == "e_bottom_zero"
        and row["periodic_cache_hit"] == "True"
        and row["periodic_cache_build_count"] == "0"
        and row["periodic_cache_fingerprint"] == "analysis-cache-fingerprint"
        and row["periodic_cache_path"] == str(cache_path)
        and row["periodic_cache_path_sha256"] == tool._sha256(cache_path)
        and row["periodic_model"] == "infinite_physical"
        and row["component"] == "numerical:physical_k0"
        and row["component_kind"] == "numerical_decomposition"
        and row["status"] == "available"
        for row in wrench_rows
    )
    assert all(
        row["torque_origin_policy"] == "moving_geometric_area_centroid"
        and float(row["torque_origin_x_m"]) == pytest.approx(1.0)
        and float(row["torque_origin_y_m"]) == pytest.approx(2.0)
        and float(row["torque_origin_z_m"]) == pytest.approx(3.0)
        for row in wrench_rows
        if row["status"] == "available"
    )
    assert all(
        float(row["model_radius_m"]) == pytest.approx(3.5e-5)
        and row["radius_source"] == "release_model_summary.assumptions.radius_m"
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["radius_relative_mismatch"]) == pytest.approx(0.0)
        and float(row["radius_relative_tolerance"]) == pytest.approx(5.0e-3)
        and row["target_geometry_representation"] == "periodic2_mesh_connected"
        for row in wrench_rows
        if row["status"] == "available"
    )
    assert any(
        row["component"] == "total_external"
        and row["component_kind"] == "physical_total"
        for row in wrench_rows
    )
    assert any(
        row["component"] == "target_periodic_images"
        and row["component_kind"] == "physical_additive"
        for row in wrench_rows
    )
    assert any(
        row["component"] == "numerical:cached_kneq0_trace_correction"
        and row["component_kind"] == "numerical_diagnostic_included"
        for row in wrench_rows
    )
    wrench_groups: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in wrench_rows:
        if row["status"] != "available":
            continue
        key = (
            row["case"],
            row["periodic_model"],
            row["resolved_batch"],
            row["mesh_id"],
        )
        wrench_groups.setdefault(key, []).append(row)
    for rows in wrench_groups.values():
        physical = {
            row["component"]
            for row in rows
            if row["component_kind"].startswith("physical_")
        }
        assert physical == set(tool.PHYSICAL_OBJECT_COMPONENTS)
        numerical = {
            row["component"].removeprefix("numerical:")
            for row in rows
            if row["component"].startswith("numerical:")
        }
        expected_numerical = (
            set(tool.CACHED_NUMERICAL_COMPONENTS)
            if rows[0]["effective_far_correction"] == "cached_kneq0"
            else {"primary_free_subtraction"}
        )
        assert numerical == expected_numerical
    with (validation_root / "analysis/object_path_curves.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        curve_rows = list(csv.DictReader(stream))
    curve_keys = [
        (
            row["case"],
            row["periodic_model"],
            row["resolved_batch"],
            row["mesh_id"],
            row["component"],
            row["displacement_m"],
        )
        for row in curve_rows
    ]
    assert len(curve_keys) == len(set(curve_keys))
    assert all(
        row["torque_origin_policy"] == "moving_geometric_area_centroid"
        and float(row["torque_origin_x_m"]) == pytest.approx(1.0)
        and float(row["torque_origin_y_m"]) == pytest.approx(2.0)
        and float(row["torque_origin_z_m"])
        == pytest.approx(3.0 + float(row["displacement_m"]))
        for row in curve_rows
    )
    curve_groups: dict[tuple[str, ...], dict[str, list[tuple[str, ...]]]] = {}
    for row in curve_rows:
        key = (
            row["case"],
            row["periodic_model"],
            row["resolved_batch"],
            row["mesh_id"],
        )
        signature = (
            row["displacement_m"],
            row["torque_origin_policy"],
            row["torque_origin_x_m"],
            row["torque_origin_y_m"],
            row["torque_origin_z_m"],
        )
        curve_groups.setdefault(key, {}).setdefault(row["component"], []).append(
            signature
        )
    for components in curve_groups.values():
        assert set(components) == set(tool.PHYSICAL_OBJECT_COMPONENTS)
        reference_signature = components["total_external"]
        assert all(values == reference_signature for values in components.values())
    assert all(
        float(row["model_radius_m"]) == pytest.approx(3.5e-5)
        and row["radius_source"] == "release_model_summary.assumptions.radius_m"
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["radius_relative_mismatch"]) == pytest.approx(0.0)
        and float(row["radius_relative_tolerance"]) == pytest.approx(5.0e-3)
        and row["target_geometry_representation"] == "periodic2_mesh_connected"
        for row in curve_rows
    )
    with (validation_root / "analysis/object_path_summary.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        path_rows = list(csv.DictReader(stream))
    assert any(
        row["mesh_id"] == "7"
        and row["resolved_batch"] == "280000"
        and row["periodic_model"] == "configured"
        and row["status"] == "available"
        and row["potential_work_available"] == "True"
        and row["numerically_qualified"] == "True"
        and row["barrier_free_from_rest"] == "True"
        and row["minimum_available_energy_J"] == "-5e-19"
        and row["endpoint_available_energy_J"] == "6e-18"
        for row in path_rows
    )
    declared_rows = [row for row in path_rows if row["status"] == "available"]
    assert declared_rows
    expected_mass = 4.0 * np.pi * 3000.0 * (3.5e-5) ** 3 / 3.0
    assert all(float(row["radius_m"]) == pytest.approx(3.5e-5) for row in declared_rows)
    assert all(
        float(row["mass_kg"]) == pytest.approx(expected_mass) for row in declared_rows
    )
    assert all(
        row["radius_source"] == "release_model_summary.assumptions.radius_m"
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["radius_relative_mismatch"]) == pytest.approx(0.0)
        and float(row["radius_relative_tolerance"]) == pytest.approx(5.0e-3)
        and row["target_geometry_representation"] == "periodic2_mesh_connected"
        and row["torque_origin_policy"] == "moving_geometric_area_centroid"
        and float(row["initial_torque_origin_x_m"]) == pytest.approx(1.0)
        and float(row["initial_torque_origin_y_m"]) == pytest.approx(2.0)
        and float(row["initial_torque_origin_z_m"]) == pytest.approx(3.0)
        and float(row["path_start_m"]) == pytest.approx(0.0)
        and float(row["path_end_m"]) == pytest.approx(7.0e-5)
        for row in declared_rows
    )
    assert any(
        row["mesh_id"] == "7"
        and row["periodic_model"] == "infinite_physical"
        and row["status"] == "available"
        and row["path_status"] == "not_converged"
        and row["potential_work_available"] == "False"
        and row["numerically_qualified"] == "False"
        and row["speed_status"] == "numerically_unqualified"
        for row in path_rows
    )
    with (validation_root / "analysis/finite_shell_convergence.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        shell_rows = list(csv.DictReader(stream))
    assert shell_rows[-1]["selected_image_layers"] == "3"
    assert shell_rows[-1]["force_tail_proxy_N"] == "0.1"
    assert shell_rows[-1]["work_tail_proxy_J"] == "0.2"
    assert shell_rows[-1]["reference_model"] == "infinite_physical"
    assert shell_rows[-1]["reference_force_error_N"] == "0.1"
    assert shell_rows[-1]["reference_work_error_J"] == "0.2"
    assert shell_rows[-1]["reference_converged"] == "True"
    assert all(
        float(row["model_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["path_start_m"]) == pytest.approx(0.0)
        and float(row["path_end_m"]) == pytest.approx(7.0e-5)
        for row in shell_rows
    )
    assert report["physics_evaluation"]["status"] == "available"
    review = (validation_root / "analysis/review_ja.md").read_text(encoding="utf-8")
    assert "mesh 6" in review
    assert "batch 280000" in review
    assert "fixed-discretization path/work/shell gate: path_work=True" in review
    assert "barrier on fixed saved discretization=True" in review
    assert "平面 oracle の誤差を sphere error bar に流用しません" in review
    assert "physical_k0 Fz" in review
    with (validation_root / "analysis/comparison_matrix.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        comparison_rows = list(csv.DictReader(stream))
    assert any(
        row["comparison_kind"] == "frozen_field_override"
        and row["metric"] == "force_z_N"
        and row["resolved_batch"] == "280000"
        and row["mesh_id"] == "6"
        for row in comparison_rows
    )
    assert any(
        row["comparison_kind"] == "frozen_field_override"
        and row["metric"] == "barrier_free_from_rest"
        for row in comparison_rows
    )
    assert any(
        row["comparison_kind"] == "frozen_field_override"
        and row["metric"] == "minimum_available_energy_J"
        for row in comparison_rows
    )

    strict_run = FakeBeach(
        archive_run / "work/latest",
        config_path=archive_run / "input/beach.toml",
    )
    _wrench, _curves, _paths, _shells, strict_physics = (
        tool._evaluate_object_physics_at_step(
            archive_run=archive_run,
            validation_root=validation_root,
            library=binary,
            run_rows=[{"case": "archived_v1_3", "status": "available"}],
            runs={"archived_v1_3": strict_run},
            step=None,
            evaluation_target_ids=(6,),
            run_shell=False,
            allow_geometry_radius_fallback=False,
            require_complete_contract=True,
        )
    )
    assert strict_physics["status"] == "available"
    assert strict_physics["successful_target_models"] == 2
    assert strict_physics["failures"] == []


def test_release_parameters_reject_invalid_vdw_config(tmp_path: Path) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "input/release_kernel_base.toml").write_text(
        '[adhesion]\nmodel = "vdw_work"\nhamaker_constant = 1.0e-19\n',
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="invalid vdw adhesion"):
        tool._release_parameters(archive, 3.5e-5)


def test_release_parameters_rejects_declared_geometry_radius_mismatch(
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    _write_release_mechanics_fixture(archive)

    within_tolerance = tool._release_parameters(
        archive,
        3.51e-5,
        allow_geometry_radius_fallback=False,
    )
    assert within_tolerance["radius_m"] == pytest.approx(3.5e-5)
    assert within_tolerance["geometry_radius_m"] == pytest.approx(3.51e-5)
    # Mass (~5e-10 kg) and work (~1e-16 J) sit at or below pytest.approx's default 1e-12 absolute
    # tolerance, so these physical values are compared relative only.
    assert within_tolerance["mass_kg"] == pytest.approx(
        4.0 * np.pi * 3000.0 * (3.5e-5) ** 3 / 3.0, rel=1.0e-12, abs=0.0
    )
    declared_coefficient = 3.0 * 0.1 * 1.0e-19 * (3.5e-5 / 2.0) / 6.0
    assert within_tolerance["adhesion_force_N"] == pytest.approx(
        declared_coefficient / (0.4e-9) ** 2, rel=1.0e-12, abs=0.0
    )
    assert within_tolerance["adhesion_work_J"] == pytest.approx(
        0.5 * declared_coefficient * (1.0 / 0.4e-9 - 1.0 / 10.0e-9), rel=1.0e-12, abs=0.0
    )

    with pytest.raises(tool.ValidationError, match="radius.*mismatch"):
        tool._release_parameters(
            archive,
            3.6e-5,
            allow_geometry_radius_fallback=False,
        )


def test_release_parameters_geometry_fallback_is_nonstrict_and_missing_only(
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "input/release_kernel_base.toml").write_text(
        '[adhesion]\nmodel = "none"\n',
        encoding="utf-8",
    )

    mechanics = tool._release_parameters(
        archive,
        3.5e-5,
        allow_geometry_radius_fallback=True,
    )

    assert mechanics["radius_m"] == pytest.approx(3.5e-5)
    assert mechanics["geometry_radius_m"] == pytest.approx(3.5e-5)
    assert mechanics["radius_source"] == "probe.vertex_bounding_radius_m_fallback"
    with pytest.raises(tool.ValidationError, match="declared dust radius is missing"):
        tool._release_parameters(
            archive,
            3.5e-5,
            allow_geometry_radius_fallback=False,
        )

    summary = archive / "analysis/local_release/release_model_summary.json"
    summary.parent.mkdir(parents=True)
    summary.write_text(
        json.dumps({"assumptions": {"radius_m": None}}),
        encoding="utf-8",
    )
    with pytest.raises(tool.ValidationError, match="declared dust radius"):
        tool._release_parameters(
            archive,
            3.5e-5,
            allow_geometry_radius_fallback=True,
        )


@pytest.mark.parametrize(
    ("representation", "invalid_radius", "message"),
    [
        ("periodic2_mesh_connected", True, "bounding radius"),
        ("periodic2_mesh_connected", np.array([3.5e-5]), "bounding radius"),
        ("periodic2_mesh_connected", None, "bounding radius is unavailable"),
        ("periodic2_mesh_connected", np.nan, "finite and positive"),
        ("periodic2_mesh_connected", 0.0, "finite and positive"),
        ("raw_saved", 3.5e-5, "representation mismatch"),
    ],
)
def test_probe_geometry_contract_rejects_invalid_radius_or_representation(
    representation: str,
    invalid_radius: object,
    message: str,
) -> None:
    tool = _load_tool()
    probe = SimpleNamespace(
        target_geometry_representation=representation,
        vertex_bounding_radius_m=invalid_radius,
    )

    with pytest.raises(tool.ValidationError, match=message):
        tool._probe_geometry_contract(probe)


def test_torque_origin_provenance_wraps_array_conversion_errors() -> None:
    tool = _load_tool()

    with pytest.raises(tool.ValidationError, match="torque origin metadata"):
        tool._torque_origin_provenance(
            {
                "torque_origin_policy": "moving_geometric_area_centroid",
                "torque_origin_m": object(),
            }
        )


def test_production_torque_origin_contract_tracks_final_displacement_grid() -> None:
    tool = _load_tool()
    probe = SimpleNamespace(geometric_area_centroid_m=np.array([1.0, 2.0, 3.0]))
    displacement = np.array([0.0, 0.1, 0.25])
    origins = np.broadcast_to(np.array([1.0, 2.0, 3.0]), (3, 3)).copy()
    origins[:, 2] += displacement

    tool._validate_production_torque_origin_contract(
        probe=probe,
        wrench_policy="moving_geometric_area_centroid",
        wrench_origin_m=origins[0],
        path_policy="moving_geometric_area_centroid",
        path_origin_m=origins,
        displacement_m=displacement,
    )
    origins[1, 0] += 1.0e-3
    with pytest.raises(tool.ValidationError, match="moving torque origins"):
        tool._validate_production_torque_origin_contract(
            probe=probe,
            wrench_policy="moving_geometric_area_centroid",
            wrench_origin_m=np.array([1.0, 2.0, 3.0]),
            path_policy="moving_geometric_area_centroid",
            path_origin_m=origins,
            displacement_m=displacement,
        )


@pytest.mark.parametrize(
    ("malformed_source", "message"),
    [
        ("adhesion_toml", "failed to read TOML"),
        ("release_summary_json", "release model summary"),
    ],
)
def test_release_parameters_rejects_malformed_assumption_files(
    tmp_path: Path,
    malformed_source: str,
    message: str,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "analysis/local_release").mkdir(parents=True)
    if malformed_source == "adhesion_toml":
        (archive / "input/release_kernel_base.toml").write_text(
            "[adhesion\nmodel = 'none'\n", encoding="utf-8"
        )
    else:
        (archive / "analysis/local_release/release_model_summary.json").write_text(
            "{invalid", encoding="utf-8"
        )

    with pytest.raises(tool.ValidationError, match=message):
        tool._release_parameters(archive, 3.5e-5)


@pytest.mark.parametrize(
    ("model", "contact", "cutoff", "message"),
    [
        ("unknown", 1.0e-9, 2.0e-9, "unsupported adhesion model"),
        ("vdw_work", 2.0e-9, 1.0e-9, "cutoff_distance"),
    ],
)
def test_release_parameters_rejects_unknown_or_nonphysical_adhesion(
    tmp_path: Path,
    model: str,
    contact: float,
    cutoff: float,
    message: str,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "input/release_kernel_base.toml").write_text(
        "\n".join(
            [
                "[adhesion]",
                f'model = "{model}"',
                "hamaker_constant = 1.0e-19",
                f"contact_distance = {contact}",
                f"cutoff_distance = {cutoff}",
                "roughness_factor = 0.5",
                "contact_count = 2.0",
                "peel_factor = 1.0",
                'contact_geometry = "sphere_sphere"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match=message):
        tool._release_parameters(archive, 3.5e-5)
