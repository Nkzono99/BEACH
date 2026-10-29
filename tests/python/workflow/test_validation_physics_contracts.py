"""Wrench, path, shell, comparison, and qualification contracts of the analysis."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from validation_support import (
    _load_tool,
    _wrench_component,
    _write_release_mechanics_fixture,
)

from beach import RunDiagnostics, RunField, RunHistory, RunMesh, RunReceipts, RunSummary

pytestmark = pytest.mark.workflow


def test_strict_wrench_component_contract_requires_complete_additive_sets() -> None:
    tool = _load_tool()
    first = _wrench_component(np.array([1.0, 0.0, 0.0]), np.zeros(3), 1.0)
    second = _wrench_component(np.array([0.0, 2.0, 0.0]), np.zeros(3), 2.0)
    third = _wrench_component(np.array([0.0, 0.0, 3.0]), np.zeros(3), 3.0)
    total = _wrench_component(np.array([1.0, 2.0, 3.0]), np.zeros(3), 6.0)
    components = {
        "other_objects_all_images": first,
        "target_periodic_images": second,
        "external_uniform": third,
        "total_external": total,
    }
    diagnostic = {
        "force_N": np.zeros(3),
        "torque_Nm": np.zeros(3),
        "potential_energy_J": 0.0,
    }
    cached_metadata = {
        name: diagnostic
        for name in (
            "periodic_kneq0",
            "physical_k0",
            "primary_free_subtraction",
            "cached_kneq0_trace_correction",
        )
    }

    def validate(
        component_values: dict[str, SimpleNamespace],
        metadata_values: dict[str, object],
        *,
        effective_far_correction: str,
        wrench_force_N: np.ndarray = total.force_N,
        wrench_torque_Nm: np.ndarray = total.torque_Nm,
    ) -> None:
        tool._validate_strict_wrench_component_contract(
            component_values,
            metadata_values,
            effective_far_correction=effective_far_correction,
            wrench_force_N=wrench_force_N,
            wrench_torque_Nm=wrench_torque_Nm,
        )

    validate(
        components,
        cached_metadata,
        effective_far_correction="cached_kneq0",
    )
    with pytest.raises(tool.ValidationError, match="physical wrench component set"):
        validate(
            {
                name: value
                for name, value in components.items()
                if name != "external_uniform"
            },
            cached_metadata,
            effective_far_correction="cached_kneq0",
        )
    with pytest.raises(
        tool.ValidationError, match="cached numerical wrench component set"
    ):
        validate(
            components,
            {
                name: value
                for name, value in cached_metadata.items()
                if name != "physical_k0"
            },
            effective_far_correction="cached_kneq0",
        )
    nonfinite_metadata = dict(cached_metadata)
    nonfinite_metadata["physical_k0"] = {
        "force_N": np.zeros(3),
        "torque_Nm": np.zeros(3),
        "potential_energy_J": np.nan,
    }
    with pytest.raises(tool.ValidationError, match="numerical.*potential energy"):
        validate(
            components,
            nonfinite_metadata,
            effective_far_correction="cached_kneq0",
        )
    with pytest.raises(tool.ValidationError, match="additive physical wrench"):
        inconsistent = dict(components)
        inconsistent["total_external"] = _wrench_component(
            np.array([1.0, 2.0, 4.0]), np.zeros(3), 6.0
        )
        validate(
            inconsistent,
            cached_metadata,
            effective_far_correction="cached_kneq0",
        )
    with pytest.raises(
        tool.ValidationError, match="finite numerical wrench component set"
    ):
        validate(
            components,
            cached_metadata,
            effective_far_correction="none",
        )
    validate(
        components,
        {"primary_free_subtraction": diagnostic},
        effective_far_correction="none",
    )
    with pytest.raises(tool.ValidationError, match="primary_free_subtraction"):
        validate(
            components,
            {},
            effective_far_correction="none",
        )
    with pytest.raises(tool.ValidationError, match="wrench total"):
        validate(
            components,
            cached_metadata,
            effective_far_correction="cached_kneq0",
            wrench_force_N=np.array([1.0, 2.0, 4.0]),
        )
    with pytest.raises(tool.ValidationError, match="wrench total"):
        validate(
            components,
            cached_metadata,
            effective_far_correction="cached_kneq0",
            wrench_torque_Nm=np.array([0.0, 0.0, 1.0]),
        )


def test_additive_identity_tolerance_scales_with_cancelling_components() -> None:
    tool = _load_tool()

    tool._require_additive_identity(
        np.array([1.0]),
        [np.array([1.0e16]), np.array([1.0]), np.array([-1.0e16])],
        label="cancelling test",
    )
    tool._require_additive_identity(
        np.array([2.0]),
        [np.array([1.0e308]), np.array([-1.0e308]), np.array([2.0])],
        label="large cancelling test",
    )


@pytest.mark.parametrize(
    "components",
    [
        [np.array([1.0e16]), np.array([1.0]), np.array([-1.0e16])],
        [np.array([1.0e308]), np.array([-1.0e308]), np.array([2.0])],
    ],
)
def test_additive_identity_rejects_false_total_under_cancellation(
    components: list[np.ndarray],
) -> None:
    tool = _load_tool()

    with pytest.raises(tool.ValidationError, match="additive components"):
        tool._require_additive_identity(
            np.array([0.0]),
            components,
            label="cancelling test",
        )


@pytest.mark.parametrize(
    ("total_charge", "potential_energy"),
    [(np.nan, 1.0), (1.0, np.inf)],
)
def test_wrench_row_rejects_nonfinite_scalar_fields(
    total_charge: float,
    potential_energy: float,
) -> None:
    tool = _load_tool()
    mechanics = {
        "radius_m": 3.5e-5,
        "radius_source": "release_model_summary.assumptions.radius_m",
        "geometry_radius_m": 3.5e-5,
        "radius_relative_mismatch": 0.0,
        "radius_relative_tolerance": 5.0e-3,
    }

    with pytest.raises(tool.ValidationError, match="finite"):
        tool._wrench_row(
            case="case",
            periodic_model="configured",
            effective_far_correction="none",
            zero_mode_policy="legacy_not_decomposed",
            lower_boundary_model="legacy_implicit",
            periodic_cache=None,
            step_selector="final",
            resolved_batch=1,
            mesh_id=6,
            mechanics=mechanics,
            target_geometry_representation="periodic2_mesh_connected",
            component="total_external",
            component_kind="physical_total",
            force=np.zeros(3),
            torque=np.zeros(3),
            torque_origin_policy="moving_geometric_area_centroid",
            torque_origin_m=np.zeros(3),
            total_charge=total_charge,
            potential_energy=potential_energy,
        )


def _strict_path_component_fixture() -> SimpleNamespace:
    first = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    second = np.array([[0.0, 2.0, 0.0], [0.0, 3.0, 0.0]])
    third = np.array([[0.0, 0.0, 3.0], [0.0, 0.0, 4.0]])
    total = first + second + third
    zeros = np.zeros_like(total)
    return SimpleNamespace(
        displacement_m=np.array([0.0, 1.0]),
        force_N=total.copy(),
        torque_Nm=zeros.copy(),
        component_force_N={
            "other_objects_all_images": first,
            "target_periodic_images": second,
            "external_uniform": third,
            "total_external": total,
        },
        component_torque_Nm={
            "other_objects_all_images": zeros,
            "target_periodic_images": zeros,
            "external_uniform": zeros,
            "total_external": zeros,
        },
    )


def test_strict_path_component_contract_accepts_complete_additive_components() -> None:
    tool = _load_tool()
    path = _strict_path_component_fixture()

    tool._validate_strict_path_component_contract(
        path,
        expected_start_m=0.0,
        expected_end_m=1.0,
    )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("missing_force", "path force component set"),
        ("missing_torque", "path torque component set"),
        ("additive_force", "additive physical path force"),
        ("additive_torque", "additive physical path torque"),
        ("total_force", "path total arrays"),
        ("total_torque", "path total arrays"),
    ],
)
def test_strict_path_component_contract_rejects_invalid_components(
    mutation: str,
    message: str,
) -> None:
    tool = _load_tool()
    path = _strict_path_component_fixture()
    if mutation == "missing_force":
        path.component_force_N.pop("external_uniform")
    elif mutation == "missing_torque":
        path.component_torque_Nm.pop("external_uniform")
    elif mutation == "additive_force":
        path.component_force_N["total_external"] = path.component_force_N[
            "total_external"
        ].copy()
        path.component_force_N["total_external"][0, 0] += 1.0
    elif mutation == "additive_torque":
        path.component_torque_Nm["total_external"] = path.component_torque_Nm[
            "total_external"
        ].copy()
        path.component_torque_Nm["total_external"][0, 0] += 1.0
    elif mutation == "total_force":
        path.force_N[0, 0] += 1.0
    else:
        path.torque_Nm[0, 0] += 1.0

    with pytest.raises(tool.ValidationError, match=message):
        tool._validate_strict_path_component_contract(
            path,
            expected_start_m=0.0,
            expected_end_m=1.0,
        )


@pytest.mark.parametrize(
    "displacement",
    [
        np.array([1.0e-6, 1.0]),
        np.array([0.0, 0.9]),
        np.array([0.0, 0.75, 0.5, 1.0]),
    ],
)
def test_strict_path_component_contract_rejects_invalid_final_grid(
    displacement: np.ndarray,
) -> None:
    tool = _load_tool()
    zeros = np.zeros((displacement.size, 3))
    path = SimpleNamespace(
        displacement_m=displacement,
        force_N=zeros,
        torque_Nm=zeros,
        component_force_N={name: zeros for name in tool.PHYSICAL_OBJECT_COMPONENTS},
        component_torque_Nm={name: zeros for name in tool.PHYSICAL_OBJECT_COMPONENTS},
    )

    with pytest.raises(tool.ValidationError, match="path displacement"):
        tool._validate_strict_path_component_contract(
            path,
            expected_start_m=0.0,
            expected_end_m=1.0,
        )


def _valid_shell_reference_contract() -> SimpleNamespace:
    return SimpleNamespace(
        image_layers=np.array([0, 1, 2, 3]),
        increment_converged=np.array([False, True, True]),
        status="converged",
        selected_image_layers=3,
        reference_model="infinite_physical",
        reference_converged=np.array([False, False, True, True]),
        reference_force_error_N=np.array([3.0, 2.0, 1.0, 0.1]),
        reference_work_error_J=np.array([6.0, 4.0, 2.0, 0.2]),
        force_increment_error_N=np.array([1.0, 0.1, 0.01]),
        work_increment_error_J=np.array([2.0, 0.2, 0.02]),
        force_tail_proxy_N=np.array([10.0, 1.0, 0.1]),
        work_tail_proxy_J=np.array([20.0, 2.0, 0.2]),
    )


def test_shell_reference_contract_accepts_selected_physical_reference() -> None:
    tool = _load_tool()
    shell = _valid_shell_reference_contract()

    tool._validate_shell_reference_contract(shell)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("layer_start", "start at zero"),
        ("layer_gap", "consecutive"),
        ("status", "status"),
        ("increment_gate", "successive increment"),
        ("reference_gate_penultimate", "successive physical reference"),
        ("reference_gate_final", "successive physical reference"),
        ("reference_negative", "reference errors.*nonnegative"),
        ("tail_negative", "increment/tail"),
        ("tail_shape", "increment/tail"),
        ("increment_integer", "boolean"),
        ("increment_nonfinite", "boolean"),
        ("reference_integer", "boolean"),
        ("reference_nonfinite", "boolean"),
        ("reference_model", "reference_model"),
        ("reference_shape", "physical reference arrays"),
        ("selected_not_final", "selected_image_layers"),
        ("nonconverged_selected", "cannot select"),
    ],
)
def test_shell_reference_contract_rejects_invalid_reference_or_selection(
    mutation: str,
    message: str,
) -> None:
    tool = _load_tool()
    shell = _valid_shell_reference_contract()
    if mutation == "layer_start":
        shell.image_layers = np.array([1, 2, 3, 4])
    elif mutation == "layer_gap":
        shell.image_layers = np.array([0, 1, 3, 3])
    elif mutation == "status":
        shell.status = "invalid"
    elif mutation == "increment_gate":
        shell.increment_converged = np.array([False, False, True])
    elif mutation == "reference_gate_penultimate":
        shell.reference_converged = np.array([False, False, False, True])
    elif mutation == "reference_gate_final":
        shell.reference_converged = np.array([False, False, True, False])
    elif mutation == "reference_negative":
        shell.reference_force_error_N = np.array([3.0, 2.0, -1.0, 0.1])
    elif mutation == "tail_negative":
        shell.force_tail_proxy_N = np.array([10.0, -1.0, 0.1])
    elif mutation == "tail_shape":
        shell.force_tail_proxy_N = np.array([10.0, 1.0])
    elif mutation == "increment_integer":
        shell.increment_converged = np.array([0, 1, 1])
    elif mutation == "increment_nonfinite":
        shell.increment_converged = np.array([False, np.nan, True])
    elif mutation == "reference_integer":
        shell.reference_converged = np.array([0, 0, 1, 1])
    elif mutation == "reference_nonfinite":
        shell.reference_converged = np.array([False, False, np.inf, True])
    elif mutation == "reference_model":
        shell.reference_model = "configured"
    elif mutation == "reference_shape":
        shell.reference_force_error_N = np.array([3.0, 2.0, 1.0])
    elif mutation == "selected_not_final":
        shell.selected_image_layers = 2
    else:
        shell.status = "not_converged"

    with pytest.raises(tool.ValidationError, match=message):
        tool._validate_shell_reference_contract(shell)


def test_object_physics_requires_both_target_meshes(tmp_path: Path) -> None:
    tool = _load_tool()

    class IncompleteRun:
        result = SimpleNamespace(
            summary=RunSummary(mesh_nelem=1, processed_particles=0, absorbed=0, escaped=0, batches=0),
            mesh=RunMesh(mesh_ids=np.array([6])),
            history=RunHistory(),
            field=RunField(),
            receipts=RunReceipts(),
            diagnostics=RunDiagnostics(),
        )

        def object_interaction_snapshot(self, **_kwargs):
            raise AssertionError("snapshot must not be built for incomplete targets")

    _wrench, _curves, _paths, _shells, physics = tool._evaluate_object_physics(
        archive_run=tmp_path / "archive",
        validation_root=tmp_path / "validation",
        library=tmp_path / "kernel.so",
        run_rows=[{"case": "archived_v1_3", "status": "available"}],
        runs={"archived_v1_3": IncompleteRun()},
    )

    assert physics["status"] == "not_evaluated"
    assert physics["successful_target_models"] == 0
    assert "required target mesh ids" in physics["failures"][0]["message"]


@pytest.mark.parametrize(
    ("metadata_override", "expected_message"),
    [
        ({"effective_far_correction": "free"}, "effective far correction mismatch"),
        (
            {"target_geometry_representation": "raw_saved"},
            "target geometry representation mismatch",
        ),
        ({"vertex_bounding_radius_m": 4.0e-5}, "metadata radius"),
        (
            {"torque_origin_m": np.array([1.0, 0.0, 0.0])},
            "metadata torque origin does not match",
        ),
    ],
)
def test_object_physics_rejects_invalid_wrench_contract_metadata(
    tmp_path: Path,
    metadata_override: dict[str, object],
    expected_message: str,
) -> None:
    tool = _load_tool()
    triangles = np.array(
        [
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]],
        ]
    )

    class FakeProbe:
        def __init__(self, periodic_model: str):
            self.periodic_model = periodic_model
            self.vertex_bounding_radius_m = 3.5e-5
            self.target_geometry_representation = "periodic2_mesh_connected"

        def wrench(self, **_kwargs):
            metadata = {
                "effective_far_correction": (
                    "none" if self.periodic_model == "configured" else "cached_kneq0"
                ),
                "target_geometry_representation": "periodic2_mesh_connected",
                "vertex_bounding_radius_m": 3.5e-5,
                "torque_origin_policy": "moving_geometric_area_centroid",
                "torque_origin_m": np.zeros(3),
            }
            metadata.update(metadata_override)
            return SimpleNamespace(
                total_charge_C=1.0e-15,
                force_N=np.zeros(3),
                torque_Nm=np.zeros(3),
                torque_origin_m=np.zeros(3),
                components={},
                numerical_metadata=metadata,
            )

        def vertical_path(self, *_args, **_kwargs):
            raise AssertionError("mismatched effective model must stop before path")

    class FakeSnapshot:
        def __init__(self, periodic_model: str):
            self.periodic_model = periodic_model

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def object_probe(self, _mesh_id: int):
            return FakeProbe(self.periodic_model)

    class FakeRun:
        result = SimpleNamespace(
            summary=RunSummary(
                batches=280000,
                mesh_nelem=0,
                processed_particles=0,
                absorbed=0,
                escaped=0,
            ),
            mesh=RunMesh(mesh_ids=np.array([6, 7]), triangles=triangles),
            history=RunHistory(),
            field=RunField(),
            receipts=RunReceipts(),
            diagnostics=RunDiagnostics(),
        )

        def object_interaction_snapshot(self, **kwargs):
            return FakeSnapshot(kwargs["periodic_model"])

    _wrench, _curves, _paths, _shells, physics = tool._evaluate_object_physics_at_step(
        archive_run=tmp_path / "archive",
        validation_root=tmp_path / "validation",
        library=tmp_path / "kernel.so",
        run_rows=[{"case": "archived_v1_3", "status": "available"}],
        runs={"archived_v1_3": FakeRun()},
        step=None,
        evaluation_target_ids=(6,),
        run_shell=False,
    )

    assert physics["status"] == "not_evaluated"
    assert physics["successful_target_models"] == 0
    assert all(
        expected_message in failure["message"] for failure in physics["failures"]
    )


def test_new_infinite_final_charge_gets_infinite_shell_reference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    archive_run = tmp_path / "archive"
    (archive_run / "input").mkdir(parents=True)
    _write_release_mechanics_fixture(archive_run)
    calls: list[tuple[object, str, int]] = []
    triangles = np.array(
        [
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]],
        ]
    )

    class FakeSnapshot:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def object_probe(self, mesh_id: int):
            return SimpleNamespace(
                mesh_id=mesh_id,
                vertex_bounding_radius_m=3.5e-5,
                target_geometry_representation="periodic2_mesh_connected",
                geometric_area_centroid_m=np.array([1.0, 2.0, 3.0]),
            )

    class FakeRun:
        result = SimpleNamespace(
            summary=RunSummary(
                batches=280000,
                mesh_nelem=0,
                processed_particles=0,
                absorbed=0,
                escaped=0,
            ),
            mesh=RunMesh(mesh_ids=np.array([6, 7]), triangles=triangles),
            history=RunHistory(),
            field=RunField(),
            receipts=RunReceipts(),
            diagnostics=RunDiagnostics(),
        )

        def object_interaction_snapshot(self, *, step, periodic_model, **_kwargs):
            calls.append((step, periodic_model, -1))
            return FakeSnapshot()

    def fake_shell(_snapshot, probe, displacement_m, **_kwargs):
        calls.append((None, "shell", probe.mesh_id))
        assert len(displacement_m) == 17
        assert displacement_m[-1] == pytest.approx(7.0e-5)
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

    monkeypatch.setattr(tool, "finite_shell_convergence", fake_shell)
    rows, failures = tool._evaluate_new_infinite_shell_reference(
        archive_run=archive_run,
        validation_root=tmp_path / "validation",
        library=tmp_path / "kernel.so",
        run=FakeRun(),
    )

    snapshot_calls = [call for call in calls if call[1] == "infinite_physical"]
    shell_calls = [call for call in calls if call[1] == "shell"]
    assert snapshot_calls == [(None, "infinite_physical", -1)]
    assert len(shell_calls) == 2
    assert {call[2] for call in shell_calls} == {6, 7}
    assert failures == []
    assert len(rows) == 8
    assert {(row["mesh_id"], row["resolved_batch"]) for row in rows} == {
        (6, 280000),
        (7, 280000),
    }
    assert all(row["periodic_model"] == "infinite_physical" for row in rows)
    assert all(row["effective_far_correction"] == "cached_kneq0" for row in rows)
    assert all(row["zero_mode_policy"] == "exclude_k0" for row in rows)
    assert all(row["lower_boundary_model"] == "e_bottom_zero" for row in rows)
    assert all(float(row["model_radius_m"]) == pytest.approx(3.5e-5) for row in rows)
    assert all(float(row["geometry_radius_m"]) == pytest.approx(3.5e-5) for row in rows)
    assert all(float(row["path_start_m"]) == pytest.approx(0.0) for row in rows)
    assert all(float(row["path_end_m"]) == pytest.approx(7.0e-5) for row in rows)

    def fake_shell_without_reference(*args, **kwargs):
        shell = fake_shell(*args, **kwargs)
        shell.reference_model = None
        return shell

    calls.clear()
    monkeypatch.setattr(tool, "finite_shell_convergence", fake_shell_without_reference)
    invalid_rows, invalid_failures = tool._evaluate_new_infinite_shell_reference(
        archive_run=archive_run,
        validation_root=tmp_path / "validation",
        library=tmp_path / "kernel.so",
        run=FakeRun(),
    )
    assert invalid_rows == []
    assert len(invalid_failures) == 2
    assert all("reference_model" in failure["message"] for failure in invalid_failures)


def test_charge_vector_comparison_preserves_spatial_difference() -> None:
    tool = _load_tool()
    left = np.array([1.0, 0.0])
    right = np.array([0.0, 1.0])

    rows = tool._charge_vector_comparison_rows(
        comparison="left_to_right",
        comparison_kind="boundary_history_response",
        sample_kind="history",
        resolved_batch=1001,
        scope="object",
        mesh_id=6,
        left_snapshot_id="left:1001",
        right_snapshot_id="right:1001",
        left=left,
        right=right,
        interpretation="test",
    )
    by_metric = {row["metric"]: row for row in rows}

    assert by_metric["charge_sum_C"]["signed_difference"] == 0.0
    assert by_metric["charge_sum_C"]["absolute_difference"] == 0.0
    assert by_metric["element_charge_l1_C"]["left_value"] == 1.0
    assert by_metric["element_charge_l1_C"]["right_value"] == 1.0
    assert by_metric["element_charge_l1_C"]["absolute_difference"] == 2.0
    assert by_metric["element_charge_l2_C"]["absolute_difference"] == pytest.approx(
        np.sqrt(2.0)
    )
    assert by_metric["element_charge_linf_C"]["absolute_difference"] == 1.0


def test_geometry_comparison_allows_declared_roundoff_not_order_change() -> None:
    tool = _load_tool()
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0e-3, 0.0, 0.0], [0.0, 1.0e-3, 0.0]],
            [[0.0, 0.0, 1.0e-3], [1.0e-3, 0.0, 1.0e-3], [0.0, 1.0e-3, 1.0e-3]],
        ]
    )
    left = SimpleNamespace(
        mesh=RunMesh(mesh_ids=np.array([6, 7]), triangles=triangles),
        history=RunHistory(),
        field=RunField(),
        receipts=RunReceipts(),
        diagnostics=RunDiagnostics(),
    )

    def compare_to(mesh_ids: np.ndarray, right_triangles: np.ndarray) -> dict:
        return tool._geometry_comparison(
            left,
            SimpleNamespace(
                mesh=RunMesh(mesh_ids=mesh_ids, triangles=right_triangles),
                history=RunHistory(),
                field=RunField(),
                receipts=RunReceipts(),
                diagnostics=RunDiagnostics(),
            ),
        )

    # Fixed offsets for the 1 mm mesh, independent of the tolerance the tool reports: 1e-18 m is
    # roundoff (a few ulp at the 1e-3 m scale); 1e-12 m (1e-9 of the mesh size) is a real change.
    roundoff = triangles.copy()
    roundoff[0, 0, 0] += 1.0e-18
    moved = triangles.copy()
    moved[0, 0, 0] += 1.0e-12
    within = compare_to(np.array([6, 7]), roundoff)
    outside = compare_to(np.array([6, 7]), moved)
    reordered = compare_to(np.array([7, 6]), triangles)

    assert within["status"] == "match"
    assert within["max_coordinate_difference_m"] == 1.0e-18
    assert outside["status"] == "coordinate_mismatch"
    assert outside["max_coordinate_difference_m"] == 1.0e-12
    assert reordered["status"] == "mesh_id_order_mismatch"


def _complete_local_qualification_fixture(tool):
    paths: list[dict[str, object]] = []
    wrenches: list[dict[str, object]] = []
    for case, model, batch, mesh_id in sorted(tool._expected_object_evaluation_keys()):
        common = {
            "case": case,
            "periodic_model": model,
            "resolved_batch": batch,
            "mesh_id": mesh_id,
            "status": "available",
        }
        paths.append({**common, "numerically_qualified": True})
        wrenches.append({**common, "component": "total_external"})
    shell_rows = [
        {
            "case": case,
            "periodic_model": "infinite_physical",
            "resolved_batch": 280000,
            "mesh_id": mesh_id,
            "status": "converged",
        }
        for case in (
            "archived_v1_3",
            "new_finite_configured",
            "new_infinite_physical",
        )
        for mesh_id in (6, 7)
    ]
    cases = {
        case: {"status": "available"}
        for case in (
            "archived_v1_3",
            "new_finite_configured",
            "new_infinite_physical",
        )
    }
    return paths, wrenches, shell_rows, cases


def test_local_qualification_reports_fixed_discretization_scope_when_qualified() -> (
    None
):
    tool = _load_tool()
    paths, wrenches, shell_rows, cases = _complete_local_qualification_fixture(tool)

    qualification = tool._local_model_numerical_qualification(
        {"status": "available"},
        wrenches,
        paths,
        shell_rows,
        cases,
    )

    assert qualification["status"] == "qualified"
    assert qualification["status_semantics"] == (
        "path_work_shell_on_fixed_saved_discretization"
    )
    assert qualification["verified_numerical_axes"] == [
        "fixed_saved_discretization_coverage",
        "path_integration",
        "work_potential_consistency",
        "force_and_barrier_decision_resolution",
        "finite_shell_to_infinite_reference",
    ]
    assert qualification["unverified_numerical_axes"] == [
        "saved_sphere_mesh_refinement",
        "source_discretization_refinement",
        "sphere_absolute_force_error",
        "sphere_absolute_torque_error",
    ]
    assert qualification["saved_sphere_mesh_refinement_status"] == "not_evaluated"
    assert qualification["source_discretization_refinement_status"] == "not_evaluated"
    assert qualification["sphere_absolute_force_error_status"] == "not_evaluated"
    assert qualification["sphere_absolute_torque_error_status"] == "not_evaluated"
    assert qualification["plane_oracle_used_as_sphere_error_bar"] is False
    assert qualification["claim_scope"] == (
        "fixed_saved_mesh_and_source_discretization; "
        "local_frozen_field_0_to_2R_path_work_shell_only; "
        "not full_discretization_or_escape_to_infinity"
    )


def test_local_qualification_rejects_duplicate_keys_masking_missing_state() -> None:
    tool = _load_tool()
    paths, wrenches, shell_rows, cases = _complete_local_qualification_fixture(tool)
    paths[-1] = dict(paths[0])
    wrenches[-1] = dict(wrenches[0])

    qualification = tool._local_model_numerical_qualification(
        {"status": "available"},
        wrenches,
        paths,
        shell_rows,
        cases,
    )

    assert qualification["status"] == "not_qualified"
    assert qualification["duplicate_path_keys"]
    assert qualification["missing_path_keys"]
    assert qualification["duplicate_wrench_keys"]
    assert qualification["missing_wrench_keys"]


def test_comparison_contract_rejects_wrong_frozen_snapshot_semantics() -> None:
    tool = _load_tool()
    snapshots = [
        {"snapshot_id": "case:history:1", "status": "available"},
        {"snapshot_id": "other:history:1", "status": "available"},
    ]
    effective = {
        "archive_version_drift": ("none", "none"),
        "frozen_field_override": ("none", "cached_kneq0"),
        "boundary_history_response_common_evaluator": (
            "cached_kneq0",
            "cached_kneq0",
        ),
        "actual_end_to_end": ("none", "cached_kneq0"),
    }
    rows = [
        {
            "comparison": kind,
            "comparison_kind": kind,
            "metric": metric,
            "left_snapshot_id": "case:history:1",
            "right_snapshot_id": "case:history:1",
            "left_effective_far_correction": pair[0],
            "right_effective_far_correction": pair[1],
            "status": "computed",
        }
        for kind, pair in effective.items()
        for metric in (
            "force_z_N",
            "endpoint_work_J",
            "minimum_available_energy_J",
            "barrier_free_from_rest",
            "endpoint_reachable_from_rest",
        )
    ]
    assert tool._comparison_artifact_contract(rows, snapshots)["status"] == "complete"
    frozen = next(
        row for row in rows if row["comparison_kind"] == "frozen_field_override"
    )
    frozen["right_snapshot_id"] = "other:history:1"

    contract = tool._comparison_artifact_contract(rows, snapshots)
    assert contract["status"] == "incomplete"
    assert "different charge snapshots" in contract["semantic_errors"][0]


def test_run_comparisons_use_canonical_comparison_kinds() -> None:
    tool = _load_tool()
    common = {
        "status": "available",
        "processed_particles": 100,
        "absorbed": 60,
        "escaped": 40,
        "total_charge_C": 1.0,
    }
    rows = tool._comparison_rows(
        [
            {"case": "archived_v1_3", **common},
            {"case": "new_finite_configured", **common},
            {"case": "new_infinite_physical", **common},
        ]
    )

    assert {row["comparison_kind"] for row in rows} == {
        "archive_version_drift",
        "actual_end_to_end",
    }
