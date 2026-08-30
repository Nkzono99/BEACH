import importlib.util
import json
import os
import shutil
import sys
from collections.abc import Mapping
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from beach import (
    FieldKernel,
    FieldKernelOptions,
    FiniteShellConvergenceResult,
    FortranRunResult,
    ObjectForcePath,
    ObjectInteractionSnapshot,
    RigidTransform,
    finite_shell_convergence,
    finite_shell_wrench,
)
from beach.fortran_results.constants import K_COULOMB
import beach.fortran_results.periodic_force_oracle as oracle_module


EPS0 = 1.0 / (4.0 * np.pi * K_COULOMB)
ROOT = Path(__file__).resolve().parents[2]
VALIDATION_TOOL_PATH = ROOT / "tools" / "periodic_object_validation.py"


def _load_validation_tool():
    spec = importlib.util.spec_from_file_location(
        "beach_periodic_object_validation_native_oracle",
        VALIDATION_TOOL_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _kernel_lib() -> Path:
    path = Path("build/libbeach_field_kernel.so")
    if not path.exists():
        pytest.skip("field kernel shared library is not built; run `make build-kernel`")
    return path


def _triangles_at(positions: np.ndarray) -> np.ndarray:
    offsets = np.array(
        [[-0.02, -0.02, 0.0], [0.04, -0.02, 0.0], [-0.02, 0.04, 0.0]]
    )
    return np.asarray(positions, dtype=float)[:, None, :] + offsets[None, :, :]


def _result(
    directory: Path,
    charges: np.ndarray,
    *,
    positions: np.ndarray | None = None,
) -> FortranRunResult:
    if positions is None:
        positions = np.array([[0.5, 0.5, 0.0], [1.5, 0.5, 0.3]])
    return FortranRunResult(
        directory=directory,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.asarray(charges, dtype=float),
        triangles=_triangles_at(positions),
        mesh_ids=np.array([1, 2]),
        field_source_model="triangle_p0",
        field_kernel_id="triangle_p0_exact_p2m_near",
    )


def _write_config(path: Path) -> None:
    path.write_text(
        """
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
box_min = [0.0, 0.0, -1.0]
box_max = [2.0, 2.0, 1.0]
tree_theta = 0.1
tree_leaf_max = 64
tree_order = 2
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
field_periodic_ewald_layers = 4
field_periodic_generation_tolerance = 1.0e-6
e0 = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )


def _write_panel_periodic_config(path: Path) -> None:
    path.write_text(
        """
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
box_min = [0.0, 0.0, -1.0]
box_max = [2.0, 2.0, 1.0]
tree_theta = 0.2
tree_leaf_max = 16
tree_order = 4
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
field_periodic_ewald_layers = 4
field_periodic_generation_tolerance = 1.0e-8
e0 = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )


def _write_free_interaction_config(
    path: Path,
) -> None:
    path.write_text(
        """
[sim]
field_solver = "fmm"
field_bc_mode = "free"
tree_theta = 0.2
tree_leaf_max = 16
tree_order = 4
box_min = [-1.0, -1.0, -1.0]
box_max = [3.0, 3.0, 1.0]
e0 = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )


def _plane_triangles(cells_per_axis: int, *, length: float = 2.0) -> np.ndarray:
    triangles: list[np.ndarray] = []
    spacing = length / cells_per_axis
    for iy in range(cells_per_axis):
        for ix in range(cells_per_axis):
            x0 = ix * spacing
            x1 = (ix + 1) * spacing
            y0 = iy * spacing
            y1 = (iy + 1) * spacing
            p00 = np.array([x0, y0, 0.0])
            p10 = np.array([x1, y0, 0.0])
            p11 = np.array([x1, y1, 0.0])
            p01 = np.array([x0, y1, 0.0])
            triangles.extend((np.array([p00, p10, p11]), np.array([p00, p11, p01])))
    return np.asarray(triangles)


def _panel_result(
    directory: Path,
    triangles: np.ndarray,
    charges: np.ndarray,
) -> FortranRunResult:
    return FortranRunResult(
        directory=directory,
        mesh_nelem=len(charges),
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.asarray(charges, dtype=float),
        triangles=np.asarray(triangles, dtype=float),
        mesh_ids=np.ones(len(charges), dtype=np.int64),
        field_source_model="triangle_p0",
        field_kernel_id="triangle_p0_exact_p2m_near",
    )


def _physical_field(
    snapshot: ObjectInteractionSnapshot,
    points: np.ndarray,
) -> np.ndarray:
    field = snapshot._periodic.eval_e(points)
    assert snapshot._zero_mode is not None
    _, zero_ez = snapshot._zero_mode.eval(points[:, 2], trace="plus")
    field[:, 2] += zero_ez
    return field


def _constant_path(force_z: float) -> ObjectForcePath:
    return ObjectForcePath(
        displacement_m=np.array([0.0, 1.0]),
        force_N=np.array([[0.0, 0.0, force_z], [0.0, 0.0, force_z]]),
        torque_Nm=np.zeros((2, 3)),
        electrostatic_work_J=np.array([0.0, force_z]),
        potential_energy_J=np.array([0.0, -force_z]),
        potential_difference_work_J=np.array([0.0, force_z]),
        numerical_metadata={
            "relative_tolerance": 1.0e-2,
            "work_absolute_tolerance_J": 0.0,
            "status_reason": "fixed_grid",
        },
        status="converged",
        work_relative_mismatch=0.0,
        work_absolute_mismatch_J=0.0,
    )


def _valid_converged_shell_record_kwargs() -> dict[str, object]:
    paths = tuple(_constant_path(1.0) for _ in range(3))
    return {
        "image_layers": np.array([0, 1, 2]),
        "symmetric_paths": paths,
        "corrected_paths": paths,
        "force_increment_error_N": np.zeros(2),
        "work_increment_error_J": np.zeros(2),
        "increment_converged": np.ones(2, dtype=bool),
        "status": "converged",
        "selected_image_layers": 2,
        "selected_path": paths[-1],
        "reference_force_error_N": np.zeros(3),
        "reference_work_error_J": np.zeros(3),
        "reference_converged": np.ones(3, dtype=bool),
        "reference_model": "infinite_physical",
    }


def _bottom_zero_closure_context(
    total_charge: float,
) -> tuple[SimpleNamespace, SimpleNamespace]:
    snapshot = SimpleNamespace(
        _options=SimpleNamespace(
            periodic2=((0, 1), (1.0, 1.0), (0.0, 0.0), 0, "none", 0.0, 4)
        ),
        _charges_C=np.array([1.0, total_charge - 1.0]),
        _triangles_m=np.zeros((2, 3, 3)),
        _centers_m=np.zeros((2, 3)),
    )
    probe = SimpleNamespace(
        _target_mask=np.array([True, False]),
        _target_points_m=np.zeros((1, 3)),
        _target_charge_weights_C=np.array([1.0]),
        _geometric_area_centroid_m=np.zeros(3),
    )
    return snapshot, probe


class _PathProbe:
    def __init__(self, path: ObjectForcePath) -> None:
        self.path = path
        self.calls = 0

    def vertical_path(self, *_args, **_kwargs) -> ObjectForcePath:
        self.calls += 1
        return self.path


def _install_fake_shell_paths(
    monkeypatch: pytest.MonkeyPatch,
    paths: list[ObjectForcePath],
) -> None:
    @contextmanager
    def fake_finite_probe(_snapshot, _probe, layer):
        yield _PathProbe(paths[layer])

    monkeypatch.setattr(oracle_module, "_finite_probe", fake_finite_probe)
    monkeypatch.setattr(
        oracle_module,
        "_correct_path",
        lambda _snapshot, _probe, path: (path, np.zeros(3)),
    )


def _canonical_periodic_sources(points: np.ndarray, length: float) -> np.ndarray:
    result = np.mod(points, length)
    for axis in (0, 1):
        ordered = np.sort(result[:, axis])
        gaps = np.diff(ordered)
        wrap_gap = ordered[0] + length - ordered[-1]
        if gaps.size and np.max(gaps) >= wrap_gap:
            index = int(np.argmax(gaps))
            origin = ordered[index] + 0.5 * gaps[index]
        else:
            origin = np.mod(ordered[-1] + 0.5 * wrap_gap, length)
        result[result[:, axis] < origin, axis] += length
    return result


def test_identity_rigid_transform_preserves_panel_center_bits_about_nonzero_origin() -> None:
    points = np.mean(_plane_triangles(4), axis=1)
    transformed = RigidTransform.identity().apply(
        points,
        origin=np.array([1.0, 1.0, 0.0]),
    )

    assert transformed.tobytes() == points.tobytes()


def test_panel_wrench_identity_transform_records_primary_self_subtraction(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    triangles = _plane_triangles(8)
    coarse_triangles = _plane_triangles(4)
    centers = np.mean(triangles, axis=1)
    coarse_centers = np.mean(coarse_triangles, axis=1)
    charges = np.zeros(triangles.shape[0])
    for center in coarse_centers:
        matches = np.flatnonzero(
            np.all(np.isclose(centers, center, rtol=0.0, atol=1.0e-14), axis=1)
        )
        assert matches.size == 1
        charges[int(matches[0])] = EPS0 * 4.0 / coarse_triangles.shape[0]
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=triangles.shape[0],
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=np.ones(triangles.shape[0], dtype=np.int64),
        field_source_model="triangle_p0",
        field_kernel_id="triangle_p0_exact_p2m_near",
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        periodic_model="configured",
        library_path=_kernel_lib(),
    ) as snapshot:
        probe = snapshot.object_probe(1)
        transformed_points = RigidTransform.identity().apply(
            probe._target_points_m,
            origin=probe._geometric_area_centroid_m,
        )
        primary_field = probe._primary.eval_e_direct(transformed_points)
        wrench = probe.wrench()

    residual = np.sum(
        probe._target_charge_weights_C[:, None] * primary_field,
        axis=0,
    )
    primary_subtraction = wrench.numerical_metadata["primary_free_subtraction"]
    assert isinstance(primary_subtraction, Mapping)
    np.testing.assert_allclose(
        primary_subtraction["force_N"],
        -residual,
        rtol=0.0,
        atol=0.0,
    )


def test_free_space_moved_single_source_has_no_primary_image_residual(
    tmp_path: Path,
) -> None:
    config = tmp_path / "triangle_p0.toml"
    _write_free_interaction_config(config)
    triangles = _plane_triangles(1)[:1]
    charge = 1.0e-12
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
        charges=np.array([charge]),
        triangles=triangles,
        mesh_ids=np.ones(1, dtype=np.int64),
        field_source_model="triangle_p0",
        field_kernel_id="triangle_p0_exact_p2m_near",
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        periodic_model="configured",
        library_path=_kernel_lib(),
    ) as snapshot:
        wrench = snapshot.object_probe(1).wrench(
            transform=RigidTransform.translation([0.0, 0.0, 0.25]),
        )

    residual = wrench.components["target_periodic_images"]
    force_scale = K_COULOMB * charge**2 / 0.25**2
    force_tolerance = 512.0 * np.finfo(np.float64).eps * force_scale
    assert np.linalg.norm(residual.force_N) <= force_tolerance
    assert residual.potential_energy_J is not None
    energy_scale = K_COULOMB * charge**2 / 0.25
    energy_tolerance = 512.0 * np.finfo(np.float64).eps * energy_scale
    assert abs(residual.potential_energy_J) <= energy_tolerance


def test_periodic_seam_split_and_unwrapped_lattice_representations_match(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    charges = np.array([1.0e-9, -0.4e-9])
    split = np.array([[0.1, 0.7, 0.0], [1.9, 0.7, 0.0]])
    unwrapped = np.array([[2.1, 0.7, 0.0], [1.9, 0.7, 0.0]])

    def result(positions: np.ndarray) -> FortranRunResult:
        value = _result(tmp_path, charges, positions=positions)
        return replace(value, mesh_ids=np.ones(2, dtype=np.int64))

    wrenches = []
    for positions in (split, unwrapped):
        with ObjectInteractionSnapshot.from_result(
            result(positions),
            step=None,
            config_path=config,
            periodic_model="configured",
            library_path=_kernel_lib(),
        ) as snapshot:
            wrenches.append(
                snapshot.object_probe(1).wrench(
                    transform=RigidTransform.translation([0.0, 0.0, 0.25]),
                )
            )

    left, right = wrenches
    np.testing.assert_allclose(left.force_N, right.force_N, rtol=2.0e-13, atol=1.0e-20)
    np.testing.assert_allclose(left.torque_Nm, right.torque_Nm, rtol=2.0e-13, atol=1.0e-20)
    for name in ("target_periodic_images", "total_external"):
        left_component = left.components[name]
        right_component = right.components[name]
        np.testing.assert_allclose(
            left_component.force_N,
            right_component.force_N,
            rtol=2.0e-13,
            atol=1.0e-20,
        )
        np.testing.assert_allclose(
            left_component.torque_Nm,
            right_component.torque_Nm,
            rtol=2.0e-13,
            atol=1.0e-20,
        )
        assert left_component.potential_energy_J == pytest.approx(
            right_component.potential_energy_J,
            rel=2.0e-13,
            abs=1.0e-20,
        )


def test_finite_shell_m1_matches_explicit_replicated_direct_sum(tmp_path: Path) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    charges = np.array([1.0e-9, 2.0e-9])
    result = _result(tmp_path, charges)
    with ObjectInteractionSnapshot.from_result(
        result, step=None, config_path=config, library_path=_kernel_lib()
    ) as snapshot:
        probe = snapshot.object_probe(1)
        shell = finite_shell_wrench(snapshot, probe, None, 1, "symmetric")
        target_points = np.array(probe._target_points_m, copy=True)
        target_weights = np.array(probe._target_charge_weights_C, copy=True)
        torque_origin = np.array(probe.geometric_area_centroid_m, copy=True)

    base = np.array([[0.5, 0.5, 0.0], [1.5, 0.5, 0.3]])
    canonical = _canonical_periodic_sources(base.copy(), 2.0)
    replicated_pos: list[np.ndarray] = []
    replicated_q: list[float] = []
    for ix in range(-1, 2):
        for iy in range(-1, 2):
            shift = np.array([2.0 * ix, 2.0 * iy, 0.0])
            for source in range(2):
                position = canonical[source] + shift
                if source == 0 and np.array_equal(position, base[0]):
                    continue
                replicated_pos.append(position)
                replicated_q.append(charges[source])
    with FieldKernel(
        _triangles_at(np.asarray(replicated_pos)),
        np.asarray(replicated_q),
        options=FieldKernelOptions(leaf_max=64),
        library_path=_kernel_lib(),
    ) as direct:
        expected_field = direct.eval_e_direct(target_points)
    expected_force_samples = target_weights[:, None] * expected_field
    expected_force = np.sum(expected_force_samples, axis=0)
    expected_torque = np.sum(
        np.cross(target_points - torque_origin[None, :], expected_force_samples),
        axis=0,
    )

    np.testing.assert_allclose(
        shell.symmetric.force_N,
        expected_force,
        rtol=2.0e-12,
        atol=1.0e-20,
    )
    np.testing.assert_allclose(
        shell.symmetric.torque_Nm,
        expected_torque,
        rtol=2.0e-12,
        atol=1.0e-20,
    )
    assert shell.selected is shell.symmetric


def test_finite_shell_exposes_analytic_bottom_zero_closure(tmp_path: Path) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    charges = np.array([1.0e-9, 2.0e-9])
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path, charges),
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        probe = snapshot.object_probe(1)
        shell = finite_shell_wrench(snapshot, probe, None, 1, "e_bottom_zero")

    expected_ez = np.sum(charges) / (2.0 * EPS0 * 4.0)
    np.testing.assert_allclose(shell.closure_shift_V_m, [0.0, 0.0, expected_ez])
    np.testing.assert_allclose(
        shell.e_bottom_zero.force_N - shell.symmetric.force_N,
        [0.0, 0.0, charges[0] * expected_ez],
        rtol=2.0e-14,
        atol=1.0e-20,
    )
    assert shell.selected is shell.e_bottom_zero


def test_finite_shell_convergence_requires_two_successive_increments(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path, np.zeros(2)),
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        probe = snapshot.object_probe(1)
        converged = finite_shell_convergence(
            snapshot,
            probe,
            np.array([0.0, 0.1]),
            max_layers=2,
        )
        unresolved = finite_shell_convergence(
            snapshot,
            probe,
            np.array([0.0, 0.1]),
            max_layers=1,
            relative_tolerance=0.0,
            force_floor_N=0.0,
            work_floor_J=0.0,
        )

    assert converged.status == "converged"
    assert converged.selected_image_layers == 2
    assert converged.selected_path is converged.corrected_paths[-1]
    assert unresolved.status == "not_converged"
    assert unresolved.selected_image_layers is None
    assert unresolved.selected_path is None


def test_finite_shell_convergence_record_freezes_path_sequences() -> None:
    path = ObjectForcePath(
        displacement_m=np.array([0.0, 1.0]),
        force_N=np.zeros((2, 3)),
        torque_Nm=np.zeros((2, 3)),
        electrostatic_work_J=np.zeros(2),
    )
    symmetric = [path]
    corrected = [path]

    result = FiniteShellConvergenceResult(
        image_layers=np.array([0]),
        symmetric_paths=symmetric,  # type: ignore[arg-type]
        corrected_paths=corrected,  # type: ignore[arg-type]
        force_increment_error_N=np.empty(0),
        work_increment_error_J=np.empty(0),
        increment_converged=np.empty(0, dtype=bool),
        status="not_converged",
        selected_image_layers=None,
        selected_path=None,
    )
    symmetric.clear()
    corrected.clear()

    assert result.symmetric_paths == (path,)
    assert result.corrected_paths == (path,)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        (
            "increment_converged",
            np.array([False, True]),
            "successive combined",
        ),
        (
            "increment_converged",
            np.array([True, False]),
            "successive combined",
        ),
        (
            "reference_converged",
            np.array([False, False, True]),
            "physical reference",
        ),
        (
            "reference_converged",
            np.array([False, True, False]),
            "physical reference",
        ),
        ("increment_converged", np.array([False, np.nan]), "boolean"),
        (
            "reference_converged",
            np.array([False, np.inf, True]),
            "boolean",
        ),
    ],
    ids=[
        "increment-penultimate-fails",
        "increment-final-fails",
        "reference-penultimate-fails",
        "reference-final-fails",
        "increment-nonboolean",
        "reference-nonboolean",
    ],
)
def test_converged_shell_record_rejects_invalid_gate_contract(
    field: str,
    value: np.ndarray,
    message: str,
) -> None:
    kwargs = _valid_converged_shell_record_kwargs()
    kwargs[field] = value

    with pytest.raises(ValueError, match=message):
        FiniteShellConvergenceResult(**kwargs)


def test_infinite_shell_selection_requires_physical_reference_agreement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference = _constant_path(1.02)
    probe = _PathProbe(reference)
    _install_fake_shell_paths(
        monkeypatch,
        [_constant_path(1.0), _constant_path(1.005), _constant_path(1.01)],
    )

    result = finite_shell_convergence(
        SimpleNamespace(periodic_model="infinite_physical"),
        probe,  # type: ignore[arg-type]
        np.array([0.0, 1.0]),
        max_layers=2,
        relative_tolerance=1.0e-2,
        force_floor_N=0.0,
        work_floor_J=0.0,
    )

    assert probe.calls == 1
    np.testing.assert_array_equal(result.increment_converged, [False, True])
    np.testing.assert_array_equal(result.reference_converged, [False, False, True])
    assert np.all(result.reference_force_error_N > 0.0)
    assert np.all(result.reference_work_error_J > 0.0)
    assert result.status == "not_converged"
    assert result.selected_image_layers is None
    assert result.selected_path is None


def test_configured_shell_selection_uses_layer_scaled_tail_proxy(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    probe = _PathProbe(_constant_path(99.0))
    _install_fake_shell_paths(
        monkeypatch,
        [_constant_path(1.0), _constant_path(1.009), _constant_path(1.018)],
    )

    result = finite_shell_convergence(
        SimpleNamespace(periodic_model="configured"),
        probe,  # type: ignore[arg-type]
        np.array([0.0, 1.0]),
        max_layers=2,
        relative_tolerance=1.0e-2,
        force_floor_N=0.0,
        work_floor_J=0.0,
    )

    assert probe.calls == 0
    np.testing.assert_allclose(result.force_tail_proxy_N, [0.009, 0.018])
    np.testing.assert_allclose(result.work_tail_proxy_J, [0.009, 0.018])
    np.testing.assert_array_equal(result.increment_converged, [True, False])
    assert result.status == "not_converged"
    assert result.selected_path is None


def test_bottom_zero_closure_rechecks_work_potential_status_with_corrected_scale() -> None:
    relative_tolerance = 2.0e-2
    symmetric = ObjectForcePath(
        displacement_m=np.array([0.0, 1.0]),
        force_N=np.array([[0.0, 0.0, 100.0], [0.0, 0.0, 100.0]]),
        torque_Nm=np.zeros((2, 3)),
        electrostatic_work_J=np.array([0.0, 100.0]),
        potential_energy_J=np.array([0.0, -99.0]),
        potential_difference_work_J=np.array([0.0, 99.0]),
        numerical_metadata={
            "relative_tolerance": relative_tolerance,
            "work_absolute_tolerance_J": 0.0,
            "status_reason": "fixed_grid",
        },
        status="converged",
        work_relative_mismatch=1.0e-2,
        work_absolute_mismatch_J=1.0,
    )
    total_charge = -199.0 * EPS0
    snapshot, probe = _bottom_zero_closure_context(total_charge)

    corrected, _ = oracle_module._correct_path(snapshot, probe, symmetric)

    assert symmetric.status == "converged"
    assert corrected.status == "not_converged"
    assert corrected.numerical_metadata["status_reason"] == (
        "work_potential_mismatch"
    )
    assert corrected.work_absolute_mismatch_J == pytest.approx(1.0)
    assert corrected.work_relative_mismatch > relative_tolerance


def test_bottom_zero_closure_does_not_hide_source_refinement_failure() -> None:
    symmetric = ObjectForcePath(
        displacement_m=np.array([0.0, 1.0]),
        force_N=np.zeros((2, 3)),
        torque_Nm=np.zeros((2, 3)),
        electrostatic_work_J=np.zeros(2),
        potential_energy_J=np.zeros(2),
        potential_difference_work_J=np.zeros(2),
        numerical_metadata={
            "relative_tolerance": 1.0e-2,
            "work_absolute_tolerance_J": 1.0e-18,
            "status_reason": "max_refinement_reached",
        },
        status="not_converged",
        work_relative_mismatch=0.0,
        work_absolute_mismatch_J=0.0,
    )
    snapshot, probe = _bottom_zero_closure_context(0.0)

    corrected, _ = oracle_module._correct_path(snapshot, probe, symmetric)

    assert corrected.status == "not_converged"
    assert corrected.numerical_metadata["status_reason"] == (
        "max_refinement_reached"
    )


@pytest.mark.skipif(
    os.environ.get("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS") != "1",
    reason="physical infinite comparison requires opt-in periodic cache generation",
)
def test_non_neutral_local_shell_convergence_is_rejected_by_physical_reference(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    charges = np.array([1.0e-9, 2.0e-9])
    displacement = np.linspace(0.0, 0.1, 9)
    relative_tolerance = 1.0e-2
    force_floor = 1.0e-16
    work_floor = 1.0e-18
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path, charges),
        step=None,
        config_path=config,
        periodic_model="infinite_physical",
        cache_dir=Path(".beach_cache/periodic2-task7-oracle"),
        library_path=_kernel_lib(),
    ) as snapshot:
        probe = snapshot.object_probe(1)
        infinite = probe.vertical_path(
            displacement,
            adaptive=False,
            relative_tolerance=relative_tolerance,
            force_absolute_tolerance_N=force_floor,
            work_absolute_tolerance_J=work_floor,
        )
        shells = finite_shell_convergence(
            snapshot,
            probe,
            displacement,
            max_layers=2,
            relative_tolerance=relative_tolerance,
            force_floor_N=force_floor,
            work_floor_J=work_floor,
        )

    assert infinite.status == "converged"
    assert np.any(shells.force_increment_error_N > force_floor)
    assert np.any(shells.work_increment_error_J > work_floor)
    assert shells.reference_force_error_N is not None
    assert shells.reference_work_error_J is not None
    assert shells.reference_converged is not None
    force_error = float(
        np.max(
            np.linalg.norm(shells.corrected_paths[-1].force_N - infinite.force_N, axis=1)
        )
    )
    force_scale = max(
        float(np.max(np.linalg.norm(shells.corrected_paths[-1].force_N, axis=1))),
        float(np.max(np.linalg.norm(infinite.force_N, axis=1))),
    )
    work_error = float(
        np.max(
            np.abs(
                shells.corrected_paths[-1].electrostatic_work_J
                - infinite.electrostatic_work_J
            )
        )
    )
    work_scale = max(
        float(np.max(np.abs(shells.corrected_paths[-1].electrostatic_work_J))),
        float(np.max(np.abs(infinite.electrostatic_work_J))),
    )
    assert shells.reference_force_error_N[-1] == pytest.approx(force_error)
    assert shells.reference_work_error_J[-1] == pytest.approx(work_error)
    assert force_error > force_floor + relative_tolerance * force_scale
    assert work_error > work_floor + relative_tolerance * work_scale
    assert shells.reference_converged[-1] is np.False_
    assert shells.status == "not_converged"
    assert shells.selected_image_layers is None
    assert shells.selected_path is None


@pytest.mark.skipif(
    os.environ.get("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS") != "1",
    reason="physical infinite comparison requires opt-in periodic cache generation",
)
def test_non_neutral_near_field_shell_converges_to_physical_infinite_reference(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    displacement = np.linspace(0.0, 0.1, 9)
    with ObjectInteractionSnapshot.from_result(
        _result(
            tmp_path,
            np.array([1.0e-9, 2.0e-9]),
            positions=np.array([[0.5, 0.5, 0.0], [0.51, 0.5, -0.01]]),
        ),
        step=None,
        config_path=config,
        periodic_model="infinite_physical",
        cache_dir=Path(".beach_cache/periodic2-task7-oracle"),
        library_path=_kernel_lib(),
    ) as snapshot:
        shells = finite_shell_convergence(
            snapshot,
            snapshot.object_probe(1),
            displacement,
            max_layers=4,
            relative_tolerance=1.0e-2,
            force_floor_N=1.0e-16,
            work_floor_J=1.0e-18,
        )

    assert shells.status == "converged"
    assert shells.selected_image_layers is not None
    assert shells.selected_image_layers >= 2
    assert shells.selected_path is not None
    assert shells.reference_converged is not None
    assert shells.reference_converged[-1]
    np.testing.assert_array_equal(shells.increment_converged[-2:], [True, True])


@pytest.mark.skipif(
    os.environ.get("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS") != "1"
    or not os.environ.get("BEACH_FIELD_KERNEL_LIB"),
    reason=(
        "native validation-tool oracles require opt-in cache generation and "
        "BEACH_FIELD_KERNEL_LIB"
    ),
)
def test_validation_tool_native_plane_oracles_match_receipt_contract(
    tmp_path: Path,
) -> None:
    tool = _load_validation_tool()
    library_value = os.environ.get("BEACH_FIELD_KERNEL_LIB")
    assert library_value is not None
    provided_library = Path(library_value).expanduser().resolve()
    assert provided_library.is_file(), (
        f"native field-kernel library does not exist: {provided_library}"
    )
    library = tmp_path / "lib" / provided_library.name
    library.parent.mkdir(parents=True)
    shutil.copy2(provided_library, library)
    library.chmod(0o555)
    build_origin = tool._library_build_info(library)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "validation_root": str(tmp_path.resolve()),
                "analysis_library": {
                    "provided_path": str(provided_library),
                    "staged_path": str(library),
                    "sha256": tool._sha256(library),
                },
                "build_origin": build_origin,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    oracle_dir = tmp_path / "provenance" / "oracles"
    oracle_dir.mkdir(parents=True)
    triangle_config = oracle_dir / "periodic_plane.toml"
    tool._oracle_panel_config(triangle_config)
    cache_dir = tmp_path / "cache" / "oracles"
    cache_dir.mkdir(parents=True)

    kernel_oracles = {
        "triangle_p0": tool._run_periodic_plane_kernel_oracle(
            root=tmp_path,
            config_path=triangle_config,
            cache_dir=cache_dir,
            library_path=library,
        ),
    }
    for label, result in kernel_oracles.items():
        tool._verify_periodic_oracle_metrics(result, label=label)
        evaluations = result["cache_evaluations"]
        assert tuple(row["label"] for row in evaluations) == (
            tool.ORACLE_CACHE_EVALUATION_LABELS
        )
        identities = result["cache_identities"]
        assert result["cache_diagnostics"] == identities["4"]
        assert identities["4"]["fingerprint"] != identities["8"]["fingerprint"]
        assert identities["4"]["path"] != identities["8"]["path"]
        assert all(row["hit"] is True for row in evaluations)
        assert all(row["build_count"] == 0 for row in evaluations)
        cosine = result["neutral_cosine_plane"]
        assert cosine["sample_abs_z_m"] == [0.25, 0.5]
        assert cosine["expected_decay_ratio"] == pytest.approx(
            np.exp(-np.pi * 0.25)
        )
        assert cosine["decay_ratio_relative_tolerance"] == 0.18
        for row in cosine["errors"]:
            assert row["field_decay_ratio_relative_error"] == pytest.approx(
                abs(
                    row["field_decay_ratio"]
                    - cosine["expected_decay_ratio"]
                )
                / cosine["expected_decay_ratio"]
            )
            assert row["potential_decay_ratio_relative_error"] == pytest.approx(
                abs(
                    row["potential_decay_ratio"]
                    - cosine["expected_decay_ratio"]
                )
                / cosine["expected_decay_ratio"]
            )
        fine_cosine = cosine["errors"][1]
        assert fine_cosine["field_decay_ratio_relative_error"] <= 0.18
        assert fine_cosine["potential_decay_ratio_relative_error"] <= 0.18
        for group, group_labels in tool.ORACLE_CACHE_EVALUATION_GROUPS.items():
            canonical = identities[group]
            grouped = [row for row in evaluations if row["label"] in group_labels]
            assert {row["label"] for row in grouped} == set(group_labels)
            assert all(
                (row["fingerprint"], row["path"], row["sha256"])
                == (
                    canonical["fingerprint"],
                    canonical["path"],
                    canonical["sha256"],
                )
                for row in grouped
            )

    receipt = {
        "receipt_schema_version": 1,
        "oracle_schema_version": 3,
        "status": "qualified",
        "manifest_sha256": tool._sha256(manifest),
        "library": str(library),
        "library_sha256": tool._sha256(library),
        "config": str(triangle_config),
        "config_sha256": tool._sha256(triangle_config),
        "kernel_configs": {
            "triangle_p0": {
                "path": str(triangle_config),
                "sha256": tool._sha256(triangle_config),
            },
        },
        "cache_dir": str(cache_dir),
        "cache_files": tool._output_inventory(cache_dir),
        "execution_job_id": "pytest-native",
        "kernel_oracles": kernel_oracles,
        "library_build_origin": build_origin,
        "verified_at": "pytest-native",
    }
    receipt_path = oracle_dir / "periodic_plane.json"
    tool._publish_execution_receipt(receipt_path, receipt)
    verified = tool._verify_periodic_oracle_receipt(
        tmp_path,
        library,
        expected_job_id="pytest-native",
    )

    assert verified["status"] == "qualified"
    assert set(verified["kernel_oracles"]) == {"triangle_p0"}


@pytest.mark.skipif(
    os.environ.get("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS") != "1",
    reason="analytic periodic panel plane requires opt-in cache generation",
)
def test_uniform_triangle_plane_matches_bottom_zero_jump_pv_and_object_wrench(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_panel_periodic_config(config)
    area_xy = 4.0
    surface_charge_density = EPS0
    triangles = _plane_triangles(4)
    triangle_area = area_xy / triangles.shape[0]
    charges = np.full(triangles.shape[0], surface_charge_density * triangle_area)
    total_charge = float(np.sum(charges))
    cache_dir = Path(".beach_cache/periodic2-task7-oracle")
    points = np.array(
        [
            [0.37, 0.61, -0.25],
            [0.37, 0.61, 0.0],
            [0.37, 0.61, 0.25],
        ]
    )
    with ObjectInteractionSnapshot.from_result(
        _panel_result(tmp_path, triangles, charges),
        step=None,
        config_path=config,
        periodic_model="infinite_physical",
        cache_dir=cache_dir,
        library_path=_kernel_lib(),
    ) as snapshot:
        field = _physical_field(snapshot, points)
        order3 = snapshot.object_probe(1, quadrature_order=3).wrench()
        order7 = snapshot.object_probe(1, quadrature_order=7).wrench()

    fmm_panel_relative_tolerance = 1.2e-1
    expected_ez = np.array([0.0, 0.5, 1.0])
    assert field[0, 2] == pytest.approx(
        expected_ez[0],
        abs=fmm_panel_relative_tolerance,
    )
    np.testing.assert_allclose(
        field[1:, 2],
        expected_ez[1:],
        rtol=fmm_panel_relative_tolerance,
        atol=0.0,
    )
    np.testing.assert_allclose(
        field[:, :2],
        0.0,
        atol=fmm_panel_relative_tolerance,
    )
    expected_force = total_charge * surface_charge_density / (2.0 * EPS0)
    for wrench in (order3, order7):
        assert wrench.force_N[2] == pytest.approx(
            expected_force,
            rel=fmm_panel_relative_tolerance,
        )
        assert np.linalg.norm(wrench.force_N[:2]) <= (
            fmm_panel_relative_tolerance * expected_force
        )
        assert np.linalg.norm(wrench.torque_Nm) <= (
            fmm_panel_relative_tolerance * expected_force * 2.0
        )
    assert np.linalg.norm(order3.force_N - order7.force_N) <= (
        1.0e-2 * abs(expected_force)
    )


@pytest.mark.skipif(
    os.environ.get("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS") != "1",
    reason="analytic periodic cosine plane requires opt-in cache generation",
)
def test_neutral_cosine_triangle_plane_converges_to_nonzero_mode_solution(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_panel_periodic_config(config)
    length = 2.0
    wave_number = 2.0 * np.pi / length
    sigma_amplitude = 2.0 * EPS0
    points = np.array(
        [
            [0.20, 0.73, 0.25],
            [0.55, 0.73, 0.25],
            [1.10, 0.73, 0.25],
            [1.65, 0.73, 0.25],
        ]
    )
    decay = np.exp(-wave_number * points[:, 2])
    phase = wave_number * points[:, 0]
    expected_field = np.zeros_like(points)
    expected_field[:, 0] = (
        sigma_amplitude / (2.0 * EPS0) * np.sin(phase) * decay
    )
    expected_field[:, 2] = (
        sigma_amplitude / (2.0 * EPS0) * np.cos(phase) * decay
    )
    expected_potential = (
        sigma_amplitude
        / (2.0 * EPS0 * wave_number)
        * np.cos(phase)
        * decay
    )
    errors: list[tuple[float, float]] = []
    for cells_per_axis in (4, 8):
        triangles = _plane_triangles(cells_per_axis, length=length)
        centers = np.mean(triangles, axis=1)
        triangle_area = length**2 / triangles.shape[0]
        charges = (
            sigma_amplitude
            * np.cos(wave_number * centers[:, 0])
            * triangle_area
        )
        with ObjectInteractionSnapshot.from_result(
            _panel_result(tmp_path, triangles, charges),
            step=None,
            config_path=config,
            periodic_model="infinite_physical",
            cache_dir=Path(".beach_cache/periodic2-task7-oracle"),
            library_path=_kernel_lib(),
        ) as snapshot:
            field = snapshot._periodic.eval_e(points)
            potential = snapshot._periodic.eval_phi(points)
        field_error = float(
            np.linalg.norm(field - expected_field) / np.linalg.norm(expected_field)
        )
        potential_error = float(
            np.linalg.norm(potential - expected_potential)
            / np.linalg.norm(expected_potential)
        )
        errors.append((field_error, potential_error))

    assert errors[1][0] < errors[0][0]
    assert errors[1][1] < errors[0][1]
    assert errors[1][0] <= 8.0e-2
    assert errors[1][1] <= 8.0e-2


def test_finite_shell_closure_path_force_and_potential_work_agree(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    charges = np.array([1.0e-9, 2.0e-9])
    displacement = np.linspace(0.0, 0.1, 9)
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path, charges),
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        result = finite_shell_convergence(
            snapshot,
            snapshot.object_probe(1),
            displacement,
            max_layers=0,
        )

    symmetric = result.symmetric_paths[0]
    corrected = result.corrected_paths[0]
    expected_ez = np.sum(charges) / (2.0 * EPS0 * 4.0)
    expected_work = charges[0] * expected_ez * displacement
    np.testing.assert_allclose(
        corrected.electrostatic_work_J - symmetric.electrostatic_work_J,
        expected_work,
        rtol=2.0e-14,
        atol=1.0e-20,
    )
    np.testing.assert_allclose(
        corrected.potential_difference_work_J
        - symmetric.potential_difference_work_J,
        expected_work,
        rtol=2.0e-14,
        atol=1.0e-20,
    )
