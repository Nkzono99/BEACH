from pathlib import Path

import numpy as np
import pytest

import beach.fortran_results.kernel as kernel_module
import beach.fortran_results.object_interaction as interaction_module
from beach import (
    calc_object_forces_kernel,
    FieldKernel,
    FieldKernelDiagnostics,
    FieldKernelError,
    FieldKernelOptions,
    FortranRunResult,
    ObjectInteractionSnapshot,
    RigidTransform,
)
from beach.fortran_results.constants import K_COULOMB


def _kernel_lib() -> Path:
    path = Path("build/libbeach_field_kernel.so")
    if not path.exists():
        pytest.skip("field kernel shared library is not built; run `make build-kernel`")
    return path


def _triangles_at(positions: np.ndarray) -> np.ndarray:
    positions = np.asarray(positions, dtype=float)
    offsets = np.array(
        [
            [-0.02, -0.02, 0.0],
            [0.04, -0.02, 0.0],
            [-0.02, 0.04, 0.0],
        ]
    )
    return positions[:, None, :] + offsets[None, :, :]


def _result(
    directory: Path,
    positions: np.ndarray,
    charges: np.ndarray,
    mesh_ids: np.ndarray,
    *,
    source_model: str = "point",
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
        triangles=_triangles_at(positions),
        mesh_ids=np.asarray(mesh_ids, dtype=np.int64),
        field_source_model=source_model,
    )


def _write_free_config(path: Path) -> None:
    path.write_text(
        """
[sim]
field_bc_mode = "free"
softening = 0.0
box_min = [-2.0, -2.0, -2.0]
box_max = [2.0, 2.0, 2.0]
e0 = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )


def _write_periodic_config(
    path: Path,
    *,
    far_correction: str,
    outer_model: str | None = None,
) -> None:
    outer = "" if outer_model is None else f'\n[outer_plasma]\nmodel = "{outer_model}"\n'
    path.write_text(
        (
            """
[sim]
field_bc_mode = "periodic2"
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
box_min = [0.0, 0.0, 0.0]
box_max = [4.0, 4.0, 4.0]
field_periodic_image_layers = 1
field_periodic_ewald_layers = 4
field_periodic_generation_tolerance = 1.0e-8
e0 = [7.0, 8.0, 9.0]
""".strip()
            + f'\nfield_periodic_far_correction = "{far_correction}"\n'
            + outer
        ),
        encoding="utf-8",
    )


def test_free_single_object_excludes_only_its_primary_source(tmp_path: Path) -> None:
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    result = _result(
        tmp_path,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([1.0e-9]),
        np.array([1]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        wrench = snapshot.object_probe(1).wrench()

    np.testing.assert_allclose(wrench.force_N, 0.0, atol=1.0e-20)
    np.testing.assert_allclose(wrench.torque_Nm, 0.0, atol=1.0e-20)
    assert tuple(wrench.components) == (
        "other_objects_all_images",
        "target_periodic_images",
        "external_uniform",
        "total_external",
    )


def test_free_two_object_force_is_coulomb_action_reaction(tmp_path: Path) -> None:
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    charges = np.array([1.0e-9, 2.0e-9])
    result = _result(tmp_path, positions, charges, np.array([1, 2]))

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        first = snapshot.object_probe(1).wrench(torque_origin="origin")
        second = snapshot.object_probe(2).wrench(torque_origin="origin")
        shifted_origin = np.array([0.25, 0.5, 0.75])
        first_shifted = snapshot.object_probe(1).wrench(
            torque_origin=shifted_origin
        )

    delta = positions[0] - positions[1]
    expected = K_COULOMB * charges[0] * charges[1] * delta / np.linalg.norm(delta) ** 3
    np.testing.assert_allclose(first.force_N, expected, rtol=2.0e-14, atol=1.0e-20)
    np.testing.assert_allclose(second.force_N, -expected, rtol=2.0e-14, atol=1.0e-20)
    np.testing.assert_allclose(first.torque_Nm, np.cross(positions[0], expected), atol=1.0e-20)
    np.testing.assert_allclose(
        first_shifted.torque_Nm,
        first.torque_Nm - np.cross(shifted_origin, first.force_N),
        atol=1.0e-20,
    )


class _FakeFieldKernel:
    instances: list["_FakeFieldKernel"] = []
    fail_charge_sum: float | None = None

    def __init__(
        self,
        source_positions: np.ndarray,
        source_charges: np.ndarray,
        *,
        options,
        library_path=None,
        source_triangles=None,
    ) -> None:
        del library_path, source_triangles
        self.source_positions = np.array(source_positions, copy=True)
        self.current = np.array(source_charges, copy=True)
        self.options = options
        self.closed = False
        self.eval_calls = 0
        self.update_history = [self.current.copy()]
        self.is_periodic = options.periodic2 is not None
        self.fail_update_sum: float | None = None
        self.fail_close_once = False
        type(self).instances.append(self)

    def update_charges(self, charges: np.ndarray) -> None:
        if self.fail_update_sum is not None and float(np.sum(charges)) == self.fail_update_sum:
            self.fail_update_sum = None
            raise RuntimeError("injected update failure")
        self.current = np.array(charges, copy=True)
        self.update_history.append(self.current.copy())

    def eval_e(self, target_points: np.ndarray) -> np.ndarray:
        self.eval_calls += 1
        total = float(np.sum(self.current))
        if self.fail_charge_sum is not None and total == self.fail_charge_sum:
            type(self).fail_charge_sum = None
            raise RuntimeError("injected periodic evaluation failure")
        return np.tile([total, 2.0 * total, 3.0 * total], (len(target_points), 1))

    def eval_phi(self, target_points: np.ndarray) -> np.ndarray:
        return np.full(len(target_points), float(np.sum(self.current)))

    def eval_e_direct(self, target_points: np.ndarray) -> np.ndarray:
        total = float(np.sum(self.current))
        return np.tile([4.0 * total, 0.0, 0.0], (len(target_points), 1))

    def eval_phi_direct(self, target_points: np.ndarray) -> np.ndarray:
        return np.full(len(target_points), 4.0 * float(np.sum(self.current)))

    def diagnostics(self) -> FieldKernelDiagnostics:
        return FieldKernelDiagnostics(False, 1, "fake", Path("fake.bin"))

    def close(self) -> None:
        if self.fail_close_once:
            self.fail_close_once = False
            raise RuntimeError("injected close failure")
        self.closed = True


class _FakeZeroMode:
    instances: list["_FakeZeroMode"] = []
    forbidden = False
    plus_delta = 0.0

    def __init__(
        self,
        source_heights_m: np.ndarray,
        source_charges_C: np.ndarray,
        area_xy_m2: float,
        **_kwargs,
    ) -> None:
        if type(self).forbidden:
            raise AssertionError("finite configured model must not construct a zero mode")
        self.source_heights_m = np.array(source_heights_m, copy=True)
        self.current = np.array(source_charges_C, copy=True)
        self.area_xy_m2 = area_xy_m2
        self.closed = False
        self.update_history = [self.current.copy()]
        self.trace_history: list[str] = []
        type(self).instances.append(self)

    def update_charges(self, charges: np.ndarray) -> None:
        self.current = np.array(charges, copy=True)
        self.update_history.append(self.current.copy())

    def eval(self, z_m: np.ndarray, trace: str = "principal_value"):
        assert trace in {"principal_value", "plus"}
        self.trace_history.append(trace)
        total = float(np.sum(self.current))
        trace_delta = type(self).plus_delta * total if trace == "plus" else 0.0
        return (
            np.full(len(z_m), 5.0 * total),
            np.full(len(z_m), 6.0 * total + trace_delta),
        )

    def close(self) -> None:
        self.closed = True


@pytest.fixture
def fake_native(monkeypatch: pytest.MonkeyPatch):
    _FakeFieldKernel.instances.clear()
    _FakeFieldKernel.fail_charge_sum = None
    _FakeZeroMode.instances.clear()
    _FakeZeroMode.forbidden = False
    _FakeZeroMode.plus_delta = 0.0
    monkeypatch.setattr(interaction_module, "FieldKernel", _FakeFieldKernel)
    monkeypatch.setattr(interaction_module, "PeriodicZeroMode", _FakeZeroMode)
    return _FakeFieldKernel, _FakeZeroMode


def test_cached_configured_and_override_compose_p_plus_z_minus_primary(
    tmp_path: Path,
    fake_native,
) -> None:
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="cached_kneq0")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )

    wrenches = []
    for model in ("configured", "infinite_physical"):
        with ObjectInteractionSnapshot.from_result(
            result,
            step=None,
            config_path=config,
            periodic_model=model,
        ) as snapshot:
            wrenches.append(snapshot.object_probe(1).wrench(components=True))

    for wrench in wrenches:
        np.testing.assert_allclose(
            wrench.components["other_objects_all_images"].force_N,
            [2.0, 4.0, 18.0],
        )
        np.testing.assert_allclose(
            wrench.components["target_periodic_images"].force_N,
            [-3.0, 2.0, 9.0],
        )
        np.testing.assert_allclose(
            wrench.components["external_uniform"].force_N,
            [7.0, 8.0, 9.0],
        )
        np.testing.assert_allclose(wrench.force_N, [6.0, 14.0, 36.0])
        assert wrench.components["other_objects_all_images"].potential_energy_J == pytest.approx(12.0)
        assert wrench.components["target_periodic_images"].potential_energy_J == pytest.approx(2.0)
        assert wrench.components["external_uniform"].potential_energy_J == pytest.approx(-24.0)
        assert wrench.components["total_external"].potential_energy_J == pytest.approx(-10.0)
        np.testing.assert_allclose(
            sum(
                wrench.components[name].force_N
                for name in (
                    "other_objects_all_images",
                    "target_periodic_images",
                    "external_uniform",
                )
            ),
            wrench.force_N,
        )
        assert set(wrench.numerical_metadata) >= {
            "periodic_kneq0",
            "physical_k0",
            "primary_free_subtraction",
        }
        np.testing.assert_allclose(
            wrench.numerical_metadata["periodic_kneq0"]["force_N"],  # type: ignore[index]
            [3.0, 6.0, 9.0],
        )
        np.testing.assert_allclose(
            wrench.numerical_metadata["physical_k0"]["force_N"],  # type: ignore[index]
            [0.0, 0.0, 18.0],
        )
        np.testing.assert_allclose(
            wrench.numerical_metadata["primary_free_subtraction"][  # type: ignore[index]
                "force_N"
            ],
            [-4.0, 0.0, 0.0],
        )
    np.testing.assert_allclose(wrenches[0].force_N, wrenches[1].force_N)
    periodic_instances = [
        instance for instance in _FakeFieldKernel.instances if instance.is_periodic
    ]
    assert all(instance.options.external_e0 == (0.0, 0.0, 0.0) for instance in periodic_instances)


def test_cached_mechanical_trace_restores_native_plus_subtraction(
    tmp_path: Path,
    fake_native,
) -> None:
    _, fake_zero = fake_native
    fake_zero.plus_delta = 2.0
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="cached_kneq0")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    ) as snapshot:
        wrench = snapshot.object_probe(1).wrench()

    np.testing.assert_allclose(
        wrench.components["other_objects_all_images"].force_N,
        [2.0, 4.0, 22.0],
    )
    np.testing.assert_allclose(
        wrench.components["target_periodic_images"].force_N,
        [-3.0, 2.0, 11.0],
    )
    np.testing.assert_allclose(wrench.force_N, [6.0, 14.0, 42.0])
    np.testing.assert_allclose(
        wrench.numerical_metadata["periodic_kneq0"]["force_N"],  # type: ignore[index]
        [3.0, 6.0, 15.0],
    )
    np.testing.assert_allclose(
        wrench.numerical_metadata["physical_k0"]["force_N"],  # type: ignore[index]
        [0.0, 0.0, 18.0],
    )
    np.testing.assert_allclose(
        wrench.numerical_metadata["cached_kneq0_trace_correction"][  # type: ignore[index]
            "force_N"
        ],
        [0.0, 0.0, 6.0],
    )
    assert set(_FakeZeroMode.instances[0].trace_history) == {
        "plus",
        "principal_value",
    }


def test_configured_finite_shell_does_not_construct_or_add_zero_mode(
    tmp_path: Path,
    fake_native,
) -> None:
    _, fake_zero = fake_native
    fake_zero.forbidden = True
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="none")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    ) as snapshot:
        wrench = snapshot.object_probe(1).wrench()

    np.testing.assert_allclose(
        wrench.components["other_objects_all_images"].force_N,
        [2.0, 4.0, 6.0],
    )
    np.testing.assert_allclose(
        wrench.components["target_periodic_images"].force_N,
        [-3.0, 2.0, 3.0],
    )


def test_uniform_potential_uses_simulator_global_origin_gauge(
    tmp_path: Path,
    fake_native,
) -> None:
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="none")
    text = config.read_text(encoding="utf-8")
    text = text.replace(
        "box_min = [0.0, 0.0, 0.0]\nbox_max = [4.0, 4.0, 4.0]",
        "box_min = [-1.0, -2.0, -3.0]\nbox_max = [3.0, 2.0, 1.0]",
    )
    config.write_text(text, encoding="utf-8")
    result = _result(
        tmp_path,
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    ) as snapshot:
        uniform = snapshot.object_probe(1).wrench().components["external_uniform"]

    assert uniform.potential_energy_J == pytest.approx(0.0)


def test_neutral_target_has_zero_vertical_whole_object_force_in_symmetric_shell(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="none")
    result = _result(
        tmp_path,
        np.array([[1.5, 2.0, 1.0], [2.5, 2.0, 1.0]]),
        np.array([1.0e-9, -1.0e-9]),
        np.array([1, 1]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        wrench = snapshot.object_probe(1).wrench()

    assert wrench.force_N[2] == pytest.approx(0.0, abs=1.0e-13)


def test_failed_component_evaluation_restores_full_source_state(
    tmp_path: Path,
    fake_native,
) -> None:
    fake_field, _ = fake_native
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="cached_kneq0")
    charges = np.array([1.0, 2.0])
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        charges,
        np.array([1, 2]),
    )
    snapshot = ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    )
    fake_field.fail_charge_sum = 2.0

    with pytest.raises(RuntimeError, match="injected"):
        snapshot.object_probe(1).wrench()

    periodic = next(instance for instance in fake_field.instances if instance.is_periodic)
    zero = _FakeZeroMode.instances[0]
    np.testing.assert_array_equal(periodic.current, charges)
    np.testing.assert_array_equal(zero.current, charges)
    clean = snapshot.object_probe(1).wrench()
    np.testing.assert_allclose(clean.force_N, [6.0, 14.0, 36.0])
    snapshot.close()


def test_snapshot_defensively_copies_saved_source_arrays(
    tmp_path: Path,
    fake_native,
) -> None:
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="none")
    charges = np.array([1.0, 2.0])
    expected_charges = charges.copy()
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        charges,
        np.array([1, 2]),
    )
    snapshot = ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    )
    result.charges[:] = 99.0
    result.triangles[:] = 99.0  # type: ignore[index]

    np.testing.assert_array_equal(snapshot.source_charges_C, expected_charges)
    np.testing.assert_allclose(snapshot.source_positions_m[0], [1.0, 1.0, 1.0])
    with pytest.raises(ValueError):
        snapshot.source_charges_C[0] = 0.0
    snapshot.close()


def test_active_outer_model_rejected_before_native_construction(
    tmp_path: Path,
    fake_native,
) -> None:
    fake_field, _ = fake_native
    config = tmp_path / "beach.toml"
    _write_periodic_config(
        config,
        far_correction="cached_kneq0",
        outer_model="unified_linear_response",
    )
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0]]),
        np.array([1.0]),
        np.array([1]),
    )

    with pytest.raises(ValueError, match="outer_plasma"):
        ObjectInteractionSnapshot.from_result(result, config_path=config)

    assert fake_field.instances == []


def test_missing_full_config_fails_closed_before_native_construction(
    tmp_path: Path,
    fake_native,
) -> None:
    fake_field, _ = fake_native
    result = _result(
        tmp_path,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([1.0]),
        np.array([1]),
    )

    with pytest.raises(ValueError, match="beach.toml"):
        ObjectInteractionSnapshot.from_result(result, step=None)

    assert fake_field.instances == []


def test_incomplete_box_config_fails_before_native_construction(
    tmp_path: Path,
    fake_native,
) -> None:
    fake_field, _ = fake_native
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    config.write_text(
        config.read_text(encoding="utf-8").replace(
            "box_max = [2.0, 2.0, 2.0]\n", ""
        ),
        encoding="utf-8",
    )
    result = _result(
        tmp_path,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([1.0]),
        np.array([1]),
    )

    with pytest.raises(ValueError, match="box_min.*box_max"):
        ObjectInteractionSnapshot.from_result(result, step=None, config_path=config)

    assert fake_field.instances == []


def test_probe_prevalidates_transformed_points_before_native_calls(
    tmp_path: Path,
    fake_native,
) -> None:
    fake_field, _ = fake_native
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="none")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0]]),
        np.array([1.0]),
        np.array([1]),
    )
    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    ) as snapshot:
        probe = snapshot.object_probe(1)
        periodic = next(instance for instance in fake_field.instances if instance.is_periodic)
        before = periodic.eval_calls
        with pytest.raises(ValueError, match="box"):
            probe.wrench(transform=RigidTransform.translation([0.0, 0.0, 4.0]))
        assert periodic.eval_calls == before


def test_restore_failure_poisons_snapshot_and_blocks_reuse(
    tmp_path: Path,
    fake_native,
) -> None:
    fake_field, _ = fake_native
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="cached_kneq0")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )
    snapshot = ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    )
    probe = snapshot.object_probe(1)
    periodic = next(instance for instance in fake_field.instances if instance.is_periodic)
    periodic.fail_update_sum = 3.0

    with pytest.raises(FieldKernelError, match="restore"):
        probe.wrench()
    with pytest.raises(FieldKernelError, match="unavailable"):
        probe.wrench()

    snapshot.close()


def test_failed_close_can_be_retried_without_reusing_partial_snapshot(
    tmp_path: Path,
    fake_native,
) -> None:
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="none")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0]]),
        np.array([1.0]),
        np.array([1]),
    )
    snapshot = ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    )
    probe = snapshot.object_probe(1)
    probe._primary.fail_close_once = True

    with pytest.raises(FieldKernelError, match="close"):
        snapshot.close()
    with pytest.raises(FieldKernelError, match="unavailable"):
        probe.wrench()
    snapshot.close()

    assert probe._primary.closed
    with pytest.raises(FieldKernelError, match="closed"):
        probe.wrench()


def test_numpy_wrench_aggregation_matches_native_force_on_charges() -> None:
    source_positions = np.array([[0.0, 0.0, 0.0], [0.5, -0.2, 0.3]])
    source_charges = np.array([1.0e-9, -0.5e-9])
    target_positions = np.array([[0.3, 0.4, 0.5], [-0.2, 0.1, 0.7]])
    target_charges = np.array([2.0e-9, -3.0e-9])
    origin = np.array([0.1, -0.1, 0.2])

    with FieldKernel(
        source_positions,
        source_charges,
        options=FieldKernelOptions(softening=0.05),
        library_path=_kernel_lib(),
    ) as kernel:
        field = kernel.eval_e(target_positions)
        potential = kernel.eval_phi(target_positions)
        native_force, native_torque = kernel.force_on_charges(
            target_positions,
            target_charges,
            origin=origin,
        )
    numpy_wrench = interaction_module._aggregate_wrench(
        target_positions,
        target_charges,
        field,
        potential,
        origin,
    )

    np.testing.assert_allclose(numpy_wrench.force_N, native_force, rtol=1.0e-14)
    np.testing.assert_allclose(numpy_wrench.torque_Nm, native_torque, rtol=1.0e-14)


class _LegacyCaptureKernel:
    instances: list["_LegacyCaptureKernel"] = []

    def __init__(self, _positions, charges, **_kwargs) -> None:
        self.updates = [np.array(charges, copy=True)]
        type(self).instances.append(self)

    def update_charges(self, charges: np.ndarray) -> None:
        self.updates.append(np.array(charges, copy=True))

    def force_on_charges(self, _positions, _charges, *, origin):
        del origin
        return np.zeros(3), np.zeros(3)

    def close(self) -> None:
        pass

    def __enter__(self):
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()


def test_legacy_object_force_helper_still_zeros_the_target_periodic_lattice(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _LegacyCaptureKernel.instances.clear()
    monkeypatch.setattr(kernel_module, "FieldKernel", _LegacyCaptureKernel)
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="none")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )

    calc_object_forces_kernel(
        result,
        step=None,
        target_mesh_ids=1,
        config_path=config,
    )

    legacy = _LegacyCaptureKernel.instances[0]
    np.testing.assert_array_equal(legacy.updates[-1], [0.0, 2.0])


def test_repeated_snapshot_and_probe_contexts_close_every_native_handle(
    tmp_path: Path,
    fake_native,
) -> None:
    fake_field, fake_zero = fake_native
    config = tmp_path / "beach.toml"
    _write_periodic_config(config, far_correction="cached_kneq0")
    result = _result(
        tmp_path,
        np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 2.0]]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )

    for _ in range(5):
        with ObjectInteractionSnapshot.from_result(
            result,
            step=None,
            config_path=config,
        ) as snapshot:
            with snapshot.object_probe(1) as probe:
                probe.wrench()

    assert all(instance.closed for instance in fake_field.instances)
    assert all(instance.closed for instance in fake_zero.instances)


def test_triangle_source_is_never_silently_downgraded_to_centroids(
    tmp_path: Path,
    fake_native,
) -> None:
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    result = _result(
        tmp_path,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([1.0]),
        np.array([1]),
        source_model="triangle_p0",
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    ) as snapshot:
        automatic = snapshot.object_probe(1).wrench()
        compatibility = snapshot.object_probe(
            1,
            target_integration="centroid_compatibility",
        ).wrench()

    assert automatic.numerical_metadata["target_integration"] == "gauss_duffy"
    assert automatic.numerical_metadata["quadrature_order"] == 7
    assert compatibility.numerical_metadata["target_integration"] == (
        "centroid_compatibility"
    )


def test_unknown_periodic_model_and_self_policy_fail_closed(tmp_path: Path) -> None:
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    result = _result(
        tmp_path,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([1.0]),
        np.array([1]),
    )

    with pytest.raises(ValueError, match="periodic_model"):
        ObjectInteractionSnapshot.from_result(
            result,
            config_path=config,
            periodic_model="finite",
        )
    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        with pytest.raises(ValueError, match="self_policy"):
            snapshot.object_probe(1, self_policy="exclude_target_lattice")
