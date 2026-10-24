from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from beach import (
    AdhesionProfile,
    ObjectForcePath,
    ObjectWrench,
    RigidTransform,
    WrenchComponent,
)


def _path(force_z: np.ndarray, displacement: np.ndarray | None = None) -> ObjectForcePath:
    h = np.asarray(displacement if displacement is not None else [0.0, 1.0, 2.0])
    force = np.zeros((h.size, 3))
    force[:, 2] = force_z
    return ObjectForcePath.from_samples(h, force, np.zeros_like(force))


def test_release_integrates_gravity_adhesion_and_dissipation() -> None:
    path = _path(np.array([2.0, 2.0, 2.0]))
    result = path.evaluate_release(
        mass_kg=1.0,
        gravity_m_s2=1.0,
        adhesion=AdhesionProfile.finite_range_constant(0.25, 0.5),
        dissipation_work_J=np.array([0.0, 0.25, 0.5]),
    )

    np.testing.assert_allclose(result.electrostatic_work_J, [0.0, 2.0, 4.0])
    np.testing.assert_allclose(result.adhesion_work_J, [0.0, 0.125, 0.125])
    np.testing.assert_allclose(result.dissipation_work_J, [0.0, 0.25, 0.5])
    np.testing.assert_allclose(result.available_energy_J, [0.0, 0.625, 1.375])
    assert result.barrier_free_from_rest
    assert result.first_inaccessible_displacement_m is None
    assert result.instantaneous_force_margin_N == pytest.approx(0.75)
    assert result.numerically_qualified


@pytest.mark.parametrize("scale", [1.0, 1.0e-18], ids=["unit", "physical"])
def test_release_detects_force_curve_barrier_at_relevant_scales(scale: float) -> None:
    path = _path(scale * np.array([0.0, 4.0, 4.0]))
    result = path.evaluate_release(
        mass_kg=1.0,
        gravity_m_s2=scale,
        adhesion=AdhesionProfile.none(),
        energy_tolerance_J=0.0,
    )

    assert result.available_energy_J[-1] > 0.0
    assert result.minimum_available_energy_J == pytest.approx(
        -0.125 * scale,
        rel=1.0e-12,
        abs=1.0e-14 * scale,
    )
    assert not result.barrier_free_from_rest
    assert result.first_inaccessible_displacement_m == pytest.approx(0.0, abs=1.0e-12)
    assert result.endpoint_positive
    assert not result.endpoint_reachable_from_rest
    assert result.maximum_reachable_speed_m_s == pytest.approx(0.0)


def test_release_checks_adhesion_breakpoints_between_force_samples() -> None:
    path = _path(np.array([2.0, 2.0]), np.array([0.0, 2.0]))
    profiles = (
        (AdhesionProfile.finite_range_constant(3.0, 0.5), -0.5),
        (
            AdhesionProfile.tabulated(
                np.array([0.0, 0.5, 1.0, 2.0]),
                np.array([4.0, 4.0, 0.0, 0.0]),
            ),
            -1.25,
        ),
    )

    for adhesion, minimum_energy in profiles:
        result = path.evaluate_release(
            mass_kg=1.0,
            gravity_m_s2=0.0,
            adhesion=adhesion,
            energy_tolerance_J=0.0,
        )

        assert result.available_energy_J[-1] > 0.0
        assert result.minimum_available_energy_J == pytest.approx(minimum_energy)
        assert not result.barrier_free_from_rest


def test_translation_partition_changes_speed_but_not_barrier() -> None:
    path = _path(np.array([2.0, 2.0, 2.0]))
    full = path.evaluate_release(1.0, 1.0, AdhesionProfile.none())
    half = path.evaluate_release(
        1.0, 1.0, AdhesionProfile.none(), eta_translation=0.5
    )

    np.testing.assert_allclose(full.available_energy_J, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(half.available_energy_J, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(full.speed_m_s, [0.0, np.sqrt(2.0), 2.0])
    np.testing.assert_allclose(half.speed_m_s, [0.0, 1.0, np.sqrt(2.0)])
    assert full.barrier_free_from_rest
    assert half.barrier_free_from_rest


def test_release_qualification_follows_source_path_status() -> None:
    path = replace(_path(np.array([1.0, 1.0, 1.0])), status="not_converged")

    result = path.evaluate_release(1.0, 0.0, AdhesionProfile.none())

    assert result.source_path_status == "not_converged"
    assert not result.numerically_qualified


def test_adhesion_tables_are_validated() -> None:
    invalid_tables = (
        ([0.0], [1.0], "at least two"),
        ([0.1, 1.0], [1.0, 0.0], "start at zero"),
        ([0.0, 1.0, 1.0], [1.0, 1.0, 0.0], "strictly increasing"),
        ([0.0, 1.0], [1.0, -1.0], "non-negative"),
    )

    for displacement, force, message in invalid_tables:
        with pytest.raises(ValueError, match=message):
            AdhesionProfile.tabulated(displacement, force)


def test_dissipation_work_is_validated() -> None:
    invalid_work = (
        (np.array([0.0, 1.0]), "shape"),
        (np.array([1.0, 1.0, 1.0]), "start at zero"),
        (np.array([0.0, 2.0, 1.0]), "non-decreasing"),
    )

    for dissipation, message in invalid_work:
        with pytest.raises(ValueError, match=message):
            _path(np.array([2.0, 2.0, 2.0])).evaluate_release(
                1.0,
                0.0,
                AdhesionProfile.none(),
                dissipation_work_J=dissipation,
            )


def test_records_defensively_copy_and_freeze_arrays_and_mappings() -> None:
    force = np.array([1.0, 2.0, 3.0])
    component = WrenchComponent(force_N=force, torque_Nm=np.zeros(3))
    wrench = ObjectWrench(
        mesh_id=2,
        step=4,
        total_charge_C=5.0,
        force_N=force,
        torque_Nm=np.zeros(3),
        torque_origin_m=np.zeros(3),
        components={"total_external": component},
        numerical_metadata={"method": "test"},
    )
    path = _path(np.array([2.0, 2.0, 2.0]))
    result = path.evaluate_release(1.0, 0.0, AdhesionProfile.none())
    force[:] = 9.0

    np.testing.assert_array_equal(component.force_N, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(wrench.force_N, [1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        component.force_N[0] = 4.0
    with pytest.raises(ValueError):
        wrench.force_N[0] = 4.0
    with pytest.raises(TypeError):
        wrench.components["other"] = component
    with pytest.raises(ValueError):
        result.available_energy_J[0] = 1.0
    with pytest.raises(FrozenInstanceError):
        component.potential_energy_J = 2.0


def test_wrench_defensively_freezes_transform_and_nested_metadata() -> None:
    transform = RigidTransform.translation([1.0, 2.0, 3.0])
    nested_array = np.array([4.0, 5.0])
    large_integer = np.array([2**60 + 1], dtype=np.int64)
    complex_value = np.array([1.0 + 2.0j], dtype=np.complex128)
    nested_list = [nested_array]
    mutable_bytes = bytearray(b"abc")
    wrench = ObjectWrench(
        mesh_id=1,
        step=None,
        total_charge_C=0.0,
        force_N=np.zeros(3),
        torque_Nm=np.zeros(3),
        torque_origin_m=np.zeros(3),
        transform=transform,
        numerical_metadata={
            "nested": {"values": nested_list},
            "large_integer": large_integer,
            "complex_value": complex_value,
            "labels": {"a", "b"},
            "bytes": mutable_bytes,
        },
    )
    transform.translation_m[:] = 9.0
    nested_array[:] = 9.0
    large_integer[:] = 0
    complex_value[:] = 0.0
    nested_list.append(np.array([6.0]))
    mutable_bytes[:] = b"xyz"

    assert isinstance(wrench.transform, RigidTransform)
    np.testing.assert_array_equal(wrench.transform.translation_m, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(
        wrench.numerical_metadata["nested"]["values"][0],  # type: ignore[index]
        [4.0, 5.0],
    )
    assert wrench.numerical_metadata["labels"] == frozenset({"a", "b"})
    assert wrench.numerical_metadata["bytes"] == b"abc"
    assert wrench.numerical_metadata["large_integer"].dtype == np.dtype(  # type: ignore[union-attr]
        np.int64
    )
    assert wrench.numerical_metadata["large_integer"][0] == 2**60 + 1  # type: ignore[index]
    assert wrench.numerical_metadata["complex_value"][0] == 1.0 + 2.0j  # type: ignore[index]
    with pytest.raises(ValueError):
        wrench.transform.translation_m[0] = 0.0
    with pytest.raises(ValueError):
        wrench.numerical_metadata["nested"]["values"][0][0] = 0.0  # type: ignore[index]


def test_wrench_rejects_unfreezable_metadata_objects() -> None:
    invalid_values = (
        (object(), "metadata value"),
        (np.array([object()], dtype=object), "object dtype"),
    )

    for value, message in invalid_values:
        with pytest.raises(TypeError, match=message):
            ObjectWrench(
                mesh_id=1,
                step=None,
                total_charge_C=0.0,
                force_N=np.zeros(3),
                torque_Nm=np.zeros(3),
                torque_origin_m=np.zeros(3),
                numerical_metadata={"bad": value},
            )
