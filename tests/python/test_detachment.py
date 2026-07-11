from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from beach import (
    AdhesionProfile,
    ObjectForcePath,
    ObjectWrench,
    WrenchComponent,
)


def _path(force_z: np.ndarray, displacement: np.ndarray | None = None) -> ObjectForcePath:
    h = np.asarray(displacement if displacement is not None else [0.0, 1.0, 2.0])
    force = np.zeros((h.size, 3))
    force[:, 2] = force_z
    return ObjectForcePath.from_samples(h, force, np.zeros_like(force))


def test_release_integrates_work_and_finite_range_resistance() -> None:
    path = _path(np.array([2.0, 2.0, 2.0]))
    result = path.evaluate_release(
        mass_kg=1.0,
        gravity_m_s2=1.0,
        adhesion=AdhesionProfile.finite_range_constant(0.25, 0.5),
    )

    np.testing.assert_allclose(result.electrostatic_work_J, [0.0, 2.0, 4.0])
    np.testing.assert_allclose(result.adhesion_work_J, [0.0, 0.125, 0.125])
    np.testing.assert_allclose(result.available_energy_J, [0.0, 0.875, 1.875])
    assert result.barrier_free_from_rest
    assert result.first_inaccessible_displacement_m is None
    assert result.endpoint_reachable_from_rest


def test_release_detects_force_curve_barrier_between_samples() -> None:
    path = _path(np.array([0.0, 4.0, 4.0]))
    result = path.evaluate_release(
        mass_kg=1.0,
        gravity_m_s2=1.0,
        adhesion=AdhesionProfile.none(),
        energy_tolerance_J=0.0,
    )

    assert result.available_energy_J[-1] > 0.0
    assert result.minimum_available_energy_J == pytest.approx(-0.125)
    assert not result.barrier_free_from_rest
    assert result.first_inaccessible_displacement_m == pytest.approx(0.0, abs=1.0e-12)
    assert result.maximum_reachable_speed_m_s == pytest.approx(0.0)


def test_continuous_barrier_detection_is_invariant_at_physical_energy_scales() -> None:
    path = _path(np.array([0.0, 4.0e-18, 4.0e-18]))
    result = path.evaluate_release(
        mass_kg=1.0,
        gravity_m_s2=1.0e-18,
        adhesion=AdhesionProfile.none(),
        energy_tolerance_J=0.0,
    )

    assert result.minimum_available_energy_J == pytest.approx(
        -1.25e-19, rel=1.0e-12, abs=1.0e-32
    )
    assert not result.barrier_free_from_rest


def test_release_includes_finite_range_endpoint_in_continuous_barrier_check() -> None:
    path = _path(np.array([2.0, 2.0]), np.array([0.0, 2.0]))
    result = path.evaluate_release(
        mass_kg=1.0,
        gravity_m_s2=0.0,
        adhesion=AdhesionProfile.finite_range_constant(3.0, 0.5),
        energy_tolerance_J=0.0,
    )

    assert result.available_energy_J[-1] > 0.0
    assert result.minimum_available_energy_J == pytest.approx(-0.5)
    assert not result.barrier_free_from_rest


def test_release_includes_tabulated_adhesion_knots_between_force_samples() -> None:
    path = _path(np.array([2.0, 2.0]), np.array([0.0, 2.0]))
    adhesion = AdhesionProfile.tabulated(
        np.array([0.0, 0.5, 1.0, 2.0]),
        np.array([4.0, 4.0, 0.0, 0.0]),
    )
    result = path.evaluate_release(
        mass_kg=1.0,
        gravity_m_s2=0.0,
        adhesion=adhesion,
        energy_tolerance_J=0.0,
    )

    assert result.available_energy_J[-1] > 0.0
    assert result.minimum_available_energy_J == pytest.approx(-1.25)
    assert not result.barrier_free_from_rest


def test_translation_partition_changes_speed_but_not_barrier() -> None:
    path = _path(np.array([2.0, 2.0, 2.0]))
    full = path.evaluate_release(1.0, 1.0, AdhesionProfile.none())
    half = path.evaluate_release(
        1.0, 1.0, AdhesionProfile.none(), eta_translation=0.5
    )

    np.testing.assert_allclose(half.available_energy_J, full.available_energy_J)
    np.testing.assert_allclose(half.speed_m_s, np.sqrt(0.5) * full.speed_m_s)
    assert half.barrier_free_from_rest is full.barrier_free_from_rest


@pytest.mark.parametrize(
    "displacement, force, message",
    [
        ([0.0], [1.0], "at least two"),
        ([0.1, 1.0], [1.0, 0.0], "start at zero"),
        ([0.0, 1.0, 1.0], [1.0, 1.0, 0.0], "strictly increasing"),
        ([0.0, 1.0], [1.0, -1.0], "non-negative"),
    ],
)
def test_adhesion_tables_are_validated(displacement, force, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        AdhesionProfile.tabulated(displacement, force)


@pytest.mark.parametrize(
    "dissipation, message",
    [
        (np.array([0.0, 1.0]), "shape"),
        (np.array([1.0, 1.0, 1.0]), "start at zero"),
        (np.array([0.0, 2.0, 1.0]), "non-decreasing"),
    ],
)
def test_dissipation_work_is_validated(dissipation: np.ndarray, message: str) -> None:
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

    np.testing.assert_array_equal(wrench.force_N, [1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        wrench.force_N[0] = 4.0
    with pytest.raises(TypeError):
        wrench.components["other"] = component
    with pytest.raises(ValueError):
        result.available_energy_J[0] = 1.0
    with pytest.raises(FrozenInstanceError):
        component.potential_energy_J = 2.0
