from pathlib import Path

import numpy as np
import pytest

import beach.fortran_results.object_interaction as interaction_module
from beach import (
    FieldKernelError,
    FortranRunResult,
    ObjectInteractionSnapshot,
    RigidTransform,
)
from beach.fortran_results.panel_quadrature import panel_target_quadrature


def _kernel_lib() -> Path:
    path = Path("build/libbeach_field_kernel.so")
    if not path.exists():
        pytest.skip("field kernel shared library is not built; run `make build-kernel`")
    return path


def _result(
    directory: Path,
    triangles: np.ndarray,
    charges: np.ndarray,
    mesh_ids: np.ndarray,
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
        mesh_ids=np.asarray(mesh_ids, dtype=np.int64),
        field_source_model="triangle_p0",
        field_kernel_id="triangle_p0_exact_p2m_near",
    )


def _write_free_config(
    path: Path,
    *,
    e0: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> None:
    path.write_text(
        f"""
[sim]
field_bc_mode = "free"
softening = 0.0
tree_theta = 0.3
tree_leaf_max = 8
box_min = [-10.0, -10.0, -10.0]
box_max = [10.0, 10.0, 10.0]
e0 = [{e0[0]}, {e0[1]}, {e0[2]}]
""".strip()
        + "\n",
        encoding="utf-8",
    )


def _unit_triangle(z: float, *, offset_x: float = 0.0) -> np.ndarray:
    return np.array(
        [
            [offset_x + 0.0, 0.0, z],
            [offset_x + 1.0, 0.0, z],
            [offset_x + 0.0, 1.0, z],
        ],
        dtype=float,
    )


@pytest.mark.parametrize("order", [3, 7])
def test_gauss_duffy_conserves_each_panel_charge(order: int) -> None:
    triangles = np.array([_unit_triangle(0.0), 2.0 * _unit_triangle(1.0)])
    charges = np.array([2.5e-9, -1.25e-9])

    points, charge_weights, element_index = panel_target_quadrature(
        triangles,
        charges,
        order,
    )

    assert points.shape == (2 * order * order, 3)
    assert charge_weights.shape == (2 * order * order,)
    assert element_index.shape == (2 * order * order,)
    for index, charge in enumerate(charges):
        assert np.sum(charge_weights[element_index == index]) == pytest.approx(
            charge,
            rel=1.0e-15,
            abs=1.0e-30,
        )


@pytest.mark.parametrize("order", [3.9, "3", None])
def test_gauss_duffy_rejects_non_integer_orders(order: object) -> None:
    with pytest.raises(ValueError, match="order"):
        panel_target_quadrature(
            np.array([_unit_triangle(0.0)]),
            np.array([1.0e-9]),
            order,  # type: ignore[arg-type]
        )


def test_single_triangle_primary_exclusion_leaves_only_uniform_wrench(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    e0 = np.array([2.0, -3.0, 4.0])
    _write_free_config(config, e0=tuple(e0))
    charge = 2.0e-9
    result = _result(
        tmp_path,
        np.array([_unit_triangle(0.0)]),
        np.array([charge]),
        np.array([1]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        wrench = snapshot.object_probe(1).wrench()

    np.testing.assert_allclose(wrench.force_N, charge * e0, rtol=2.0e-12, atol=1.0e-20)
    np.testing.assert_allclose(wrench.torque_Nm, 0.0, atol=1.0e-20)
    np.testing.assert_allclose(
        wrench.components["target_periodic_images"].force_N,
        0.0,
        atol=1.0e-20,
    )
    assert wrench.numerical_metadata["target_integration"] == "gauss_duffy"
    assert wrench.numerical_metadata["quadrature_order"] == 7


def test_orders_three_and_seven_converge_for_separated_panels(tmp_path: Path) -> None:
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    result = _result(
        tmp_path,
        np.array([_unit_triangle(0.0), _unit_triangle(6.5)]),
        np.array([1.0e-9, 2.0e-9]),
        np.array([1, 2]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        order3 = snapshot.object_probe(2, quadrature_order=3).wrench()
        order7 = snapshot.object_probe(2, quadrature_order=7).wrench()

    assert np.linalg.norm(order3.force_N - order7.force_N) <= (
        1.0e-2 * np.linalg.norm(order7.force_N) + 1.0e-20
    )
    assert np.linalg.norm(order3.torque_Nm - order7.torque_Nm) <= (
        1.0e-2 * np.linalg.norm(order7.torque_Nm) + 1.0e-20
    )


def test_near_surface_centroid_compatibility_is_labeled_and_differs(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    result = _result(
        tmp_path,
        np.array([_unit_triangle(0.0), _unit_triangle(0.05)]),
        np.array([1.0e-9, 1.0e-9]),
        np.array([1, 2]),
    )

    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
        library_path=_kernel_lib(),
    ) as snapshot:
        centroid = snapshot.object_probe(
            2,
            target_integration="centroid_compatibility",
        ).wrench()
        quadrature = snapshot.object_probe(2, quadrature_order=7).wrench()

    assert centroid.numerical_metadata["target_integration"] == "centroid_compatibility"
    assert quadrature.numerical_metadata["target_integration"] == "gauss_duffy"
    difference = np.linalg.norm(centroid.force_N - quadrature.force_N)
    assert difference > 1.0e-3 * np.linalg.norm(quadrature.force_N)


class _CovariantFieldKernel:
    mode = "radial"

    def __init__(self, _positions, charges, **_kwargs) -> None:
        self.current = np.array(charges, copy=True)
        self.closed = False

    def update_charges(self, charges: np.ndarray) -> None:
        self.current = np.array(charges, copy=True)

    def eval_e(self, points: np.ndarray) -> np.ndarray:
        if float(np.sum(self.current)) != 2.0:
            return np.zeros_like(points)
        if type(self).mode == "nonlinear_radial":
            radius2 = np.sum(points * points, axis=1)
            return radius2[:, None] * points
        if type(self).mode == "radial":
            return np.array(points, copy=True)
        return np.broadcast_to(np.array([1.0, -2.0, 0.5]), points.shape).copy()

    def eval_phi(self, points: np.ndarray) -> np.ndarray:
        if float(np.sum(self.current)) != 2.0:
            return np.zeros(len(points))
        if type(self).mode == "nonlinear_radial":
            radius2 = np.sum(points * points, axis=1)
            return -0.25 * radius2 * radius2
        if type(self).mode == "radial":
            return -0.5 * np.sum(points * points, axis=1)
        return -(points @ np.array([1.0, -2.0, 0.5]))

    def eval_e_direct(self, points: np.ndarray) -> np.ndarray:
        return np.zeros_like(points)

    def eval_phi_direct(self, points: np.ndarray) -> np.ndarray:
        return np.zeros(len(points))

    def diagnostics(self):
        raise FieldKernelError("not available")

    def close(self) -> None:
        self.closed = True


def test_triangle_wrench_obeys_rigid_rotation_and_torque_origin_covariance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(interaction_module, "FieldKernel", _CovariantFieldKernel)
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    result = _result(
        tmp_path,
        np.array([_unit_triangle(0.0, offset_x=1.0), _unit_triangle(-2.0)]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )
    rotation = RigidTransform.from_axis_angle(
        [0.0, 0.0, 1.0],
        angle_deg=90.0,
    )
    matrix = rotation.rotation

    _CovariantFieldKernel.mode = "nonlinear_radial"
    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    ) as snapshot:
        probe = snapshot.object_probe(1)
        base = probe.wrench(transform_origin="origin")
        rotated = probe.wrench(
            transform=rotation,
            transform_origin="origin",
        )
        origin_a = np.array([0.2, -0.1, 0.3])
        origin_b = np.array([-0.4, 0.5, -0.2])
        about_a = probe.wrench(torque_origin=origin_a)
        about_b = probe.wrench(torque_origin=origin_b)

    assert np.linalg.norm(base.torque_Nm) > 1.0e-3
    np.testing.assert_allclose(rotated.force_N, matrix @ base.force_N, atol=1.0e-13)
    np.testing.assert_allclose(rotated.torque_Nm, matrix @ base.torque_Nm, atol=1.0e-13)
    np.testing.assert_allclose(
        about_b.torque_Nm,
        about_a.torque_Nm - np.cross(origin_b - origin_a, about_a.force_N),
        atol=1.0e-13,
    )


def test_triangle_wrench_translation_moves_geometric_torque_origin(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(interaction_module, "FieldKernel", _CovariantFieldKernel)
    config = tmp_path / "beach.toml"
    _write_free_config(config)
    result = _result(
        tmp_path,
        np.array([_unit_triangle(0.0), _unit_triangle(-2.0)]),
        np.array([1.0, 2.0]),
        np.array([1, 2]),
    )
    translation = np.array([0.4, -0.3, 0.2])
    transform = RigidTransform.translation(translation)

    _CovariantFieldKernel.mode = "uniform"
    with ObjectInteractionSnapshot.from_result(
        result,
        step=None,
        config_path=config,
    ) as snapshot:
        probe = snapshot.object_probe(1)
        base = probe.wrench()
        moved = probe.wrench(transform=transform)

    np.testing.assert_allclose(moved.force_N, base.force_N, atol=1.0e-13)
    np.testing.assert_allclose(moved.torque_Nm, base.torque_Nm, atol=1.0e-13)
    np.testing.assert_allclose(
        moved.torque_origin_m,
        base.torque_origin_m + translation,
        atol=1.0e-13,
    )
