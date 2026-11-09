"""Periodic pair torques against an independent finite-image surface integral."""

from pathlib import Path

import numpy as np
import pytest

from beach import K_COULOMB, RunResult, calc_coulomb
from beach.analysis.field_kernel import field_kernel_from_result
from beach.io.types import RunField, RunMesh, RunSummary


def _panel_samples(
    vertices: np.ndarray, charge: float
) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(16)
    u, v = np.meshgrid((nodes + 1.0) / 2.0, (nodes + 1.0) / 2.0, indexing="ij")
    wu, wv = np.meshgrid(weights / 2.0, weights / 2.0, indexing="ij")
    a, b, c = vertices
    points = a + u[..., None] * (b - a) + ((1.0 - u) * v)[..., None] * (c - a)
    return points.reshape(-1, 3), (2.0 * charge * (1.0 - u) * wu * wv).ravel()


def _image_wrench(target, source, origin):
    points, weights = target
    source_points, source_weights = source
    forces = np.zeros_like(points)
    for nx in range(-1, 2):
        for ny in range(-1, 2):
            d = (
                points[:, None, :]
                - source_points[None, :, :]
                - [2.0 * nx, 3.0 * ny, 0.0]
            )
            forces += (
                K_COULOMB
                * weights[:, None]
                * np.sum(
                    source_weights[None, :, None]
                    * d
                    / np.linalg.norm(d, axis=2)[..., None] ** 3,
                    axis=1,
                )
            )
    return forces.sum(axis=0), np.cross(points - origin, forces).sum(axis=0)


@pytest.mark.parametrize(
    "borrow_kernel", [False, True], ids=["source-only", "borrowed"]
)
@pytest.mark.parametrize("torque_origin", ["origin", "target_center", "source_center"])
def test_periodic_coulomb_evaluates_both_torques(
    tmp_path: Path,
    field_kernel_library: Path,
    borrow_kernel: bool,
    torque_origin: str,
) -> None:
    (tmp_path / "beach.toml").write_text(
        """
[domain]
box_min = [0.0, 0.0, -1.0]
box_max = [2.0, 3.0, 1.0]
periodic_axes = ["x", "y"]
[fields]
boundary = "periodic2"
[fields.external]
electric_v_m = [0.3, -0.4, 0.2]
[fields.solver]
method = "fmm"
[fields.solver.tree]
theta = 0.01
leaf_max = 64
[fields.periodic]
backend = "finite_images"
image_layers = 1
""",
        encoding="utf-8",
    )
    centers = np.array([[0.2, 0.3, 0.1], [0.8, 0.7, 0.4]])
    triangles = centers[:, None, :] + [
        [-0.02, -0.01, 0.0],
        [0.04, -0.01, 0.0],
        [-0.02, 0.02, 0.0],
    ]
    charges = np.array([1.0e-9, -2.0e-9])
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
        mesh=RunMesh(triangles=triangles, mesh_ids=np.array([1, 2])),
        field=RunField(
            source_model="triangle_p0", kernel_id="triangle_p0_exact_p2m_near"
        ),
    )
    origin = {
        "origin": np.zeros(3),
        "target_center": centers[0],
        "source_center": centers[1],
    }[torque_origin]
    panels = [
        _panel_samples(triangle, charge) for triangle, charge in zip(triangles, charges)
    ]
    expected_a = _image_wrench(panels[0], panels[1], origin)
    expected_b = _image_wrench(panels[1], panels[0], origin)
    # This geometry has non-central image forces: the two cell torques do not cancel.
    assert np.linalg.norm(expected_a[1] + expected_b[1]) > 1.0e-10
    if borrow_kernel:
        with field_kernel_from_result(
            result, step=None, library_path=field_kernel_library
        ) as kernel:
            # The pair interaction must exclude a prescribed external field on both groups.
            interaction = calc_coulomb(
                result, 1, 2, step=None, torque_origin=torque_origin, kernel=kernel
            )
            np.testing.assert_array_equal(kernel.source_charges, charges)
    else:
        interaction = calc_coulomb(
            result,
            1,
            2,
            step=None,
            torque_origin=torque_origin,
            library_path=field_kernel_library,
        )
    for actual, expected in (
        (interaction.force_on_a_N, expected_a[0]),
        (interaction.force_on_b_N, expected_b[0]),
        (interaction.torque_on_a_Nm, expected_a[1]),
        (interaction.torque_on_b_Nm, expected_b[1]),
    ):
        np.testing.assert_allclose(actual, expected, rtol=2.0e-9, atol=1.0e-20)
