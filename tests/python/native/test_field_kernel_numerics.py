"""FieldKernel numerics against the real native kernel (the native gate)."""

from pathlib import Path

import numpy as np
import pytest

from beach import (
    BeachScene,
    FieldKernel,
    FieldKernelOptions,
    RunResult,
    calc_object_forces_kernel,
)
from beach.analysis.field_kernel import field_kernel_from_result
from beach.analysis.panel_quadrature import panel_target_quadrature
from beach.cli.main import main as beachx_main
from beach.io.constants import K_COULOMB
from beach.io.types import RunField, RunMesh, RunSummary


def _uniform_triangle_field(
    vertices: np.ndarray, charge_C: float, target: np.ndarray
) -> np.ndarray:
    """Field of a uniformly charged triangle by tensor Gauss-Legendre on the Duffy square."""

    nodes, weights = np.polynomial.legendre.leggauss(64)
    s = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    a, b, c = vertices
    twice_area = float(np.linalg.norm(np.cross(b - a, c - a)))
    # x(s, t) = a + s (b - a) + s t (c - b) has Jacobian s * twice_area on [0, 1]^2.
    points = (
        a[None, None, :]
        + s[:, None, None] * (b - a)[None, None, :]
        + (s[:, None] * s[None, :])[:, :, None] * (c - b)[None, None, :]
    )
    weight = (w[:, None] * w[None, :]) * s[:, None] * twice_area
    separation = target[None, None, :] - points
    distance = np.linalg.norm(separation, axis=2)
    density = charge_C / (0.5 * twice_area)
    integrand = separation / distance[:, :, None] ** 3
    return K_COULOMB * density * np.sum(weight[:, :, None] * integrand, axis=(0, 1))


@pytest.mark.parametrize("edge_m", [1.0e-6, 1.0e-9])
def test_native_kernel_accepts_small_and_translated_panels(
    field_kernel_library: Path, edge_m: float
) -> None:
    lib = field_kernel_library
    local = edge_m * np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    target_local = edge_m * np.array([[0.3, 0.2, 0.7]])
    fields = []
    for shift in (np.zeros(3), np.array([10.0, -10.0, 10.0])):
        with FieldKernel((shift + local)[None], np.array([1.0e-18]), library_path=lib) as kernel:
            fields.append(kernel.eval_e_direct(shift + target_local))
    # Translating by 10 m keeps ~eps*10/edge relative resolution of the local geometry.
    np.testing.assert_allclose(fields[1], fields[0], rtol=1.0e3 * 10.0 * 2.3e-16 / edge_m)

    # A sliver just above the shared 128 eps lref^2 threshold passes both checks.
    height = 400.0 * np.finfo(np.float64).eps
    sliver = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, height, 0.0]]])
    with FieldKernel(sliver, np.array([1.0e-12]), library_path=lib):
        pass


def test_field_kernel_exact_panel_direct_scales_with_charge_refresh(
    field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    centers = np.array(
        [[-0.6, -0.2, 0.0], [0.0, 0.4, 0.1], [0.7, -0.3, -0.1]],
        dtype=float,
    )
    triangles = np.stack(
        (
            centers + [-0.05, -0.05, 0.0],
            centers + [0.05, -0.05, 0.0],
            centers + [-0.05, 0.05, 0.0],
        ),
        axis=1,
    )
    source_q = np.array([1.0e-9, -1.2e-9, 0.8e-9])
    targets = np.array([[0.13, -0.08, 0.4], [-0.4, 0.2, 0.3]], dtype=float)
    with FieldKernel(triangles, source_q, library_path=lib) as kernel:
        ordinary_e = kernel.eval_e(targets)
        direct_e = kernel.eval_e_direct(targets)
        direct_phi = kernel.eval_phi_direct(targets)
        kernel.update_charges(-0.25 * source_q)
        refreshed_e = kernel.eval_e_direct(targets)
        refreshed_phi = kernel.eval_phi_direct(targets)

    # Mixed-sign panels and off-plane targets: the ordinary evaluation equals the exact panel sum.
    np.testing.assert_allclose(ordinary_e, direct_e, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(refreshed_e, -0.25 * direct_e, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(
        refreshed_phi, -0.25 * direct_phi, rtol=2.0e-14, atol=1.0e-14
    )


def test_field_kernel_fmm_eval_follows_charge_refresh(
    field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    rng = np.random.default_rng(3)
    corners = rng.uniform(-1.0, 1.0, size=(40, 1, 3))
    triangles = corners + np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.02]])
    source_q = rng.standard_normal(40) * 1.0e-9
    refreshed_q = rng.standard_normal(40) * 1.0e-9
    targets = rng.uniform(-1.5, 1.5, size=(6, 3))
    options = FieldKernelOptions(theta=0.5, leaf_max=4, order=4)
    with FieldKernel(triangles, refreshed_q, options=options, library_path=lib) as fresh:
        expected_e = fresh.eval_e(targets)
        expected_phi = fresh.eval_phi(targets)
        expected_force = fresh.force_on_charges(targets, source_q[:6])

    # FMM after FMM, and FMM after exact-direct use, must both see the refreshed charges.
    with FieldKernel(triangles, source_q, options=options, library_path=lib) as kernel:
        stale_e = kernel.eval_e(targets)
        kernel.update_charges(refreshed_q)
        np.testing.assert_array_equal(kernel.eval_e(targets), expected_e)
        np.testing.assert_array_equal(kernel.eval_phi(targets), expected_phi)
    with FieldKernel(triangles, source_q, options=options, library_path=lib) as kernel:
        kernel.eval_e_direct(targets)
        kernel.update_charges(refreshed_q)
        force, torque = kernel.force_on_charges(targets, source_q[:6])
        np.testing.assert_array_equal(force, expected_force[0])
        np.testing.assert_array_equal(torque, expected_force[1])
        np.testing.assert_array_equal(kernel.eval_e(targets), expected_e)
    assert not np.array_equal(stale_e, expected_e)


def test_field_kernel_from_result_uses_triangle_geometry(
    tmp_path: Path, field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
        dtype=float,
    )
    charges = np.array([2.5e-9], dtype=float)
    result = RunResult(
        directory=tmp_path,
        charges=charges,
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
        field=RunField(
            source_model="triangle_p0", kernel_id="triangle_p0_exact_p2m_near"
        ),
    )
    target = np.array([[0.25, 0.2, 0.4]], dtype=float)

    with field_kernel_from_result(result, step=None, library_path=lib) as panel_kernel:
        panel_field = panel_kernel.eval_e(target)
        direct_field = panel_kernel.eval_e_direct(target)

    np.testing.assert_allclose(panel_field, direct_field, rtol=2.0e-14, atol=1.0e-14)
    # Both paths share the factory geometry, so the saved triangle is also checked against an
    # independent quadrature of the uniformly charged surface.
    expected = _uniform_triangle_field(triangles[0], charges[0], target[0])
    np.testing.assert_allclose(panel_field[0], expected, rtol=1.0e-9)


def test_panel_direct_and_ordinary_use_same_principal_value_at_target_quadrature(
    field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
        dtype=float,
    )
    charges = np.array([2.0e-9])
    points, _, _ = panel_target_quadrature(triangles, charges, order=7)

    with FieldKernel(
        triangles,
        charges,
        library_path=lib,
    ) as kernel:
        ordinary_e = kernel.eval_e(points)
        direct_e = kernel.eval_e_direct(points)
        ordinary_phi = kernel.eval_phi(points)
        direct_phi = kernel.eval_phi_direct(points)

    np.testing.assert_allclose(ordinary_e, direct_e, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(ordinary_phi, direct_phi, rtol=2.0e-14, atol=1.0e-14)


def test_field_kernel_periodic2_wraps_batched_equivalent_targets(
    field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [0.75, 0.0, 0.0], [0.0, 0.75, 0.0]]],
        dtype=float,
    )
    charges = np.array([2.0e-9], dtype=float)
    targets = np.array(
        [
            [0.25, 0.25, 1.0],
            [1.25, -0.75, 1.0],
        ],
        dtype=float,
    )
    options = FieldKernelOptions(
        periodic2=((0, 1), (1.0, 1.0), (0.0, 0.0), 1, "none", 0.0, 4),
        box_min=(0.0, 0.0, -1.0),
        box_max=(1.0, 1.0, 1.0),
    )

    with FieldKernel(
        triangles,
        charges,
        options=options,
        library_path=lib,
    ) as kernel:
        potential = kernel.eval_phi(targets)
        field = kernel.eval_e(targets)

    np.testing.assert_allclose(potential[0], potential[1], rtol=1.0e-14)
    np.testing.assert_allclose(field[0], field[1], rtol=1.0e-14, atol=1.0e-14)


def test_field_kernel_adds_uniform_external_e0(field_kernel_library: Path) -> None:
    lib = field_kernel_library
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0, 0.2, 0.0]]],
        dtype=float,
    )
    source_q = np.array([0.0], dtype=float)
    target = np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0]], dtype=float)
    target_q = np.array([2.0, 3.0], dtype=float)

    with FieldKernel(
        triangles,
        source_q,
        options=FieldKernelOptions(external_e0=(1.0, 0.0, 0.0)),
        library_path=lib,
    ) as kernel:
        field = kernel.eval_e(target)
        direct_field = kernel.eval_e_direct(target)
        force, torque = kernel.force_on_charges(target, target_q, origin=(0.0, 0.0, 0.0))

    np.testing.assert_allclose(field, np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    np.testing.assert_array_equal(direct_field, np.zeros((2, 3)))
    np.testing.assert_allclose(force, np.array([5.0, 0.0, 0.0]))
    np.testing.assert_allclose(torque, np.array([0.0, 0.0, -8.0]))


def test_calc_object_forces_kernel_excludes_self_sources(
    tmp_path: Path, field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    (tmp_path / "beach.toml").write_text(
        '[sim]\nfield_bc_mode = "free"\n',
        encoding="utf-8",
    )
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
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
        mesh=RunMesh(triangles=triangles, mesh_ids=np.array([1, 2], dtype=np.int64)),
    )

    records = calc_object_forces_kernel(result, step=None, library_path=lib)

    assert [record.mesh_id for record in records] == [1, 2]
    assert records[0].force_N[0] < 0.0
    assert records[1].force_N[0] > 0.0
    np.testing.assert_allclose(
        records[0].force_N, -records[1].force_N, rtol=1.0e-12, atol=1.0e-20
    )
    assert records[0].total_charge_C == pytest.approx(charges[0])


def test_scene_move_updates_kernel_force_geometry(
    tmp_path: Path, field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    (tmp_path / "beach.toml").write_text(
        '[sim]\nfield_bc_mode = "free"\n',
        encoding="utf-8",
    )
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
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
        mesh=RunMesh(triangles=triangles, mesh_ids=np.array([1, 2], dtype=np.int64)),
    )

    scene = BeachScene.from_result(result, step=None)
    moved = scene.move(2, by=[1.0, 0.0, 0.0])
    original_records = scene.calc_object_forces_kernel(
        target_mesh_ids=1,
        library_path=lib,
    )
    reference_records = calc_object_forces_kernel(
        result,
        step=None,
        target_mesh_ids=1,
        library_path=lib,
    )
    records = moved.calc_object_forces_kernel(
        target_mesh_ids=1,
        library_path=lib,
    )

    assert records[0].force_N[0] < 0.0
    assert abs(records[0].force_N[0]) < abs(original_records[0].force_N[0])
    np.testing.assert_allclose(
        original_records[0].force_N,
        reference_records[0].force_N,
        rtol=1.0e-12,
        atol=1.0e-24,
    )
    np.testing.assert_allclose(scene.centers[1], triangles[1].mean(axis=0))


def test_kernel_forces_cli_writes_csv(tmp_path: Path, field_kernel_library: Path) -> None:
    lib = field_kernel_library
    out = tmp_path / "run"
    out.mkdir()
    config_path = tmp_path / "beach.toml"
    config_path.write_text("[sim]\ntree_leaf_max = 8\n", encoding="utf-8")
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=0",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-9\n2,2.0e-9\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,0.0,0.3,0.0,0.0,0.0,0.3,1.0e-9,1\n"
        "2,1.0,0.0,0.0,1.0,0.3,0.0,1.0,0.0,0.3,2.0e-9,2\n",
        encoding="utf-8",
    )
    save_csv = tmp_path / "forces.csv"

    beachx_main(
        [
            "kernel-forces",
            str(out),
            "--config",
            str(config_path),
            "--library",
            str(lib),
            "--save",
            str(save_csv),
        ]
    )

    text = save_csv.read_text(encoding="utf-8")
    assert "mesh_id,step,total_charge_C" in text
    assert "1,-1,1e-09" in text
