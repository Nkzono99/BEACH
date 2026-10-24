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
from beach.io.types import RunField, RunMesh, RunSummary


def test_field_kernel_eval_e_matches_exact_panel_direct(field_kernel_library: Path) -> None:
    lib = field_kernel_library
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0, 0.2, 0.0]],
            [[1.0, 0.0, 0.0], [1.2, 0.0, 0.0], [1.0, 0.2, 0.0]],
        ],
        dtype=float,
    )
    source_q = np.array([1.0e-9, -2.0e-9], dtype=float)
    target = np.array([[0.0, 1.0, 0.3]], dtype=float)

    with FieldKernel(triangles, source_q, library_path=lib) as kernel:
        field = kernel.eval_e(target)
        exact = kernel.eval_e_direct(target)

    np.testing.assert_allclose(field, exact, rtol=2.0e-14, atol=1.0e-14)


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
        direct_e = kernel.eval_e_direct(targets)
        direct_phi = kernel.eval_phi_direct(targets)
        kernel.update_charges(-0.25 * source_q)
        refreshed_e = kernel.eval_e_direct(targets)
        refreshed_phi = kernel.eval_phi_direct(targets)

    np.testing.assert_allclose(refreshed_e, -0.25 * direct_e, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(
        refreshed_phi, -0.25 * direct_phi, rtol=2.0e-14, atol=1.0e-14
    )


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
