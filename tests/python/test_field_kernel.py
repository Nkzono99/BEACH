from pathlib import Path

import numpy as np
import pytest

from beach.cli.kernel_forces import main as kernel_forces_main
from beach import (
    BeachScene,
    FieldKernel,
    FieldKernelOptions,
    FortranRunResult,
    calc_object_forces_kernel,
)
from beach.fortran_results.constants import K_COULOMB


def _kernel_lib() -> Path:
    path = Path("build/libbeach_field_kernel.so")
    if not path.exists():
        pytest.skip("field kernel shared library is not built; run `make build-kernel`")
    return path


def test_field_kernel_eval_e_matches_direct_sum() -> None:
    lib = _kernel_lib()
    source_pos = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    source_q = np.array([1.0e-9, -2.0e-9], dtype=float)
    target = np.array([[0.0, 1.0, 0.0]], dtype=float)

    with FieldKernel(source_pos, source_q, library_path=lib) as kernel:
        field = kernel.eval_e(target)

    delta = target[0] - source_pos
    expected = K_COULOMB * np.sum(
        source_q[:, None] * delta / (np.linalg.norm(delta, axis=1)[:, None] ** 3),
        axis=0,
    )
    np.testing.assert_allclose(field[0], expected, rtol=1.0e-14, atol=1.0e-14)


def test_field_kernel_adds_uniform_external_e0() -> None:
    lib = _kernel_lib()
    source_pos = np.array([[0.0, 0.0, 0.0]], dtype=float)
    source_q = np.array([0.0], dtype=float)
    target = np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0]], dtype=float)
    target_q = np.array([2.0, 3.0], dtype=float)

    with FieldKernel(
        source_pos,
        source_q,
        options=FieldKernelOptions(external_e0=(1.0, 0.0, 0.0)),
        library_path=lib,
    ) as kernel:
        field = kernel.eval_e(target)
        force, torque = kernel.force_on_charges(target, target_q, origin=(0.0, 0.0, 0.0))

    np.testing.assert_allclose(field, np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    np.testing.assert_allclose(force, np.array([5.0, 0.0, 0.0]))
    np.testing.assert_allclose(torque, np.array([0.0, 0.0, -8.0]))


def test_calc_object_forces_kernel_excludes_self_sources(tmp_path: Path) -> None:
    lib = _kernel_lib()
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    records = calc_object_forces_kernel(result, step=None, library_path=lib)

    centers = triangles.mean(axis=1)
    delta_12 = centers[0] - centers[1]
    expected_1 = K_COULOMB * charges[0] * charges[1] * delta_12 / (
        np.linalg.norm(delta_12) ** 3
    )
    assert [record.mesh_id for record in records] == [1, 2]
    np.testing.assert_allclose(records[0].force_N, expected_1, rtol=1.0e-14)
    np.testing.assert_allclose(records[1].force_N, -expected_1, rtol=1.0e-14)
    assert records[0].total_charge_C == pytest.approx(charges[0])


def test_calc_object_forces_kernel_uses_explicit_config_softening(tmp_path: Path) -> None:
    lib = _kernel_lib()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    config_path = tmp_path / "case.toml"
    config_path.write_text("[sim]\nsoftening = 0.5\ntree_leaf_max = 8\n", encoding="utf-8")
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
    result = FortranRunResult(
        directory=run_dir,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    records = calc_object_forces_kernel(
        result,
        step=None,
        library_path=lib,
        config_path=config_path,
    )

    centers = triangles.mean(axis=1)
    delta_12 = centers[0] - centers[1]
    soft2 = 0.5**2
    expected_1 = K_COULOMB * charges[0] * charges[1] * delta_12 / (
        (np.dot(delta_12, delta_12) + soft2) ** 1.5
    )
    np.testing.assert_allclose(records[0].force_N, expected_1, rtol=1.0e-14)


def test_scene_move_updates_kernel_force_geometry(tmp_path: Path) -> None:
    lib = _kernel_lib()
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    scene = BeachScene.from_result(result, step=None)
    moved = scene.move(2, by=[1.0, 0.0, 0.0])
    records = moved.calc_object_forces_kernel(
        target_mesh_ids=1,
        library_path=lib,
    )

    delta_12 = moved.centers[0] - moved.centers[1]
    expected_1 = K_COULOMB * charges[0] * charges[1] * delta_12 / (
        np.linalg.norm(delta_12) ** 3
    )
    np.testing.assert_allclose(records[0].force_N, expected_1, rtol=1.0e-14)
    np.testing.assert_allclose(scene.centers[1], triangles[1].mean(axis=0))


def test_scene_rotate_keeps_charges_attached_to_mesh_elements(tmp_path: Path) -> None:
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
            [[3.0, 0.0, 0.0], [3.0, 0.3, 0.0], [3.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9, 3.0e-9], dtype=float)
    mesh_ids = np.array([1, 2, 2], dtype=np.int64)
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=3,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=mesh_ids,
    )

    scene = BeachScene.from_result(result, step=None)
    rotated = scene.rotate(2, axis=[0.0, 0.0, 1.0], angle_deg=90.0)

    mask = mesh_ids == 2
    origin = scene.centers[mask].mean(axis=0)
    angle = np.deg2rad(90.0)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    expected_centers = (scene.centers[mask] - origin) @ rotation.T + origin
    np.testing.assert_allclose(rotated.centers[mask], expected_centers)
    np.testing.assert_allclose(rotated.charges, scene.charges)


def test_kernel_forces_cli_writes_csv(tmp_path: Path) -> None:
    lib = _kernel_lib()
    out = tmp_path / "run"
    out.mkdir()
    config_path = tmp_path / "case.toml"
    config_path.write_text("[sim]\nsoftening = 0.0\n", encoding="utf-8")
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=0",
                "last_rel_change=0.0",
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

    kernel_forces_main(
        [
            str(out),
            "--config",
            str(config_path),
            "--library",
            str(lib),
            "--save-csv",
            str(save_csv),
        ]
    )

    text = save_csv.read_text(encoding="utf-8")
    assert "mesh_id,step,total_charge_C" in text
    assert "1,-1,1e-09" in text
