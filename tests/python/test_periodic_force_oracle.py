from pathlib import Path

import numpy as np
import pytest

from beach import (
    FieldKernel,
    FieldKernelOptions,
    FortranRunResult,
    ObjectInteractionSnapshot,
    finite_shell_convergence,
    finite_shell_wrench,
)
from beach.fortran_results.constants import K_COULOMB


EPS0 = 1.0 / (4.0 * np.pi * K_COULOMB)


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


def _result(directory: Path, charges: np.ndarray) -> FortranRunResult:
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
        field_source_model="point",
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
softening = 0.05
tree_theta = 0.1
tree_leaf_max = 64
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
field_periodic_ewald_layers = 4
e0 = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
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
        np.asarray(replicated_pos),
        np.asarray(replicated_q),
        options=FieldKernelOptions(softening=0.05, leaf_max=64),
        library_path=_kernel_lib(),
    ) as direct:
        expected_field = direct.eval_e_direct(base[[0]])[0]

    np.testing.assert_allclose(
        shell.symmetric.force_N,
        charges[0] * expected_field,
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
