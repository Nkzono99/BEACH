
from beach.io.types import RunField, RunMesh, RunSummary
from pathlib import Path

import numpy as np
import pytest

import beach.analysis.object_interaction as interaction_module
from beach import (
    FieldKernelError,
    RunResult,
    ObjectInteractionSnapshot,
)


def _triangles_at(positions: np.ndarray) -> np.ndarray:
    offsets = np.array(
        [[-0.02, -0.02, 0.0], [0.04, -0.02, 0.0], [-0.02, 0.04, 0.0]]
    )
    return np.asarray(positions, dtype=float)[:, None, :] + offsets[None, :, :]


def _result(
    directory: Path,
    *,
    triangles: np.ndarray | None = None,
    source_model: str = "triangle_p0",
) -> RunResult:
    if triangles is None:
        triangles = _triangles_at(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0]]))
    return RunResult(
        directory=directory,
        charges=np.array([1.0, 2.0]),
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
        mesh=RunMesh(
            triangles=np.asarray(triangles, dtype=float), mesh_ids=np.array([1, 2])
        ),
        field=RunField(source_model=source_model),
    )


def _write_config(path: Path, *, box_max_z: float = 4.0) -> None:
    path.write_text(
        f"""
[sim]
field_bc_mode = "free"
box_min = [-2.0, -2.0, -2.0]
box_max = [2.0, 2.0, {box_max_z}]
e0 = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )


_PATH_L = 1.0e-3
_PATH_F0 = 1.0e-12


class _AnalyticPathKernel:
    mode = "inverse_square"
    eval_calls = 0

    def __init__(self, _positions, charges, **_kwargs) -> None:
        self.current = np.array(charges, copy=True)

    def update_charges(self, charges: np.ndarray) -> None:
        self.current = np.array(charges, copy=True)

    @classmethod
    def _field_potential(cls, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        z = points[:, 2]
        field = np.zeros_like(points)
        if cls.mode == "inverse_square":
            field[:, 2] = 1.0 / (1.0 + z) ** 2
            potential = 1.0 / (1.0 + z)
        elif cls.mode == "small_nonlinear_crossing":
            scale = 1.0e-15
            field[:, 2] = scale * (z * z + 0.6 * z - 0.4)
            potential = -scale * (z**3 / 3.0 + 0.3 * z * z - 0.4 * z)
        elif cls.mode == "sharp":
            field[:, 2] = 1.0 / (0.01 + z) ** 2
            potential = 1.0 / (0.01 + z)
        elif cls.mode == "inconsistent_constant":
            field[:, 2] = 1.0
            potential = np.zeros_like(z)
        elif cls.mode in ("quadratic_conservative", "quadratic_nonconservative"):
            # F = F0 (1 + 0.04 s^2) with s = z / L on a unit target charge; the
            # conservative U = -F0 L (s + 0.04 s^3 / 3) gives F = -dU/dz exactly.
            s = z / _PATH_L
            cubic = 0.04 if cls.mode == "quadratic_conservative" else 0.2
            field[:, 2] = _PATH_F0 * (1.0 + 0.04 * s * s)
            potential = -_PATH_F0 * _PATH_L * (s + cubic * s**3 / 3.0)
        else:
            raise AssertionError("unknown analytic path mode")
        return field, potential

    def eval_e(self, points: np.ndarray) -> np.ndarray:
        type(self).eval_calls += 1
        if float(np.sum(self.current)) != 2.0:
            return np.zeros_like(points)
        return self._field_potential(points)[0]

    def eval_phi(self, points: np.ndarray) -> np.ndarray:
        type(self).eval_calls += 1
        if float(np.sum(self.current)) != 2.0:
            return np.zeros(len(points))
        return self._field_potential(points)[1]

    def eval_e_direct(self, points: np.ndarray) -> np.ndarray:
        type(self).eval_calls += 1
        return np.zeros_like(points)

    def eval_phi_direct(self, points: np.ndarray) -> np.ndarray:
        type(self).eval_calls += 1
        return np.zeros(len(points))

    def diagnostics(self):
        raise FieldKernelError("not cached")

    def close(self) -> None:
        pass


@pytest.fixture
def analytic_kernel(monkeypatch: pytest.MonkeyPatch):
    _AnalyticPathKernel.eval_calls = 0
    _AnalyticPathKernel.mode = "inverse_square"
    monkeypatch.setattr(interaction_module, "FieldKernel", _AnalyticPathKernel)
    return _AnalyticPathKernel


def test_vertical_path_keeps_sources_frozen_and_matches_analytic_force_and_work(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path), step=None, config_path=config
    ) as snapshot:
        probe = snapshot.object_probe(1)
        source_positions = np.array(snapshot.source_positions_m, copy=True)
        source_charges = np.array(snapshot.source_charges_C, copy=True)
        initial = probe.wrench()
        path = probe.vertical_path(np.linspace(0.0, 1.0, 65), adaptive=True)
        np.testing.assert_array_equal(snapshot.source_positions_m, source_positions)
        np.testing.assert_array_equal(snapshot.source_charges_C, source_charges)

    expected_force_z = 1.0 / (1.0 + path.displacement_m) ** 2
    expected_potential = 1.0 / (1.0 + path.displacement_m)
    expected_work = path.displacement_m / (1.0 + path.displacement_m)
    assert path.potential_energy_J is not None
    assert path.potential_difference_work_J is not None
    np.testing.assert_allclose(
        path.force_N,
        np.column_stack(
            (
                np.zeros(path.displacement_m.size),
                np.zeros(path.displacement_m.size),
                expected_force_z,
            )
        )
    )
    np.testing.assert_allclose(path.potential_energy_J, expected_potential)
    np.testing.assert_allclose(path.potential_difference_work_J, expected_work)

    np.testing.assert_allclose(path.force_N[0], initial.force_N)
    np.testing.assert_allclose(path.torque_Nm[0], initial.torque_Nm)
    np.testing.assert_allclose(
        path.electrostatic_work_J,
        expected_work,
        rtol=5.0e-3,
        atol=1.0e-18,
    )
    assert path.status == "converged"


def test_vertical_path_records_moving_and_fixed_torque_origins(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    displacement = np.array([0.0, 0.25, 0.5])
    explicit_origin = np.array([0.3, -0.2, 0.1])

    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path), step=None, config_path=config
    ) as snapshot:
        probe = snapshot.object_probe(1)
        centroid = np.array(probe.geometric_area_centroid_m, copy=True)
        moving = probe.vertical_path(displacement, adaptive=False)
        fixed_origin = probe.vertical_path(
            displacement,
            adaptive=False,
            torque_origin="origin",
        )
        fixed_explicit = probe.vertical_path(
            displacement,
            adaptive=False,
            torque_origin=explicit_origin,
        )

    moving_origins = np.asarray(moving.numerical_metadata["torque_origin_m"])
    expected_moving = np.broadcast_to(centroid, (displacement.size, 3)).copy()
    expected_moving[:, 2] += displacement
    assert moving.numerical_metadata["torque_origin_policy"] == (
        "moving_geometric_area_centroid"
    )
    np.testing.assert_allclose(moving_origins, expected_moving)
    assert fixed_origin.numerical_metadata["torque_origin_policy"] == "fixed_origin"
    np.testing.assert_allclose(
        fixed_origin.numerical_metadata["torque_origin_m"],
        np.zeros((displacement.size, 3)),
    )
    assert fixed_explicit.numerical_metadata["torque_origin_policy"] == (
        "fixed_explicit"
    )
    np.testing.assert_allclose(
        fixed_explicit.numerical_metadata["torque_origin_m"],
        np.broadcast_to(explicit_origin, (displacement.size, 3)),
    )


def test_vertical_path_converges_after_refinement_and_reports_limit(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    analytic_kernel.mode = "sharp"
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path), step=None, config_path=config
    ) as snapshot:
        probe = snapshot.object_probe(1)
        refined = probe.vertical_path(
            np.array([0.0, 1.0]),
            relative_tolerance=1.0e-3,
            force_absolute_tolerance_N=0.0,
            work_absolute_tolerance_J=0.0,
            max_refinement=12,
        )
        failed = probe.vertical_path(
            np.array([0.0, 1.0]),
            relative_tolerance=1.0e-12,
            force_absolute_tolerance_N=0.0,
            work_absolute_tolerance_J=0.0,
            max_refinement=0,
        )

    assert refined.displacement_m.size > 2
    assert refined.refinement_count > 0
    assert refined.status == "converged"
    # The refined trapezoid work must reach the closed form phi(0) - phi(z) of the sharp field to the
    # requested relative tolerance; the unrefined two-point trapezoid is about 50 times too large.
    z = refined.displacement_m
    np.testing.assert_allclose(refined.force_N[:, 2], 1.0 / (0.01 + z) ** 2, rtol=1.0e-12)
    np.testing.assert_allclose(
        refined.electrostatic_work_J[1:],
        1.0 / 0.01 - 1.0 / (0.01 + z[1:]),
        rtol=1.0e-3,
    )
    assert failed.status == "not_converged"
    assert failed.numerical_metadata["status_reason"] == "max_refinement_reached"


def test_vertical_path_zero_crossing_uses_absolute_tolerance(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    analytic_kernel.mode = "small_nonlinear_crossing"
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path), step=None, config_path=config
    ) as snapshot:
        path = snapshot.object_probe(1).vertical_path(
            np.array([0.0, 1.0]),
            relative_tolerance=0.0,
            force_absolute_tolerance_N=1.0e-14,
            work_absolute_tolerance_J=1.0e-14,
            max_refinement=1,
        )

    assert path.status == "converged"
    assert path.refinement_count == 0
    assert path.force_N[0, 2] < 0.0 < path.force_N[-1, 2]
    assert path.numerical_metadata["max_force_refinement_error_N"] > 0.0
    assert path.numerical_metadata["max_work_refinement_error_J"] > 0.0


def test_vertical_path_reports_potential_defect_separately_from_quadrature_error(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config)
    analytic_kernel.mode = "inconsistent_constant"
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path), step=None, config_path=config
    ) as snapshot:
        path = snapshot.object_probe(1).vertical_path(
            np.array([0.0, 1.0]),
            relative_tolerance=0.0,
            force_absolute_tolerance_N=0.0,
            work_absolute_tolerance_J=0.0,
            max_refinement=0,
        )

    assert path.status == "not_converged"
    assert path.numerical_metadata["status_reason"] == "work_potential_mismatch"
    assert path.numerical_metadata["max_work_refinement_error_J"] == 0.0
    assert path.numerical_metadata["max_work_potential_error_J"] > 0.0


def _quadratic_path(
    tmp_path: Path,
    analytic_kernel,
    mode: str,
    max_refinement: int,
    grid_in_L: tuple[float, ...] = (0.0, 1.0),
):
    config = tmp_path / "beach.toml"
    _write_config(config)
    analytic_kernel.mode = mode
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path), step=None, config_path=config
    ) as snapshot:
        return snapshot.object_probe(1).vertical_path(
            _PATH_L * np.array(grid_in_L), max_refinement=max_refinement
        )


def test_vertical_path_refines_until_the_returned_work_matches_the_potential(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    path = _quadratic_path(tmp_path, analytic_kernel, "quadratic_conservative", 8)

    # On [0, L] the trapezoid work F0 L (1 + 0.02) misses the potential difference
    # F0 L (1 + 0.04 / 3) by 6.7e-18 J, above the 6.1e-18 J allowance, while each
    # half interval of the midpoint test passes; one midpoint brings it to 1.7e-18 J.
    assert path.status == "converged"
    assert path.numerical_metadata["status_reason"] == "tolerances_satisfied"
    np.testing.assert_allclose(path.displacement_m, [0.0, 0.5 * _PATH_L, _PATH_L])
    assert path.refinement_count == 1
    assert path.numerical_metadata["inserted_point_count"] == 1
    assert path.numerical_metadata["refinement_rounds"] == 1

    s = path.displacement_m / _PATH_L
    exact_work = _PATH_F0 * _PATH_L * (s + 0.04 * s**3 / 3.0)
    force_z = path.force_N[:, 2]
    recomputed = np.concatenate(
        ([0.0], np.cumsum(0.5 * np.diff(path.displacement_m) * (force_z[1:] + force_z[:-1])))
    )
    np.testing.assert_allclose(path.electrostatic_work_J, recomputed, rtol=1.0e-14, atol=0.0)
    np.testing.assert_allclose(path.potential_difference_work_J, exact_work, rtol=1.0e-12)
    mismatch = float(np.max(np.abs(recomputed - exact_work)))
    assert path.work_absolute_mismatch_J == pytest.approx(mismatch, rel=1.0e-6)
    assert mismatch <= 1.0e-18 + 5.0e-3 * float(np.max(np.abs(exact_work)))


def test_vertical_path_refines_intervals_whose_errors_add_up_along_the_path(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    path = _quadratic_path(
        tmp_path, analytic_kernel, "quadratic_conservative", 8, (0.0, 1.0, 2.0)
    )

    # Each unit interval misses its potential difference by F0 L 0.08 / 12 = 6.7e-18 J
    # and passes the local tests, but the two add to 1.33e-17 J against an allowance
    # of 1.15e-17 J on [0, 2L]; halving both brings the sum to 3.3e-18 J.
    assert path.status == "converged"
    np.testing.assert_allclose(path.displacement_m / _PATH_L, [0.0, 0.5, 1.0, 1.5, 2.0])
    assert path.refinement_count == 2
    assert path.numerical_metadata["refinement_rounds"] == 1


def test_vertical_path_without_refinement_budget_reports_the_mismatch(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    path = _quadratic_path(tmp_path, analytic_kernel, "quadratic_conservative", 0)

    assert path.status == "not_converged"
    assert path.numerical_metadata["status_reason"] == "work_potential_mismatch"
    np.testing.assert_array_equal(path.displacement_m, [0.0, _PATH_L])
    assert path.refinement_count == 0
    assert path.numerical_metadata["refinement_rounds"] == 0
    assert path.work_absolute_mismatch_J == pytest.approx(
        _PATH_F0 * _PATH_L * (0.02 - 0.04 / 3.0), rel=1.0e-6
    )


def test_vertical_path_keeps_a_nonconservative_potential_not_converged(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    path = _quadratic_path(tmp_path, analytic_kernel, "quadratic_nonconservative", 8)

    assert path.status == "not_converged"
    assert path.numerical_metadata["status_reason"] == "work_potential_mismatch"
    assert path.numerical_metadata["refinement_rounds"] == 8
    assert path.refinement_count == path.displacement_m.size - 2
    # The force integrates to F0 L (1 + 0.04 / 3) but U differs by F0 L (0.2 - 0.04) / 3.
    assert path.work_absolute_mismatch_J == pytest.approx(
        _PATH_F0 * _PATH_L * 0.16 / 3.0, rel=1.0e-3
    )


def test_vertical_path_checks_triangle_vertices_before_native_evaluation(
    tmp_path: Path,
    analytic_kernel,
) -> None:
    config = tmp_path / "beach.toml"
    _write_config(config, box_max_z=1.0)
    triangles = np.array(
        [
            [[0.0, 0.0, 0.99], [0.2, 0.0, 0.0], [0.0, 0.2, 0.0]],
            [[0.5, 0.0, 0.0], [0.7, 0.0, 0.0], [0.5, 0.2, 0.0]],
        ]
    )
    with ObjectInteractionSnapshot.from_result(
        _result(tmp_path, triangles=triangles, source_model="triangle_p0"),
        step=None,
        config_path=config,
    ) as snapshot:
        probe = snapshot.object_probe(1)
        before = analytic_kernel.eval_calls
        with pytest.raises(ValueError, match="box"):
            probe.vertical_path(np.array([0.0, 0.02]))
        assert analytic_kernel.eval_calls == before
