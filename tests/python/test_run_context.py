from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import beach.fortran_results.potential as potential_module
from beach.fortran_results.context import RunContext, resolve_result
from beach.fortran_results.history import FortranChargeHistory
from beach.fortran_results.potential import compute_potential_mesh
from beach.fortran_results.types import FortranRunResult


def _result(directory: Path) -> FortranRunResult:
    return FortranRunResult(
        directory=directory,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([2.0e-9]),
        triangles=np.array(
            [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
        ),
        mesh_ids=np.array([4]),
    )


class _BeachLike:
    def __init__(self, result: FortranRunResult, config_path: Path | None) -> None:
        self.result = result
        self.config_path = config_path


def _write_config(path: Path, theta: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"[sim]\ntree_theta = {theta}\n", encoding="utf-8")


def test_run_context_preserves_result_protocol_and_errors(tmp_path: Path) -> None:
    result = _result(tmp_path)
    beach_like = _BeachLike(result, None)

    assert resolve_result(result) is result
    assert resolve_result(beach_like) is result
    with pytest.raises(TypeError, match="result must be FortranRunResult or Beach\\."):
        resolve_result(object())


def test_run_context_config_precedence_and_lazy_cache(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "runs" / "case" / "output"
    output_dir.mkdir(parents=True)
    auto_path = output_dir / "beach.toml"
    parent_path = output_dir.parent / "beach.toml"
    grandparent_path = output_dir.parent.parent / "beach.toml"
    inherited_path = tmp_path / "inherited.toml"
    explicit_path = tmp_path / "explicit.toml"
    _write_config(auto_path, 1.0)
    _write_config(parent_path, 1.1)
    _write_config(grandparent_path, 1.2)
    _write_config(inherited_path, 2.0)
    _write_config(explicit_path, 30.0)
    beach_like = _BeachLike(_result(output_dir), inherited_path)
    context = RunContext.from_value(beach_like, config_path=explicit_path)

    _write_config(explicit_path, 3.0)
    assert context.output_dir == output_dir
    assert context.config_path == explicit_path
    assert context.sim == {"tree_theta": 3.0}
    _write_config(explicit_path, 30.0)
    assert context.config == {"sim": {"tree_theta": 3.0}}

    inherited = RunContext.from_value(beach_like)
    assert inherited.config_path == inherited_path
    assert inherited.sim == {"tree_theta": 2.0}

    auto = RunContext.from_value(_result(output_dir))
    assert auto.config_path == auto_path
    assert auto.sim == {"tree_theta": 1.0}

    auto_path.unlink()
    parent = RunContext.from_value(_result(output_dir))
    assert parent.config_path == parent_path
    assert parent.sim == {"tree_theta": 1.1}

    parent_path.unlink()
    grandparent = RunContext.from_value(_result(output_dir))
    assert grandparent.config_path == grandparent_path
    assert grandparent.sim == {"tree_theta": 1.2}

    missing = RunContext.from_value(
        _result(output_dir),
        config_path=tmp_path / "missing.toml",
    )
    with pytest.raises(ValueError, match='config file is not found: ".*missing.toml"'):
        _ = missing.config


def test_run_context_selection_and_missing_data_contracts(tmp_path: Path) -> None:
    result = _result(tmp_path)
    context = RunContext.from_value(result)

    np.testing.assert_allclose(context.charges_at(None), np.array([2.0e-9]))
    np.testing.assert_allclose(context.charges_at(-1), np.array([2.0e-9]))
    np.testing.assert_array_equal(context.mesh_ids_or_default(), np.array([4]))
    np.testing.assert_allclose(
        context.require_triangles(),
        np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]),
    )
    with pytest.raises(ValueError, match="charge_history.csv is required"):
        context.charges_at(2)
    with pytest.raises(ValueError, match="charge_history.csv is not found or empty"):
        context.require_history()

    history = FortranChargeHistory.from_arrays(
        tmp_path / "charge_history.csv",
        mesh_nelem=1,
        history=np.array([[1.0e-9, 3.0e-9]]),
        processed_particles_by_batch=np.array([10, 20]),
        rel_change_by_batch=np.array([0.2, 0.1]),
        batch_indices=np.array([2, 5]),
    )
    historical = RunContext.from_value(replace(result, history=history))
    np.testing.assert_allclose(historical.charges_at(2), np.array([1.0e-9]))
    np.testing.assert_allclose(historical.charges_at(-1), np.array([3.0e-9]))
    np.testing.assert_array_equal(
        historical.require_history().batch_indices,
        np.array([2, 5]),
    )

    missing_mesh_data = RunContext.from_value(
        replace(result, mesh_ids=None, triangles=None)
    )
    np.testing.assert_array_equal(
        missing_mesh_data.mesh_ids_or_default(),
        np.ones(1, dtype=np.int64),
    )
    with pytest.raises(ValueError, match="mesh_triangles.csv is not found"):
        missing_mesh_data.require_triangles()


def test_potential_call_inherits_config_for_kernel_and_reference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    config_path = tmp_path / "external.toml"
    config_path.write_text(
        """
[sim]
field_bc_mode = "free"
tree_theta = 0.4

[particles]
[[particles.species]]
inject_face = "x_low"
pos_low = [0.0, 0.0, 0.0]
pos_high = [0.0, 2.0, 2.0]
""".lstrip(),
        encoding="utf-8",
    )
    beach_like = _BeachLike(_result(output_dir), config_path)
    theta_seen: list[float] = []

    class FakeKernel:
        def __init__(self, _triangles, _charges, **kwargs) -> None:
            theta_seen.append(float(kwargs["options"].theta))

        def __enter__(self) -> FakeKernel:
            return self

        def __exit__(self, *_args: object) -> None:
            pass

        def eval_phi(self, points: np.ndarray) -> np.ndarray:
            samples = np.asarray(points, dtype=float)
            return samples @ np.array([1.0, 2.0, 3.0])

    monkeypatch.setattr(potential_module, "FieldKernel", FakeKernel)

    potential = compute_potential_mesh(
        beach_like,
        reference_point="species1_injection_center",
    )

    np.testing.assert_allclose(potential, np.array([-4.0]))
    assert theta_seen
    assert theta_seen[-1] == pytest.approx(0.4)
