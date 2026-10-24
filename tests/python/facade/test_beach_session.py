"""Session contract of ``Beach``: one context and one kernel cache per run."""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

import beach
import beach.analysis.field_kernel as field_kernel_module
import beach.analysis.potential as potential_module
from beach import Beach, BeachScene, ObjectInteractionSnapshot


class _RecordingKernel:
    built: list["_RecordingKernel"] = []

    def __init__(self, triangles, charges, *, options, **_kwargs) -> None:
        self.options = options
        self.charges = np.asarray(charges, dtype=float).copy()
        self.closed = False
        _RecordingKernel.built.append(self)

    def close(self) -> None:
        self.closed = True

    def __enter__(self) -> "_RecordingKernel":
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    @property
    def source_charges(self) -> np.ndarray:
        return self.charges.copy()

    def update_charges(self, charges: np.ndarray) -> None:
        self.charges = np.asarray(charges, dtype=float).copy()

    def eval_phi(self, points: np.ndarray) -> np.ndarray:
        return np.zeros(np.asarray(points).shape[0])

    def eval_e(self, points: np.ndarray) -> np.ndarray:
        return np.zeros((np.asarray(points).shape[0], 3))

    def force_on_charges(self, points, charges, *, origin):
        return np.zeros(3), np.zeros(3)


def _write_run(out: Path) -> Path:
    out.mkdir(parents=True)
    (out / "summary.txt").write_text(
        "mesh_nelem=1\nprocessed_particles=1\nabsorbed=1\nescaped=0\n"
        "batches=2\nlast_rel_change=0.0\nfield_source_model=triangle_p0\n",
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,2.0e-9\n", encoding="utf-8")
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0,0,0,1,0,0,0,1,0,2.0e-9,1\n",
        encoding="utf-8",
    )
    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        "1,1,0.0,1,1.0e-9\n"
        "2,1,0.0,1,2.0e-9\n",
        encoding="utf-8",
    )
    (out / "mesh_potential.csv").write_text(
        "elem_idx,potential_V\n1,5.0\n", encoding="utf-8"
    )
    (out / "beach.toml").write_text(
        '[sim]\ntree_theta = 0.3\nfield_bc_mode = "free"\n', encoding="utf-8"
    )
    return out


@pytest.fixture(autouse=True)
def _recording_kernel(monkeypatch: pytest.MonkeyPatch) -> None:
    _RecordingKernel.built = []
    monkeypatch.setattr(field_kernel_module, "FieldKernel", _RecordingKernel)
    monkeypatch.setattr(potential_module, "FieldKernel", _RecordingKernel)


def test_kernel_is_cached_per_plan_and_closed_with_the_session(tmp_path: Path) -> None:
    out = _write_run(tmp_path / "run")

    with Beach(out) as run:
        latest = run.kernel()
        assert run.kernel(step=-1) is latest
        # The latest history snapshot equals charges.csv, so both share one plan.
        assert run.kernel(step=None) is latest
        assert run.kernel(step=1) is not latest
        assert run.kernel(theta=0.7) is not latest
        assert latest.options.theta == pytest.approx(0.3)

    assert len(_RecordingKernel.built) == 3
    assert all(kernel.closed for kernel in _RecordingKernel.built)


def test_reload_releases_kernels_and_rebuilds_on_demand(tmp_path: Path) -> None:
    out = _write_run(tmp_path / "run")
    run = Beach(out)
    first = run.kernel()

    run.reload()

    assert first.closed
    assert run.kernel() is not first


def test_context_owns_the_session_config(tmp_path: Path) -> None:
    out = _write_run(tmp_path / "run")
    explicit = tmp_path / "explicit.toml"
    explicit.write_text("[sim]\ntree_theta = 0.9\n", encoding="utf-8")

    run = Beach(out, config_path=explicit)

    assert run.context.config_path == explicit
    assert run.context.result is run.result
    assert run.kernel().options.theta == pytest.approx(0.9)


def test_methods_share_the_session_kernel_and_keep_its_charges(tmp_path: Path) -> None:
    out = _write_run(tmp_path / "run")

    with Beach(out) as run:
        np.testing.assert_allclose(run.compute_potential(), np.array([5.0]))
        assert _RecordingKernel.built == []
        run.compute_potential(source="kernel")
        run.compute_electric_field(np.array([[0.0, 0.0, 1.0]]))
        run.calc_object_forces_kernel()
        run.compute_potential_points(np.array([[0.0, 0.0, 1.0]]), library_path=tmp_path)

    assert len(_RecordingKernel.built) == 2
    shared, own_library = _RecordingKernel.built
    np.testing.assert_allclose(shared.charges, np.array([2.0e-9]))
    assert shared.closed and own_library.closed


_FORWARDED = {
    "calc_coulomb": beach.calc_coulomb,
    "analyze_coulomb_mobility": beach.analyze_coulomb_mobility,
    "object_interaction_snapshot": ObjectInteractionSnapshot.from_result,
    "calc_object_forces_kernel": beach.calc_object_forces_kernel,
    "scene": BeachScene.from_result,
    "plot_mesh": beach.plot_charge_mesh,
    "plot_bar": beach.plot_charges,
    "compute_potential": beach.compute_potential_mesh,
    "compute_potential_points": beach.compute_potential_points,
    "compute_potential_slices": beach.compute_potential_slices,
    "plot_potential": beach.plot_potential_mesh,
    "plot_potential_slices": beach.plot_potential_slices,
    "plot_mesh_source_boxplot": beach.plot_mesh_source_boxplot,
    "plot_coulomb_force_matrix": beach.plot_coulomb_force_matrix,
    "animate_mesh": beach.animate_history_mesh,
    "compute_electric_field": beach.compute_electric_field_points,
    "trace_field_lines": beach.trace_field_lines,
    "plot_field_lines": beach.plot_field_lines_3d,
}


@pytest.mark.parametrize("name", sorted(_FORWARDED))
def test_methods_show_the_module_signature_and_docstring(name: str) -> None:
    function = _FORWARDED[name]
    method = getattr(Beach, name)

    function_parameters = list(inspect.signature(function).parameters.values())
    method_parameters = list(inspect.signature(method).parameters.values())
    assert method_parameters[0].name == "self"
    assert method_parameters[1:] == function_parameters[1:]
    assert method.__name__ == name
    assert inspect.getdoc(method).startswith(inspect.getdoc(function))


def test_forwarded_analysis_uses_the_requested_history_step(tmp_path: Path) -> None:
    out = _write_run(tmp_path / "run")
    with Beach(out) as run:
        records = run.calc_object_forces_kernel(step=1)
        assert records[0].total_charge_C == 1.0e-9
        history_kernel = run.kernel(step=1)
        assert _RecordingKernel.built == [history_kernel]
        np.testing.assert_array_equal(history_kernel.source_charges, [1.0e-9])


def test_session_reference_iterator_is_consumed_once(tmp_path: Path) -> None:
    out = _write_run(tmp_path / "run")
    with Beach(out) as run:
        actual = run.compute_potential(reference_point=(v for v in (0.0, 0.0, 2.0)))
        np.testing.assert_array_equal(actual, [0.0])
        with pytest.raises(ValueError, match='source="saved" cannot be combined'):
            run.compute_potential(source="saved", reference_point=(v for v in (0.0, 0.0, 2.0)))


def test_equivalent_resolved_options_and_history_charges_share_a_kernel(tmp_path: Path) -> None:
    out = _write_run(tmp_path / "run")
    # Final charges differ from the latest history snapshot.
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,3.0e-9\n", encoding="utf-8")
    with Beach(out) as run:
        latest = run.kernel()
        assert run.kernel(step=2, theta=0.3, leaf_max=16, order=4) is latest
        assert run.kernel(step=None) is not latest


def test_public_names_still_resolve() -> None:
    assert beach.Beach.__name__ == "Beach"
    assert "load_run" in dir(beach)
    with pytest.raises(AttributeError):
        beach.no_such_name  # noqa: B018
