"""Session contract of ``Beach``: one context and one kernel cache per run."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import beach.analysis.field_kernel as field_kernel_module
from beach import Beach


class _RecordingKernel:
    built: list["_RecordingKernel"] = []

    def __init__(self, triangles, charges, *, options, **_kwargs) -> None:
        self.options = options
        self.charges = np.asarray(charges, dtype=float).copy()
        self.closed = False
        _RecordingKernel.built.append(self)

    def close(self) -> None:
        self.closed = True


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
    (out / "beach.toml").write_text("[sim]\ntree_theta = 0.3\n", encoding="utf-8")
    return out


@pytest.fixture(autouse=True)
def _recording_kernel(monkeypatch: pytest.MonkeyPatch) -> None:
    _RecordingKernel.built = []
    monkeypatch.setattr(field_kernel_module, "FieldKernel", _RecordingKernel)


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
