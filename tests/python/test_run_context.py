from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest

import beach.fortran_results.context as context_module
from beach.fortran_results.context import RunContext, resolve_result
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
        triangles=np.zeros((1, 3, 3)),
        mesh_ids=np.array([4]),
    )


class _BeachLike:
    def __init__(self, result: FortranRunResult, config_path: Path | None) -> None:
        self.result = result
        self.config_path = config_path


def _write_config(path: Path, softening: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"[sim]\nsoftening = {softening}\n", encoding="utf-8")


def test_run_context_preserves_result_protocol_and_errors(tmp_path: Path) -> None:
    result = _result(tmp_path)
    beach_like = _BeachLike(result, None)

    assert resolve_result(result) is result
    assert resolve_result(beach_like) is result
    with pytest.raises(TypeError, match="result must be FortranRunResult or Beach\\."):
        resolve_result(object())


def test_run_context_config_precedence_and_lazy_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "runs" / "case" / "output"
    output_dir.mkdir(parents=True)
    auto_path = output_dir / "beach.toml"
    inherited_path = tmp_path / "inherited.toml"
    explicit_path = tmp_path / "explicit.toml"
    _write_config(auto_path, 1.0)
    _write_config(inherited_path, 2.0)
    _write_config(explicit_path, 3.0)

    calls: list[Path] = []
    original_loader = context_module.load_toml

    def counting_loader(path: Path) -> dict[str, object]:
        calls.append(path)
        return original_loader(path)

    monkeypatch.setattr(context_module, "load_toml", counting_loader)
    beach_like = _BeachLike(_result(output_dir), inherited_path)
    context = RunContext.from_value(beach_like, config_path=explicit_path)

    assert calls == []
    assert context.output_dir == output_dir
    assert context.config_path == explicit_path
    assert context.sim == {"softening": 3.0}
    assert context.config is context.config
    assert context.sim is context.sim
    assert calls == [explicit_path]

    inherited = RunContext.from_value(beach_like)
    assert inherited.config_path == inherited_path
    assert inherited.sim == {"softening": 2.0}

    auto = RunContext.from_value(_result(output_dir))
    assert auto.config_path == auto_path
    assert auto.sim == {"softening": 1.0}


def test_run_context_missing_config_and_selection_contracts(tmp_path: Path) -> None:
    result = _result(tmp_path)
    context = RunContext.from_value(result, config_path=tmp_path / "missing.toml")

    with pytest.raises(ValueError, match='config file is not found: ".*missing.toml"'):
        _ = context.config
    assert np.array_equal(context.charges_at(None), result.charges)
    assert np.array_equal(context.charges_at(-1), result.charges)
    assert np.array_equal(context.mesh_ids_or_default(), np.array([4]))
    assert context.require_triangles() is result.triangles
    with pytest.raises(FrozenInstanceError):
        context.result = result  # type: ignore[misc]


def test_potential_call_inherits_beach_config_and_parses_it_once(
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
softening = 2.0

[particles]
[[particles.species]]
inject_face = "x_low"
pos_low = [0.0, 0.0, 0.0]
pos_high = [0.0, 2.0, 2.0]
""".lstrip(),
        encoding="utf-8",
    )
    beach_like = _BeachLike(_result(output_dir), config_path)
    calls: list[Path] = []
    original_loader = context_module.load_toml

    def counting_loader(path: Path) -> dict[str, object]:
        calls.append(path)
        return original_loader(path)

    monkeypatch.setattr(context_module, "load_toml", counting_loader)

    potential = compute_potential_mesh(
        beach_like,
        reference_point="species1_injection_center",
    )

    assert potential.shape == (1,)
    assert np.all(np.isfinite(potential))
    assert calls == [config_path]
