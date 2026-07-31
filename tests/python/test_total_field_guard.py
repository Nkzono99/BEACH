from pathlib import Path

import numpy as np
import pytest

from beach.fortran_results.context import RunContext
from beach.fortran_results.kernel import (
    FieldKernelOptions,
    _require_total_field_reconstruction,
)
from beach.fortran_results.periodic import Periodic2Config
from beach.fortran_results.potential import _potential_history
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
        charges=np.array([0.0]),
        triangles=np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]),
    )


def test_total_field_guard_requires_config(tmp_path: Path) -> None:
    context = RunContext.from_value(_result(tmp_path))

    with pytest.raises(ValueError, match="requires the run's beach.toml"):
        _require_total_field_reconstruction(
            context,
            FieldKernelOptions(),
            operation="test operation",
        )


def test_total_field_guard_rejects_cached_kneq0_component(tmp_path: Path) -> None:
    (tmp_path / "beach.toml").write_text(
        '[external_boundary.field]\nmodel = "none"\n',
        encoding="utf-8",
    )
    context = RunContext.from_value(_result(tmp_path))
    periodic2 = Periodic2Config(
        axes=(0, 1),
        lengths=(1.0, 1.0),
        origins=(0.0, 0.0),
        image_layers=1,
        far_correction="cached_kneq0",
        ewald_alpha=0.0,
        ewald_layers=4,
    )

    with pytest.raises(ValueError, match="only the k!=0 component"):
        _require_total_field_reconstruction(
            context,
            FieldKernelOptions(periodic2=periodic2),
            operation="test operation",
        )


def test_total_field_guard_accepts_configured_complete_kernel(tmp_path: Path) -> None:
    (tmp_path / "beach.toml").write_text(
        '[external_boundary.field]\nmodel = "none"\n',
        encoding="utf-8",
    )
    context = RunContext.from_value(_result(tmp_path))

    _require_total_field_reconstruction(
        context,
        FieldKernelOptions(),
        operation="test operation",
    )


def test_total_field_guard_recommends_none_for_explicit_historical_oracle(
    tmp_path: Path,
) -> None:
    config = tmp_path / "historical.toml"
    config.write_text(
        '[sim]\nfield_periodic_far_correction = "m2l_root_oracle"\n',
        encoding="utf-8",
    )
    context = RunContext.from_value(_result(tmp_path), config_path=config)

    with pytest.raises(ValueError, match='removed; use "none"'):
        _require_total_field_reconstruction(
            context,
            FieldKernelOptions(),
            operation="test operation",
        )


def test_potential_history_rejects_cached_component_without_loading_kernel() -> None:
    periodic2 = Periodic2Config(
        axes=(0, 1),
        lengths=(1.0, 1.0),
        origins=(0.0, 0.0),
        image_layers=1,
        far_correction="cached_kneq0",
        ewald_alpha=0.0,
        ewald_layers=4,
    )

    with pytest.raises(ValueError, match="only the k!=0 component"):
        _potential_history(
            np.array([[1.0]]),
            np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]),
            field_options=FieldKernelOptions(periodic2=periodic2),
        )
