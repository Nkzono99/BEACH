"""examples/analyze_periodic_object_detachment.py: input validation and shell summary."""

from __future__ import annotations

import sys

import numpy as np
import pytest

import examples.analyze_periodic_object_detachment as example_module
from beach import (
    FiniteShellConvergenceResult,
    ObjectForcePath,
)


def _simple_path(force_z: float) -> ObjectForcePath:
    h = np.array([0.0, 1.0])
    force = np.array([[0.0, 0.0, force_z], [0.0, 0.0, force_z]])
    torque = np.zeros((2, 3))
    potential = np.array([force_z, 0.0])
    return ObjectForcePath.from_samples(h, force, torque, potential)


def test_example_shell_summary_names_closure_raw_corrected_and_errors() -> None:
    symmetric = tuple(_simple_path(value) for value in (1.0, 1.1, 1.11))
    corrected = tuple(_simple_path(value) for value in (2.0, 2.1, 2.11))
    shells = FiniteShellConvergenceResult(
        image_layers=np.array([0, 1, 2]),
        symmetric_paths=symmetric,
        corrected_paths=corrected,
        force_increment_error_N=np.array([0.1, 0.01]),
        work_increment_error_J=np.array([0.1, 0.01]),
        increment_converged=np.array([True, True]),
        status="converged",
        selected_image_layers=2,
        selected_path=corrected[-1],
        reference_force_error_N=np.array([0.2, 0.1, 0.01]),
        reference_work_error_J=np.array([0.2, 0.1, 0.01]),
        reference_converged=np.array([False, True, True]),
        reference_model="infinite_physical",
    )

    summary = example_module._finite_shell_summary(shells)

    assert summary["selected_closure"] == "e_bottom_zero"
    assert summary["raw_symmetric"][0]["image_layers"] == 0
    assert summary["corrected_e_bottom_zero"][-1]["image_layers"] == 2
    assert summary["force_increment_error_N"] == [0.1, 0.01]
    assert summary["work_increment_error_J"] == [0.1, 0.01]
    assert summary["force_tail_proxy_N"] == [0.1, 0.02]
    assert summary["work_tail_proxy_J"] == [0.1, 0.02]
    assert summary["reference_model"] == "infinite_physical"
    assert summary["reference_converged"] == [False, True, True]


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        pytest.param(
            ["--eta-translation", "1.1"],
            "eta_translation must be finite and between 0 and 1",
            id="eta-above-one",
        ),
        pytest.param(
            ["--gravity-m-s2", "-1"],
            "gravity_m_s2 must be finite and non-negative",
            id="negative-gravity",
        ),
        pytest.param(
            ["--shell-max-layers", "-1"],
            "shell_max_layers must be non-negative",
            id="negative-shell-count",
        ),
    ],
)
def test_example_validates_inputs_before_constructing_beach(
    monkeypatch: pytest.MonkeyPatch,
    extra: list[str],
    message: str,
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "analyze-periodic-object-detachment",
            "unused-output",
            "--config",
            "unused.toml",
            "--target-mesh-id",
            "1",
            "--z-max-m",
            "1",
            "--z-points",
            "3",
            "--mass-kg",
            "1",
            *extra,
        ],
    )
    monkeypatch.setattr(
        example_module,
        "Beach",
        lambda *_args, **_kwargs: pytest.fail("Beach constructed before validation"),
    )

    with pytest.raises(ValueError, match=message):
        example_module.main()
