from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import beach.cli.object_detachment as cli_module
import examples.analyze_periodic_object_detachment as example_module
from beach import (
    FiniteShellConvergenceResult,
    ObjectForcePath,
    ObjectWrench,
    RigidTransform,
    WrenchComponent,
)
from beach.cli.main import main as beachx_main


_ARTIFACT_NAMES = (
    "instantaneous_wrench.csv",
    "path.csv",
    "summary.json",
    "report.md",
)
_PATH_COMPONENTS = (
    "other_objects_all_images",
    "target_periodic_images",
    "external_uniform",
    "total_external",
)


def _required_args(tmp_path: Path) -> list[str]:
    return [
        "object-detachment",
        str(tmp_path / "run"),
        "--config",
        str(tmp_path / "beach.toml"),
        "--target-mesh-id",
        "6",
        "--periodic-model",
        "configured",
        "--z-max-m",
        "1.0",
        "--z-points",
        "3",
        "--mass-kg",
        "1.0",
    ]


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        (["--mass-kg", "0"], "--mass-kg must be > 0"),
        (["--z-points", "1"], "--z-points must be >= 2"),
        (["--z-max-m", "0"], "--z-max-m must be > 0"),
        (["--gravity-m-s2", "-1"], "--gravity-m-s2 must be >= 0"),
        (["--eta-translation", "1.1"], "--eta-translation must be between 0 and 1"),
    ],
)
def test_object_detachment_cli_validates_numeric_contracts(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    extra: list[str],
    message: str,
) -> None:
    args = _required_args(tmp_path)
    for index in range(0, len(extra), 2):
        option = extra[index]
        value = extra[index + 1]
        if option in args:
            args[args.index(option) + 1] = value
        else:
            args.extend([option, value])

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(args)

    assert excinfo.value.code == 2
    assert message in capsys.readouterr().err


def test_object_detachment_cli_requires_complete_adhesion_model(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    args = _required_args(tmp_path) + ["--adhesion-force-n", "1.0"]

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(args)

    assert excinfo.value.code == 2
    assert "must be supplied together" in capsys.readouterr().err


def test_object_detachment_cli_rejects_invalid_torque_origin(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    args = _required_args(tmp_path) + ["--torque-origin", "0,not-a-number,1"]

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(args)

    assert excinfo.value.code == 2
    assert "--torque-origin must be" in capsys.readouterr().err


class _FakeProbe:
    last_wrench_kwargs: dict[str, object] = {}
    last_path_kwargs: dict[str, object] = {}

    def wrench(self, **kwargs) -> ObjectWrench:
        type(self).last_wrench_kwargs = dict(kwargs)
        components = {
            "other_objects_all_images": WrenchComponent(
                force_N=np.array([1.0, 2.0, 3.0]),
                torque_Nm=np.array([0.1, 0.2, 0.3]),
                potential_energy_J=4.0,
            ),
            "target_periodic_images": WrenchComponent(
                force_N=np.array([0.5, 0.0, 1.0]),
                torque_Nm=np.array([0.0, 0.1, 0.0]),
                potential_energy_J=2.0,
            ),
            "external_uniform": WrenchComponent(
                force_N=np.zeros(3),
                torque_Nm=np.zeros(3),
                potential_energy_J=0.0,
            ),
            "total_external": WrenchComponent(
                force_N=np.array([1.5, 2.0, 4.0]),
                torque_Nm=np.array([0.1, 0.3, 0.3]),
                potential_energy_J=6.0,
            ),
        }
        return ObjectWrench(
            mesh_id=6,
            step=-1,
            total_charge_C=2.0e-9,
            force_N=np.array([1.5, 2.0, 4.0]),
            torque_Nm=np.array([0.1, 0.3, 0.3]),
            torque_origin_m=np.array([0.0, 0.0, 0.0]),
            transform=RigidTransform.identity(),
            transform_origin_m=np.zeros(3),
            components=components,
            numerical_metadata={
                "periodic_model": "configured",
                "effective_far_correction": "none",
                "target_integration": "point_centroid",
                "quadrature_order": None,
                "periodic_cache": {
                    "hit": None,
                    "build_count": 0,
                    "fingerprint": None,
                    "path": None,
                },
                "periodic_kneq0": {
                    "force_N": np.array([0.5, 1.0, 2.0]),
                    "torque_Nm": np.array([0.01, 0.02, 0.03]),
                    "potential_energy_J": 3.0,
                },
                "physical_k0": {
                    "force_N": np.array([0.0, 0.0, 0.5]),
                    "torque_Nm": np.zeros(3),
                    "potential_energy_J": 1.0,
                },
                "primary_free_subtraction": {
                    "force_N": np.array([-0.5, 0.0, 0.0]),
                    "torque_Nm": np.zeros(3),
                    "potential_energy_J": -1.0,
                },
                "cached_kneq0_trace_correction": {
                    "force_N": np.array([0.0, 0.0, 0.25]),
                    "torque_Nm": np.zeros(3),
                    "potential_energy_J": 0.0,
                },
            },
        )

    def vertical_path(self, displacement_m: np.ndarray, **kwargs) -> ObjectForcePath:
        type(self).last_path_kwargs = dict(kwargs)
        h = np.asarray(displacement_m)
        force = np.column_stack((np.zeros(h.size), np.zeros(h.size), 4.0 - h))
        torque = np.zeros((h.size, 3))
        base = ObjectForcePath.from_samples(
            h,
            force,
            torque,
            potential_energy_J=None,
        )
        component_force = {
            "other_objects_all_images": np.column_stack(
                (np.zeros(h.size), np.zeros(h.size), 3.0 - h)
            ),
            "target_periodic_images": np.column_stack(
                (np.zeros(h.size), np.zeros(h.size), np.ones(h.size))
            ),
            "external_uniform": np.zeros((h.size, 3)),
            "total_external": force,
        }
        component_torque = {
            name: np.zeros((h.size, 3)) for name in _PATH_COMPONENTS
        }
        return ObjectForcePath(
            displacement_m=base.displacement_m,
            force_N=base.force_N,
            torque_Nm=base.torque_Nm,
            electrostatic_work_J=base.electrostatic_work_J,
            component_force_N=component_force,
            component_torque_Nm=component_torque,
        )


class _FakeSnapshot:
    def __init__(self, step: int | None = None) -> None:
        self.step = step

    def __enter__(self):
        return self

    def __exit__(self, *_args) -> None:
        pass

    def object_probe(self, mesh_id: int):
        assert mesh_id == 6
        return _FakeProbe()


class _FakeHistory:
    def __init__(self, batch_indices: list[int]) -> None:
        self.batch_indices = np.asarray(batch_indices, dtype=np.int64)

    @property
    def has_data(self) -> bool:
        return self.batch_indices.size > 0


class _FakeResult:
    def __init__(
        self,
        *,
        batches: int = 123,
        history: _FakeHistory | None = None,
    ) -> None:
        self.batches = batches
        self.history = history


def test_object_detachment_cli_writes_deterministic_artifacts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    config = tmp_path / "beach.toml"
    config.write_text("[sim]\n", encoding="utf-8")
    monkeypatch.setattr(
        cli_module,
        "load_fortran_result",
        lambda _path: _FakeResult(),
    )
    captured: dict[str, object] = {}

    def fake_snapshot(*_args, **kwargs):
        captured.update(kwargs)
        return _FakeSnapshot(step=kwargs["step"])

    monkeypatch.setattr(
        cli_module.ObjectInteractionSnapshot,
        "from_result",
        fake_snapshot,
    )
    artifact_dir = tmp_path / "analysis"
    args = _required_args(tmp_path)
    args.extend(
        [
            "--gravity-m-s2",
            "0",
            "--output-dir",
            str(artifact_dir),
        ]
    )

    beachx_main(args)

    assert {path.name for path in artifact_dir.iterdir()} == set(_ARTIFACT_NAMES)
    wrench_header = (artifact_dir / "instantaneous_wrench.csv").read_text(
        encoding="utf-8"
    ).splitlines()[0]
    assert wrench_header == (
        "component,component_kind,force_x_N,force_y_N,force_z_N,"
        "torque_x_Nm,torque_y_Nm,torque_z_Nm,potential_energy_J"
    )
    path_header = (artifact_dir / "path.csv").read_text(
        encoding="utf-8"
    ).splitlines()[0]
    assert "potential_difference_work_J" in path_header
    assert "speed_m_s" in path_header
    for component in _PATH_COMPONENTS:
        assert f"{component}_force_z_N" in path_header
        assert f"{component}_torque_z_Nm" in path_header
    expected_component_fields = [
        f"{component}_{kind}_{axis}_{unit}"
        for component in _PATH_COMPONENTS
        for kind, unit in (("force", "N"), ("torque", "Nm"))
        for axis in ("x", "y", "z")
    ]
    assert path_header.split(",")[-len(expected_component_fields) :] == (
        expected_component_fields
    )
    with (artifact_dir / "instantaneous_wrench.csv").open(
        encoding="utf-8",
        newline="",
    ) as stream:
        rows = {row["component"]: row for row in csv.DictReader(stream)}
    assert rows["total_external"]["component_kind"] == "physical_total"
    assert rows["periodic_kneq0"]["component_kind"] == "numerical_decomposition"
    assert rows["cached_kneq0_trace_correction"]["component_kind"] == (
        "numerical_diagnostic_included"
    )
    physical_fz = sum(
        float(rows[name]["force_z_N"])
        for name in (
            "other_objects_all_images",
            "target_periodic_images",
            "external_uniform",
        )
    )
    assert physical_fz == float(rows["total_external"]["force_z_N"])
    assert physical_fz + float(
        rows["cached_kneq0_trace_correction"]["force_z_N"]
    ) != float(rows["total_external"]["force_z_N"])
    summary_text = (artifact_dir / "summary.json").read_text(encoding="utf-8")
    assert "NaN" not in summary_text
    summary = json.loads(summary_text)
    assert summary["schema_version"] == 1
    assert summary["physics_policy"]["self_policy"] == (
        "exclude_primary_keep_images"
    )
    assert summary["physics_policy"]["surface_trace"] == "principal_value"
    assert summary["path"]["potential_difference_available"] is False
    assert summary["path"]["potential_difference_unavailable_reason"]
    assert summary["release"]["numerically_qualified"] is False
    assert summary["release"]["mechanics_record_numerically_qualified"] is True
    assert summary["inputs"]["requested_step"] == "final"
    assert summary["inputs"]["resolved_snapshot_step"] is None
    assert summary["inputs"]["charge_source"] == "charges.csv"
    assert summary["inputs"]["resolved_charge_batch"] == 123
    assert captured["step"] is None
    assert "Numerically qualified: `False`" in (
        artifact_dir / "report.md"
    ).read_text(encoding="utf-8")


def test_object_detachment_cli_forwards_model_step_and_path_options(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "run").mkdir()
    (tmp_path / "beach.toml").write_text("[sim]\n", encoding="utf-8")
    captured: dict[str, object] = {}

    def fake_snapshot(*_args, **kwargs):
        captured.update(kwargs)
        return _FakeSnapshot(step=kwargs["step"])

    monkeypatch.setattr(
        cli_module,
        "load_fortran_result",
        lambda _path: _FakeResult(),
    )
    monkeypatch.setattr(
        cli_module.ObjectInteractionSnapshot,
        "from_result",
        fake_snapshot,
    )
    artifact_dir = tmp_path / "analysis"
    args = _required_args(tmp_path)
    args[args.index("configured")] = "infinite-physical"
    args.extend(
        [
            "--step",
            "17",
            "--fixed-grid",
            "--relative-tolerance",
            "0.02",
            "--force-absolute-tolerance-n",
            "3e-12",
            "--work-absolute-tolerance-j",
            "4e-18",
            "--max-refinement",
            "5",
            "--torque-origin",
            "1,2,3",
            "--output-dir",
            str(artifact_dir),
        ]
    )

    beachx_main(args)

    assert captured["periodic_model"] == "infinite_physical"
    assert captured["step"] == 17
    np.testing.assert_allclose(
        _FakeProbe.last_wrench_kwargs["torque_origin"],
        [1.0, 2.0, 3.0],
    )
    path_kwargs = dict(_FakeProbe.last_path_kwargs)
    np.testing.assert_allclose(path_kwargs.pop("torque_origin"), [1.0, 2.0, 3.0])
    assert path_kwargs == {
        "adaptive": False,
        "relative_tolerance": 0.02,
        "force_absolute_tolerance_N": 3.0e-12,
        "work_absolute_tolerance_J": 4.0e-18,
        "max_refinement": 5,
        "components": True,
    }
    summary = json.loads((artifact_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["inputs"]["requested_step"] == 17
    assert summary["inputs"]["resolved_snapshot_step"] == 17
    assert summary["inputs"]["charge_source"] == "charge_history.csv"
    assert summary["inputs"]["resolved_charge_batch"] == 17


def test_object_detachment_cli_rejects_ambiguous_step_word(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    args = _required_args(tmp_path) + ["--step", "latest-history"]

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(args)

    assert excinfo.value.code == 2
    assert "final or an integer batch" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("history", "expected_source", "expected_batch"),
    [
        (_FakeHistory([5, 9]), "charge_history.csv", 9),
        (None, "charges.csv", 123),
    ],
)
def test_object_detachment_cli_resolves_minus_one_charge_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    history: _FakeHistory | None,
    expected_source: str,
    expected_batch: int,
) -> None:
    (tmp_path / "run").mkdir()
    (tmp_path / "beach.toml").write_text("[sim]\n", encoding="utf-8")
    monkeypatch.setattr(
        cli_module,
        "load_fortran_result",
        lambda _path: _FakeResult(history=history),
    )
    monkeypatch.setattr(
        cli_module.ObjectInteractionSnapshot,
        "from_result",
        lambda *_args, **kwargs: _FakeSnapshot(step=kwargs["step"]),
    )
    artifact_dir = tmp_path / "analysis"
    args = _required_args(tmp_path) + [
        "--step",
        "-1",
        "--output-dir",
        str(artifact_dir),
    ]

    beachx_main(args)

    summary = json.loads((artifact_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["inputs"]["requested_step"] == -1
    assert summary["inputs"]["resolved_snapshot_step"] == -1
    assert summary["inputs"]["charge_source"] == expected_source
    assert summary["inputs"]["resolved_charge_batch"] == expected_batch


def _install_fake_evaluator(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "run").mkdir(exist_ok=True)
    (tmp_path / "beach.toml").write_text("[sim]\n", encoding="utf-8")
    monkeypatch.setattr(
        cli_module,
        "load_fortran_result",
        lambda _path: _FakeResult(),
    )
    monkeypatch.setattr(
        cli_module.ObjectInteractionSnapshot,
        "from_result",
        lambda *_args, **kwargs: _FakeSnapshot(step=kwargs["step"]),
    )


def test_object_detachment_strict_json_failure_preserves_existing_artifacts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_evaluator(tmp_path, monkeypatch)
    artifact_dir = tmp_path / "analysis"
    artifact_dir.mkdir()
    for name in _ARTIFACT_NAMES:
        (artifact_dir / name).write_text(f"old:{name}\n", encoding="utf-8")
    monkeypatch.setattr(
        cli_module,
        "_build_summary",
        lambda **_kwargs: {"nonfinite": float("nan")},
    )
    args = _required_args(tmp_path) + ["--output-dir", str(artifact_dir)]

    with pytest.raises(ValueError, match="non-finite"):
        beachx_main(args)

    for name in _ARTIFACT_NAMES:
        assert (artifact_dir / name).read_text(encoding="utf-8") == f"old:{name}\n"
    assert {path.name for path in artifact_dir.iterdir()} == set(_ARTIFACT_NAMES)


def test_object_detachment_replace_failure_rolls_back_artifact_set(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_evaluator(tmp_path, monkeypatch)
    artifact_dir = tmp_path / "analysis"
    artifact_dir.mkdir()
    for name in _ARTIFACT_NAMES:
        (artifact_dir / name).write_text(f"old:{name}\n", encoding="utf-8")
    real_replace = os.replace
    replace_count = 0

    def fail_second_replace(source, destination):
        nonlocal replace_count
        replace_count += 1
        if replace_count == 2:
            raise OSError("injected replace failure")
        return real_replace(source, destination)

    monkeypatch.setattr(os, "replace", fail_second_replace)
    args = _required_args(tmp_path) + ["--output-dir", str(artifact_dir)]

    with pytest.raises(OSError, match="injected replace failure"):
        beachx_main(args)

    for name in _ARTIFACT_NAMES:
        assert (artifact_dir / name).read_text(encoding="utf-8") == f"old:{name}\n"
    assert {path.name for path in artifact_dir.iterdir()} == set(_ARTIFACT_NAMES)


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
    "overrides",
    [
        {"eta_translation": 1.1},
        {"gravity_m_s2": -1.0},
        {"shell_max_layers": -1},
    ],
)
def test_example_validates_inputs_before_constructing_beach(
    monkeypatch: pytest.MonkeyPatch,
    overrides: dict[str, object],
) -> None:
    values = {
        "z_max_m": 1.0,
        "z_points": 3,
        "mass_kg": 1.0,
        "gravity_m_s2": 0.0,
        "eta_translation": 1.0,
        "adhesion_force_n": 0.0,
        "adhesion_range_m": 0.0,
        "shell_max_layers": None,
        "periodic_model": "infinite-physical",
    }
    values.update(overrides)
    monkeypatch.setattr(
        example_module,
        "Beach",
        lambda *_args, **_kwargs: pytest.fail("Beach constructed before validation"),
    )

    with pytest.raises(ValueError):
        example_module._validate_args(SimpleNamespace(**values))


def _kernel_lib() -> Path:
    path = Path("build/libbeach_field_kernel.so")
    if not path.exists():
        pytest.skip("field kernel shared library is not built; run `make build-kernel`")
    return path


def test_object_detachment_cli_native_two_object_smoke(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "summary.txt").write_text(
        """
mesh_nelem=2
processed_particles=0
absorbed=0
escaped=0
batches=1
last_rel_change=0.0
field_source_model=triangle_p0
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (run_dir / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-9\n2,2.0e-9\n",
        encoding="utf-8",
    )
    (run_dir / "mesh_triangles.csv").write_text(
        """elem_idx,mesh_id,v0_x_m,v0_y_m,v0_z_m,v1_x_m,v1_y_m,v1_z_m,v2_x_m,v2_y_m,v2_z_m
1,1,-0.02,-0.02,0.0,0.04,-0.02,0.0,-0.02,0.04,0.0
2,2,0.98,-0.02,0.0,1.04,-0.02,0.0,0.98,0.04,0.0
""",
        encoding="utf-8",
    )
    config = tmp_path / "beach.toml"
    config.write_text(
        """
[sim]
field_bc_mode = "free"
box_min = [-2.0, -2.0, -2.0]
box_max = [2.0, 2.0, 2.0]
e0 = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )
    artifact_dir = tmp_path / "native-analysis"

    beachx_main(
        [
            "object-detachment",
            str(run_dir),
            "--config",
            str(config),
            "--target-mesh-id",
            "1",
            "--periodic-model",
            "configured",
            "--z-max-m",
            "0.1",
            "--z-points",
            "5",
            "--mass-kg",
            "1.0",
            "--gravity-m-s2",
            "0",
            "--library",
            str(_kernel_lib()),
            "--output-dir",
            str(artifact_dir),
        ]
    )

    summary = json.loads((artifact_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["instantaneous_wrench"]["force_N"][0] < 0.0
    assert summary["path"]["point_count"] >= 5
