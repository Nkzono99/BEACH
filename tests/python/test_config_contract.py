"""Check that ``beachx lint`` returns the Fortran preflight verdict.

The JSON Schema owns key names, types, and ranges; the Fortran preflight owns
combinations and derived values. ``make test-config-contract`` supplies the
freshly built executable. Ordinary Python-only runs skip these cases instead of
selecting a stale binary.
"""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
import subprocess

import pytest
import tomli_w

from beach.cli.lint import run_lint
from beach.config import default_config
from beach.config._layout import to_runtime_layout
from beach.config._toml import load_toml_file


ROOT = Path(__file__).resolve().parents[2]


def _periodic2_cell() -> dict:
    config = default_config()
    config["domain"]["periodic_axes"] = ["x", "y"]
    config["field_boundary"]["mode"] = "periodic2"
    config["sim"].update(
        field_solver="fmm", field_periodic_image_layers=1,
        field_periodic_far_correction="none",
    )
    return config


def _cases() -> list:
    cases = [pytest.param(default_config(), True, id="tutorial")]
    for name in (
        "periodic2_zhao_fixed_current",
        "periodic2_zhao_no_photo_fixed_current",
        "periodic2_zhao_outflow_refresh",
    ):
        config = load_toml_file(ROOT / "examples" / f"{name}.toml")
        cases.append(pytest.param(config, True, id=name))
    for table_path in (
        ("sim",), ("domain",), ("mesh",), ("output",), ("particles",),
        ("field_boundary",), ("particle_boundary",),
        ("particles", "species", 0), ("mesh", "templates", 0),
    ):
        config = default_config()
        table = config
        for component in table_path:
            table = table[component]
        table["unrecognized_key"] = 1
        label = ".".join(str(component) for component in table_path)
        cases.append(pytest.param(config, False, id=f"unknown-{label}"))
    for key, value, accepted in (
        ("dt", -1.0, False),
        ("dt", float("nan"), False),
        ("dt", float("inf"), False),
        ("dt", True, False),
        ("dtt", 1.0e-8, False),
        ("batch_count", 0, False),
        ("max_step", 0, False),
        ("rng_seed", 2**32 + 123, False),
        ("rng_seed", -(2**32) + 123, False),
        ("tree_theta", -1.0, True),  # inactive for direct evaluation
        ("field_solver", "DIRECT", True),
        ("e0_phi_z_deg", 90.0, False),  # requires e0_abs
    ):
        config = default_config()
        config["sim"][key] = value
        cases.append(pytest.param(config, accepted, id=f"{key}-{value}"))

    config = default_config()
    config["sim"].update(field_solver="treecode", tree_theta=-1.0)
    cases.append(pytest.param(config, False, id="active-tree-theta"))

    config = default_config()
    config["sim"].update(e0=[0.0, 0.0, -1.0], e0_abs=1.0)
    cases.append(pytest.param(config, False, id="e0-vector-and-magnitude"))

    config = default_config()
    config["output"].update(resume=True, write_files=False)
    cases.append(pytest.param(config, False, id="resume-without-file-output"))

    config = default_config()
    config["mesh"]["templates"][0]["kind"] = "unknown"
    cases.append(pytest.param(config, False, id="unsupported-template-kind"))

    config = default_config()
    config["mesh"]["mode"] = "auto"
    cases.append(pytest.param(config, False, id="auto-mesh-without-surface-side"))

    for name_length, accepted in ((64, True), (65, False)):
        config = default_config()
        config["mesh"]["groups"] = {"g" * name_length: {"scale": 1.0}}
        cases.append(pytest.param(config, accepted, id=f"mesh-group-name-length-{name_length}"))

    config = default_config()
    config["sim"].update(batch_duration=1.0, batch_duration_step=1.0)
    cases.append(pytest.param(config, False, id="exclusive-batch-duration"))

    for dt, steps, label in ((1.0e308, 1.0e308, "overflow"), (1.0e-300, 1.0e-300, "underflow")):
        config = default_config()
        config["sim"].update(dt=dt, batch_duration_step=steps)
        cases.append(pytest.param(config, False, id=f"resolved-batch-duration-{label}"))

    config = default_config()
    config["sim"]["b0"] = [0.0, float("nan"), 0.0]
    cases.append(pytest.param(config, False, id="nonfinite-vector"))

    config = default_config()
    config["particles"]["species"][0]["species_key"] = "s" * 65
    cases.append(pytest.param(config, False, id="overlong-species-key"))

    config = default_config()
    config["domain"] = {"box_origin": [0.0, 0.0, 0.0], "box_size": [1.0, 1.0, 1.0]}
    cases.append(pytest.param(config, True, id="domain-authoring"))

    for species_override in (False, True):
        config = default_config()
        del config["domain"]
        config["particle_boundary"] = {}
        if species_override:
            config["particles"]["species"][0]["boundary"] = {"z_high": "reflect"}
        else:
            config["particle_boundary"]["z_high"] = "reflect"
        cases.append(pytest.param(config, False, id=f"boundary-without-domain-{species_override}"))

    config = _periodic2_cell()
    config["particle_boundary"]["x_low"] = "reflect"
    cases.append(pytest.param(config, False, id="global-override-of-periodic-face"))
    config = _periodic2_cell()
    config["particles"]["species"][0]["boundary"] = {"x_high": "open"}
    cases.append(pytest.param(config, False, id="species-override-of-periodic-face"))

    config = default_config()
    config["particles"]["species"][0]["q_particle"] = 0.0
    cases.append(pytest.param(config, False, id="zero-charge"))
    for enabled in (False, True):
        config = default_config()
        species = copy.deepcopy(config["particles"]["species"][0])
        species.update(species_key="species_1", enabled=enabled)
        config["particles"]["species"].append(species)
        cases.append(pytest.param(config, False, id=f"duplicate-species-key-{enabled}"))

    config = default_config()
    config["particles"]["species"][0]["enabled"] = False
    cases.append(pytest.param(config, False, id="disabled-only-workload"))

    config = default_config()
    disabled = copy.deepcopy(config["particles"]["species"][0])
    disabled.update(enabled=False, q_particle=0.0)
    config["particles"]["species"].append(disabled)
    cases.append(pytest.param(config, True, id="active-source-with-disabled-zero-charge"))

    reservoir = to_runtime_layout(load_toml_file(ROOT / "examples/beach.toml"))
    first_species = reservoir["particles"]["species"][0]
    first_species.pop("target_macro_particles_per_batch")
    first_species["w_particle"] = 1.0e6
    cases.append(pytest.param(reservoir, True, id="reservoir-fixed-weight"))
    for density in (float("nan"), float("inf")):
        config = copy.deepcopy(reservoir)
        config["particles"]["species"][0]["number_density_cm3"] = density
        cases.append(pytest.param(config, False, id=f"reservoir-density-{density}"))
    config = copy.deepcopy(reservoir)
    config["particles"]["species"][0]["number_density_m3"] = 1.0e6
    cases.append(pytest.param(config, False, id="reservoir-two-density-forms"))
    config = copy.deepcopy(reservoir)
    config["particles"]["species"][0]["boundary"] = {"z_high": "reflect"}
    cases.append(pytest.param(config, False, id="reservoir-inflow-through-reflecting-face"))

    for target, accepted in ((-1.0e-9, True), (1.0e-9, False)):
        config = copy.deepcopy(reservoir)
        config["particles"]["species"][0].update(
            surface_charge_closure="fixed_current", target_absorbed_current_a=target,
        )
        cases.append(pytest.param(config, accepted, id=f"fixed-current-absorbed-{target}"))
    config = copy.deepcopy(reservoir)
    config["particles"]["species"][0]["target_absorbed_current_a"] = -1.0e-9
    cases.append(pytest.param(config, False, id="target-current-with-explicit-closure"))

    periodic = load_toml_file(ROOT / "examples/periodic2_closed_photoelectron.toml")
    periodic["periodic2"] = {"lower_boundary_model": "symmetric_vacuum"}
    cases.append(pytest.param(copy.deepcopy(periodic), False, id="partial-periodic-wrong-solver"))
    periodic["sim"]["field_solver"] = "direct"
    cases.append(pytest.param(copy.deepcopy(periodic), True, id="partial-periodic-reference"))
    config = copy.deepcopy(periodic)
    config["sim"]["field_periodic_image_layers"] = -1
    cases.append(pytest.param(config, False, id="negative-periodic-image-layers"))
    periodic["sim"].update(field_solver="fmm", field_periodic_far_correction="cached_kneq0")
    periodic["periodic2"] = {
        "nonzero_mode_backend": "cached_kneq0",
        "max_nonzero_mode_potential_step": 1.0e-2,
    }
    cases.append(pytest.param(copy.deepcopy(periodic), True, id="partial-periodic-cached"))
    config = copy.deepcopy(periodic)
    config["particles"]["species"][-1]["boundary"]["z_high"] = "open"
    cases.append(pytest.param(config, False, id="neutral-return-without-reflecting-face"))

    config = _periodic2_cell()
    config["sim"].update(batch_duration=1.0e-6, field_periodic_far_correction="cached_kneq0")
    config["periodic2"] = {
        "nonzero_mode_backend": "cached_kneq0",
        "zero_mode_policy": "exclude_k0",
        "lower_boundary_model": "e_bottom_zero",
        "max_nonzero_mode_potential_step": 1.0e-2,
    }
    cases.append(pytest.param(config, False, id="adaptive-step-with-volume-seed"))
    return cases


def _write_config(config: dict, directory: Path) -> Path:
    config = copy.deepcopy(config)
    config["output"]["dir"] = str(directory / "simulation-output")
    path = directory / "beach.toml"
    path.write_text(tomli_w.dumps(config), encoding="utf-8")
    return path


@pytest.mark.parametrize("config,accepted", _cases())
def test_lint_returns_the_preflight_verdict(
    config: dict, accepted: bool, native_checker: Path, tmp_path: Path,
) -> None:
    path = _write_config(config, tmp_path)
    args = argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=native_checker)
    if accepted:
        run_lint(args)
    else:
        with pytest.raises(SystemExit):
            run_lint(args)
    assert not (tmp_path / "simulation-output").exists(), "configuration checks must not start a simulation"


def test_lint_reports_the_fortran_diagnostic(native_checker: Path, tmp_path: Path) -> None:
    config = _periodic2_cell()
    config["particle_boundary"]["x_low"] = "reflect"
    path = _write_config(config, tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=native_checker))
    message = str(excinfo.value)
    assert "preflight check failed" in message
    assert "cannot override a periodic domain face" in message
    assert "Backtrace" not in message


def test_lint_requires_the_beach_executable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = _write_config(default_config(), tmp_path)
    monkeypatch.delenv("BEACH_CONFIG_CHECK_EXE", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(SystemExit, match="beach executable not found"):
        run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=None))


def test_native_check_requires_one_path(native_checker: Path, tmp_path: Path) -> None:
    result = subprocess.run(
        [str(native_checker), "--check-config"], cwd=tmp_path,
        text=True, capture_output=True, timeout=20,
    )
    assert result.returncode != 0
    assert "usage: beach --check-config beach.toml" in result.stderr


@pytest.mark.parametrize("grouped", [False, True], ids=["flat", "grouped"])
@pytest.mark.parametrize("removed", ["reservoir_face", "plane_source", "source_normal"])
def test_removed_particle_source_inputs_are_rejected(
    grouped: bool, removed: str, native_checker: Path, tmp_path: Path,
) -> None:
    from beach.config._layout import to_grouped_layout

    config = default_config()
    if grouped:
        config = to_grouped_layout(config)
        source = config["particles"]["species"][0]["source"]
        source["normal" if removed == "source_normal" else "mode"] = (
            [0.0, 0.0, -1.0] if removed == "source_normal" else removed
        )
        config["output"]["dir"] = str(tmp_path / "simulation-output")
        path = tmp_path / "beach.toml"
        path.write_text(tomli_w.dumps(config), encoding="utf-8")
    else:
        species = config["particles"]["species"][0]
        species["source_normal" if removed == "source_normal" else "source_mode"] = (
            [0.0, 0.0, -1.0] if removed == "source_normal" else removed
        )
        path = _write_config(config, tmp_path)
    with pytest.raises(SystemExit):
        run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=native_checker))
    native = subprocess.run(
        [str(native_checker), "--check-config", str(path)], cwd=tmp_path,
        text=True, capture_output=True, timeout=20,
    )
    assert native.returncode != 0
    assert not (tmp_path / "simulation-output").exists()


@pytest.mark.parametrize("grouped", [False, True], ids=["flat", "grouped"])
@pytest.mark.parametrize("value", ["escape", "potential_barrier"])
def test_removed_ordinary_open_key_is_rejected(
    grouped: bool, value: str, native_checker: Path, tmp_path: Path,
) -> None:
    from beach.config._layout import to_grouped_layout

    config = default_config()
    if grouped:
        config = to_grouped_layout(config)
        key = "open_model"
        config["particles"]["boundary"][key] = value
    else:
        key = "ordinary_open_model"
        config["particle_boundary"][key] = value
    path = tmp_path / "beach.toml"
    path.write_text(tomli_w.dumps(config), encoding="utf-8")

    with pytest.raises(SystemExit):
        run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=native_checker))
    native = subprocess.run(
        [str(native_checker), "--check-config", str(path)], cwd=tmp_path,
        text=True, capture_output=True, timeout=20,
    )
    assert native.returncode != 0
    assert key in native.stdout + native.stderr
