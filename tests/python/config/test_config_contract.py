"""Check that ``beachx lint`` returns the schema and Fortran preflight verdicts.

The JSON Schema owns key names, types, and ranges; the Fortran preflight owns
combinations and derived values. Schema cases run without a native executable.
``make test-config-contract`` supplies the freshly built executable for the
preflight cases; ordinary Python-only runs skip them instead of selecting a stale
binary.
"""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
import subprocess

import pytest
import tomli_w

from beach.cli.lint import run as run_lint
from beach.config import default_config, to_grouped_layout, to_runtime_layout
from beach.config._toml import load_toml_file


ROOT = Path(__file__).resolve().parents[3]


def _periodic2_cell() -> dict:
    config = default_config()
    config["domain"]["periodic_axes"] = ["x", "y"]
    config["field_boundary"]["mode"] = "periodic2"
    config["sim"].update(
        field_solver="fmm", field_periodic_image_layers=1,
        field_periodic_far_correction="none",
    )
    return config


def _schema_cases() -> list:
    # Key names, types, and ranges belong to the JSON Schema (and the finite-number
    # check of the Python loader), so lint rejects these before any native check.
    cases = []
    config = to_grouped_layout(default_config())
    config["particles"]["tracking"]["dtt"] = 1.0e-8
    cases.append(pytest.param(config, "('dtt' was unexpected)", id="dtt-1e-08"))
    for key, value, diagnostic in (
        ("dt", -1.0, "schema error at particles.tracking.dt_s"),
        ("dt", float("nan"), "config error at particles.tracking.dt_s: number must be finite"),
        ("dt", float("inf"), "config error at particles.tracking.dt_s: number must be finite"),
        ("dt", True, "schema error at particles.tracking.dt_s"),
        ("batch_count", 0, "schema error at run.batch_count"),
        ("max_step", 0, "schema error at particles.tracking.max_steps_per_particle"),
        ("e0_phi_z_deg", 90.0, "schema error at fields.external"),  # requires e0_abs
    ):
        config = default_config()
        config["sim"][key] = value
        cases.append(pytest.param(config, diagnostic, id=f"{key}-{value}"))

    config = default_config()
    config["sim"].update(e0=[0.0, 0.0, -1.0], e0_abs=1.0)
    cases.append(pytest.param(
        config, "'required': ['electric_v_m', 'electric_magnitude_v_m']",
        id="e0-vector-and-magnitude",
    ))

    config = default_config()
    config["output"].update(resume=True, write_files=False)
    cases.append(pytest.param(config, "schema error at output.enabled", id="resume-without-file-output"))

    config = default_config()
    config["mesh"]["parts"][0]["kind"] = "unknown"
    cases.append(pytest.param(config, "schema error at mesh.parts[0].kind", id="unsupported-template-kind"))

    config = default_config()
    config["mesh"]["groups"] = {"g" * 65: {"scale": 1.0}}
    cases.append(pytest.param(config, "schema error at mesh.groups", id="mesh-group-name-length-65"))

    config = default_config()
    config["sim"].update(batch_duration=1.0, batch_duration_step=1.0)
    cases.append(pytest.param(
        config, "'required': ['duration_s', 'duration_steps']", id="exclusive-batch-duration",
    ))

    reservoir = to_runtime_layout(load_toml_file(ROOT / "examples/beach.toml"))
    reservoir["particles"]["species"][0]["number_density_m3"] = 1.0e6
    cases.append(pytest.param(
        reservoir, "'required': ['number_density_cm3', 'number_density_m3']",
        id="reservoir-two-density-forms",
    ))
    config = default_config()
    config["particles"]["species"][0]["target_macro_particles_per_batch"] = 10
    cases.append(pytest.param(
        config, "'required': ['target_macro_particles_per_batch']",
        id="volume-seed-with-reservoir-target",
    ))
    return cases


def _preflight_cases() -> list:
    # Bundled examples (including the tutorial default) are accepted by
    # test_config_layout.test_grouped_examples_pass_native_config_check. Rules whose
    # diagnostics tests/fortran/test_app_config_parser.f90 already checks in L1 are
    # not repeated here; test_lint_reports_the_fortran_diagnostic compares the relay
    # with the native stderr. A diagnostic of None means the configuration is valid.
    cases = []
    for key, value in (
        ("tree_theta", -1.0),  # inactive for direct evaluation
        ("field_solver", "DIRECT"),
    ):
        config = default_config()
        config["sim"][key] = value
        cases.append(pytest.param(config, None, id=f"{key}-{value}"))

    config = default_config()
    config["sim"].update(field_solver="treecode", tree_theta=-1.0)
    cases.append(pytest.param(
        config, "fields.solver.tree.theta must satisfy 0 < theta <= 1", id="active-tree-theta",
    ))

    config = default_config()
    config["mesh"]["groups"] = {"g" * 64: {"scale": 1.0}}
    cases.append(pytest.param(config, None, id="mesh-group-name-length-64"))

    for dt, steps, label in ((1.0e308, 1.0e308, "overflow"), (1.0e-300, 1.0e-300, "underflow")):
        config = default_config()
        del config["sim"]["batch_duration"]
        config["sim"].update(dt=dt, batch_duration_step=steps)
        cases.append(pytest.param(
            config, "run.batch.duration_steps produced invalid run.batch.duration_s",
            id=f"resolved-batch-duration-{label}",
        ))

    config = default_config()
    config["domain"] = {"box_origin": [0.0, 0.0, 0.0], "box_size": [1.0, 1.0, 1.0]}
    cases.append(pytest.param(config, None, id="domain-authoring"))

    for species_override, diagnostic in (
        (False, "[particles.boundary] requires a finite [domain]"),
        (True, "particles.species[1].boundary requires a finite [domain]"),
    ):
        config = default_config()
        del config["domain"]
        config["particle_boundary"] = {}
        if species_override:
            config["particles"]["species"][0]["boundary"] = {"z_high": "reflect"}
        else:
            config["particle_boundary"]["z_high"] = "reflect"
        cases.append(pytest.param(config, diagnostic, id=f"boundary-without-domain-{species_override}"))

    config = _periodic2_cell()
    config["particles"]["species"][0]["boundary"] = {"x_high": "open"}
    cases.append(pytest.param(
        config, "particles.species[1].boundary cannot override a periodic domain face",
        id="species-override-of-periodic-face",
    ))

    config = default_config()
    config["particles"]["species"][0]["q_particle"] = 0.0
    cases.append(pytest.param(
        config, "particles.species[1].charge_c must be finite and non-zero", id="zero-charge",
    ))
    config = default_config()
    species = copy.deepcopy(config["particles"]["species"][0])
    species.update(species_key="species_1", enabled=False)
    config["particles"]["species"].append(species)
    cases.append(pytest.param(
        config, "particles.species[2].species_key values must be unique",
        id="duplicate-disabled-species-key",
    ))

    config = default_config()
    config["particles"]["species"][0]["enabled"] = False
    cases.append(pytest.param(
        config, "At least one enabled [[particles.species]] entry", id="disabled-only-workload",
    ))

    config = default_config()
    disabled = copy.deepcopy(config["particles"]["species"][0])
    disabled.update(enabled=False, q_particle=0.0)
    config["particles"]["species"].append(disabled)
    cases.append(pytest.param(config, None, id="active-source-with-disabled-zero-charge"))

    reservoir = to_runtime_layout(load_toml_file(ROOT / "examples/beach.toml"))
    first_species = reservoir["particles"]["species"][0]
    first_species.pop("target_macro_particles_per_batch")
    first_species["w_particle"] = 1.0e6
    cases.append(pytest.param(reservoir, None, id="reservoir-fixed-weight"))
    config = copy.deepcopy(reservoir)
    config["domain"]["periodic_axes"] = ["z"]
    cases.append(pytest.param(
        config, "particles.boundary.z_low cannot override a periodic domain face",
        id="reservoir-inflow-through-periodic-face",
    ))

    config = copy.deepcopy(reservoir)
    config["particles"]["species"][0].update(
        surface_charge_closure="fixed_current", target_absorbed_current_a=-1.0e-9,
    )
    cases.append(pytest.param(config, None, id="fixed-current-absorbed"))

    periodic = to_runtime_layout(load_toml_file(ROOT / "examples/closed_photoelectron.toml"))
    periodic["periodic2"] = {"lower_boundary_model": "symmetric_vacuum"}
    cases.append(pytest.param(
        copy.deepcopy(periodic), "panel_spectral_reference requires field backend direct",
        id="partial-periodic-wrong-solver",
    ))
    periodic["sim"]["field_solver"] = "direct"
    cases.append(pytest.param(copy.deepcopy(periodic), None, id="partial-periodic-reference"))
    config = copy.deepcopy(periodic)
    config["sim"]["field_periodic_image_layers"] = -1
    cases.append(pytest.param(
        config, "fields.periodic.image_layers must be >= 0", id="negative-periodic-image-layers",
    ))
    periodic["sim"].update(field_solver="fmm", field_periodic_far_correction="cached_kneq0")
    periodic["periodic2"] = {
        "nonzero_mode_backend": "cached_kneq0",
        "max_nonzero_mode_potential_step": 1.0e-2,
    }
    cases.append(pytest.param(copy.deepcopy(periodic), None, id="partial-periodic-cached"))
    config = copy.deepcopy(periodic)
    config["particles"]["species"][-1]["boundary"]["z_high"] = "open"
    cases.append(pytest.param(
        config, "neutral_return requires a negative photo_raycast species",
        id="neutral-return-without-reflecting-face",
    ))

    config = _periodic2_cell()
    config["sim"].update(batch_duration=1.0e-6, field_periodic_far_correction="cached_kneq0")
    config["periodic2"] = {
        "nonzero_mode_backend": "cached_kneq0",
        "zero_mode_policy": "exclude_k0",
        "lower_boundary_model": "e_bottom_zero",
        "max_nonzero_mode_potential_step": 1.0e-2,
    }
    cases.append(pytest.param(
        config, "max_nonzero_mode_potential_step_v requires time-scaled reservoir or photo_raycast",
        id="adaptive-step-with-volume-seed",
    ))
    for label, mesh in (
        ("absent", {}), ("empty", {"parts": []}),
        ("disabled", {"parts": [{"kind": "plane", "enabled": False}]}),
    ):
        config = default_config()
        config["mesh"] = mesh
        cases.append(pytest.param(config, None, id=f"no-enabled-mesh-{label}"))
    return cases


def _assert_schema_rejection(path: Path, diagnostic: str, beach: Path) -> None:
    with pytest.raises(SystemExit) as caught:
        run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=beach))
    message = str(caught.value)
    assert message.startswith(("schema validation failed:", "config error at ")), message
    assert diagnostic in message


def _write_config(config: dict, directory: Path) -> Path:
    config = copy.deepcopy(config)
    config["output"]["dir"] = str(directory / "simulation-output")
    path = directory / "beach.toml"
    path.write_text(tomli_w.dumps(to_grouped_layout(config)), encoding="utf-8")
    return path


@pytest.mark.parametrize("config,diagnostic", _schema_cases())
def test_lint_rejects_schema_violations_before_the_preflight(
    config: dict, diagnostic: str, tmp_path: Path,
) -> None:
    # A missing executable proves the native preflight was never consulted.
    _assert_schema_rejection(_write_config(config, tmp_path), diagnostic, tmp_path / "missing-beach")


@pytest.mark.parametrize("config,diagnostic", _preflight_cases())
def test_lint_returns_the_preflight_verdict(
    config: dict, diagnostic: str | None, native_checker: Path, tmp_path: Path,
) -> None:
    path = _write_config(config, tmp_path)
    args = argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=native_checker)
    if diagnostic is None:
        run_lint(args)
    else:
        with pytest.raises(SystemExit) as caught:
            run_lint(args)
        message = str(caught.value)
        assert message.startswith(f"preflight check failed: {path}"), message
        assert diagnostic in message
    assert not (tmp_path / "simulation-output").exists(), "configuration checks must not start a simulation"


def _periodic_face_override() -> dict:
    config = _periodic2_cell()
    config["particle_boundary"]["x_low"] = "reflect"
    return config


def _wrong_sign_absorbed_target() -> dict:
    config = to_runtime_layout(load_toml_file(ROOT / "examples/beach.toml"))
    config["particles"]["species"][0].update(
        surface_charge_closure="fixed_current", target_absorbed_current_a=1.0e-9,
    )
    return config


def _no_volume_macro_particles() -> dict:
    config = default_config()
    config["particles"]["species"][0]["npcls_per_step"] = 0
    return config


@pytest.mark.parametrize("make_config, diagnostic", [
    (_periodic_face_override, "particles.boundary.x_low cannot override a periodic domain face"),
    (_wrong_sign_absorbed_target, "same sign as particles.species[1].charge_c"),
    (_no_volume_macro_particles, "sampling.volume_macro_particles_per_batch > 0"),
])
def test_lint_reports_the_fortran_diagnostic(
    make_config, diagnostic: str, native_checker: Path, tmp_path: Path,
) -> None:
    path = _write_config(make_config(), tmp_path)
    native = subprocess.run(
        [str(native_checker), "--check-config", str(path)], cwd=tmp_path,
        text=True, capture_output=True, timeout=20,
    )
    assert native.returncode != 0
    with pytest.raises(SystemExit) as excinfo:
        run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=native_checker))
    lines = str(excinfo.value).splitlines()
    assert lines[:2] == [f"preflight check failed: {path}", f"beach={native_checker}"]
    relayed = lines[2:]
    assert any(diagnostic in line for line in relayed)
    assert set(relayed) <= set(native.stderr.splitlines())
    assert not any("Backtrace" in line for line in relayed)


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


@pytest.mark.parametrize("removed, schema_diagnostic, native_diagnostic", [
    ("reservoir_face", "schema error at particles.species[0].source.mode",
     "particles.species[1].source.mode"),
    ("plane_source", "schema error at particles.species[0].source.mode",
     "particles.species[1].source.mode"),
    ("source_normal", "('normal' was unexpected)", "particles.species.source.normal"),
])
def test_removed_particle_source_inputs_are_rejected(
    removed: str, schema_diagnostic: str, native_diagnostic: str,
    native_checker: Path, tmp_path: Path,
) -> None:
    config = to_grouped_layout(default_config())
    config["particles"]["species"][0]["source"]["normal" if removed == "source_normal" else "mode"] = (
        [0.0, 0.0, -1.0] if removed == "source_normal" else removed
    )
    _assert_rejected_by_lint_and_native(
        config, schema_diagnostic, native_diagnostic, native_checker, tmp_path,
    )


@pytest.mark.parametrize("value", ["escape", "potential_barrier"])
def test_removed_ordinary_open_key_is_rejected(
    value: str, native_checker: Path, tmp_path: Path,
) -> None:
    config = to_grouped_layout(default_config())
    config["particles"]["boundary"]["open_model"] = value
    _assert_rejected_by_lint_and_native(
        config, "('open_model' was unexpected)", "particles.boundary.open_model",
        native_checker, tmp_path,
    )


def _assert_rejected_by_lint_and_native(
    config, schema_diagnostic, native_diagnostic, native_checker, tmp_path,
):
    """Removed inputs fail the schema in lint and the native loader when beach runs directly."""
    config["output"]["dir"] = str(tmp_path / "simulation-output")
    path = tmp_path / "beach.toml"
    path.write_text(tomli_w.dumps(config), encoding="utf-8")
    _assert_schema_rejection(path, schema_diagnostic, native_checker)
    native = subprocess.run(
        [str(native_checker), "--check-config", str(path)], cwd=tmp_path,
        text=True, capture_output=True, timeout=20,
    )
    assert native.returncode != 0
    assert native_diagnostic in native.stderr
    assert not (tmp_path / "simulation-output").exists()


@pytest.mark.parametrize("charge", [0.0, "invalid"])
def test_native_species_diagnostic_uses_public_key_and_index(
    charge: float | str, native_checker: Path, tmp_path: Path,
) -> None:
    config = default_config()
    second = copy.deepcopy(config["particles"]["species"][0])
    second["q_particle"] = charge
    config["particles"]["species"].append(second)
    path = _write_config(config, tmp_path)
    native = subprocess.run(
        [str(native_checker), "--check-config", str(path)], cwd=tmp_path,
        text=True, capture_output=True, timeout=20,
    )
    assert native.returncode != 0
    assert "particles.species[2].charge_c" in native.stderr
    assert "q_particle" not in native.stderr


def _photo_raycast_tutorial() -> dict:
    config = default_config()
    config["particles"]["species"] = [{
        "source_mode": "photo_raycast",
        "q_particle": -1.602176634e-19,
        "m_particle": 9.1093837139e-31,
        "temperature_ev": 2.0,
        "emit_current_density_a_m2": 1.0e-4,
        "rays_per_batch": 100,
        "inject_face": "z_high",
        "pos_low": [0.0, 0.0, 0.07],
        "pos_high": [0.04, 0.04, 0.07],
    }]
    return config


def _without_domain(config: dict) -> dict:
    del config["domain"]
    config["particle_boundary"] = {}
    return config


@pytest.mark.parametrize("edit, diagnostic", [
    (lambda config: config, None),
    (_without_domain, "photo_raycast requires a finite [domain]"),
    (lambda config: config["sim"].update(batch_duration=0.0) or config,
     "photo_raycast requires run.batch.duration_s > 0"),
    (lambda config: config["sim"].update(batch_count=0) or config, "run.batch_count must be > 0"),
    (lambda config: config["sim"].update(batch_duration=-1.0e-8) or config, "run.batch.duration_s must be >= 0"),
    (lambda config: config["particles"]["species"][0].update(ray_direction=[0.0, 0.0, 1.0]) or config,
     "ray_direction must point inward"),
])
def test_photo_raycast_inputs_are_owned_by_preflight(edit, diagnostic, native_checker, tmp_path):
    path = _write_config(edit(_photo_raycast_tutorial()), tmp_path)
    native = subprocess.run([str(native_checker), "--check-config", str(path)], cwd=tmp_path,
                            text=True, capture_output=True, timeout=20)
    if diagnostic is None:
        assert native.returncode == 0, native.stdout + native.stderr
    else:
        assert native.returncode != 0
        assert diagnostic in native.stdout + native.stderr


# group, path, and kind rules: test_app_config_parser mesh_parts_public_contract.
@pytest.mark.parametrize("part, diagnostic", [
    ({"kind": "obj", "path": "unused.obj", "placement_mode": "absolute"}, "placement_mode"),
    ({"kind": "obj", "path": "unused.obj", "center": [0, 0, 0]}, "center"),
    ({"kind": "plane", "scale": 1.0}, "scale"),
    ({"kind": "plane", "rotation": [0, 0, 0]}, "rotation"),
    ({"kind": "sphere", "tessellation": "uniform"}, "tessellation"),
    ({"kind": "disk", "tessellation": "icosphere"}, "tessellation"),
    ({"kind": "plane", "tessellation": "uniform"}, "tessellation"),
])
def test_mesh_part_combinations_are_owned_by_preflight(part, diagnostic, native_checker, tmp_path):
    config = default_config()
    part = {**part, "surface_side": "normal_plus"}
    if part.get("kind") != "obj":
        part["center"] = [0, 0, 0]
    config["mesh"]["parts"] = [part]
    path = _write_config(config, tmp_path)
    native = subprocess.run([str(native_checker), "--check-config", str(path)],
                            text=True, capture_output=True, timeout=20)
    assert native.returncode != 0
    assert f"mesh.parts[1].{diagnostic}" in native.stdout + native.stderr
    with pytest.raises(SystemExit) as caught:
        run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20, beach=native_checker))
    assert "mesh.parts" in str(caught.value)
