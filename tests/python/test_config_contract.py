"""Run the same authored TOML through Python and the actual Fortran preflight.

``make test-config-contract`` supplies the freshly built executable. Ordinary
Python-only runs skip these integration tests instead of selecting a stale binary.
"""

from __future__ import annotations

import argparse
import copy
import os
from pathlib import Path
import subprocess

import pytest

from beach.cli.config import run_validate
from beach.cli.lint import run_lint
from beach.config import ConfigError, default_config, dump_beach_toml, load_config_file
from beach.config._toml import load_toml_file


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def native_checker() -> Path:
    executable = os.environ.get("BEACH_CONFIG_CHECK_EXE")
    if not executable:
        pytest.skip("run make test-config-contract for Python/Fortran parity")
    path = Path(executable).resolve()
    assert path.is_file(), f"missing Fortran checker: {path}"
    return path


def _cases() -> list:
    cases = [pytest.param(default_config(), True, id="tutorial")]
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
    ):
        config = default_config()
        config["sim"][key] = value
        cases.append(pytest.param(config, accepted, id=f"{key}-{value}"))

    config = default_config()
    config["sim"].update(field_solver="treecode", tree_theta=-1.0)
    cases.append(pytest.param(config, False, id="active-tree-theta"))

    config = default_config()
    config["output"].update(resume=True, write_files=False)
    cases.append(pytest.param(config, False, id="resume-without-file-output"))

    config = default_config()
    config["mesh"]["templates"][0]["kind"] = "unknown"
    cases.append(pytest.param(config, False, id="unsupported-template-kind"))

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
        config["particle_boundary"] = {"ordinary_open_model": "escape"}
        if species_override:
            config["particles"]["species"][0]["boundary"] = {"z_high": "reflect"}
        else:
            config["particle_boundary"]["z_high"] = "reflect"
        cases.append(pytest.param(config, False, id=f"boundary-without-domain-{species_override}"))

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

    reservoir = load_toml_file(ROOT / "examples/beach.toml")
    first_species = reservoir["particles"]["species"][0]
    first_species.pop("target_macro_particles_per_batch")
    first_species["w_particle"] = 1.0e6
    cases.append(pytest.param(reservoir, True, id="reservoir-fixed-weight"))
    for density in (float("nan"), float("inf")):
        config = copy.deepcopy(reservoir)
        config["particles"]["species"][0]["number_density_cm3"] = density
        cases.append(pytest.param(config, False, id=f"reservoir-density-{density}"))

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
    cases.append(pytest.param(periodic, True, id="partial-periodic-cached"))
    return cases


@pytest.mark.parametrize("config,accepted", _cases())
def test_same_config_in_all_entry_points(
    config: dict, accepted: bool, native_checker: Path, tmp_path: Path,
) -> None:
    config = copy.deepcopy(config)
    output = tmp_path / "simulation-output"
    config["output"]["dir"] = str(output)
    path = tmp_path / "beach.toml"
    path.write_text(dump_beach_toml(config), encoding="utf-8")

    for check in (
        lambda: load_config_file(path),
        lambda: run_validate(argparse.Namespace(config_path=path)),
        lambda: run_lint(argparse.Namespace(config_path=path, schema=None, max_errors=20)),
    ):
        if accepted:
            check()
        else:
            with pytest.raises((ConfigError, SystemExit)):
                check()

    result = subprocess.run(
        [str(native_checker), "--check-config", str(path)],
        cwd=tmp_path, text=True, capture_output=True, timeout=20,
    )
    assert (result.returncode == 0) == accepted, result.stdout + result.stderr
    if accepted:
        assert "status=ok" in result.stdout
    assert not output.exists(), "configuration checks must not start a simulation"


def test_native_check_requires_one_path(native_checker: Path, tmp_path: Path) -> None:
    result = subprocess.run(
        [str(native_checker), "--check-config"], cwd=tmp_path,
        text=True, capture_output=True, timeout=20,
    )
    assert result.returncode != 0
    assert "usage: beach --check-config beach.toml" in result.stderr
