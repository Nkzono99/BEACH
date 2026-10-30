"""beachx config init/diff and beachx lint."""

from __future__ import annotations

from pathlib import Path

import pytest
from config_fixtures import (
    _write_base_config,
)

from beach.cli.main import main as beachx_main
from beach.cli.lint import resolve_native_checker, run_native_check
from beach.config import (
    load_config_file,
)


@pytest.mark.parametrize("use_environment", [False, True])
def test_lint_runs_relative_checker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, use_environment: bool,
) -> None:
    monkeypatch.chdir(tmp_path)
    checker = tmp_path / "local-checker"
    checker.write_text('#!/bin/sh\n[ "$1" = --check-config ] && printf "status=ok\\n"\n')
    checker.chmod(0o755)
    if use_environment:
        monkeypatch.setenv("BEACH_CONFIG_CHECK_EXE", "./local-checker")
    selected = resolve_native_checker(None if use_environment else Path("./local-checker"))
    run_native_check(selected, Path("input.toml"))


def _write_high_level_authoring_config(path: Path) -> None:
    path.write_text(
        """
[run]
batch_count = 2
rng_seed = 12345

[run.batch]
duration_steps = 10.0

[domain]
box_origin = [
    1.0,
    2.0,
    3.0,
]
box_size = [
    2.0,
    4.0,
    6.0,
]
periodic_axes = [
    "x",
    "y",
]

[mesh]

[mesh.groups.cavity_unit]
placement_mode = "box_anchor"
anchor = "box_center"
scale_from = "box_x"
scale_factor = 0.5

[[mesh.parts]]
group = "cavity_unit"
kind = "sphere"
surface_side = "outward_closed"
radius = 0.2
n_lon = 8
n_lat = 4
center_local = [
    0.0,
    0.0,
    0.0,
]

[particles.tracking]
dt_s = 2e-08
max_steps_per_particle = 100

[particles.boundary]
z_low = "open"
z_high = "open"

[[particles.species]]
charge_c = -1.602176634e-19
mass_kg = 9.10938356e-31

[particles.species.source]
mode = "photo_raycast"
emit_current_density_a_m2 = 0.001
inject_face = "z_high"
region_mode = "face_fraction"
uv_low = [
    0.25,
    0.5,
]
uv_high = [
    0.75,
    1.0,
]

[particles.species.sampling]
rays_per_batch = 8

[particles.species.distribution]
temperature_k = 0.0
drift_velocity_m_s = [
    0.0,
    0.0,
    -1.0,
]

[fields]
boundary = "periodic2"

[fields.solver]
method = "fmm"

[output]
enabled = true
dir = "outputs/latest"

[output.history]
stride_batches = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )


def test_config_cli_init_and_diff(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.chdir(tmp_path)

    beachx_main(["config", "init"])
    init_streams = capsys.readouterr()
    assert "saved=beach.toml" in init_streams.out
    assert (tmp_path / "beach.toml").exists()
    initialized_text = (tmp_path / "beach.toml").read_text(encoding="utf-8")

    with pytest.raises(SystemExit, match="config file already exists"):
        beachx_main(["config", "init"])
    assert (tmp_path / "beach.toml").read_text(encoding="utf-8") == initialized_text

    initialized = load_config_file(tmp_path / "beach.toml")
    assert initialized["domain"]["periodic_axes"] == ["x", "y"]
    assert initialized["field_boundary"]["mode"] == "periodic2"
    assert initialized["sim"]["batch_count"] == 10000
    species = initialized["particles"]["species"]
    assert [item["species_key"] for item in species] == ["electron", "proton"]
    assert all(item["boundary_inflow"] == {"z_high": "reservoir"} for item in species)

    modified = tmp_path / "modified.toml"
    modified.write_text(
        (tmp_path / "beach.toml")
        .read_text(encoding="utf-8")
        .replace("batch_count = 10000", "batch_count = 10001"),
        encoding="utf-8",
    )
    beachx_main(["config", "diff", "beach.toml", str(modified)])
    diff_streams = capsys.readouterr()
    assert "sim.batch_count: 10000 -> 10001" in diff_streams.out


def test_lint_cli_accepts_valid_config(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    native_checker: Path,
) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)

    beachx_main(["lint", str(config_path), "--beach", str(native_checker)])
    streams = capsys.readouterr()

    assert f"config={config_path}" in streams.out
    assert "schema=package:beach.config/schemas/beach.schema.json" in streams.out
    assert "checks=toml,schema,preflight" in streams.out
    assert "status=ok" in streams.out


def test_lint_cli_accepts_output_restart_from(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    native_checker: Path,
) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "stride_batches = 1",
        'stride_batches = 1\n[run.restart]\nfrom = "outputs/parent"',
    )
    config_path.write_text(text, encoding="utf-8")

    beachx_main(["lint", str(config_path), "--beach", str(native_checker)])
    streams = capsys.readouterr()

    assert "checks=toml,schema,preflight" in streams.out
    assert "status=ok" in streams.out


def test_lint_cli_accepts_high_level_authoring_config(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    native_checker: Path,
) -> None:
    config_path = tmp_path / "beach.toml"
    _write_high_level_authoring_config(config_path)

    beachx_main(["lint", str(config_path), "--beach", str(native_checker)])
    streams = capsys.readouterr()

    assert "checks=toml,schema,preflight" in streams.out
    assert "status=ok" in streams.out


def test_lint_cli_reports_toml_parse_error(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    config_path.write_text("[sim\n", encoding="utf-8")

    with pytest.raises(SystemExit, match="TOML parse error"):
        beachx_main(["lint", str(config_path)])


def test_lint_custom_schema_adds_constraints_without_replacing_beach_contract(
    tmp_path: Path, native_checker: Path,
) -> None:
    config_path = tmp_path / "beach.toml"
    schema_path = tmp_path / "extra.schema.json"
    _write_base_config(config_path)
    schema_path.write_text("{}", encoding="utf-8")
    beachx_main(["lint", str(config_path), "--schema", str(schema_path), "--beach", str(native_checker)])

    schema_path.write_text(
        '{"properties": {"run": {"properties": {"batch_count": {"maximum": 1}}}}}',
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="batch_count.*maximum"):
        beachx_main(["lint", str(config_path), "--schema", str(schema_path)])

    schema_path.write_text("{}", encoding="utf-8")
    config_path.write_text(
        config_path.read_text(encoding="utf-8").replace("batch_count = 2", "batch_count = true"),
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="batch_count.*integer"):
        beachx_main(["lint", str(config_path), "--schema", str(schema_path)])


def test_lint_cli_reports_authoring_schema_error(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_high_level_authoring_config(config_path)
    config_path.write_text(
        config_path.read_text(encoding="utf-8").replace(
            "scale_factor", "scale_factorr"
        ),
        encoding="utf-8",
    )

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["lint", str(config_path)])

    message = str(excinfo.value)
    assert "schema validation failed" in message
    assert "schema phase=authoring" in message
    assert "schema error at mesh.groups.cavity_unit" in message


def test_lint_cli_reports_schema_error(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "enabled = true",
        'enabled = "yes"',
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["lint", str(config_path)])

    message = str(excinfo.value)
    assert "schema validation failed" in message
    assert "schema error at output.enabled" in message
