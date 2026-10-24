from __future__ import annotations

import copy
from pathlib import Path

import pytest

from beach.config._layout import PATHS, to_grouped_layout

from beach.cli.main import main as beachx_main
from beach.config import (
    ConfigError,
    ConfigValidationError,
    default_config,
    load_config_file,
    normalize_config_document,
)


def _write_base_config(path: Path, *, field_bc_mode: str = "periodic2") -> None:
    path.write_text(
        f"""
[run]
batch_count = 2
rng_seed = 12345

[run.batch]
duration_steps = 10.0

[domain]
box_min = [
    0.0,
    0.0,
    0.0,
]
box_max = [
    1.0,
    1.0,
    10.0,
]
periodic_axes = [
    "x",
    "y",
]

[mesh]

[[mesh.parts]]
kind = "plane"
enabled = true
surface_side = "normal_plus"
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [
    0.5,
    0.5,
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
mode = "volume_seed"

[particles.species.sampling]
volume_macro_particles_per_batch = 10

[fields]
boundary = "{field_bc_mode}"

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


def _default_periodic2_config() -> dict[str, object]:
    config = default_config()
    config["domain"]["periodic_axes"] = ["x", "y"]
    config["field_boundary"]["mode"] = "periodic2"
    config["sim"].update(
        {
            "field_solver": "fmm",
            "field_periodic_image_layers": 1,
            "field_periodic_far_correction": "none",
        }
    )
    return config


def test_load_config_file_accepts_direct_beach_toml(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)

    result = load_config_file(config_path)

    assert result["field_boundary"]["mode"] == "periodic2"
    assert result["particles"]["species"][0]["npcls_per_step"] == 10
    assert result["mesh"]["parts"][0]["kind"] == "plane"


def test_declared_enum_normalization_preserves_leading_spaces() -> None:
    config = default_config()
    config["sim"]["field_solver"] = "DIRECT "
    assert normalize_config_document(to_grouped_layout(config))["sim"]["field_solver"] == "direct"
    config["sim"]["field_solver"] = " DIRECT"
    with pytest.raises(ConfigValidationError, match="method"):
        normalize_config_document(to_grouped_layout(config))


def test_default_config_uses_free_space_without_periodic_options() -> None:
    config = default_config()

    assert config["domain"]["periodic_axes"] == []
    assert config["field_boundary"]["mode"] == "free"
    assert "field_periodic_far_correction" not in config["sim"]
    assert "softening" not in config["sim"]
    assert "field" not in config


@pytest.mark.parametrize(
    "key",
    ["softening", "multiple_box_events_soft_discard_count_limit"],
)
def test_config_rejects_removed_sim_keys(key: str) -> None:
    config = default_config()
    config = to_grouped_layout(config)
    config["particles"]["tracking"][key] = 1.0

    with pytest.raises(ConfigValidationError, match=rf"Additional properties.*{key}"):
        normalize_config_document(to_grouped_layout(config))


def test_default_config_matches_official_tutorial_case() -> None:
    assert normalize_config_document(to_grouped_layout(default_config())) == load_config_file(
        Path("examples/tutorial_electron_beam.toml")
    )


def test_normalization_preserves_case_in_paths_and_species_identifiers() -> None:
    config = default_config()
    config["sim"]["field_solver"] = "DIRECT"
    config["output"]["dir"] = "Outputs/MixedCase"
    config["particles"]["species"][0]["species_key"] = "ElectronA"

    normalized = normalize_config_document(to_grouped_layout(config))

    assert normalized["sim"]["field_solver"] == "direct"
    assert config["sim"]["field_solver"] == "DIRECT"
    assert normalized["output"]["dir"] == "Outputs/MixedCase"
    assert normalized["particles"]["species"][0]["species_key"] == "ElectronA"


def test_separated_boundary_tables_preserve_the_local_source_contract() -> None:
    runtime = default_config()
    authoring = copy.deepcopy(runtime)
    authoring["reservoir"] = {"inflow_model": "infinity_barrier", "phi_infty": 1.0}

    normalized = normalize_config_document(to_grouped_layout(authoring))
    assert normalized["particle_boundary"] == authoring["particle_boundary"]
    assert normalized["reservoir"] == authoring["reservoir"]
    assert normalized["sim"] == runtime["sim"]


def test_load_config_file_resolves_high_level_notation(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    config_path.write_text(
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

[[mesh.parts]]
kind = "plane"
enabled = true
surface_side = "normal_plus"
size_mode = "box_fraction"
size_frac = [
    0.5,
    0.25,
]
nx = 20
ny = 20
placement_mode = "box_anchor"
anchor = "box_center"

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

    result = load_config_file(config_path)

    assert result["domain"]["box_min"] == [1.0, 2.0, 3.0]
    assert result["domain"]["box_max"] == [3.0, 6.0, 9.0]
    species = result["particles"]["species"][0]
    assert species["pos_low"] == [1.5, 4.0, 9.0]
    assert species["pos_high"] == [2.5, 6.0, 9.0]
    template = result["mesh"]["parts"][0]
    assert template["center"] == [2.0, 4.0, 6.0]
    assert template["size_x"] == pytest.approx(1.0)
    assert template["size_y"] == pytest.approx(1.0)


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


def test_load_config_file_rejects_unknown_top_level_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "unknown_keys.toml"
    config_path.write_text(
        """
schema_version = 1
use_presets = ["sim/periodic2_fmm"]

[run]
rng_seed = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigError, match="Additional properties"):
        load_config_file(config_path)


def test_config_rejects_unimplemented_dielectric_inputs() -> None:
    config = default_config()
    config["mesh"]["parts"][0]["surface_model"] = "dielectric"
    with pytest.raises(ConfigValidationError, match="surface_model.*dielectric"):
        normalize_config_document(to_grouped_layout(config))

    config = default_config()
    config["mesh"]["parts"][0]["epsilon_r"] = 3.9
    with pytest.raises(ConfigValidationError, match="epsilon_r"):
        normalize_config_document(to_grouped_layout(config))


@pytest.mark.parametrize(
    ("key", "value"),
    [
        # Range bounds and their activation: test_config_schema soft-discard bounds.
        ("multiple_box_events_soft_discard_count_grace", True),
        ("multiple_box_events_soft_discard_fraction_limit", True),
        ("multiple_box_events_soft_discard_fraction_limit", float("inf")),
        ("multiple_box_events_soft_discard_fraction_limit", float("nan")),
        ("multiple_box_events_soft_discard_abs_charge_limit", True),
        ("multiple_box_events_soft_discard_abs_charge_limit", float("nan")),
    ],
)
def test_active_soft_discard_controls_reject_invalid_values(
    key: str, value: object,
) -> None:
    config = default_config()
    config["sim"]["multiple_box_events_policy"] = "soft_discard"
    config["sim"][key] = value

    with pytest.raises(ConfigValidationError, match=PATHS["sim." + key].split(".")[-1]):
        normalize_config_document(to_grouped_layout(config))


def test_load_config_file_rejects_nonfinite_template_scalar(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "size_x = 1.0", "size_x = inf"
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(ConfigValidationError, match="mesh.parts\\[0\\].size_x"):
        load_config_file(config_path)


def test_load_config_file_rejects_removed_photo_escape_model(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    config_path.write_text(
        config_path.read_text(encoding="utf-8").replace(
            "volume_macro_particles_per_batch = 10",
            'volume_macro_particles_per_batch = 10\nphoto_escape_model = "boltzmann_cutoff"',
        ),
        encoding="utf-8",
    )

    with pytest.raises(ConfigValidationError, match="photo_escape_model"):
        load_config_file(config_path)


def test_zhao_stationary_accepts_a_nondrifting_electron_reservoir() -> None:
    root = Path(__file__).resolve().parents[2]
    config = load_config_file(root / "examples/zero_current.toml")
    config["particles"]["species"][0]["drift_velocity"] = [0.0, 0.0, 0.0]
    normalized = normalize_config_document(to_grouped_layout(config))
    assert normalized["particles"]["species"][0]["drift_velocity"] == [0.0, 0.0, 0.0]


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
    assert initialized["domain"]["periodic_axes"] == []
    assert initialized["sim"]["field_solver"] == "direct"
    assert initialized["field_boundary"]["mode"] == "free"
    assert initialized["sim"]["batch_count"] == 150
    assert initialized["sim"]["batch_duration"] == 1.0e-8
    assert len(initialized["particles"]["species"]) == 1
    species = initialized["particles"]["species"][0]
    assert species["source_mode"] == "volume_seed"
    assert species["npcls_per_step"] == 200
    assert species["w_particle"] == 4000.0
    assert species["pos_low"] == [0.008, 0.008, 0.06]
    assert species["pos_high"] == [0.032, 0.032, 0.06]
    assert species["drift_velocity"] == [0.0, 0.0, -1.3262e6]

    modified = tmp_path / "modified.toml"
    modified.write_text(
        (tmp_path / "beach.toml")
        .read_text(encoding="utf-8")
        .replace("batch_count = 150", "batch_count = 151"),
        encoding="utf-8",
    )
    beachx_main(["config", "diff", "beach.toml", str(modified)])
    diff_streams = capsys.readouterr()
    assert "sim.batch_count: 150 -> 151" in diff_streams.out


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


def test_config_cli_rejects_removed_migrate_command(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        beachx_main(["config", "migrate", "old.toml", "new.toml"])
    assert exc.value.code == 2
    assert "invalid choice: 'migrate'" in capsys.readouterr().err
