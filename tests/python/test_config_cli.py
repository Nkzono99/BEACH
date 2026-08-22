from __future__ import annotations

import copy
from pathlib import Path

import pytest

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
[sim]
dt = 2.0e-8
batch_duration_step = 10.0
batch_count = 2
max_step = 100
rng_seed = 12345
field_solver = "fmm"

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "{field_bc_mode}"

[particle_boundary]
z_low = "open"
z_high = "open"

[particles]
[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 10

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_side = "normal_plus"
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.0]

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )


def test_load_config_file_accepts_direct_beach_toml(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)

    result = load_config_file(config_path)

    assert result["field_boundary"]["mode"] == "periodic2"
    assert result["particles"]["species"][0]["npcls_per_step"] == 10
    assert result["mesh"]["templates"][0]["kind"] == "plane"


def test_default_config_uses_no_periodic_far_correction() -> None:
    config = default_config()

    assert config["sim"]["field_periodic_far_correction"] == "none"
    assert "softening" not in config["sim"]
    assert "field" not in config


def test_config_rejects_removed_point_source_options() -> None:
    config = default_config()
    config["sim"]["softening"] = 1.0e-6
    with pytest.raises(ConfigValidationError, match="removed sim key.*softening"):
        normalize_config_document(config)

    config = default_config()
    config["field"] = {"element_kernel": "point"}
    with pytest.raises(ConfigError, match="unsupported top-level key.*field"):
        normalize_config_document(config)


def test_explicit_auto_mesh_requires_obj_surface_side() -> None:
    config = default_config()
    config["mesh"]["mode"] = "auto"

    with pytest.raises(ConfigValidationError, match="mesh.surface_side"):
        normalize_config_document(config)

    config["mesh"]["surface_side"] = "normal_plus"
    assert normalize_config_document(config)["mesh"]["mode"] == "auto"


def test_default_config_matches_official_tutorial_case() -> None:
    assert normalize_config_document(default_config()) == load_config_file(
        Path("examples/tutorial_insulator.toml")
    )


def test_periodic2_accepts_symmetric_vacuum_and_rejects_unknown_lower_model() -> None:
    config = load_config_file(Path("examples/periodic2_closed_photoelectron.toml"))
    config["periodic2"] = {}
    config["periodic2"]["lower_boundary_model"] = "symmetric_vacuum"

    normalized = normalize_config_document(config)

    assert normalized["periodic2"]["lower_boundary_model"] == "symmetric_vacuum"
    config["periodic2"]["lower_boundary_model"] = "unknown"
    with pytest.raises(ConfigValidationError, match="lower_boundary_model"):
        normalize_config_document(config)


def test_adaptive_nonzero_mode_requires_cached_time_scaled_sources() -> None:
    config = load_config_file(Path("examples/periodic2_closed_photoelectron.toml"))
    config["sim"]["field_solver"] = "fmm"
    config["sim"]["field_periodic_far_correction"] = "cached_kneq0"
    config["sim"]["field_periodic_image_layers"] = 1
    config["periodic2"] = {}
    config["periodic2"]["nonzero_mode_backend"] = "cached_kneq0"
    config["periodic2"]["max_nonzero_mode_potential_step"] = 1.0e-2

    normalized = normalize_config_document(config)
    assert normalized["periodic2"]["max_nonzero_mode_potential_step"] == 1.0e-2

    fixed_weight = copy.deepcopy(config)
    fixed_weight["particles"]["species"][0]["w_particle"] = 1.0
    fixed_weight["particles"]["species"][0].pop(
        "target_macro_particles_per_batch", None
    )
    with pytest.raises(ConfigValidationError, match="target_macro_particles_per_batch"):
        normalize_config_document(fixed_weight)

    spectral = copy.deepcopy(config)
    spectral["sim"]["field_solver"] = "direct"
    spectral["sim"]["field_periodic_far_correction"] = "none"
    spectral["periodic2"]["nonzero_mode_backend"] = "panel_spectral_reference"
    spectral["periodic2"]["zero_mode_policy"] = "exclude_k0"
    spectral["periodic2"]["lower_boundary_model"] = "e_bottom_zero"
    with pytest.raises(ConfigValidationError, match="cached_kneq0"):
        normalize_config_document(spectral)

    volume = default_config()
    volume["sim"]["batch_duration"] = 1.0e-6
    volume["sim"]["field_periodic_far_correction"] = "cached_kneq0"
    volume["periodic2"] = {
        "nonzero_mode_backend": "cached_kneq0",
        "zero_mode_policy": "exclude_k0",
        "lower_boundary_model": "e_bottom_zero",
        "max_nonzero_mode_potential_step": 1.0e-2,
    }
    with pytest.raises(ConfigValidationError, match="time-scaled"):
        normalize_config_document(volume)


def test_separated_boundary_tables_preserve_the_local_source_contract() -> None:
    runtime = default_config()
    authoring = copy.deepcopy(runtime)
    authoring["particle_boundary"]["ordinary_open_model"] = "potential_barrier"
    authoring["reservoir"] = {"inflow_model": "infinity_barrier", "phi_infty": 1.0}

    normalized = normalize_config_document(authoring)
    assert normalized["particle_boundary"] == authoring["particle_boundary"]
    assert normalized["reservoir"] == authoring["reservoir"]
    assert normalized["sim"] == runtime["sim"]


def test_particle_boundary_overrides_resolve_after_global_defaults() -> None:
    config = load_config_file(Path("examples/periodic2_closed_photoelectron.toml"))
    photoelectron = config["particles"]["species"][-1]
    photoelectron["boundary"]["z_high"] = "inherit"
    photoelectron.pop("q_particle")
    for inflow_species in config["particles"]["species"][:-1]:
        inflow_species.setdefault("boundary", {})["z_high"] = "open"
    config["particle_boundary"]["z_high"] = "reflect"

    normalized = normalize_config_document(config)
    assert normalized["particle_boundary"]["z_high"] == "reflect"
    assert normalized["particles"]["species"][-1]["boundary"]["z_high"] == "inherit"

    config["particle_boundary"]["z_high"] = "redistributed_reflect"
    normalized = normalize_config_document(config)
    assert normalized["particle_boundary"]["z_high"] == "redistributed_reflect"

    config["particle_boundary"]["z_high"] = "open"
    with pytest.raises(ConfigValidationError, match="reflecting action on inject_face"):
        normalize_config_document(config)


def test_particle_boundary_cannot_override_periodic_topology() -> None:
    config = default_config()
    config["particle_boundary"]["x_low"] = "reflect"
    with pytest.raises(ConfigValidationError, match="periodic domain face"):
        normalize_config_document(config)

    config = default_config()
    config["particles"]["species"][0]["boundary"] = {"x_high": "open"}
    with pytest.raises(ConfigValidationError, match="periodic domain face"):
        normalize_config_document(config)


def test_boundary_inflow_accepts_nonperiodic_faces_and_rejects_periodic_faces() -> None:
    config = default_config()
    config["sim"]["batch_duration"] = 1.0e-6
    species = config["particles"]["species"][0]
    species["number_density_m3"] = 1.0e6
    species["temperature_k"] = 2.0e4
    species["boundary_inflow"] = {
        "z_low": "reservoir",
        "z_high": "reservoir",
    }

    normalized = normalize_config_document(config)
    assert normalized["particles"]["species"][0]["boundary_inflow"] == {
        "z_low": "reservoir",
        "z_high": "reservoir",
    }

    reflecting = copy.deepcopy(config)
    reflecting["particles"]["species"][0]["boundary"] = {"z_high": "reflect"}
    with pytest.raises(ConfigValidationError, match="requires an open"):
        normalize_config_document(reflecting)

    config["particles"]["species"][0]["boundary_inflow"]["x_low"] = "reservoir"
    with pytest.raises(ConfigValidationError, match="periodic domain face"):
        normalize_config_document(config)


def test_boundary_inflow_requires_reservoir_physics_and_volume_source_mode() -> None:
    config = default_config()
    config["sim"]["batch_duration"] = 1.0e-6
    species = config["particles"]["species"][0]
    species["boundary_inflow"] = {"z_high": "reservoir"}

    with pytest.raises(ConfigValidationError, match="number_density"):
        normalize_config_document(config)

    species["number_density_cm3"] = 5.0
    photo = copy.deepcopy(config)
    photo_species = photo["particles"]["species"][0]
    photo_species["source_mode"] = "photo_raycast"
    photo_species["emit_current_density_a_m2"] = 1.0
    photo_species["rays_per_batch"] = 1
    photo_species["inject_face"] = "z_high"
    with pytest.raises(ConfigValidationError, match="cannot combine"):
        normalize_config_document(photo)

    species["source_mode"] = "plane_source"
    species["source_normal"] = [0.0, 0.0, -1.0]
    species["pos_low"] = [0.1, 0.1, 0.5]
    species["pos_high"] = [0.9, 0.9, 0.5]
    species.pop("npcls_per_step")
    with pytest.raises(ConfigValidationError, match="cannot combine"):
        normalize_config_document(config)


def test_plane_source_requires_internal_rectangle_and_axis_normal() -> None:
    config = default_config()
    config["sim"]["batch_duration"] = 1.0e-6
    species = config["particles"]["species"][0]
    species["source_mode"] = "plane_source"
    species.pop("npcls_per_step")
    species["number_density_m3"] = 1.0e6
    species["temperature_k"] = 2.0e4
    species["pos_low"] = [0.1, 0.1, 0.5]
    species["pos_high"] = [0.9, 0.9, 0.5]
    species["source_normal"] = [0.0, 0.0, -1.0]

    normalized = normalize_config_document(config)
    assert normalized["particles"]["species"][0]["source_normal"] == [
        0.0,
        0.0,
        -1.0,
    ]

    species["pos_low"][2] = 0.0
    species["pos_high"][2] = 0.0
    with pytest.raises(ConfigValidationError, match="strictly inside"):
        normalize_config_document(config)

    species["pos_low"][2] = 0.5
    species["pos_high"][2] = 0.5
    species["source_normal"] = [0.0, 0.5, -1.0]
    with pytest.raises(ConfigValidationError, match="non-zero axis-aligned"):
        normalize_config_document(config)

    species["source_normal"] = [0.0, 0.0, -2.0]
    assert normalize_config_document(config)["particles"]["species"][0][
        "source_normal"
    ] == [0.0, 0.0, -2.0]


def test_runtime_validator_rejects_removed_external_boundary_contract() -> None:
    config = default_config()
    config["external_boundary"] = {
        "field": {"model": "none"},
        "particles": {"mode": "local_source"},
    }

    with pytest.raises(ConfigError, match="unsupported top-level key"):
        normalize_config_document(config)


def test_load_config_file_resolves_high_level_notation(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    config_path.write_text(
        """
[sim]
dt = 2.0e-8
batch_duration_step = 10.0
batch_count = 2
max_step = 100
rng_seed = 12345
field_solver = "fmm"

[domain]
box_origin = [1.0, 2.0, 3.0]
box_size = [2.0, 4.0, 6.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"

[particles]
[[particles.species]]
source_mode = "reservoir_face"
number_density_m3 = 1.0
temperature_k = 0.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
w_particle = 1.0
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.25, 0.5]
uv_high = [0.75, 1.0]
drift_velocity = [0.0, 0.0, -1.0]

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_side = "normal_plus"
size_mode = "box_fraction"
size_frac = [0.5, 0.25]
nx = 20
ny = 20
placement_mode = "box_anchor"
anchor = "box_center"

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
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
    template = result["mesh"]["templates"][0]
    assert template["center"] == [2.0, 4.0, 6.0]
    assert template["size_x"] == pytest.approx(1.0)
    assert template["size_y"] == pytest.approx(1.0)


def _write_high_level_authoring_config(path: Path) -> None:
    path.write_text(
        """
[sim]
dt = 2.0e-8
batch_duration_step = 10.0
batch_count = 2
max_step = 100
rng_seed = 12345
field_solver = "fmm"

[domain]
box_origin = [1.0, 2.0, 3.0]
box_size = [2.0, 4.0, 6.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"

[particles]
[[particles.species]]
source_mode = "reservoir_face"
number_density_m3 = 1.0
temperature_k = 0.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
w_particle = 1.0
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.25, 0.5]
uv_high = [0.75, 1.0]
drift_velocity = [0.0, 0.0, -1.0]

[mesh]
mode = "template"

[mesh.groups.cavity_unit]
placement_mode = "box_anchor"
anchor = "box_center"
scale_from = "box_x"
scale_factor = 0.5

[[mesh.templates]]
group = "cavity_unit"
kind = "sphere"
surface_side = "outward_closed"
radius = 0.2
n_lon = 8
n_lat = 4
center_local = [0.0, 0.0, 0.0]

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )


def test_load_config_file_rejects_legacy_composition_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "legacy_composition.toml"
    config_path.write_text(
        """
schema_version = 1
use_presets = ["sim/periodic2_fmm"]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigError, match="reserved top-level key"):
        load_config_file(config_path)


def test_load_config_file_rejects_conductor_with_periodic2(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8")
    text = text.replace(
        "center = [0.5, 0.5, 0.0]",
        'center = [0.5, 0.5, 0.0]\nsurface_model = "conductor"',
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(ConfigValidationError, match='surface_model="conductor"'):
        load_config_file(config_path)


def test_config_rejects_unimplemented_dielectric_inputs() -> None:
    config = default_config()
    config["mesh"]["templates"][0]["surface_model"] = "dielectric"
    with pytest.raises(ConfigValidationError, match="dielectric.*not implemented"):
        normalize_config_document(config)

    config = default_config()
    config["mesh"]["templates"][0]["epsilon_r"] = 3.9
    with pytest.raises(ConfigValidationError, match="epsilon_r was removed"):
        normalize_config_document(config)


def test_inactive_sim_controls_do_not_enforce_backend_specific_bounds() -> None:
    config = default_config()
    config["sim"].update(
        {
            "field_solver": "direct",
            "field_periodic_far_correction": "none",
            "field_periodic_generation_tolerance": -1.0,
            "field_periodic_cache_dir": "",
            "tree_theta": -1.0,
            "tree_leaf_max": 0,
            "tree_min_nelem": 0,
            "multiple_box_events_policy": "abort",
            "multiple_box_events_soft_discard_count_limit": 0,
            "multiple_box_events_soft_discard_abs_charge_limit": -1.0,
            "raycast_max_bounce": 0,
        }
    )
    config["field_boundary"]["mode"] = "free"

    normalized = normalize_config_document(config)

    assert normalized["sim"]["field_solver"] == "direct"


def test_load_config_file_rejects_nonfinite_template_scalar(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "size_x = 1.0", "size_x = inf"
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(ConfigValidationError, match="mesh.templates\\[1\\].size_x"):
        load_config_file(config_path)


def test_load_config_file_rejects_removed_photo_escape_model(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    config_path.write_text(
        """
[sim]
batch_count = 1
batch_duration = 1.0e-6

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]
periodic_axes = []

[particles]
[[particles.species]]
source_mode = "photo_raycast"
emit_current_density_a_m2 = 1.0e-3
rays_per_batch = 10
deposit_opposite_charge_on_emit = true
photo_escape_model = "boltzmann_cutoff"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
temperature_ev = 1.5
inject_face = "z_high"
pos_low = [0.0, 0.0, 1.0]
pos_high = [1.0, 1.0, 1.0]

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_side = "normal_plus"
size_x = 1.0
size_y = 1.0
nx = 1
ny = 1
center = [0.5, 0.5, 0.0]

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigValidationError, match="photo_escape_model"):
        load_config_file(config_path)


def test_fixed_current_requires_signed_independent_channels() -> None:
    config = default_config()
    species = config["particles"]["species"][0]
    species["surface_charge_closure"] = "fixed_current"
    species["target_absorbed_current_a"] = -2.0e-6
    normalized = normalize_config_document(config)
    assert normalized["particles"]["species"][0][
        "target_absorbed_current_a"
    ] == pytest.approx(-2.0e-6)

    wrong_sign = copy.deepcopy(config)
    wrong_sign["particles"]["species"][0]["target_absorbed_current_a"] = 2.0e-6
    with pytest.raises(ConfigValidationError, match="same sign as q_particle"):
        normalize_config_document(wrong_sign)

    net_only = copy.deepcopy(config)
    net_only["particles"]["species"][0].pop("target_absorbed_current_a")
    with pytest.raises(ConfigValidationError, match="at least one target current"):
        normalize_config_document(net_only)


def test_zhao_stationary_model_supplies_fixed_current_targets() -> None:
    root = Path(__file__).resolve().parents[2]
    normalized = load_config_file(root / "examples/periodic2_zhao_fixed_current.toml")
    assert normalized["surface_current_model"]["model"] == "zhao_stationary"
    assert all(
        item["surface_charge_closure"] == "fixed_current"
        for item in normalized["particles"]["species"]
    )

    hot_ions = copy.deepcopy(normalized)
    hot_ions["particles"]["species"][1]["temperature_ev"] = 2.0
    with pytest.raises(ConfigValidationError, match="cold ions"):
        normalize_config_document(hot_ions)

    malformed_drift = copy.deepcopy(normalized)
    malformed_drift["particles"]["species"][0]["drift_velocity"] = [0.0, 0.0]
    with pytest.raises(ConfigValidationError, match="inward z-high"):
        normalize_config_document(malformed_drift)

    reflected_photoelectrons = copy.deepcopy(normalized)
    reflected_photoelectrons["particles"]["species"][2]["boundary"] = {
        "z_high": "reflect"
    }
    with pytest.raises(ConfigValidationError, match="open z-high"):
        normalize_config_document(reflected_photoelectrons)

    no_photo = load_config_file(
        root / "examples/periodic2_zhao_no_photo_fixed_current.toml"
    )
    assert no_photo["surface_current_model"]["photoelectron_source_scale"] == 0.0
    assert len(no_photo["particles"]["species"]) == 2

    explicit_no_photo_type_c = copy.deepcopy(no_photo)
    explicit_no_photo_type_c["surface_current_model"]["zhao_branch"] = "c"
    normalized_no_photo_type_c = normalize_config_document(explicit_no_photo_type_c)
    assert normalized_no_photo_type_c["surface_current_model"]["zhao_branch"] == "c"

    stale_photo_setting = copy.deepcopy(no_photo)
    stale_photo_setting["surface_current_model"]["solar_elevation_deg"] = 60.0
    with pytest.raises(ConfigValidationError, match="omitting all photoelectron"):
        normalize_config_document(stale_photo_setting)

    invalid_no_photo_branch = copy.deepcopy(no_photo)
    invalid_no_photo_branch["surface_current_model"]["zhao_branch"] = "a"
    with pytest.raises(ConfigValidationError, match="zhao_branch.*auto.*c"):
        normalize_config_document(invalid_no_photo_branch)


def test_config_cli_init_validate_and_diff(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.chdir(tmp_path)

    beachx_main(["config", "init"])
    init_streams = capsys.readouterr()
    assert "saved=beach.toml" in init_streams.out
    assert (tmp_path / "beach.toml").exists()

    beachx_main(["config", "validate"])
    validate_streams = capsys.readouterr()
    assert "status=ok" in validate_streams.out

    initialized = load_config_file(tmp_path / "beach.toml")
    assert initialized["domain"]["periodic_axes"] == ["x", "y"]
    assert initialized["sim"]["field_solver"] == "fmm"
    assert initialized["field_boundary"]["mode"] == "periodic2"
    assert initialized["sim"]["field_periodic_image_layers"] == 1
    assert initialized["sim"]["field_periodic_far_correction"] == "none"
    assert len(initialized["particles"]["species"]) == 1
    species = initialized["particles"]["species"][0]
    assert species["source_mode"] == "volume_seed"
    assert species["npcls_per_step"] == 1
    assert species["pos_low"] == [0.5, 0.5, 0.8]
    assert species["pos_high"] == [0.5, 0.5, 0.8]
    assert species["drift_velocity"] == [0.0, 0.0, -1.0e6]

    modified = tmp_path / "modified.toml"
    modified.write_text(
        (tmp_path / "beach.toml")
        .read_text(encoding="utf-8")
        .replace("batch_count = 1", "batch_count = 2"),
        encoding="utf-8",
    )
    beachx_main(["config", "diff", "beach.toml", str(modified)])
    diff_streams = capsys.readouterr()
    assert "sim.batch_count: 1 -> 2" in diff_streams.out


def test_lint_cli_accepts_valid_config(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)

    beachx_main(["lint", str(config_path)])
    streams = capsys.readouterr()

    assert f"config={config_path}" in streams.out
    assert "schema=package:beach.config/schemas/beach.schema.json" in streams.out
    assert "checks=toml,schema,semantic" in streams.out
    assert "status=ok" in streams.out


def test_lint_cli_accepts_output_restart_from(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "history_stride = 1",
        'history_stride = 1\nresume = true\nrestart_from = "outputs/parent"',
    )
    config_path.write_text(text, encoding="utf-8")

    beachx_main(["lint", str(config_path)])
    streams = capsys.readouterr()

    assert "checks=toml,schema,semantic" in streams.out
    assert "status=ok" in streams.out


def test_lint_cli_accepts_high_level_authoring_config(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    config_path = tmp_path / "beach.toml"
    _write_high_level_authoring_config(config_path)

    beachx_main(["lint", str(config_path)])
    streams = capsys.readouterr()

    assert "checks=toml,schema,semantic" in streams.out
    assert "status=ok" in streams.out


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
        "write_files = true",
        'write_files = "yes"',
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["lint", str(config_path)])

    message = str(excinfo.value)
    assert "schema validation failed" in message
    assert "schema error at output.write_files" in message


def test_lint_cli_rejects_restart_from_without_resume(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "history_stride = 1",
        'history_stride = 1\nrestart_from = "outputs/parent"',
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["lint", str(config_path)])

    message = str(excinfo.value)
    assert "schema validation failed" in message
    assert "resume" in message


def test_lint_cli_reports_semantic_error(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "center = [0.5, 0.5, 0.0]",
        'center = [0.5, 0.5, 0.0]\nsurface_model = "conductor"',
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["lint", str(config_path)])

    assert 'surface_model="conductor"' in str(excinfo.value)
