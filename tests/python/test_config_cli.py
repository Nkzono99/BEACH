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
    validate_runtime_config,
)
from beach.config._toml import load_toml_file


def _write_base_config(path: Path, *, field_bc_mode: str = "periodic2") -> None:
    path.write_text(
        f"""
[sim]
dt = 2.0e-8
batch_duration_step = 10.0
batch_count = 2
max_step = 100
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "{field_bc_mode}"

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

    assert result["sim"]["field_bc_mode"] == "periodic2"
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


def test_load_config_file_accepts_kinetic_profile_return_split() -> None:
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))

    assert config["periodic2"]["nonzero_mode_backend"] == "panel_spectral_reference"
    assert config["outer_plasma"]["model"] == "kinetic_1d"
    assert config["outer_plasma"]["return_model"] == "kinetic_1d_profile_return"
    assert (
        config["coupling"]["particle_transfer_mode"]
        == "electrostatic_1d_instant_return"
    )


def test_periodic2_accepts_symmetric_vacuum_and_rejects_unknown_lower_model() -> None:
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))
    config["periodic2"]["lower_boundary_model"] = "symmetric_vacuum"

    normalized = normalize_config_document(config)

    assert normalized["periodic2"]["lower_boundary_model"] == "symmetric_vacuum"
    config["periodic2"]["lower_boundary_model"] = "unknown"
    with pytest.raises(ConfigValidationError, match="lower_boundary_model"):
        normalize_config_document(config)


def test_photoelectron_settings_reject_unknown_density_and_nonconserving_modes() -> (
    None
):
    config = load_config_file(
        Path("examples/periodic2_zhao_charge_driven_outer.toml")
    )
    config["outer_plasma"]["photoelectron_density_model"] = "statistical_return"
    with pytest.raises(ConfigValidationError, match="statistical_return"):
        normalize_config_document(config)

    config["outer_plasma"]["photoelectron_density_model"] = "none"
    photoelectron = next(
        species
        for species in config["particles"]["species"]
        if species.get("source_mode") == "photo_raycast"
    )
    photoelectron["deposit_opposite_charge_on_emit"] = False
    with pytest.raises(ConfigValidationError, match="deposit_opposite_charge_on_emit"):
        normalize_config_document(config)


def test_photoelectron_settings_reject_removed_closure_key() -> None:
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))
    config["outer_plasma"]["photoelectron_closure"] = "individual_return"

    with pytest.raises(
        ConfigValidationError, match="photoelectron_closure was removed"
    ):
        normalize_config_document(config)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("photoelectron_histogram_enabled", True),
        ("photoelectron_histogram_bins", 16),
        ("photoelectron_histogram_energy_max", 10.0),
        ("photoelectron_ambient_charge_scale", 1.0),
        ("max_photoelectron_charge_ratio", 0.1),
    ],
)
def test_runtime_rejects_removed_photoelectron_histogram_controls(
    key: str, value: object
) -> None:
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))
    config["outer_plasma"][key] = value

    with pytest.raises(ConfigValidationError, match=key):
        normalize_config_document(config)


def test_kinetic_photoelectron_density_model_requires_matching_transfer_pair() -> None:
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))
    zhao_config = load_config_file(
        Path("examples/periodic2_zhao_charge_driven_outer.toml")
    )
    photoelectron = next(
        copy.deepcopy(species)
        for species in zhao_config["particles"]["species"]
        if species.get("source_mode") == "photo_raycast"
    )
    photoelectron["pos_low"] = [0.0, 0.0, 1.0]
    photoelectron["pos_high"] = [1.0, 1.0, 1.0]
    config["particles"]["species"].append(photoelectron)
    config["outer_plasma"]["photoelectron_density_model"] = "kinetic_mean"

    normalize_config_document(config)

    config["coupling"]["particle_transfer_mode"] = "none"
    with pytest.raises(ConfigValidationError, match="particle_transfer_mode"):
        normalize_config_document(config)

    config["coupling"]["particle_transfer_mode"] = "electrostatic_1d_instant_return"
    config["outer_plasma"]["return_model"] = "none"
    with pytest.raises(ConfigValidationError, match="particle_transfer_mode"):
        normalize_config_document(config)

    config["coupling"]["particle_transfer_mode"] = "none"
    normalize_config_document(config)


def test_zhao_charge_driven_closure_constraints() -> None:
    config = load_config_file(
        Path("examples/periodic2_zhao_charge_driven_outer.toml")
    )
    config["outer_plasma"]["zhao_branch"] = "c"

    normalize_config_document(config)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["model"] = "linear_debye"
    with pytest.raises(ConfigValidationError, match="outer_plasma.model"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["photoelectron_density_model"] = "kinetic_mean"
    with pytest.raises(ConfigValidationError, match="photoelectron_density_model"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["infinity_potential"] = 1.0
    with pytest.raises(ConfigValidationError, match="infinity_potential"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["sim"]["sheath_injection_model"] = "zhao_a"
    with pytest.raises(ConfigValidationError, match="sheath_injection_model"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["sim"]["reservoir_potential_model"] = "infinity_barrier"
    with pytest.raises(ConfigValidationError, match="reservoir_potential_model"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["sim"]["sheath_reference_coordinate"] = 0.0
    with pytest.raises(ConfigValidationError, match="sheath_reference_coordinate"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["sim"]["sheath_photoelectron_ref_density_cm3"] = 0.0
    with pytest.raises(
        ConfigValidationError, match="sheath_photoelectron_ref_density_cm3"
    ):
        normalize_config_document(invalid)

    no_photo = copy.deepcopy(config)
    no_photo["outer_plasma"]["photoelectron_source_scale"] = 0.0
    no_photo["sim"]["sheath_photoelectron_ref_density_cm3"] = 0.0
    for species in no_photo["particles"]["species"]:
        if species.get("source_mode") == "photo_raycast":
            species["enabled"] = False
    normalize_config_document(no_photo)

    invalid = copy.deepcopy(no_photo)
    invalid["sim"]["sheath_photoelectron_ref_density_cm3"] = -1.0
    with pytest.raises(
        ConfigValidationError, match="sheath_photoelectron_ref_density_cm3"
    ):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["photoelectron_source_scale"] = -1.0
    with pytest.raises(ConfigValidationError, match="photoelectron_source_scale"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["kinetic_closure"] = "absorbing_maxwellian"
    invalid["outer_plasma"]["zhao_branch"] = "auto"
    invalid["outer_plasma"]["photoelectron_source_scale"] = 0.0
    with pytest.raises(ConfigValidationError, match="photoelectron_source_scale"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["kinetic_closure"] = "absorbing_maxwellian"
    with pytest.raises(ConfigValidationError, match="zhao_branch"):
        normalize_config_document(invalid)


def test_external_boundary_contract_rejects_competing_owners_and_mismatched_pairs() -> (
    None
):
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))
    config["sim"]["open_boundary_model"] = "potential_barrier"
    normalize_config_document(config)

    invalid = copy.deepcopy(config)
    invalid["sim"]["reservoir_potential_model"] = "infinity_barrier"
    with pytest.raises(ConfigValidationError, match="owns inflow"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["return_model"] = "none"
    with pytest.raises(ConfigValidationError, match="matching"):
        normalize_config_document(invalid)


def test_external_boundary_facade_lowers_to_the_runtime_contract() -> None:
    authoring = load_toml_file(Path("examples/periodic2_kinetic_outer.toml"))
    external_boundary = authoring["external_boundary"]

    runtime = copy.deepcopy(authoring)
    del runtime["external_boundary"]
    runtime["outer_plasma"] = {
        **external_boundary["field"],
        "return_model": "kinetic_1d_profile_return",
        "interface_z": runtime["sim"]["box_max"][2],
    }
    runtime["coupling"] = {
        "update_mode": "explicit",
        "particle_transfer_mode": "electrostatic_1d_instant_return",
        "field_evolution_timescale": external_boundary["particles"][
            "field_evolution_timescale"
        ],
        "outer_queue_enabled": False,
    }

    normalized = normalize_config_document(authoring)

    assert "external_boundary" not in normalized
    assert normalized == normalize_config_document(runtime)


def test_external_boundary_facade_preserves_the_simple_default_contract() -> None:
    runtime = default_config()
    authoring = copy.deepcopy(runtime)
    authoring["external_boundary"] = {
        "field": {"model": "none"},
        "particles": {"mode": "local_source"},
    }

    assert normalize_config_document(authoring) == runtime


def test_runtime_validator_rejects_unlowered_external_boundary_facade() -> None:
    config = default_config()
    config["external_boundary"] = {
        "field": {"model": "none"},
        "particles": {"mode": "local_source"},
    }

    with pytest.raises(ConfigValidationError, match="authoring-only"):
        validate_runtime_config(config)


def test_external_boundary_facade_rejects_inert_none_options() -> None:
    config = default_config()
    config["external_boundary"] = {
        "field": {"model": "none", "debye_length": 1.0},
        "particles": {"mode": "local_source"},
    }
    with pytest.raises(ConfigError, match="does not accept outer-field"):
        normalize_config_document(config)

    config["external_boundary"]["field"] = {"model": "none"}
    config["external_boundary"]["particles"]["outer_update_stride"] = 1
    with pytest.raises(ConfigError, match="does not accept.*coupling"):
        normalize_config_document(config)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("infinity_potential", 0.0),
        ("photoelectron_histogram_enabled", False),
        ("photoelectron_histogram_bins", 16),
        ("photoelectron_histogram_energy_max", 10.0),
        ("photoelectron_ambient_charge_scale", 1.0),
        ("max_photoelectron_charge_ratio", 0.1),
        ("unified_grid_points", 129),
        ("accessible_fraction_tolerance", 0.1),
        ("max_linearity_ratio", 0.25),
    ],
)
def test_external_boundary_facade_rejects_removed_field_options(
    key: str, value: object
) -> None:
    config = load_toml_file(Path("examples/periodic2_kinetic_outer.toml"))
    config["external_boundary"]["field"][key] = value

    with pytest.raises(ConfigError, match="unsupported external_boundary.field"):
        normalize_config_document(config)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("outer_orbit_dt", 1.0e-9),
        ("outer_orbit_max_steps", 100),
        ("outer_orbit_energy_tolerance", 1.0e-4),
    ],
)
def test_external_boundary_facade_rejects_removed_particle_options(
    key: str, value: object
) -> None:
    config = load_toml_file(Path("examples/periodic2_kinetic_outer.toml"))
    config["external_boundary"]["particles"][key] = value

    with pytest.raises(ConfigError, match="unsupported external_boundary.particles"):
        normalize_config_document(config)


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("outer_plasma", "unified_grid_points", 129),
        ("outer_plasma", "accessible_fraction_tolerance", 0.1),
        ("outer_plasma", "max_linearity_ratio", 0.25),
        ("coupling", "outer_orbit_dt", 1.0e-9),
        ("coupling", "outer_orbit_max_steps", 100),
        ("coupling", "outer_orbit_energy_tolerance", 1.0e-4),
    ],
)
def test_runtime_config_rejects_removed_unified_options(
    section: str, key: str, value: object
) -> None:
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))
    config[section][key] = value

    with pytest.raises(ConfigValidationError, match=f"unsupported {section} key"):
        normalize_config_document(config)


def test_external_boundary_facade_requires_specialized_field_option_owners() -> None:
    zhao = load_toml_file(Path("examples/periodic2_kinetic_outer.toml"))
    zhao["external_boundary"]["field"]["zhao_branch"] = "auto"
    with pytest.raises(ConfigError, match="require kinetic_closure"):
        normalize_config_document(zhao)

    zhao_density = load_toml_file(
        Path("examples/periodic2_zhao_charge_driven_outer.toml")
    )
    zhao_density["external_boundary"]["field"]["photoelectron_density_model"] = "none"
    with pytest.raises(ConfigError, match="requires kinetic_closure=absorbing"):
        normalize_config_document(zhao_density)


def test_ambient_linear_photo_facade_keeps_photoelectrons_tracked_only() -> None:
    config = load_toml_file(
        Path("examples/periodic2_ambient_linear_photo_outer.toml")
    )
    normalized = normalize_config_document(config)

    assert normalized["outer_plasma"]["kinetic_closure"] == (
        "ambient_linear_debye"
    )
    assert normalized["outer_plasma"].get("photoelectron_density_model", "none") == "none"
    assert normalized["outer_plasma"]["return_model"] == (
        "kinetic_1d_profile_return"
    )

    invalid = copy.deepcopy(config)
    invalid["external_boundary"]["field"]["photoelectron_density_model"] = (
        "kinetic_mean"
    )
    with pytest.raises(ConfigError, match="requires kinetic_closure=absorbing"):
        normalize_config_document(invalid)


@pytest.mark.parametrize(
    ("example", "key", "value"),
    [
        ("periodic2_zhao_transient_outer.toml", "outer_update_stride", 1),
    ],
)
def test_external_boundary_facade_rejects_mode_specific_particle_options(
    example: str,
    key: str,
    value: object,
) -> None:
    config = load_toml_file(Path("examples") / example)
    config["external_boundary"]["particles"][key] = value

    with pytest.raises(ConfigError, match="not applicable"):
        normalize_config_document(config)


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("field", "max_gap_ratio", 0.0),
        ("field", "max_local_charge_ratio", 0.0),
        ("particles", "outer_update_stride", 0),
        ("particles", "field_evolution_timescale", -1.0),
        ("particles", "max_frozen_field_ratio", 0.0),
    ],
)
def test_external_boundary_facade_rejects_invalid_control_ranges(
    section: str,
    key: str,
    value: object,
) -> None:
    config = load_toml_file(Path("examples/periodic2_kinetic_outer.toml"))
    config["external_boundary"][section][key] = value

    with pytest.raises(ConfigError, match=key):
        normalize_config_document(config)


def test_external_boundary_facade_validates_zhao_steady_start_semantics() -> None:
    valid = load_toml_file(
        Path("examples/periodic2_zhao_charge_driven_outer.toml")
    )
    normalize_config_document(valid)

    invalid_mode = copy.deepcopy(valid)
    invalid_mode["external_boundary"]["particles"]["steady_start_mode"] = "bogus"
    with pytest.raises(ConfigError, match="steady_start_mode"):
        normalize_config_document(invalid_mode)

    invalid_mesh_id = copy.deepcopy(valid)
    invalid_mesh_id["external_boundary"]["particles"]["steady_start_mesh_id"] = 0
    with pytest.raises(ConfigError, match="steady_start_mesh_id"):
        normalize_config_document(invalid_mesh_id)

    invalid_contract = load_toml_file(Path("examples/periodic2_kinetic_outer.toml"))
    invalid_contract["external_boundary"]["particles"][
        "steady_start_mode"
    ] = "zhao_floating"
    with pytest.raises(ConfigError, match="not applicable"):
        normalize_config_document(invalid_contract)

    invalid_mesh = copy.deepcopy(valid)
    invalid_mesh["mesh"]["mode"] = "obj"
    with pytest.raises(ConfigError, match="mesh.mode=template"):
        normalize_config_document(invalid_mesh)


@pytest.mark.parametrize(
    ("owner", "value"),
    [
        ("outer_plasma", {}),
        ("coupling", {}),
        ("reservoir_potential_model", "none"),
        ("open_boundary_model", "escape"),
    ],
)
def test_external_boundary_facade_rejects_raw_owner_mixing_by_presence(
    owner: str,
    value: object,
) -> None:
    config = default_config()
    config["external_boundary"] = {
        "field": {"model": "none"},
        "particles": {"mode": "local_source"},
    }
    if owner in {"outer_plasma", "coupling"}:
        config[owner] = value
    else:
        config["sim"][owner] = value

    with pytest.raises(ConfigError, match="cannot be combined"):
        normalize_config_document(config)


def test_external_boundary_facade_enforces_inflow_ownership() -> None:
    kinetic = load_toml_file(Path("examples/periodic2_kinetic_outer.toml"))
    kinetic["external_boundary"]["particles"]["inflow_model"] = "source_vdf"
    with pytest.raises(ConfigError, match="requires.*inflow_model=auto"):
        normalize_config_document(kinetic)

    removed_inflow = copy.deepcopy(kinetic)
    removed_inflow["external_boundary"]["particles"]["inflow_model"] = "legacy_sheath"
    with pytest.raises(ConfigError, match="inflow_model"):
        normalize_config_document(removed_inflow)

    removed_option = copy.deepcopy(kinetic)
    removed_option["external_boundary"]["particles"][
        "legacy_sheath_model"
    ] = "zhao_auto"
    with pytest.raises(ConfigError, match="unsupported external_boundary.particles"):
        normalize_config_document(removed_option)


def test_kinetic_outer_requires_unique_ambient_species_and_rejects_removed_gauge_key() -> (
    None
):
    config = load_config_file(Path("examples/periodic2_kinetic_outer.toml"))

    duplicate = copy.deepcopy(config)
    duplicate["particles"]["species"].append(
        copy.deepcopy(duplicate["particles"]["species"][0])
    )
    with pytest.raises(ConfigValidationError, match="exactly one enabled negative"):
        normalize_config_document(duplicate)

    shifted_gauge = copy.deepcopy(config)
    shifted_gauge["outer_plasma"]["infinity_potential"] = 1.0
    with pytest.raises(ConfigValidationError, match="infinity_potential"):
        normalize_config_document(shifted_gauge)


def test_zhao_transient_outer_queue_constraints() -> None:
    config = load_config_file(Path("examples/periodic2_zhao_transient_outer.toml"))

    assert config["coupling"]["outer_queue_enabled"] is True
    normalize_config_document(config)

    default_dt_step = copy.deepcopy(config)
    del default_dt_step["sim"]["dt"]
    del default_dt_step["sim"]["batch_duration"]
    default_dt_step["sim"]["batch_duration_step"] = 100.0
    normalize_config_document(default_dt_step)

    invalid = copy.deepcopy(config)
    invalid["coupling"]["outer_update_stride"] = 2
    with pytest.raises(ConfigValidationError, match="outer_update_stride=1"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["sim"]["batch_duration"] = 0.0
    with pytest.raises(ConfigValidationError, match="finite positive"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["sim"]["batch_duration"] = 5.0e-6
    with pytest.raises(ConfigValidationError, match="field evolution timescale"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["kinetic_closure"] = "absorbing_maxwellian"
    with pytest.raises(ConfigValidationError, match="persistent outer queue"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["zhao_branch"] = "b"
    with pytest.raises(ConfigValidationError, match="automatic Zhao"):
        normalize_config_document(invalid)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["photoelectron_source_scale"] = 0.0
    with pytest.raises(ConfigValidationError, match="photoelectron_source_scale"):
        normalize_config_document(invalid)


def test_tracked_outer_transfer_requires_opposite_charge_deposit() -> None:
    config = load_config_file(
        Path("examples/periodic2_zhao_charge_driven_outer.toml")
    )
    photo_species = next(
        species
        for species in config["particles"]["species"]
        if species.get("source_mode") == "photo_raycast"
    )
    photo_species["deposit_opposite_charge_on_emit"] = False

    with pytest.raises(ConfigValidationError, match="deposit_opposite_charge_on_emit"):
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
use_box = true
box_origin = [1.0, 2.0, 3.0]
box_size = [2.0, 4.0, 6.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"

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

    assert result["sim"]["box_min"] == [1.0, 2.0, 3.0]
    assert result["sim"]["box_max"] == [3.0, 6.0, 9.0]
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
use_box = true
box_origin = [1.0, 2.0, 3.0]
box_size = [2.0, 4.0, 6.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"

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
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]

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
    assert initialized["sim"]["bc_x_low"] == "periodic"
    assert initialized["sim"]["bc_x_high"] == "periodic"
    assert initialized["sim"]["bc_y_low"] == "periodic"
    assert initialized["sim"]["bc_y_high"] == "periodic"
    assert initialized["sim"]["field_solver"] == "fmm"
    assert initialized["sim"]["field_bc_mode"] == "periodic2"
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
