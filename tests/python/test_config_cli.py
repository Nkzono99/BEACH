from __future__ import annotations

import copy
from pathlib import Path
from types import MappingProxyType

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


@pytest.mark.parametrize("invalid", [False, True])
@pytest.mark.parametrize("read_only", [False, True])
def test_runtime_validation_preserves_caller_configuration(
    invalid: bool, read_only: bool,
) -> None:
    config = default_config()
    if invalid:
        config["particles"]["species"][0]["source_mode"] = "unknown"
    original = copy.deepcopy(config)

    def freeze_tables(value):
        if isinstance(value, dict):
            return MappingProxyType({key: freeze_tables(item) for key, item in value.items()})
        if isinstance(value, list):
            return [freeze_tables(item) for item in value]
        return value

    document = freeze_tables(config) if read_only else config
    if invalid:
        with pytest.raises(ConfigValidationError, match="source_mode"):
            validate_runtime_config(document)
    else:
        validate_runtime_config(document)
    assert config == original


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
    assert result["mesh"]["templates"][0]["kind"] == "plane"


def test_declared_enum_normalization_preserves_leading_spaces() -> None:
    config = default_config()
    config["sim"]["field_solver"] = "DIRECT "
    assert normalize_config_document(config)["sim"]["field_solver"] == "direct"
    config["sim"]["field_solver"] = " DIRECT"
    with pytest.raises(ConfigValidationError, match="field_solver"):
        normalize_config_document(config)


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
    config["sim"][key] = 1.0

    with pytest.raises(ConfigValidationError, match=rf"removed sim key.*{key}"):
        normalize_config_document(config)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("field", {"element_kernel": "point"}),
        (
            "external_boundary",
            {"field": {"model": "none"}, "particles": {"mode": "local_source"}},
        ),
    ],
)
def test_runtime_validator_rejects_removed_top_level_contracts(
    key: str,
    value: object,
) -> None:
    config = default_config()
    config[key] = value

    with pytest.raises(ConfigError, match=rf"Additional properties.*{key}"):
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


def test_normalization_preserves_case_in_paths_and_species_identifiers() -> None:
    config = default_config()
    config["sim"]["field_solver"] = "DIRECT"
    config["output"]["dir"] = "Outputs/MixedCase"
    config["particles"]["species"][0]["species_key"] = "ElectronA"

    normalized = normalize_config_document(config)

    assert normalized["sim"]["field_solver"] == "direct"
    assert config["sim"]["field_solver"] == "DIRECT"
    assert normalized["output"]["dir"] == "Outputs/MixedCase"
    assert normalized["particles"]["species"][0]["species_key"] == "ElectronA"


def test_periodic2_accepts_symmetric_vacuum_and_rejects_unknown_lower_model() -> None:
    config = load_config_file(Path("examples/periodic2_closed_photoelectron.toml"))
    config["sim"]["field_solver"] = "direct"
    config["periodic2"] = {}
    config["periodic2"]["lower_boundary_model"] = "symmetric_vacuum"

    normalized = normalize_config_document(config)

    assert normalized["periodic2"]["lower_boundary_model"] == "symmetric_vacuum"
    assert "nonzero_mode_backend" not in normalized["periodic2"]
    wrong_solver = copy.deepcopy(config)
    wrong_solver["sim"]["field_solver"] = "fmm"
    with pytest.raises(ConfigValidationError, match="panel_spectral_reference.*direct"):
        normalize_config_document(wrong_solver)
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
    volume["sim"]["field_solver"] = "fmm"
    volume["domain"]["periodic_axes"] = ["x", "y"]
    volume["field_boundary"]["mode"] = "periodic2"
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
    config = _default_periodic2_config()
    config["particle_boundary"]["x_low"] = "reflect"
    with pytest.raises(ConfigValidationError, match="periodic domain face"):
        normalize_config_document(config)

    config = _default_periodic2_config()
    config["particles"]["species"][0]["boundary"] = {"x_high": "open"}
    with pytest.raises(ConfigValidationError, match="periodic domain face"):
        normalize_config_document(config)


@pytest.mark.parametrize("species_override", [False, True])
def test_particle_boundary_overrides_require_domain(species_override: bool) -> None:
    config = default_config()
    del config["domain"]
    config["particle_boundary"] = {"ordinary_open_model": "escape"}
    normalize_config_document(config)
    if species_override:
        config["particles"]["species"][0]["boundary"] = {"z_high": "reflect"}
    else:
        config["particle_boundary"]["z_high"] = "reflect"
    with pytest.raises(ConfigValidationError, match="boundary.*finite.*domain"):
        normalize_config_document(config)


def test_periodic_image_layers_bounds_apply_only_to_periodic_fields() -> None:
    config = default_config()
    config["sim"]["field_periodic_image_layers"] = -1
    normalize_config_document(config)
    config = _default_periodic2_config()
    config["sim"]["field_periodic_image_layers"] = 0
    normalize_config_document(config)
    config["sim"]["field_periodic_image_layers"] = -1
    with pytest.raises(ConfigValidationError, match="field_periodic_image_layers.*>= 0"):
        normalize_config_document(config)


def test_enabled_particle_species_require_nonzero_charge() -> None:
    config = default_config()
    config["particles"]["species"][0]["q_particle"] = 0.0
    with pytest.raises(ConfigValidationError, match="q_particle.*non-zero"):
        normalize_config_document(config)
    disabled = copy.deepcopy(config["particles"]["species"][0])
    disabled["enabled"] = False
    config = default_config()
    config["particles"]["species"].append(disabled)
    normalize_config_document(config)


@pytest.mark.parametrize("enabled", [False, True])
def test_particle_species_keys_cannot_duplicate_generated_defaults(enabled: bool) -> None:
    config = default_config()
    species = copy.deepcopy(config["particles"]["species"][0])
    species.update(species_key="species_1", enabled=enabled)
    config["particles"]["species"].append(species)
    with pytest.raises(ConfigValidationError, match="species_key.*unique"):
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

    config = _default_periodic2_config()
    config["sim"]["batch_duration"] = 1.0e-6
    periodic_species = config["particles"]["species"][0]
    periodic_species["number_density_m3"] = 1.0e6
    periodic_species["temperature_k"] = 2.0e4
    periodic_species["boundary_inflow"] = {"x_low": "reservoir"}
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
    with pytest.raises(ConfigValidationError, match="surface_model.*dielectric"):
        normalize_config_document(config)

    config = default_config()
    config["mesh"]["templates"][0]["epsilon_r"] = 3.9
    with pytest.raises(ConfigValidationError, match="epsilon_r"):
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
            "multiple_box_events_soft_discard_count_grace": -1,
            "multiple_box_events_soft_discard_fraction_limit": -1.0,
            "multiple_box_events_soft_discard_abs_charge_limit": -1.0,
            "raycast_max_bounce": 0,
        }
    )
    config["field_boundary"]["mode"] = "free"

    normalized = normalize_config_document(config)

    assert normalized["sim"]["field_solver"] == "direct"


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("multiple_box_events_soft_discard_count_grace", True),
        ("multiple_box_events_soft_discard_count_grace", -1),
        ("multiple_box_events_soft_discard_count_grace", 1.5),
        ("multiple_box_events_soft_discard_fraction_limit", True),
        ("multiple_box_events_soft_discard_fraction_limit", 0.0),
        ("multiple_box_events_soft_discard_fraction_limit", 1.000001),
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

    with pytest.raises(ConfigValidationError, match=key):
        normalize_config_document(config)


def test_active_soft_discard_controls_accept_range_boundaries() -> None:
    config = default_config()
    config["sim"].update(
        {
            "multiple_box_events_policy": "soft_discard",
            "multiple_box_events_soft_discard_count_grace": 0,
            "multiple_box_events_soft_discard_fraction_limit": 1.0,
        }
    )

    normalized = normalize_config_document(config)

    assert normalized["sim"]["multiple_box_events_soft_discard_count_grace"] == 0
    assert normalized["sim"]["multiple_box_events_soft_discard_fraction_limit"] == 1.0


def test_load_config_file_rejects_nonfinite_template_scalar(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace(
        "size_x = 1.0", "size_x = inf"
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(ConfigValidationError, match="mesh.templates\\[0\\].size_x"):
        load_config_file(config_path)


def test_load_config_file_rejects_removed_photo_escape_model(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    config_path.write_text(
        config_path.read_text(encoding="utf-8").replace(
            "npcls_per_step = 10",
            'npcls_per_step = 10\nphoto_escape_model = "boltzmann_cutoff"',
        ),
        encoding="utf-8",
    )

    with pytest.raises(ConfigValidationError, match="photo_escape_model"):
        load_config_file(config_path)


def test_fixed_current_requires_signed_independent_channels() -> None:
    config = default_config()
    config["sim"]["batch_duration"] = 1.0e-2
    species = config["particles"]["species"][0]
    species["surface_charge_closure"] = "fixed_current"
    species["target_absorbed_current_a"] = -2.0e-6
    normalized = normalize_config_document(config)
    assert normalized["particles"]["species"][0][
        "target_absorbed_current_a"
    ] == pytest.approx(-2.0e-6)

    stepped_duration = copy.deepcopy(config)
    stepped_duration["sim"].pop("batch_duration")
    stepped_duration["sim"]["batch_duration_step"] = 10.0
    normalize_config_document(stepped_duration)

    wrong_sign = copy.deepcopy(config)
    wrong_sign["particles"]["species"][0]["target_absorbed_current_a"] = 2.0e-6
    with pytest.raises(ConfigValidationError, match="same sign as q_particle"):
        normalize_config_document(wrong_sign)

    net_only = copy.deepcopy(config)
    net_only["particles"]["species"][0].pop("target_absorbed_current_a")
    with pytest.raises(ConfigValidationError, match="at least one target current"):
        normalize_config_document(net_only)

    zero_duration = copy.deepcopy(config)
    zero_duration["sim"]["batch_duration"] = 0.0
    with pytest.raises(ConfigValidationError, match="batch_duration must be > 0"):
        normalize_config_document(zero_duration)

    disabled = default_config()
    disabled["particles"]["species"].append(
        {
            "enabled": False,
            "surface_charge_closure": "fixed_current",
            "target_absorbed_current_a": -2.0e-6,
        }
    )
    normalize_config_document(disabled)

    photo_config = load_config_file(
        Path("examples/periodic2_closed_photoelectron.toml")
    )
    photoelectron = photo_config["particles"]["species"][-1]
    photoelectron["surface_charge_closure"] = "fixed_current"
    photoelectron["target_emission_current_a"] = 3.0e-6
    normalized_photo = normalize_config_document(photo_config)
    assert normalized_photo["particles"]["species"][-1][
        "target_emission_current_a"
    ] == pytest.approx(3.0e-6)

    wrong_emission_sign = copy.deepcopy(photo_config)
    wrong_emission_sign["particles"]["species"][-1][
        "target_emission_current_a"
    ] = -3.0e-6
    with pytest.raises(ConfigValidationError, match="opposite to q_particle"):
        normalize_config_document(wrong_emission_sign)


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
    with pytest.raises(ConfigValidationError, match="drift_velocity.*too short"):
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


def test_zhao_outflow_refresh_requires_a_photoelectron_split_periodic_cell() -> None:
    root = Path(__file__).resolve().parents[2]
    refresh = load_config_file(root / "examples/periodic2_zhao_outflow_refresh.toml")
    assert refresh["surface_current_model"]["outflow_refresh_batches"] == 2

    for value in (-1, 1.5, True):
        invalid = copy.deepcopy(refresh)
        invalid["surface_current_model"]["outflow_refresh_batches"] = value
        with pytest.raises(ConfigValidationError, match="outflow_refresh_batches"):
            normalize_config_document(invalid)

    without_split = copy.deepcopy(refresh)
    without_split.pop("periodic2")
    with pytest.raises(ConfigValidationError, match="split-zero-mode"):
        normalize_config_document(without_split)

    no_photo = load_config_file(
        root / "examples/periodic2_zhao_no_photo_fixed_current.toml"
    )
    no_photo["surface_current_model"]["outflow_refresh_batches"] = 1
    with pytest.raises(ConfigValidationError, match="photoelectron_source_scale > 0"):
        normalize_config_document(no_photo)

    fixed_root = copy.deepcopy(refresh)
    fixed_root["surface_current_model"]["outflow_refresh_batches"] = 0
    fixed_root.pop("periodic2")
    assert normalize_config_document(fixed_root)["surface_current_model"]["model"] == (
        "zhao_stationary"
    )


@pytest.mark.parametrize(
    ("role_index", "key", "value", "message"),
    [
        (0, "q_particle", -2.0 * 1.602176634e-19, "singly charged"),
        (2, "m_particle", 2.0 * 9.1093837139e-31, "matching ambient-electron"),
        (1, "m_particle", -1.0, "m_particle.*minimum"),
        (0, "temperature_ev", 0.0, "electron temperature must be positive"),
        (2, "temperature_ev", 0.0, "photoelectron temperature must be positive"),
        (1, "temperature_ev", 2.0, "cold ions"),
        (0, "drift_velocity", [0.0, 0.0, 1.0e3], "inward z-high"),
        (1, "drift_velocity", [0.0, 0.0, 0.0], "inward z-high"),
        (0, "drift_velocity", [float("nan"), 0.0, -4.0e5], "drift_velocity.*finite"),
        (1, "number_density_cm3", 0.0, "finite and > 0|positive number density"),
        (0, "number_density_cm3", 5.0, "share the solar-wind number density"),
    ],
)
def test_zhao_stationary_rejects_unsupported_species_contract(
    role_index: int, key: str, value: object, message: str
) -> None:
    root = Path(__file__).resolve().parents[2]
    config = load_config_file(root / "examples/periodic2_zhao_fixed_current.toml")
    config["particles"]["species"][role_index][key] = value

    with pytest.raises(ConfigValidationError, match=message):
        normalize_config_document(config)


def test_zhao_stationary_accepts_a_nondrifting_electron_reservoir() -> None:
    root = Path(__file__).resolve().parents[2]
    config = load_config_file(root / "examples/periodic2_zhao_fixed_current.toml")
    config["particles"]["species"][0]["drift_velocity"] = [0.0, 0.0, 0.0]
    normalized = normalize_config_document(config)
    assert normalized["particles"]["species"][0]["drift_velocity"] == [0.0, 0.0, 0.0]


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
    initialized_text = (tmp_path / "beach.toml").read_text(encoding="utf-8")

    with pytest.raises(SystemExit, match="config file already exists"):
        beachx_main(["config", "init"])
    assert (tmp_path / "beach.toml").read_text(encoding="utf-8") == initialized_text

    beachx_main(["config", "validate"])
    validate_streams = capsys.readouterr()
    assert "status=ok" in validate_streams.out

    initialized = load_config_file(tmp_path / "beach.toml")
    assert initialized["domain"]["periodic_axes"] == []
    assert initialized["sim"]["field_solver"] == "direct"
    assert initialized["field_boundary"]["mode"] == "free"
    assert initialized["sim"]["batch_count"] == 20
    assert len(initialized["particles"]["species"]) == 1
    species = initialized["particles"]["species"][0]
    assert species["source_mode"] == "volume_seed"
    assert species["npcls_per_step"] == 200
    assert species["w_particle"] == 2.0e5
    assert species["pos_low"] == [0.15, 0.15, 0.8]
    assert species["pos_high"] == [0.85, 0.85, 0.8]
    assert species["drift_velocity"] == [0.0, 0.0, -1.0e6]

    modified = tmp_path / "modified.toml"
    modified.write_text(
        (tmp_path / "beach.toml")
        .read_text(encoding="utf-8")
        .replace("batch_count = 20", "batch_count = 21"),
        encoding="utf-8",
    )
    beachx_main(["config", "diff", "beach.toml", str(modified)])
    diff_streams = capsys.readouterr()
    assert "sim.batch_count: 20 -> 21" in diff_streams.out


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


def test_lint_cli_reports_toml_parse_error(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    config_path.write_text("[sim\n", encoding="utf-8")

    with pytest.raises(SystemExit, match="TOML parse error"):
        beachx_main(["lint", str(config_path)])


def test_lint_custom_schema_adds_constraints_without_replacing_beach_contract(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "beach.toml"
    schema_path = tmp_path / "extra.schema.json"
    _write_base_config(config_path)
    schema_path.write_text("{}", encoding="utf-8")
    beachx_main(["lint", str(config_path), "--schema", str(schema_path)])

    schema_path.write_text(
        '{"properties": {"sim": {"properties": {"batch_count": {"maximum": 1}}}}}',
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
