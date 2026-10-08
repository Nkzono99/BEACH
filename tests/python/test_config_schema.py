from __future__ import annotations

from beach.config._layout import to_runtime_layout

import copy
from pathlib import Path

import pytest

from beach.config._toml import load_toml_file
from beach.config.schema import load_schema, schema_errors


ROOT = Path(__file__).resolve().parents[2]
SCHEMA_CANONICAL = ROOT / "schemas/beach.schema.json"
SCHEMA_COPIES = (
    ROOT / "beach/config/schemas/beach.schema.json",
    ROOT / "plugins/beach-context/references/schemas/beach.schema.json",
)


def _minimal_valid_config() -> dict[str, object]:
    return {"particles": {"species": [{}]}}


def _config_with_sim(**sim: object) -> dict[str, object]:
    config = _minimal_valid_config()
    config["sim"] = sim
    return config


def test_schema_distribution_copies_match_canonical() -> None:
    expected = SCHEMA_CANONICAL.read_bytes()

    for path in SCHEMA_COPIES:
        assert path.read_bytes() == expected, (
            f"{path.relative_to(ROOT)} must match schemas/beach.schema.json byte-for-byte"
        )


def test_schema_distinguishes_fortran_integer_storage_from_real_values() -> None:
    schema, _ = load_schema()
    for value in (-(2**31), 2**31 - 1):
        assert schema_errors(_config_with_sim(rng_seed=value), schema) == []
    for value in (-(2**31) - 1, 2**31, 1.0, True):
        assert schema_errors(_config_with_sim(rng_seed=value), schema)
    config = _minimal_valid_config()
    config["particles"]["species"][0]["number_density_m3"] = 10**12
    assert schema_errors(config, schema) == []


def test_tree_parameter_schema_uses_runtime_omission_semantics() -> None:
    schema, _ = load_schema()
    sim = schema["$defs"]["sim"]["properties"]

    for key in ("tree_theta", "tree_leaf_max"):
        assert "default" not in sim[key]
        assert "element-count heuristic" in sim[key]["description"]
    assert "direct evaluation to FMM" in sim["tree_min_nelem"]["description"]


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("field", {"element_kernel": "point"}),
        ("outer_plasma", {"model": "kinetic_1d"}),
        ("coupling", {"update_mode": "explicit"}),
        ("external_boundary", {}),
    ],
    ids=("field", "outer-plasma", "coupling", "external-boundary"),
)
def test_schema_rejects_removed_top_level_contracts(key: str, value: object) -> None:
    schema, _ = load_schema()
    config = _minimal_valid_config()

    assert schema_errors(config, schema) == []
    config[key] = value
    assert schema_errors(config, schema)


def test_schema_rejects_removed_sim_softening() -> None:
    schema, _ = load_schema()
    config = _config_with_sim(softening=1.0e-6)

    assert schema_errors(_config_with_sim(), schema) == []
    assert schema_errors(config, schema)


def test_schema_enforces_soft_discard_bounds_only_when_active() -> None:
    schema, _ = load_schema()
    inactive_values = {
        "multiple_box_events_soft_discard_count_grace": -1,
        "multiple_box_events_soft_discard_fraction_limit": -1.0,
        "multiple_box_events_soft_discard_abs_charge_limit": -1.0,
    }

    assert schema_errors(_config_with_sim(**inactive_values), schema) == []
    assert (
        schema_errors(
            _config_with_sim(
                multiple_box_events_policy="abort",
                **inactive_values,
            ),
            schema,
        )
        == []
    )

    active_values = {
        "multiple_box_events_policy": "soft_discard",
        "multiple_box_events_soft_discard_count_grace": 0,
        "multiple_box_events_soft_discard_fraction_limit": 1.0,
        "multiple_box_events_soft_discard_abs_charge_limit": 1.0e-30,
    }
    assert schema_errors(_config_with_sim(**active_values), schema) == []
    assert schema_errors(
        _config_with_sim(multiple_box_events_soft_discard_count_limit=1000),
        schema,
    )

    invalid_active_values = (
        ("multiple_box_events_soft_discard_count_grace", -1),
        ("multiple_box_events_soft_discard_count_grace", 0.5),
        ("multiple_box_events_soft_discard_fraction_limit", 0.0),
        ("multiple_box_events_soft_discard_fraction_limit", 1.000001),
        ("multiple_box_events_soft_discard_fraction_limit", "1e-6"),
        ("multiple_box_events_soft_discard_abs_charge_limit", 0.0),
        ("multiple_box_events_soft_discard_abs_charge_limit", "1e-12"),
    )
    for key, value in invalid_active_values:
        assert schema_errors(
            _config_with_sim(
                multiple_box_events_policy="soft_discard",
                **{key: value},
            ),
            schema,
        )


def test_schema_accepts_panel_spectral_reference() -> None:
    schema, _ = load_schema()
    config = _minimal_valid_config()
    config["periodic2"] = {
        "nonzero_mode_backend": "panel_spectral_reference",
        "zero_mode_policy": "exclude_k0",
        "lower_boundary_model": "e_bottom_zero",
    }

    assert schema_errors(config, schema) == []


def test_schema_accepts_species_reflection_actions_and_rejects_unknown_value() -> None:
    schema, _ = load_schema()
    config = to_runtime_layout(load_toml_file(ROOT / "examples/beach.toml"))
    config["particles"]["species"][0]["boundary"] = {"z_high": "reflect"}

    assert schema_errors(config, schema) == []

    redistributed = copy.deepcopy(config)
    redistributed["particles"]["species"][0]["boundary"][
        "z_high"
    ] = "redistributed_reflect"
    assert schema_errors(redistributed, schema) == []

    invalid = copy.deepcopy(config)
    invalid["particles"]["species"][0]["boundary"]["z_high"] = "unknown"
    assert schema_errors(invalid, schema)


def test_schema_accepts_boundary_inflow_and_plane_source_contracts() -> None:
    schema, _ = load_schema()
    config = {
        "sim": {"batch_duration": 1.0},
        "domain": {
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "npcls_per_step": 0,
                    "number_density_m3": 1.0,
                    "temperature_k": 0.0,
                    "w_particle": 1.0,
                    "boundary_inflow": {"z_high": "reservoir"},
                }
            ]
        },
    }
    assert schema_errors(config, schema) == []

    invalid_value = copy.deepcopy(config)
    invalid_value["particles"]["species"][0]["boundary_inflow"]["z_high"] = "open"
    assert schema_errors(invalid_value, schema)

    plane = copy.deepcopy(config)
    species = plane["particles"]["species"][0]
    species["source_mode"] = "plane_source"
    species.pop("npcls_per_step")
    species.pop("boundary_inflow")
    species["pos_low"] = [0.1, 0.1, 0.5]
    species["pos_high"] = [0.9, 0.9, 0.5]
    species["source_normal"] = [0.0, 0.0, -1.0]
    assert schema_errors(plane, schema) == []

    mixed = copy.deepcopy(plane)
    mixed["particles"]["species"][0]["boundary_inflow"] = {
        "z_high": "reservoir"
    }
    assert schema_errors(mixed, schema)

    photo_mixed = copy.deepcopy(config)
    photo_species = photo_mixed["particles"]["species"][0]
    photo_species["source_mode"] = "photo_raycast"
    photo_species.pop("npcls_per_step")
    photo_species["inject_face"] = "z_high"
    photo_species["emit_current_density_a_m2"] = 1.0
    photo_species["rays_per_batch"] = 1
    assert schema_errors(photo_mixed, schema)


def test_schema_constrains_neutral_return_to_closed_negative_photoelectrons() -> None:
    schema, _ = load_schema()
    config = load_toml_file(ROOT / "examples/periodic2_closed_photoelectron.toml")

    assert schema_errors(config, schema) == []

    species = config["particles"]["species"][-1]
    invalid_enum = copy.deepcopy(config)
    invalid_enum["particles"]["species"][-1]["surface_charge_closure"] = "unknown"
    assert schema_errors(invalid_enum, schema)

    nonphoto = copy.deepcopy(config)
    nonphoto["particles"]["species"][-1]["source_mode"] = "reservoir_face"
    assert schema_errors(nonphoto, schema)

    positive = copy.deepcopy(config)
    positive["particles"]["species"][-1]["q_particle"] = abs(species["q_particle"])
    assert schema_errors(positive, schema)

    no_countercharge = copy.deepcopy(config)
    no_countercharge["particles"]["species"][-1]["deposit_opposite_charge_on_emit"] = (
        False
    )
    assert schema_errors(no_countercharge, schema)

    no_reflect = copy.deepcopy(config)
    no_reflect["particles"]["species"][-1]["boundary"]["z_high"] = "inherit"
    assert schema_errors(no_reflect, schema) == []


def test_schema_constrains_fixed_surface_current_channels() -> None:
    schema, _ = load_schema()
    config = load_toml_file(ROOT / "examples/periodic2_closed_photoelectron.toml")
    species = config["particles"]["species"][-1]
    species["surface_charge_closure"] = "fixed_current"
    species["target_absorbed_current_a"] = -2.0e-6
    species["target_emission_current_a"] = 3.0e-6
    assert schema_errors(config, schema) == []

    no_target = copy.deepcopy(config)
    no_target["particles"]["species"][-1].pop("target_absorbed_current_a")
    no_target["particles"]["species"][-1].pop("target_emission_current_a")
    # Cross-species automatic models may supply both targets; semantic
    # validation rejects a target-less fixed_current when no model does so.
    assert schema_errors(no_target, schema) == []

    implicit = copy.deepcopy(config)
    implicit["particles"]["species"][-1].pop("surface_charge_closure")
    assert schema_errors(implicit, schema)

    nonphoto_emission = copy.deepcopy(config)
    nonphoto_emission["particles"]["species"][-1]["source_mode"] = "volume_seed"
    assert schema_errors(nonphoto_emission, schema)


def test_schema_accepts_zhao_stationary_surface_current_model() -> None:
    schema, _ = load_schema()
    config = load_toml_file(ROOT / "examples/periodic2_zhao_fixed_current.toml")
    assert schema_errors(config, schema) == []

    disabled = copy.deepcopy(config)
    disabled["surface_current_model"] = {"model": "none"}
    assert schema_errors(disabled, schema) == []

    disabled_with_model_key = copy.deepcopy(disabled)
    disabled_with_model_key["surface_current_model"]["zhao_branch"] = "a"
    assert schema_errors(disabled_with_model_key, schema)

    refresh = load_toml_file(ROOT / "examples/periodic2_zhao_outflow_refresh.toml")
    assert schema_errors(refresh, schema) == []
    disabled_with_refresh = copy.deepcopy(disabled)
    disabled_with_refresh["surface_current_model"]["outflow_refresh_batches"] = 1
    assert schema_errors(disabled_with_refresh, schema)
    negative_refresh = copy.deepcopy(refresh)
    negative_refresh["surface_current_model"]["outflow_refresh_batches"] = -1
    assert schema_errors(negative_refresh, schema)

    no_photo = load_toml_file(
        ROOT / "examples/periodic2_zhao_no_photo_fixed_current.toml"
    )
    assert schema_errors(no_photo, schema) == []

    explicit_no_photo_type_c = copy.deepcopy(no_photo)
    explicit_no_photo_type_c["surface_current_model"]["zhao_branch"] = "c"
    assert schema_errors(explicit_no_photo_type_c, schema) == []

    stale_photo_setting = copy.deepcopy(no_photo)
    stale_photo_setting["surface_current_model"]["solar_elevation_deg"] = 60.0
    assert schema_errors(stale_photo_setting, schema)

    invalid_no_photo_branch = copy.deepcopy(no_photo)
    invalid_no_photo_branch["surface_current_model"]["zhao_branch"] = "a"
    assert schema_errors(invalid_no_photo_branch, schema)


def test_schema_rejects_removed_matching_plane_contract() -> None:
    schema, _ = load_schema()
    zhao = load_toml_file(ROOT / "examples/periodic2_zhao_fixed_current.toml")
    for key, value in (
        ("response_backend", "zhao_online"),
        ("response_table_path", "outer-response.csv"),
        ("implicit_zero_mode", True),
        ("coupling_rtol", 1.0e-4),
    ):
        removed_key = copy.deepcopy(zhao)
        removed_key["surface_current_model"][key] = value
        assert schema_errors(removed_key, schema), key

    removed_model = copy.deepcopy(zhao)
    removed_model["surface_current_model"]["model"] = "matching_plane_quasistatic"
    assert schema_errors(removed_model, schema)


def test_schema_constrains_upper_panel_fourier_retry_backend() -> None:
    schema, _ = load_schema()
    config = load_toml_file(ROOT / "examples/periodic2_zhao_outflow_refresh.toml")
    config["sim"]["multiple_box_events_retry_backend"] = (
        "upper_panel_fourier"
    )

    assert schema_errors(config, schema) == []

    wrong_nonzero_backend = copy.deepcopy(config)
    wrong_nonzero_backend["periodic2"]["nonzero_mode_backend"] = (
        "panel_spectral_reference"
    )
    assert schema_errors(wrong_nonzero_backend, schema)

    wrong_solver = copy.deepcopy(config)
    wrong_solver["sim"]["field_solver"] = "direct"
    assert schema_errors(wrong_solver, schema)

    wrong_field_boundary = copy.deepcopy(config)
    wrong_field_boundary["field_boundary"]["mode"] = "free"
    assert schema_errors(wrong_field_boundary, schema)

    legacy_cached = copy.deepcopy(config)
    legacy_cached.pop("surface_current_model")
    legacy_cached.pop("periodic2")
    assert schema_errors(legacy_cached, schema) == []


def test_schema_requires_surface_side_only_for_enabled_templates() -> None:
    schema, _ = load_schema()
    disabled = to_runtime_layout(load_toml_file(ROOT / "examples/beach.toml"))
    disabled["mesh"]["templates"].append({"enabled": False})
    enabled = copy.deepcopy(disabled)
    enabled["mesh"]["templates"][-1]["enabled"] = True
    implicit_enabled = copy.deepcopy(disabled)
    implicit_enabled["mesh"]["templates"][-1].pop("enabled")

    assert schema_errors(disabled, schema) == []
    assert schema_errors(enabled, schema)
    assert schema_errors(implicit_enabled, schema)


def test_schema_requires_obj_surface_side_for_explicit_auto_or_obj_mode() -> None:
    schema, _ = load_schema()
    base = to_runtime_layout(load_toml_file(ROOT / "examples/tutorial_insulator.toml"))
    auto = copy.deepcopy(base)
    auto["mesh"]["mode"] = "auto"
    auto["mesh"].pop("surface_side", None)
    obj = copy.deepcopy(auto)
    obj["mesh"]["mode"] = "obj"
    implicit_template = copy.deepcopy(base)
    implicit_template["mesh"].pop("mode")

    assert schema_errors(auto, schema)
    assert schema_errors(obj, schema)
    auto["mesh"]["surface_side"] = "normal_plus"
    assert schema_errors(auto, schema) == []
    assert schema_errors(implicit_template, schema) == []
    assert schema_errors(base, schema) == []
