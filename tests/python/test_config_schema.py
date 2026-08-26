from __future__ import annotations

import copy
from pathlib import Path

from beach.config._toml import load_toml_file
from beach.config.schema import load_schema, schema_definition_property_names
from beach.config.schema import schema_errors


ROOT = Path(__file__).resolve().parents[2]
SCHEMA_CANONICAL = ROOT / "schemas/beach.schema.json"
SCHEMA_COPIES = (
    ROOT / "beach/config/schemas/beach.schema.json",
    ROOT / "plugins/beach-context/references/schemas/beach.schema.json",
)


def test_schema_distribution_copies_match_canonical() -> None:
    expected = SCHEMA_CANONICAL.read_bytes()

    for path in SCHEMA_COPIES:
        assert path.read_bytes() == expected, (
            f"{path.relative_to(ROOT)} must match schemas/beach.schema.json byte-for-byte"
        )


def test_packaged_schema_exposes_workload_property_contracts() -> None:
    schema, label = load_schema()

    assert label == "package:beach.config/schemas/beach.schema.json"
    assert schema_definition_property_names("sim") == frozenset(
        schema["$defs"]["sim"]["properties"]
    )
    assert schema_definition_property_names("species") == frozenset(
        schema["$defs"]["species"]["properties"]
    )


def test_schema_rejects_removed_point_kernel_configuration() -> None:
    schema, _ = load_schema()

    assert "field" not in schema["properties"]
    assert "softening" not in schema["$defs"]["sim"]["properties"]
    assert schema_errors({"field": {"element_kernel": "point"}}, schema)
    assert schema_errors({"sim": {"softening": 1.0e-6}}, schema)


def test_schema_rejects_outer_coupling_and_keeps_panel_reference() -> None:
    schema, _ = load_schema()

    assert "outer_plasma" not in schema["properties"]
    assert "coupling" not in schema["properties"]
    assert schema_errors({"outer_plasma": {"model": "kinetic_1d"}}, schema)
    assert schema_errors({"coupling": {"update_mode": "explicit"}}, schema)
    assert "panel_spectral_reference" in (
        schema["properties"]["periodic2"]["properties"]["nonzero_mode_backend"][
            "enum"
        ]
    )


def test_schema_accepts_species_reflection_actions_and_rejects_unknown_value() -> None:
    schema, _ = load_schema()
    config = load_toml_file(ROOT / "examples/beach.toml")
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

    assert schema_errors({"external_boundary": {}}, schema)


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
    disabled_with_model_key["surface_current_model"]["coupling_rtol"] = 1.0e-4
    assert schema_errors(disabled_with_model_key, schema)

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


def test_schema_accepts_matching_plane_and_rejects_model_key_mixing() -> None:
    schema, _ = load_schema()
    matching = load_toml_file(
        ROOT / "tests/fortran/matching_plane_quasistatic.toml"
    )

    assert schema_errors(matching, schema) == []

    with_atol = copy.deepcopy(matching)
    with_atol["surface_current_model"]["coupling_atol"] = [0.0, 0.05, 0.0, 0.0]
    assert schema_errors(with_atol, schema) == []

    invalid_atol = copy.deepcopy(matching)
    invalid_atol["surface_current_model"]["coupling_atol"] = [0.0, -0.05, 0.0, 0.0]
    assert schema_errors(invalid_atol, schema)

    explicit_table = copy.deepcopy(matching)
    explicit_table["surface_current_model"]["response_backend"] = "table"
    assert schema_errors(explicit_table, schema) == []

    missing_response = copy.deepcopy(matching)
    missing_response["surface_current_model"].pop("response_table_path")
    assert schema_errors(missing_response, schema)

    zhao_key = copy.deepcopy(matching)
    zhao_key["surface_current_model"]["zhao_branch"] = "auto"
    assert schema_errors(zhao_key, schema)

    zhao_online = copy.deepcopy(matching)
    zhao_online["surface_current_model"].pop("response_table_path")
    zhao_online["surface_current_model"]["response_backend"] = "zhao_online"
    assert schema_errors(zhao_online, schema) == []

    zhao_online["surface_current_model"]["zhao_branch"] = "b"
    assert schema_errors(zhao_online, schema) == []

    online_with_table = copy.deepcopy(zhao_online)
    online_with_table["surface_current_model"]["response_table_path"] = (
        "outer-response.csv"
    )
    assert schema_errors(online_with_table, schema)

    invalid_backend = copy.deepcopy(zhao_online)
    invalid_backend["surface_current_model"]["response_backend"] = "unknown"
    assert schema_errors(invalid_backend, schema)

    reference_area = copy.deepcopy(matching)
    reference_area["surface_current_model"]["reference_area_m2"] = 1.0
    assert schema_errors(reference_area, schema)

    no_split_zero_mode = copy.deepcopy(matching)
    no_split_zero_mode.pop("periodic2")
    assert schema_errors(no_split_zero_mode, schema)

    nonzero_field = copy.deepcopy(matching)
    nonzero_field["sim"]["b0"] = [0.0, 0.0, 1.0]
    assert schema_errors(nonzero_field, schema)

    generic_open_barrier = copy.deepcopy(matching)
    generic_open_barrier["particle_boundary"]["ordinary_open_model"] = (
        "potential_barrier"
    )
    assert schema_errors(generic_open_barrier, schema)

    zhao = load_toml_file(ROOT / "examples/periodic2_zhao_fixed_current.toml")
    zhao["surface_current_model"]["response_table_path"] = "outer-response.csv"
    assert schema_errors(zhao, schema)

    zhao_backend = load_toml_file(
        ROOT / "examples/periodic2_zhao_fixed_current.toml"
    )
    zhao_backend["surface_current_model"]["response_backend"] = "zhao_online"
    assert schema_errors(zhao_backend, schema)


def test_schema_requires_surface_side_only_for_enabled_templates() -> None:
    schema, _ = load_schema()
    disabled = load_toml_file(ROOT / "examples/beach.toml")
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
    base = load_toml_file(ROOT / "examples/tutorial_insulator.toml")
    auto = copy.deepcopy(base)
    auto["mesh"]["mode"] = "auto"
    auto["mesh"].pop("surface_side", None)
    obj = copy.deepcopy(auto)
    obj["mesh"]["mode"] = "obj"

    assert schema["$defs"]["mesh"]["properties"]["mode"]["default"] == "template"
    assert schema_errors(auto, schema)
    assert schema_errors(obj, schema)
    auto["mesh"]["surface_side"] = "normal_plus"
    assert schema_errors(auto, schema) == []
    assert schema_errors(base, schema) == []
