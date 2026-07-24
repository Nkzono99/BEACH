from __future__ import annotations

import copy
from pathlib import Path

from beach.config import load_config_file
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


def test_schema_rejects_removed_outer_boundary_and_sheath_keys() -> None:
    schema, _ = load_schema()
    config = load_config_file(ROOT / "examples/periodic2_kinetic_outer.toml")
    removed_outer_keys = {
        "infinity_potential": 0.0,
        "photoelectron_histogram_enabled": False,
        "photoelectron_histogram_bins": 16,
        "photoelectron_histogram_energy_max": 10.0,
        "photoelectron_ambient_charge_scale": 1.0,
        "max_photoelectron_charge_ratio": 0.1,
    }
    removed_sim_keys = {
        "sheath_injection_model": "none",
        "sheath_reference_coordinate": 0.0,
    }

    outer_properties = schema["properties"]["outer_plasma"]["properties"]
    facade_field_properties = schema["$defs"]["externalBoundaryField"]["properties"]
    sim_properties = schema["$defs"]["sim"]["properties"]
    assert not (removed_outer_keys.keys() & outer_properties.keys())
    assert not (removed_outer_keys.keys() & facade_field_properties.keys())
    assert not (removed_sim_keys.keys() & sim_properties.keys())

    for key, value in removed_outer_keys.items():
        invalid = copy.deepcopy(config)
        invalid["outer_plasma"][key] = value
        assert schema_errors(invalid, schema)

    authoring = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")
    for key, value in removed_outer_keys.items():
        invalid = copy.deepcopy(authoring)
        invalid["external_boundary"]["field"][key] = value
        assert schema_errors(invalid, schema)

    for key, value in removed_sim_keys.items():
        invalid = copy.deepcopy(config)
        invalid["sim"][key] = value
        assert schema_errors(invalid, schema)


def test_schema_accepts_external_boundary_facade_and_rejects_raw_owner_mix() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")

    assert schema_errors(authoring, schema) == []

    for owner, value in (
        ("outer_plasma", {}),
        ("coupling", {}),
    ):
        mixed = copy.deepcopy(authoring)
        mixed[owner] = value
        assert schema_errors(mixed, schema)

    for selector, value in (
        ("reservoir_potential_model", "none"),
        ("open_boundary_model", "escape"),
    ):
        mixed = copy.deepcopy(authoring)
        mixed["sim"][selector] = value
        assert schema_errors(mixed, schema)


def test_schema_accepts_linearized_photoelectron_mean_only_with_ambient_linear_closure() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(
        ROOT / "examples/periodic2_ambient_linear_photo_outer.toml"
    )

    assert schema_errors(authoring, schema) == []

    wrong_closure = copy.deepcopy(authoring)
    wrong_closure["external_boundary"]["field"]["kinetic_closure"] = (
        "absorbing_maxwellian"
    )
    assert schema_errors(wrong_closure, schema)

    wrong_density = copy.deepcopy(authoring)
    wrong_density["external_boundary"]["field"][
        "photoelectron_density_model"
    ] = "kinetic_mean"
    assert schema_errors(wrong_density, schema)


def test_schema_requires_explicit_box_for_an_active_external_field() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")

    missing_interface_box = copy.deepcopy(authoring)
    del missing_interface_box["sim"]["box_max"]
    assert schema_errors(missing_interface_box, schema)

    high_level_box = copy.deepcopy(authoring)
    box_min = high_level_box["sim"].pop("box_min")
    box_max = high_level_box["sim"].pop("box_max")
    high_level_box["sim"]["box_origin"] = box_min
    high_level_box["sim"]["box_size"] = [
        box_max[index] - box_min[index] for index in range(3)
    ]
    assert schema_errors(high_level_box, schema) == []

    del high_level_box["sim"]["box_size"]
    assert schema_errors(high_level_box, schema)


def test_schema_enforces_external_boundary_required_models_and_queue_shape() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(ROOT / "examples/periodic2_zhao_transient_outer.toml")

    assert schema_errors(authoring, schema) == []

    missing_field_model = copy.deepcopy(authoring)
    del missing_field_model["external_boundary"]["field"]["model"]
    assert schema_errors(missing_field_model, schema)

    missing_particle_mode = copy.deepcopy(authoring)
    del missing_particle_mode["external_boundary"]["particles"]["mode"]
    assert schema_errors(missing_particle_mode, schema)

    invalid_queue = copy.deepcopy(authoring)
    invalid_queue["external_boundary"]["field"]["kinetic_closure"] = (
        "absorbing_maxwellian"
    )
    assert schema_errors(invalid_queue, schema)

    invalid_queue = copy.deepcopy(authoring)
    invalid_queue["external_boundary"]["particles"]["inflow_model"] = "source_vdf"
    assert schema_errors(invalid_queue, schema)

    for field_key, invalid_value in (
        ("model", "linear_debye"),
        ("zhao_branch", "b"),
        ("photoelectron_source_scale", 0.0),
    ):
        invalid_queue = copy.deepcopy(authoring)
        invalid_queue["external_boundary"]["field"][field_key] = invalid_value
        assert schema_errors(invalid_queue, schema)

    for particle_key, invalid_value in (
        ("steady_start_mode", "zhao_floating"),
        ("outer_update_stride", 2),
        ("field_evolution_timescale", 0.0),
        ("max_frozen_field_ratio", 0.0),
    ):
        invalid_queue = copy.deepcopy(authoring)
        invalid_queue["external_boundary"]["particles"][particle_key] = invalid_value
        assert schema_errors(invalid_queue, schema)

    missing_queue_timescale = copy.deepcopy(authoring)
    del missing_queue_timescale["external_boundary"]["particles"][
        "field_evolution_timescale"
    ]
    assert schema_errors(missing_queue_timescale, schema)

    missing_queue_batch_width = copy.deepcopy(authoring)
    missing_queue_batch_width["sim"].pop("batch_duration", None)
    missing_queue_batch_width["sim"].pop("batch_duration_step", None)
    assert schema_errors(missing_queue_batch_width, schema)


def test_schema_enforces_facade_model_mode_requirements() -> None:
    schema, _ = load_schema()
    kinetic = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")

    inactive = copy.deepcopy(kinetic)
    inactive["external_boundary"]["field"] = {"model": "none"}
    assert schema_errors(inactive, schema)

    competing_inflow = copy.deepcopy(kinetic)
    competing_inflow["external_boundary"]["particles"]["inflow_model"] = (
        "infinity_barrier"
    )
    assert schema_errors(competing_inflow, schema)

    missing_scale = copy.deepcopy(kinetic)
    del missing_scale["external_boundary"]["field"]["debye_length"]
    assert schema_errors(missing_scale, schema)

    missing_scale = copy.deepcopy(kinetic)
    del missing_scale["external_boundary"]["field"]["thermal_voltage"]
    assert schema_errors(missing_scale, schema)


def test_schema_rejects_removed_unified_field_vocabulary() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")

    for key, value in (
        ("max_linearity_ratio", 0.25),
        ("unified_grid_points", 129),
        ("accessible_fraction_tolerance", 0.1),
    ):
        invalid = copy.deepcopy(authoring)
        invalid["external_boundary"]["field"][key] = value
        assert schema_errors(invalid, schema)

    invalid = copy.deepcopy(authoring)
    invalid["external_boundary"]["field"]["model"] = "unified_linear_response"
    assert schema_errors(invalid, schema)

    runtime = load_config_file(ROOT / "examples/periodic2_kinetic_outer.toml")
    runtime["outer_plasma"]["model"] = "unified_linear_response"
    assert schema_errors(runtime, schema)


def test_schema_requires_facade_specialized_field_option_owners() -> None:
    schema, _ = load_schema()

    zhao = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")
    zhao["external_boundary"]["field"]["zhao_branch"] = "auto"
    assert schema_errors(zhao, schema)

    zhao_density = load_toml_file(
        ROOT / "examples/periodic2_zhao_charge_driven_outer.toml"
    )
    zhao_density["external_boundary"]["field"][
        "photoelectron_density_model"
    ] = "none"
    assert schema_errors(zhao_density, schema)


def test_schema_rejects_facade_mode_specific_particle_options() -> None:
    schema, _ = load_schema()

    for example, key, value in (
        ("periodic2_kinetic_outer.toml", "outer_orbit_dt", 1.0e-9),
        ("periodic2_zhao_transient_outer.toml", "outer_orbit_max_steps", 100),
        ("periodic2_kinetic_outer.toml", "outer_orbit_energy_tolerance", 1.0e-4),
        ("periodic2_zhao_transient_outer.toml", "outer_update_stride", 1),
    ):
        authoring = load_toml_file(ROOT / "examples" / example)
        authoring["external_boundary"]["particles"][key] = value
        assert schema_errors(authoring, schema)


def test_schema_enforces_facade_control_ranges() -> None:
    schema, _ = load_schema()

    for section, key, value in (
        ("field", "max_gap_ratio", 0.0),
        ("field", "max_local_charge_ratio", 0.0),
        ("particles", "outer_update_stride", 0),
        ("particles", "field_evolution_timescale", -1.0),
        ("particles", "max_frozen_field_ratio", 0.0),
    ):
        authoring = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")
        authoring["external_boundary"][section][key] = value
        assert schema_errors(authoring, schema)


def test_schema_enforces_facade_zhao_closure_and_zero_source_exception() -> None:
    schema, _ = load_schema()
    no_photo = load_toml_file(
        ROOT / "examples/periodic2_zhao_no_photo_outer.toml"
    )
    no_photo["sim"]["sheath_photoelectron_ref_density_cm3"] = 0.0

    assert (
        no_photo["external_boundary"]["field"]["photoelectron_source_scale"] == 0.0
    )
    assert schema_errors(no_photo, schema) == []

    negative_density = copy.deepcopy(no_photo)
    negative_density["sim"]["sheath_photoelectron_ref_density_cm3"] = -1.0
    assert schema_errors(negative_density, schema)

    photoemitting = load_toml_file(
        ROOT / "examples/periodic2_zhao_charge_driven_outer.toml"
    )
    photoemitting["sim"]["sheath_photoelectron_ref_density_cm3"] = 0.0
    assert schema_errors(photoemitting, schema)

    invalid_model = copy.deepcopy(no_photo)
    invalid_model["external_boundary"]["field"]["model"] = "linear_debye"
    assert schema_errors(invalid_model, schema)

    invalid_density = copy.deepcopy(no_photo)
    invalid_density["external_boundary"]["field"][
        "photoelectron_density_model"
    ] = "kinetic_mean"
    assert schema_errors(invalid_density, schema)

    local_source = copy.deepcopy(no_photo)
    local_source["external_boundary"]["particles"]["mode"] = "local_source"
    local_source["external_boundary"]["particles"]["inflow_model"] = "source_vdf"
    local_source["external_boundary"]["particles"].pop(
        "field_evolution_timescale", None
    )
    assert schema_errors(local_source, schema) == []

    competing_inflow = copy.deepcopy(local_source)
    competing_inflow["external_boundary"]["particles"][
        "inflow_model"
    ] = "infinity_barrier"
    assert schema_errors(competing_inflow, schema)


def test_schema_rejects_removed_legacy_sheath_vocabulary() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")
    authoring["external_boundary"]["particles"]["inflow_model"] = "legacy_sheath"

    assert schema_errors(authoring, schema)

    authoring["external_boundary"]["particles"]["inflow_model"] = "source_vdf"
    authoring["external_boundary"]["particles"]["legacy_sheath_model"] = "zhao_auto"
    assert schema_errors(authoring, schema)


def test_schema_rejects_inert_external_boundary_none_options() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(ROOT / "examples/periodic2_kinetic_outer.toml")
    authoring["external_boundary"] = {
        "field": {"model": "none"},
        "particles": {"mode": "local_source"},
    }

    assert schema_errors(authoring, schema) == []

    no_explicit_box = copy.deepcopy(authoring)
    no_explicit_box["sim"].pop("box_min", None)
    no_explicit_box["sim"].pop("box_max", None)
    assert schema_errors(no_explicit_box, schema) == []

    inert_field = copy.deepcopy(authoring)
    inert_field["external_boundary"]["field"]["debye_length"] = 1.0
    assert schema_errors(inert_field, schema)

    inert_coupling = copy.deepcopy(authoring)
    inert_coupling["external_boundary"]["particles"]["outer_update_stride"] = 1
    assert schema_errors(inert_coupling, schema)


def test_schema_enforces_transient_zhao_queue_conditionals() -> None:
    schema, _ = load_schema()
    config = load_config_file(ROOT / "examples/periodic2_zhao_transient_outer.toml")

    assert schema_errors(config, schema) == []

    invalid_branch = copy.deepcopy(config)
    invalid_branch["outer_plasma"]["zhao_branch"] = "b"
    assert any(
        "outer_plasma.zhao_branch" in error
        for error in schema_errors(invalid_branch, schema)
    )

    missing_timescale = copy.deepcopy(config)
    del missing_timescale["coupling"]["field_evolution_timescale"]
    assert any(
        "coupling" in error for error in schema_errors(missing_timescale, schema)
    )

    missing_batch_width = copy.deepcopy(config)
    del missing_batch_width["sim"]["batch_duration"]
    missing_batch_width["sim"].pop("batch_duration_step", None)
    assert any("sim" in error for error in schema_errors(missing_batch_width, schema))


def test_schema_enforces_facade_zhao_floating_steady_start_contract() -> None:
    schema, _ = load_schema()
    authoring = load_toml_file(
        ROOT / "examples/periodic2_zhao_charge_driven_outer.toml"
    )

    assert schema_errors(authoring, schema) == []

    invalid_mode = copy.deepcopy(authoring)
    invalid_mode["external_boundary"]["particles"]["mode"] = "local_source"
    assert schema_errors(invalid_mode, schema)

    invalid_model = copy.deepcopy(authoring)
    invalid_model["external_boundary"]["field"]["model"] = "linear_debye"
    assert schema_errors(invalid_model, schema)

    invalid_mesh = copy.deepcopy(authoring)
    invalid_mesh["mesh"]["mode"] = "obj"
    assert schema_errors(invalid_mesh, schema)

    invalid_mesh_id = copy.deepcopy(authoring)
    invalid_mesh_id["external_boundary"]["particles"]["steady_start_mesh_id"] = 0
    assert schema_errors(invalid_mesh_id, schema)

    inactive_steady_start = load_toml_file(
        ROOT / "examples/periodic2_kinetic_outer.toml"
    )
    inactive_steady_start["external_boundary"]["particles"][
        "steady_start_mode"
    ] = "none"
    assert schema_errors(inactive_steady_start, schema)


def test_schema_enforces_zhao_floating_steady_start_contract() -> None:
    schema, _ = load_schema()
    config = load_config_file(ROOT / "examples/periodic2_zhao_charge_driven_outer.toml")
    config["coupling"]["steady_start_mode"] = "zhao_floating"
    config["coupling"]["steady_start_mesh_id"] = 1

    assert schema_errors(config, schema) == []

    resumed = copy.deepcopy(config)
    resumed["output"]["resume"] = True
    assert schema_errors(resumed, schema) == []

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["model"] = "linear_debye"
    assert schema_errors(invalid, schema)

    invalid = copy.deepcopy(config)
    invalid["outer_plasma"]["return_model"] = "none"
    assert schema_errors(invalid, schema)

    invalid = copy.deepcopy(config)
    invalid["coupling"]["particle_transfer_mode"] = "none"
    assert schema_errors(invalid, schema)

    invalid = copy.deepcopy(config)
    invalid["coupling"]["outer_queue_enabled"] = True
    assert schema_errors(invalid, schema)

    invalid = copy.deepcopy(config)
    invalid["coupling"]["steady_start_mesh_id"] = 0
    assert schema_errors(invalid, schema)

    invalid = copy.deepcopy(config)
    invalid["mesh"]["mode"] = "obj"
    assert schema_errors(invalid, schema)


def test_schema_rejects_removed_outer_plasma_infinity_potential() -> None:
    schema, _ = load_schema()
    kinetic = load_config_file(ROOT / "examples/periodic2_kinetic_outer.toml")
    kinetic["outer_plasma"]["infinity_potential"] = 1.0
    assert schema_errors(kinetic, schema)
