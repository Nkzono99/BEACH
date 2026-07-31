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
