from __future__ import annotations

import copy
from pathlib import Path

from beach.config import load_config_file
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
