from __future__ import annotations

import copy
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from beach.config import load_config_file

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 compatibility
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def _schema() -> dict[str, object]:
    return json.loads(_read("schemas/beach.schema.json"))


def test_parameter_section_inventory_covers_top_level_tables() -> None:
    expected_tables = (
        "sim",
        "particles",
        "mesh",
        "field",
        "periodic2",
        "outer_plasma",
        "coupling",
        "output",
    )

    for path, end_heading in (
        ("docs/Parameters.md", "## パラメータ詳細リファレンス"),
        ("docs/Parameters.en.md", "## Detailed Parameter Reference"),
    ):
        inventory = _read(path).split(end_heading, maxsplit=1)[0]
        for table in expected_tables:
            assert f"`[{table}]`" in inventory, (path, table)


def test_direct_periodic2_split_reference_matches_schema_runtime_and_docs() -> None:
    schema = _schema()
    validator = Draft202012Validator(schema)

    for path in (
        "examples/periodic2_linear_outer_reference.toml",
        "examples/periodic2_unified_linear_response.toml",
    ):
        document = tomllib.loads(_read(path))
        errors = sorted(validator.iter_errors(document), key=lambda error: list(error.path))
        assert not errors, (path, [error.message for error in errors])
        load_config_file(ROOT / path)

    description = schema["$defs"]["sim"]["properties"]["field_bc_mode"]["description"]
    assert "direct triangle_p0 panel_spectral_reference split model" in description

    for path in ("docs/FieldSolvers.md", "docs/FieldSolvers.en.md"):
        text = _read(path)
        assert "panel_spectral_reference" in text
        assert "exclude_k0" in text
        assert "triangle_p0" in text

    assert "periodic2 は fmm 必須" not in _read("docs/agent-user-guide.md")
    assert "periodic2 requires fmm" not in _read("docs/agent-user-guide.en.md")


def test_photoelectron_histogram_schema_requires_explicit_return_model() -> None:
    validator = Draft202012Validator(_schema())
    document = tomllib.loads(_read("examples/periodic2_photoelectron_return.toml"))

    errors = list(validator.iter_errors(document))
    assert not errors, [error.message for error in errors]

    del document["outer_plasma"]["return_model"]
    errors = list(validator.iter_errors(document))
    assert any(
        error.message == "'return_model' is a required property" for error in errors
    )


def test_photoelectron_schema_matches_density_transfer_and_deposit_contracts() -> None:
    validator = Draft202012Validator(_schema())
    kinetic = tomllib.loads(_read("examples/periodic2_kinetic_outer.toml"))
    kinetic["outer_plasma"]["photoelectron_density_model"] = "kinetic_mean"

    assert not list(validator.iter_errors(kinetic))

    mismatched = copy.deepcopy(kinetic)
    mismatched["coupling"]["particle_transfer_mode"] = "none"
    assert list(validator.iter_errors(mismatched))

    mismatched = copy.deepcopy(kinetic)
    mismatched["outer_plasma"]["return_model"] = "none"
    assert list(validator.iter_errors(mismatched))

    untracked = copy.deepcopy(kinetic)
    untracked["outer_plasma"]["return_model"] = "none"
    untracked["coupling"]["particle_transfer_mode"] = "none"
    assert not list(validator.iter_errors(untracked))

    tracked_photoelectron = tomllib.loads(
        _read("examples/periodic2_photoelectron_return.toml")
    )
    tracked_photoelectron["outer_plasma"]["photoelectron_histogram_enabled"] = False
    tracked_photoelectron["particles"]["species"][0][
        "deposit_opposite_charge_on_emit"
    ] = False
    assert list(validator.iter_errors(tracked_photoelectron))


def test_coupling_reference_covers_every_schema_key() -> None:
    coupling = _schema()["properties"]["coupling"]["properties"]
    for path in ("docs/Parameters.md", "docs/Parameters.en.md"):
        text = _read(path)
        for key in coupling:
            assert f"`{key}`" in text, (path, key)


def test_output_manifest_matches_implementation_and_bilingual_docs() -> None:
    manifest = json.loads(_read("schemas/beach.output-manifest.json"))
    files = manifest["files"]
    names = [entry["name"] for entry in files]
    assert len(names) == len(set(names))

    docs = (
        _read("docs/OutputGuide.md"),
        _read("docs/OutputGuide.en.md"),
        _read("docs/Parameters.md"),
        _read("docs/Parameters.en.md"),
    )
    for entry in files:
        name = entry["name"]
        assert "output.write_files=true" in entry["condition"], name
        producer = _read(entry["producer"])
        assert name in producer, (name, entry["producer"])
        for text in docs:
            assert name in text, name

        mpi_name = entry.get("mpi_name")
        if mpi_name:
            for text in docs:
                assert mpi_name in text, mpi_name

        consumer_path = entry.get("consumer")
        if consumer_path:
            assert name in _read(consumer_path), (name, consumer_path)

    for path in (
        "docs/OutputGuide.md",
        "docs/OutputGuide.en.md",
        "docs/Parameters.md",
        "docs/Parameters.en.md",
        "docs/agent-user-guide.md",
        "docs/agent-user-guide.en.md",
    ):
        text = _read(path)
        assert "`macro_residuals*.csv`" not in text, path
