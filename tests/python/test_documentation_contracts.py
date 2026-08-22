from __future__ import annotations

import json
from pathlib import Path

from jsonschema import Draft7Validator

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
        "domain",
        "field_boundary",
        "particle_boundary",
        "reservoir",
        "particles",
        "mesh",
        "periodic2",
        "output",
    )

    for path, end_heading in (
        ("docs/Parameters.md", "## パラメータ詳細リファレンス"),
        ("docs/Parameters.en.md", "## Detailed Parameter Reference"),
    ):
        inventory = _read(path).split(end_heading, maxsplit=1)[0]
        for table in expected_tables:
            assert f"`[{table}]`" in inventory, (path, table)


def test_parameter_reference_preserves_schema_coverage_and_toml_hierarchy() -> None:
    schema = _schema()
    documented_objects = {
        "sim": schema["$defs"]["sim"]["properties"],
        "species": schema["$defs"]["species"]["properties"],
        "mesh": schema["$defs"]["mesh"]["properties"],
        "mesh.groups": schema["$defs"]["meshGroup"]["properties"],
        "mesh.templates": schema["$defs"]["template"]["properties"],
        "periodic2": schema["properties"]["periodic2"]["properties"],
        "domain": schema["$defs"]["domain"]["properties"],
        "field_boundary": schema["$defs"]["fieldBoundary"]["properties"],
        "particle_boundary": schema["$defs"]["particleBoundary"]["properties"],
        "reservoir": schema["$defs"]["reservoir"]["properties"],
        "species.boundary": schema["$defs"]["speciesParticleBoundary"]["properties"],
        "output": schema["$defs"]["output"]["properties"],
    }
    expected_headings = {
        "docs/Parameters.md": (
            "### `[sim]`:",
            "### `[domain]`:",
            "### `[field_boundary]`:",
            "### `[particle_boundary]`:",
            "### `[reservoir]`:",
            "### `[periodic2]`:",
            "### `[[particles.species]]`:",
            "### `[mesh]`:",
            "#### `[[mesh.templates]]`:",
            "### 要素 source の固定規則",
            "### `[output]`:",
        ),
        "docs/Parameters.en.md": (
            "### `[sim]`:",
            "### `[domain]`:",
            "### `[field_boundary]`:",
            "### `[particle_boundary]`:",
            "### `[reservoir]`:",
            "### `[periodic2]`:",
            "### `[[particles.species]]`:",
            "### `[mesh]`:",
            "#### `[[mesh.templates]]`:",
            "### Fixed element-source rules",
            "### `[output]`:",
        ),
    }
    structural_markers = {
        ("mesh", "groups"): "`[mesh.groups.<name>]`",
        ("mesh", "templates"): "`[[mesh.templates]]`",
    }

    for path, headings in expected_headings.items():
        text = _read(path)
        for object_name, properties in documented_objects.items():
            for key in properties:
                structural_marker = structural_markers.get((object_name, key))
                if structural_marker is not None:
                    assert structural_marker in text, (path, object_name, key)
                    continue
                markers = (
                    f"`{key}`",
                    f".{key}`",
                    f"`{key}=",
                    f".{key}=",
                )
                assert any(marker in text for marker in markers), (
                    path,
                    object_name,
                    key,
                )
        for heading in headings:
            assert heading in text, (path, heading)
        assert "├── [particles]" in text
        assert "│   └── [[particles.species]]" in text
        assert "│   ├── [mesh.groups.<name>]" in text
        assert "│   └── [[mesh.templates]]" in text
        assert "├── [domain]" in text
        assert "├── [field_boundary]" in text
        assert "├── [particle_boundary]" in text
        assert "├── [reservoir]" in text
        assert "│       └── [particles.species.boundary]" in text


def test_parameter_editor_schema_uses_only_github_raw_url() -> None:
    directive = (
        "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/"
        "schemas/beach.schema.json"
    )
    for path in ("docs/Parameters.md", "docs/Parameters.en.md"):
        text = _read(path)
        assert text.count("#:schema") == 1, path
        assert directive in text
        assert "#:schema ../schemas/beach.schema.json" not in text


def test_parameter_reference_prose_paragraphs_remain_scannable() -> None:
    excluded_prefixes = ("#", "|", "```", "- ")

    for path in ("docs/Parameters.md", "docs/Parameters.en.md"):
        for block in _read(path).split("\n\n"):
            paragraph = block.strip()
            if not paragraph or paragraph.startswith(excluded_prefixes):
                continue
            if paragraph[0].isdigit() and paragraph[1:3] in {". ", ") "}:
                continue
            normalized = " ".join(paragraph.splitlines())
            assert len(normalized) <= 450, (path, normalized[:120])


def test_direct_periodic2_split_reference_matches_schema_runtime_and_docs() -> None:
    schema = _schema()
    validator = Draft7Validator(schema)

    for path in ("examples/periodic2_closed_photoelectron.toml",):
        document = tomllib.loads(_read(path))
        errors = sorted(
            validator.iter_errors(document), key=lambda error: list(error.path)
        )
        assert not errors, (path, [error.message for error in errors])
        load_config_file(ROOT / path)

    description = schema["$defs"]["fieldBoundary"]["properties"]["mode"]["description"]
    assert "direct triangle_p0 panel_spectral_reference split model" in description

    for path in ("docs/FieldSolvers.md", "docs/FieldSolvers.en.md"):
        text = _read(path)
        assert "panel_spectral_reference" in text
        assert "exclude_k0" in text
        assert "triangle_p0" in text

    assert "periodic2 は fmm 必須" not in _read("docs/agent-user-guide.md")
    assert "periodic2 requires fmm" not in _read("docs/agent-user-guide.en.md")


def test_output_manifest_matches_implementation_and_bilingual_docs() -> None:
    manifest = json.loads(_read("schemas/beach.output-manifest.json"))
    files = manifest["files"]
    names = [entry["name"] for entry in files]
    assert len(names) == len(set(names))
    matching_history = next(
        entry for entry in files if entry["name"] == "matching_plane_history.csv"
    )
    assert "matching_plane_quasistatic" in matching_history["condition"]
    assert matching_history["restart_role"] == "none"

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
