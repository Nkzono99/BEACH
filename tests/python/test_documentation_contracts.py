from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def _schema() -> dict[str, object]:
    return json.loads(_read("schemas/beach.schema.json"))


def _markdown_sections(
    text: str,
    heading_prefix: str,
    *,
    include_subsections: bool = True,
) -> str:
    lines = text.splitlines()
    outside_fence: list[bool] = []
    fence_marker: str | None = None
    for line in lines:
        stripped = line.lstrip()
        if fence_marker is None and stripped.startswith(("```", "~~~")):
            fence_marker = stripped[:3]
            outside_fence.append(False)
        elif fence_marker is not None:
            outside_fence.append(False)
            if stripped.startswith(fence_marker):
                fence_marker = None
        else:
            outside_fence.append(True)

    sections: list[str] = []
    for start, line in enumerate(lines):
        if not outside_fence[start] or not line.startswith(heading_prefix):
            continue
        level = len(line) - len(line.lstrip("#"))
        end = len(lines)
        for index in range(start + 1, len(lines)):
            if not outside_fence[index]:
                continue
            candidate = lines[index]
            candidate_level = len(candidate) - len(candidate.lstrip("#"))
            if (
                not candidate.startswith("#")
                or not candidate[candidate_level:].startswith(" ")
            ):
                continue
            if not include_subsections or candidate_level <= level:
                end = index
                break
        sections.append("\n".join(lines[start:end]))
    assert sections, heading_prefix
    return "\n".join(sections)


def test_parameter_reference_covers_schema_tables_keys_and_hierarchy() -> None:
    schema = _schema()
    documented_objects = {
        "sim": schema["$defs"]["sim"]["properties"],
        "domain": schema["$defs"]["domain"]["properties"],
        "field_boundary": schema["$defs"]["fieldBoundary"]["properties"],
        "particle_boundary": schema["$defs"]["particleBoundary"]["properties"],
        "reservoir": schema["$defs"]["reservoir"]["properties"],
        "surface_current_model": schema["$defs"]["surfaceCurrentModel"][
            "properties"
        ],
        "particles": schema["$defs"]["particles"]["properties"],
        "species": schema["$defs"]["species"]["properties"],
        "species.boundary": schema["$defs"]["speciesParticleBoundary"]["properties"],
        "species.boundary_inflow": schema["$defs"]["speciesBoundaryInflow"][
            "properties"
        ],
        "mesh": schema["$defs"]["mesh"]["properties"],
        "mesh.groups": schema["$defs"]["meshGroup"]["properties"],
        "mesh.templates": schema["$defs"]["template"]["properties"],
        "periodic2": schema["properties"]["periodic2"]["properties"],
        "output": schema["$defs"]["output"]["properties"],
    }
    structural_markers = {
        ("particles", "species"): "`[[particles.species]]`",
        ("mesh", "groups"): "`[mesh.groups.<name>]`",
        ("mesh", "templates"): "`[[mesh.templates]]`",
    }
    object_headings = {
        "sim": "### `[sim]`",
        "domain": "### `[domain]`",
        "field_boundary": "### `[field_boundary]`",
        "particle_boundary": "### `[particle_boundary]`",
        "reservoir": "### `[reservoir]`",
        "surface_current_model": "### `[surface_current_model]`",
        "species": "### `[[particles.species]]`",
        "species.boundary": "#### `[particles.species.boundary]`",
        "species.boundary_inflow": "#### `[particles.species.boundary_inflow]`",
        "mesh": "### `[mesh]`",
        "mesh.templates": "#### `[[mesh.templates]]`",
        "periodic2": "### `[periodic2]`",
        "output": "### `[output]`",
    }
    hierarchy_markers = (
        "├── [sim]",
        "├── [domain]",
        "├── [field_boundary]",
        "├── [particle_boundary]",
        "├── [reservoir]",
        "├── [surface_current_model]",
        "├── [particles]",
        "│   └── [[particles.species]]",
        "│       ├── [particles.species.boundary_inflow]",
        "│       └── [particles.species.boundary]",
        "├── [mesh]",
        "│   ├── [mesh.groups.<name>]",
        "│   └── [[mesh.templates]]",
        "├── [periodic2]",
        "└── [output]",
    )
    for path, detail_heading, helper_heading in (
        (
            "docs/Parameters.md",
            "## パラメータ詳細リファレンス",
            "## 座標・配置の補助パラメータ",
        ),
        (
            "docs/Parameters.en.md",
            "## Detailed Parameter Reference",
            "## Coordinate and Placement Helper Parameters",
        ),
    ):
        text = _read(path)
        inventory, separator, _ = text.partition(detail_heading)
        assert separator, path
        for table in schema["properties"]:
            assert f"`[{table}]`" in inventory, (path, table)

        object_sections = {
            object_name: _markdown_sections(
                text,
                heading,
                include_subsections=object_name != "mesh",
            )
            for object_name, heading in object_headings.items()
        }
        helper_section = _markdown_sections(text, helper_heading)
        helper_row_prefixes = {
            "domain": ("| `domain.",),
            "species": ("| `inject_region_mode`", "| `uv_low`"),
            "mesh.groups": ("| group ", "| Group "),
            "mesh.templates": ("| template ", "| Template "),
        }
        for object_name, prefixes in helper_row_prefixes.items():
            rows = "\n".join(
                line
                for line in helper_section.splitlines()
                if line.startswith(prefixes)
            )
            assert rows, (path, object_name)
            object_sections[object_name] = "\n".join(
                (object_sections.get(object_name, ""), rows)
            )

        for object_name, properties in documented_objects.items():
            parameter_rows = [
                line
                for line in object_sections.get(object_name, "").splitlines()
                if line.startswith("|") and line.count("|") >= 5
            ]
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
                assert any(
                    marker in row for marker in markers for row in parameter_rows
                ), (path, object_name, key)
        for marker in hierarchy_markers:
            assert marker in text, (path, marker)


def test_parameter_reference_uses_single_remote_editor_schema_directive() -> None:
    directive = (
        "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/"
        "schemas/beach.schema.json"
    )
    for path in ("docs/Parameters.md", "docs/Parameters.en.md"):
        text = _read(path)
        assert text.count("#:schema") == 1, path
        assert directive in text, path


def test_field_solver_docs_cover_the_direct_periodic2_split() -> None:
    schema = _schema()
    assert "direct" in schema["$defs"]["sim"]["properties"]["field_solver"]["enum"]
    assert (
        "panel_spectral_reference"
        in schema["properties"]["periodic2"]["properties"][
            "nonzero_mode_backend"
        ]["enum"]
    )
    assert (
        schema["properties"]["periodic2"]["properties"]["zero_mode_policy"][
            "const"
        ]
        == "exclude_k0"
    )

    for path in ("docs/FieldSolvers.md", "docs/FieldSolvers.en.md"):
        text = _read(path)
        direct_row = next(
            line for line in text.splitlines() if line.startswith("| `direct` |")
        )
        assert "panel_spectral_reference" in direct_row, path
        assert "exclude_k0" in direct_row, path
        assert "`triangle_p0`" in text, path


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

    reference_docs = (
        _read("docs/OutputReference.md"),
        _read("docs/OutputReference.en.md"),
    )
    for entry in files:
        assert {"name", "producer", "condition", "restart_role"} <= entry.keys()
        name = entry["name"]
        assert "output.write_files=true" in entry["condition"], name
        producer = _read(entry["producer"])
        assert name in producer, (name, entry["producer"])
        for text in reference_docs:
            assert name in text, name

        mpi_name = entry.get("mpi_name")
        if mpi_name:
            for text in reference_docs:
                assert mpi_name in text, mpi_name

        consumer_path = entry.get("consumer")
        if consumer_path:
            assert name in _read(consumer_path), (name, consumer_path)
