from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def _schema() -> dict[str, object]:
    schema = json.loads(_read("schemas/beach.schema.json"))
    return schema


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
    definitions = schema["$defs"]

    def resolve(node):
        while "$ref" in node:
            node = definitions[node["$ref"].rsplit("/", 1)[-1]]
        return node

    def public_tables(node, path=""):
        node = resolve(node)
        if node.get("type") == "array":
            node = resolve(node["items"])
        if node.get("type") != "object":
            return
        if path:
            yield path, node
        for key, child in node.get("properties", {}).items():
            yield from public_tables(child, f"{path}.{key}" if path else key)
        if isinstance(node.get("additionalProperties"), dict):
            yield from public_tables(node["additionalProperties"], path + ".<name>")

    objects = dict(public_tables(definitions["groupedConfig"]))
    groups = ("run", "domain", "mesh", "particles", "fields", "sheath", "output")

    def table_marker(path):
        return f"[[{path}]]" if path in {"particles.species", "mesh.parts"} else f"[{path}]"

    for path in ("docs/Parameters.md", "docs/Parameters.en.md"):
        text = _read(path)
        detail_heading = (
            "## Detailed Parameter Reference" if path.endswith(".en.md")
            else "## パラメータ詳細リファレンス"
        )
        inventory, separator, detail = text.partition(detail_heading)
        assert separator, path
        headings = [
            match.group(1)
            for match in re.finditer(r"^### `\[([a-z]+)\]`", detail, re.MULTILINE)
        ]
        assert headings == list(groups), path
        for index, group in enumerate(groups):
            branch = "└──" if index == len(groups) - 1 else "├──"
            assert f"{branch} [{group}]" in inventory, (path, group)

        for table, node in objects.items():
            # mesh.groups is a named-table container; each entry owns the keys.
            if table == "mesh.groups":
                continue
            marker = table_marker(table)
            assert f"`{marker}`" in inventory, (path, table)
            assert "── " + marker in inventory, (path, table)
            if table.startswith("particles.species"):
                level = 3 + table.removeprefix("particles.species").count(".")
            elif table == "mesh.groups.<name>":
                level = 4
            else:
                level = 3 + table.count(".")
            content = _markdown_sections(
                detail, "#" * level + f" `{marker}`",
                include_subsections=table == "mesh.parts",
            )
            rows = [line for line in content.splitlines() if line.startswith("|")]
            # Check the key column, not mentions in another parameter's explanation.
            documented_keys = {
                key
                for row in rows if row.count("|") >= 5
                for key in re.findall(r"`([^`]+)`", row.split("|")[1])
            }
            for key in node.get("properties", {}):
                child = f"{table}.{key}"
                if child in objects:
                    child_marker = table_marker(child + ".<name>" if child == "mesh.groups" else child)
                    assert f"`{child_marker}`" in detail, (path, table, key)
                else:
                    assert key in documented_keys, (path, table, key)


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
    fields = schema["$defs"]["groupedConfig"]["properties"]["fields"]["properties"]
    assert "direct" in fields["solver"]["properties"]["method"]["enum"]
    assert "panel_spectral_reference" in fields["periodic"]["properties"]["backend"]["enum"]

    for path in ("docs/FieldSolvers.md", "docs/FieldSolvers.en.md"):
        text = _read(path)
        direct_row = next(
            line for line in text.splitlines() if line.startswith("| `direct` |")
        )
        assert "panel_spectral_reference" in direct_row, path
        assert "fields.periodic.lower_boundary_model" in direct_row, path
        assert "`triangle_p0`" in text, path


def test_output_manifest_matches_implementation_and_bilingual_docs() -> None:
    manifest = json.loads(_read("schemas/beach.output-manifest.json"))
    files = manifest["files"]
    names = [entry["name"] for entry in files]
    assert len(names) == len(set(names))
    matching_history = next(
        entry for entry in files if entry["name"] == "matching_plane_history.csv"
    )
    assert "outflow_refresh_batches>0" in matching_history["condition"]
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
