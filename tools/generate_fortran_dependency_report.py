#!/usr/bin/env python3
"""Generate a Fortran source/dependency overview for BEACH documentation."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import subprocess
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOTS = (REPO_ROOT / "src", REPO_ROOT / "app")
FORTRAN_SUFFIXES = (".f90", ".F90")

MODULE_RE = re.compile(
    r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)(?P<name>[a-z]\w*)\b",
    re.IGNORECASE,
)
SUBMODULE_RE = re.compile(
    r"^\s*submodule\s*\(\s*(?P<parent>[^)]+?)\s*\)\s*(?P<name>[a-z]\w*)\b",
    re.IGNORECASE,
)
PROGRAM_RE = re.compile(r"^\s*program\s+(?P<name>[a-z]\w*)\b", re.IGNORECASE)
USE_RE = re.compile(
    r"^\s*use(?:\s*,\s*(?:non_intrinsic|intrinsic)\s*)?(?:\s*::)?\s*(?P<name>[a-z]\w*)\b",
    re.IGNORECASE,
)

GRAPH_COLORS = (
    "#dbeafe",
    "#dcfce7",
    "#fef3c7",
    "#fce7f3",
    "#e0e7ff",
    "#f3e8ff",
    "#fae8ff",
    "#ffe4e6",
    "#e0f2fe",
    "#ecfccb",
)


@dataclass(frozen=True)
class Entity:
    """Top-level documented Fortran entity."""

    name: str
    kind: str
    path: Path
    group: str
    summary: str
    raw_parent: str | None
    use_dependencies: tuple[str, ...]
    parent_dependencies: tuple[str, ...]

    @property
    def internal_dependencies(self) -> tuple[str, ...]:
        names = list(self.parent_dependencies) + list(self.use_dependencies)
        seen: set[str] = set()
        ordered: list[str] = []
        for name in names:
            key = name.lower()
            if key in seen:
                continue
            seen.add(key)
            ordered.append(name)
        return tuple(ordered)


def iter_fortran_files() -> list[Path]:
    """Return tracked Fortran source files under src/ and app/."""

    files: list[Path] = []
    for root in SOURCE_ROOTS:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.suffix not in FORTRAN_SUFFIXES:
                continue
            files.append(path)
    return sorted(files)


def strip_inline_comment(line: str) -> str:
    """Remove trailing inline comments from a line."""

    if "!" not in line:
        return line
    if line.lstrip().startswith("!"):
        return ""
    return line.split("!", 1)[0]


def logical_lines(lines: Iterable[str]) -> list[str]:
    """Join free-form continuation lines for use-statement parsing."""

    joined: list[str] = []
    buffer = ""

    for raw_line in lines:
        code = strip_inline_comment(raw_line.rstrip("\n")).rstrip()
        if not code:
            if buffer:
                joined.append(buffer.strip())
                buffer = ""
            continue

        part = code.lstrip()
        if part.startswith("&"):
            part = part[1:].lstrip()

        if buffer:
            buffer = f"{buffer} {part}".strip()
        else:
            buffer = part

        if buffer.endswith("&"):
            buffer = buffer[:-1].rstrip()
            continue

        joined.append(buffer.strip())
        buffer = ""

    if buffer:
        joined.append(buffer.strip())

    return joined


def normalize_summary(doc_lines: list[str]) -> str:
    """Collapse doc comment lines to a single readable summary."""

    if not doc_lines:
        return ""
    text = " ".join(line.strip() for line in doc_lines if line.strip())
    text = re.sub(r"\s+", " ", text).strip()
    return text


def entity_group(path: Path) -> str:
    """Return the directory bucket label for a source file."""

    rel_dir = path.relative_to(REPO_ROOT).parent
    return rel_dir.as_posix()


def parse_entities(path: Path) -> list[dict[str, str | None]]:
    """Extract top-level modules/submodules/programs from one file."""

    pending_doc: list[str] = []
    entities: list[dict[str, str | None]] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if stripped.startswith("!>"):
            pending_doc.append(stripped[2:].strip())
            continue
        if not stripped:
            continue
        if stripped.startswith("!"):
            pending_doc = []
            continue

        module_match = MODULE_RE.match(raw_line)
        if module_match:
            entities.append(
                {
                    "kind": "module",
                    "name": module_match.group("name"),
                    "raw_parent": None,
                    "summary": normalize_summary(pending_doc),
                }
            )
            pending_doc = []
            continue

        submodule_match = SUBMODULE_RE.match(raw_line)
        if submodule_match:
            entities.append(
                {
                    "kind": "submodule",
                    "name": submodule_match.group("name"),
                    "raw_parent": submodule_match.group("parent"),
                    "summary": normalize_summary(pending_doc),
                }
            )
            pending_doc = []
            continue

        program_match = PROGRAM_RE.match(raw_line)
        if program_match:
            entities.append(
                {
                    "kind": "program",
                    "name": program_match.group("name"),
                    "raw_parent": None,
                    "summary": normalize_summary(pending_doc),
                }
            )
            pending_doc = []
            continue

        pending_doc = []

    return entities


def parse_use_statements(path: Path) -> list[str]:
    """Extract use-dependencies from a file."""

    names: list[str] = []
    lines = logical_lines(path.read_text(encoding="utf-8").splitlines())
    for line in lines:
        match = USE_RE.match(line)
        if match:
            names.append(match.group("name"))
    return names


def build_entities() -> list[Entity]:
    """Build the source inventory with internal/external dependencies."""

    per_file: list[tuple[Path, list[dict[str, str | None]], list[str]]] = []
    internal_names: dict[str, str] = {}

    for path in iter_fortran_files():
        raw_entities = parse_entities(path)
        use_names = parse_use_statements(path)
        per_file.append((path, raw_entities, use_names))
        for entity in raw_entities:
            internal_names[entity["name"].lower()] = entity["name"]  # type: ignore[index]

    entities: list[Entity] = []
    for path, raw_entities, use_names in per_file:
        if not raw_entities:
            continue

        internal_use_names = []
        for name in use_names:
            internal = internal_names.get(name.lower())
            if internal is not None:
                internal_use_names.append(internal)

        for raw in raw_entities:
            parent_dependencies: list[str] = []
            raw_parent = raw["raw_parent"]
            if raw_parent:
                for part in raw_parent.split(":"):
                    candidate = internal_names.get(part.strip().lower())
                    if candidate is not None:
                        parent_dependencies.append(candidate)

            entity_name = str(raw["name"])
            normalized_uses = tuple(
                name
                for name in deduplicate(internal_use_names)
                if name.lower() != entity_name.lower()
            )
            normalized_parents = tuple(
                name
                for name in deduplicate(parent_dependencies)
                if name.lower() != entity_name.lower()
            )

            entities.append(
                Entity(
                    name=entity_name,
                    kind=str(raw["kind"]),
                    path=path.relative_to(REPO_ROOT),
                    group=entity_group(path),
                    summary=str(raw["summary"] or ""),
                    raw_parent=str(raw_parent) if raw_parent else None,
                    use_dependencies=normalized_uses,
                    parent_dependencies=normalized_parents,
                )
            )

    return sorted(entities, key=lambda entity: (entity.group, entity.kind, entity.name.lower()))


def deduplicate(names: Iterable[str]) -> list[str]:
    """Return a stable de-duplicated name list."""

    seen: set[str] = set()
    ordered: list[str] = []
    for name in names:
        key = name.lower()
        if key in seen:
            continue
        seen.add(key)
        ordered.append(name)
    return ordered


def external_dependencies(entities: list[Entity], path: Path) -> list[str]:
    """Return file-local dependencies that are not defined in this project."""

    internal_names = {entity.name.lower() for entity in entities}
    names = []
    for name in parse_use_statements(REPO_ROOT / path):
        if name.lower() not in internal_names:
            names.append(name)
    return deduplicate(names)


def graph_dot(entities: list[Entity]) -> str:
    """Render a DOT graph for module and program dependencies."""

    group_to_entities: dict[str, list[Entity]] = defaultdict(list)
    for entity in entities:
        group_to_entities[entity.group].append(entity)

    lines = [
        "digraph beach_fortran_dependencies {",
        '  graph [rankdir="LR", overlap="false", splines="true", fontname="Helvetica", fontsize="11"];',
        '  node [fontname="Helvetica", fontsize="10", shape="box", style="rounded,filled", fillcolor="#f8fafc", color="#334155"];',
        '  edge [fontname="Helvetica", fontsize="9", color="#64748b", arrowsize="0.8"];',
    ]

    groups = sorted(group_to_entities)
    for idx, group in enumerate(groups):
        color = GRAPH_COLORS[idx % len(GRAPH_COLORS)]
        cluster_name = slugify(group)
        lines.append(f'  subgraph "cluster_{cluster_name}" {{')
        lines.append(f'    label="{dot_escape(group)}";')
        lines.append('    style="rounded,filled";')
        lines.append(f'    color="{color}";')
        lines.append('    penwidth="1.2";')
        lines.append('    fillcolor="#ffffff";')
        for entity in sorted(group_to_entities[group], key=lambda item: item.name.lower()):
            shape = {"module": "box", "submodule": "component", "program": "ellipse"}[entity.kind]
            fill = {
                "module": "#f8fafc",
                "submodule": "#eff6ff",
                "program": "#fef3c7",
            }[entity.kind]
            tooltip = dot_escape(entity.summary or entity.path.as_posix())
            lines.append(
                f'    "{entity.name}" [shape="{shape}", fillcolor="{fill}", tooltip="{tooltip}"];'
            )
        lines.append("  }")

    for entity in entities:
        for dependency in entity.use_dependencies:
            lines.append(f'  "{entity.name}" -> "{dependency}";')
        for dependency in entity.parent_dependencies:
            lines.append(
                f'  "{entity.name}" -> "{dependency}" [style="dashed", color="#0f766e", label="parent"];'
            )

    lines.append("}")
    return "\n".join(lines) + "\n"


def slugify(text: str) -> str:
    """Create a DOT-safe cluster suffix."""

    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def dot_escape(text: str) -> str:
    """Escape a string for DOT output."""

    return text.replace("\\", "\\\\").replace('"', '\\"')


def markdown_report(entities: list[Entity]) -> str:
    """Render the dependency overview page."""

    file_count = len({entity.path for entity in entities})
    module_count = sum(entity.kind == "module" for entity in entities)
    submodule_count = sum(entity.kind == "submodule" for entity in entities)
    program_count = sum(entity.kind == "program" for entity in entities)
    internal_edge_count = sum(len(entity.internal_dependencies) for entity in entities)

    incoming = Counter()
    for entity in entities:
        for dependency in entity.internal_dependencies:
            incoming[dependency] += 1

    per_group_counts = Counter(entity.group for entity in entities)
    per_group_edges = Counter()
    for entity in entities:
        per_group_edges[entity.group] += len(entity.internal_dependencies)

    lines: list[str] = []
    lines.append("title: Fortran 依存関係マップ")
    lines.append("")
    lines.append("# Fortran 依存関係マップ")
    lines.append("")
    lines.append(
        "> このページは `tools/generate_fortran_dependency_report.py` から自動生成しています。"
    )
    lines.append("")
    lines.append("## 概要")
    lines.append("")
    lines.append(f"- ソースファイル数: {file_count}")
    lines.append(f"- モジュール数: {module_count}")
    lines.append(f"- submodule 数: {submodule_count}")
    lines.append(f"- program 数: {program_count}")
    lines.append(f"- 内部依存エッジ数: {internal_edge_count}")
    lines.append("")
    lines.append("## 全体グラフ")
    lines.append("")
    lines.append("![Fortran モジュール依存グラフ](../media/fortran_module_dependencies.svg)")
    lines.append("")
    lines.append("実線は `use` 依存、破線は `submodule(parent)` の親参照を表します。")
    lines.append("")
    lines.append("## ディレクトリ別サマリ")
    lines.append("")
    lines.append("| ディレクトリ | エンティティ数 | 内部依存数 |")
    lines.append("| --- | ---: | ---: |")
    for group in sorted(per_group_counts):
        lines.append(
            f"| `{group}` | {per_group_counts[group]} | {per_group_edges[group]} |"
        )
    lines.append("")
    lines.append("## 被依存の多いモジュール")
    lines.append("")
    lines.append("| エンティティ | kind | 被依存数 |")
    lines.append("| --- | --- | ---: |")
    for entity in sorted(
        entities, key=lambda item: (-incoming[item.name], item.name.lower())
    )[:10]:
        lines.append(f"| `{entity.name}` | `{entity.kind}` | {incoming[entity.name]} |")
    lines.append("")
    lines.append("## エンティティ一覧")
    lines.append("")
    lines.append("| エンティティ | kind | パス | 内部依存 | 概要 |")
    lines.append("| --- | --- | --- | --- | --- |")
    for entity in entities:
        deps = ", ".join(f"`{name}`" for name in entity.internal_dependencies) or "-"
        summary = entity.summary or "-"
        lines.append(
            f"| `{entity.name}` | `{entity.kind}` | `{entity.path.as_posix()}` | {deps} | {escape_pipes(summary)} |"
        )
    lines.append("")
    lines.append("## 詳細")
    lines.append("")
    for entity in entities:
        lines.append(f"### `{entity.name}`")
        lines.append("")
        lines.append(f"- kind: `{entity.kind}`")
        lines.append(f"- path: `{entity.path.as_posix()}`")
        lines.append(f"- group: `{entity.group}`")
        if entity.raw_parent:
            lines.append(f"- parent: `{entity.raw_parent}`")
        internal = ", ".join(f"`{name}`" for name in entity.internal_dependencies) or "なし"
        lines.append(f"- internal dependencies: {internal}")
        external = ", ".join(
            f"`{name}`" for name in external_dependencies(entities, entity.path)
        ) or "なし"
        lines.append(f"- external dependencies: {external}")
        if entity.summary:
            lines.append(f"- summary: {entity.summary}")
        lines.append("")

    return "\n".join(lines) + "\n"


def escape_pipes(text: str) -> str:
    """Escape markdown table separators."""

    return text.replace("|", "\\|")


def render_svg(dot_text: str, svg_path: Path) -> None:
    """Render DOT text to SVG with Graphviz."""

    dot_bin = shutil.which("dot")
    if dot_bin is None:
        raise RuntimeError("Graphviz 'dot' command is required to render the SVG graph.")

    svg_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [dot_bin, "-Tsvg", "-o", str(svg_path)],
        input=dot_text.encode("utf-8"),
        check=True,
    )


def write_text(path: Path, content: str) -> None:
    """Write UTF-8 text, creating parent directories as needed."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(
        description="Generate BEACH Fortran dependency docs for FORD/GitHub Pages."
    )
    parser.add_argument(
        "--markdown",
        type=Path,
        default=REPO_ROOT / "docs" / "fortran_dependency_map.md",
        help="Path to the generated Markdown page.",
    )
    parser.add_argument(
        "--dot",
        type=Path,
        default=REPO_ROOT / "build" / "fortran_module_dependencies.dot",
        help="Path to the generated DOT file.",
    )
    parser.add_argument(
        "--svg",
        type=Path,
        default=REPO_ROOT / "docs" / "media" / "fortran_module_dependencies.svg",
        help="Path to the rendered SVG graph.",
    )
    return parser.parse_args()


def main() -> int:
    """Generate the Markdown report and Graphviz assets."""

    args = parse_args()
    entities = build_entities()
    dot_text = graph_dot(entities)
    markdown_text = markdown_report(entities)

    write_text(args.dot, dot_text)
    render_svg(dot_text, args.svg)
    write_text(args.markdown, markdown_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
