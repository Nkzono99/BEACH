from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SYNC_PATH = ROOT / "tools" / "sync_starlight_docs.py"
FENCE_RE = re.compile(r"^\s{0,3}(`{3,}|~{3,})")
INLINE_MATH_RE = re.compile(r"(?<!\\)\$(?!\$)(?:\\.|[^$\n])*?(?<!\\)\$")
UNESCAPED_DOLLAR_RE = re.compile(r"(?<!\\)\$")
TEX_COMMAND_RE = re.compile(r"\\[A-Za-z]+")


def _load_sync_module():
    spec = importlib.util.spec_from_file_location("beach_docs_sync", SYNC_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _read_doc(name: str) -> str:
    return (ROOT / "docs" / name).read_text(encoding="utf-8")


def _strip_inline_code(line: str) -> str:
    result: list[str] = []
    index = 0
    while index < len(line):
        if line[index] != "`":
            result.append(line[index])
            index += 1
            continue

        marker_end = index
        while marker_end < len(line) and line[marker_end] == "`":
            marker_end += 1
        marker = line[index:marker_end]
        closing = line.find(marker, marker_end)
        if closing < 0:
            result.append(marker)
            index = marker_end
            continue
        index = closing + len(marker)

    return "".join(result)


def test_pages_sources_use_remark_math_inline_delimiters() -> None:
    module = _load_sync_module()

    for page in module.PAGES:
        text = _read_doc(page.source)

        fence_marker: str | None = None
        in_display_math = False
        for line_number, line in enumerate(text.splitlines(), start=1):
            fence_match = FENCE_RE.match(line)
            if fence_match:
                marker = fence_match.group(1)
                if fence_marker is None:
                    fence_marker = marker
                elif (
                    marker[0] == fence_marker[0]
                    and len(marker) >= len(fence_marker)
                    and not line[fence_match.end() :].strip()
                ):
                    fence_marker = None
                continue
            if fence_marker is not None:
                continue

            if line.strip() == "$$":
                in_display_math = not in_display_math
                continue
            if in_display_math:
                continue

            prose = _strip_inline_code(line)
            prose = INLINE_MATH_RE.sub("", prose)
            location = f"{page.source}:{line_number}"
            assert r"\(" not in prose, location
            assert r"\)" not in prose, location
            assert UNESCAPED_DOLLAR_RE.search(prose) is None, location
            assert TEX_COMMAND_RE.search(prose) is None, location

        assert fence_marker is None, f"{page.source}: unclosed code fence"
        assert not in_display_math, f"{page.source}: unclosed display math"


def test_generated_pages_show_development_status_freshness_and_edit_source() -> None:
    module = _load_sync_module()

    for locale in ("root", "en"):
        page = next(
            page
            for page in module.PAGES
            if page.locale == locale and page.slug == "parameters"
        )
        target, content = module.render_page(page)
        locale_dir = "en" if locale == "en" else ""
        assert target == module.CONTENT_ROOT / locale_dir / "parameters.md"
        assert "lastUpdated:" in content
        assert f"editUrl: {module.GITHUB_EDIT_ROOT}/docs/{page.source}" in content
        description = json.dumps(page.description, ensure_ascii=False)
        assert f"description: {description}" in content
        assert f"sidebar:\n  order: {page.order}" in content
        if locale == "en":
            assert "Development documentation" in content
        else:
            assert "開発版ドキュメント" in content


def test_navigation_keeps_core_workflow_pages_and_readable_sources() -> None:
    module = _load_sync_module()
    locales_by_slug: dict[str, set[str]] = {}

    for page in module.PAGES:
        locales_by_slug.setdefault(page.slug, set()).add(page.locale)
        assert (module.DOCS_ROOT / page.source).is_file(), page.source

    for slug in {
        "installation",
        "tutorial",
        "configuration",
        "execution",
        "output-guide",
        "validation-guide",
        "troubleshooting",
    }:
        assert locales_by_slug.get(slug) == {"root", "en"}, slug


def test_render_page_rewrites_supported_links_and_fences(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_sync_module()
    docs_root = tmp_path / "docs"
    content_root = tmp_path / "content"
    docs_root.mkdir()
    (docs_root / "Fixture.md").write_text(
        """title: Fixture title

Lang: [日本語](Fixture.md) | [English](Fixture.en.md)

# Fixture title

[local](Parameters.md#sim) and [English](Parameters.en.md)
[source](../src/example.f90), [ADR](adr/0001-example.md),
![image](images/example.png), and [module](../module/example.html).
The repository path `docs/OutputGuide.md` stays plain text.

```fortran
program example
end program example
```
""",
        encoding="utf-8",
    )
    monkeypatch.setattr(module, "DOCS_ROOT", docs_root)
    monkeypatch.setattr(module, "CONTENT_ROOT", content_root)
    monkeypatch.setattr(module, "source_last_updated", lambda _source: "2026-01-02")

    page = module.Page(
        source="Fixture.md",
        locale="root",
        slug="fixture",
        order=42,
        description="fixture description",
    )
    target, content = module.render_page(page)

    assert target == content_root / "fixture.md"
    assert "Lang:" not in content
    assert "\n# Fixture title\n" not in content
    assert "[local](/BEACH/parameters.html#sim)" in content
    assert "[English](/BEACH/en/parameters.html)" in content
    assert f"[source]({module.GITHUB_BLOB_ROOT}/src/example.f90)" in content
    assert f"[ADR]({module.GITHUB_BLOB_ROOT}/docs/adr/0001-example.md)" in content
    assert "![image](/BEACH/images/example.png)" in content
    assert "[module](/BEACH/fortran/module/example.html)" in content
    assert "`docs/OutputGuide.md`" in content
    assert "```fortran-free-form" in content


def test_configuration_recipes_cover_supported_meshes_and_particle_sources() -> None:
    schema = json.loads(
        (ROOT / "schemas" / "beach.schema.json").read_text(encoding="utf-8")
    )
    mesh_kinds = schema["$defs"]["template"]["properties"]["kind"]["enum"]
    source_modes = schema["$defs"]["species"]["properties"]["source_mode"]["enum"]

    for name in ("ConfigurationRecipes.md", "ConfigurationRecipes.en.md"):
        text = _read_doc(name)
        mesh_table = text.split("| `kind` |", maxsplit=1)[1].split(
            "\n\n", maxsplit=1
        )[0]
        source_table = text.split("| `source_mode` |", maxsplit=1)[1].split(
            "\n\n", maxsplit=1
        )[0]

        for kind in mesh_kinds:
            assert f"`{kind}`" in mesh_table, (name, kind)
        for source_mode in source_modes:
            assert f"`{source_mode}`" in source_table, (name, source_mode)
