from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SYNC_PATH = ROOT / "tools" / "sync_starlight_docs.py"
FENCE_RE = re.compile(r"^\s{0,3}(`{3,}|~{3,})")
INLINE_MATH_RE = re.compile(r"(?<!\\)\$(?!\$)(?:\\.|[^$\n])*?(?<!\\)\$")
UNESCAPED_DOLLAR_RE = re.compile(r"(?<!\\)\$")
TEX_COMMAND_RE = re.compile(r"\\[A-Za-z]+")
MARKDOWN_LINK_RE = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")


def _load_sync_module():
    spec = importlib.util.spec_from_file_location("beach_docs_sync", SYNC_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _read_doc(name: str) -> str:
    return (ROOT / "docs" / name).read_text(encoding="utf-8")


def _canonical_docs() -> list[Path]:
    """Return the Markdown pages published on the site (not ADRs or the manuscript)."""
    docs = ROOT / "docs"
    excluded = {"adr", "manuscript", "superpowers"}
    return [
        path
        for path in docs.rglob("*.md")
        if path.relative_to(docs).parts[0] not in excluded
    ]


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

    page = next(page for page in module.PAGES if page.slug == "parameters")
    target, content = module.render_page(page)
    assert target == module.CONTENT_ROOT / "parameters.md"
    # Uncommitted (for example, just moved) sources have no commit date yet.
    if module.source_last_updated(page.source):
        assert "lastUpdated:" in content
    assert f"editUrl: {module.GITHUB_EDIT_ROOT}/docs/{page.source}" in content
    description = json.dumps(page.description, ensure_ascii=False)
    assert f"description: {description}" in content
    assert f"sidebar:\n  order: {page.order}" in content
    assert "開発版ドキュメント" in content


def test_navigation_keeps_core_workflow_pages_and_readable_sources() -> None:
    module = _load_sync_module()
    slugs = set()

    for page in module.PAGES:
        slugs.add(page.slug)
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
        assert slug in slugs, slug


def test_navigation_publishes_every_canonical_doc() -> None:
    """Do not leave a canonical docs page outside the generated site."""

    module = _load_sync_module()
    published = {page.source for page in module.PAGES}
    canonical = {
        path.relative_to(module.DOCS_ROOT).as_posix()
        for path in _canonical_docs()
    }

    assert published == canonical


def test_canonical_document_links_resolve_to_repository_files() -> None:
    """Keep relative links valid on GitHub and in the generated site."""

    for source in _canonical_docs():
        text = source.read_text(encoding="utf-8")
        for raw_target in MARKDOWN_LINK_RE.findall(text):
            target = raw_target.strip("<>").split("#", maxsplit=1)[0]
            if not target or "://" in target or target.startswith(("/", "mailto:")):
                continue
            assert not target.endswith(".html"), (source.name, raw_target)
            assert (source.parent / target).exists(), (source.name, raw_target)


def test_render_page_rewrites_supported_links_and_fences(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_sync_module()
    docs_root = tmp_path / "docs"
    content_root = tmp_path / "content"
    (docs_root / "guide").mkdir(parents=True)
    (docs_root / "guide" / "Fixture.md").write_text(
        """title: Fixture title

# Fixture title

[local](../reference/Parameters.md#run) and [sibling](Other.md)
[SPEC](../../SPEC.md#fixture), [source](../../src/example.f90), [ADR](../adr/0001-example.md),
![image](../images/example.png) and <img src="../media/example.gif">.
The repository path `docs/OutputGuide.md` stays plain text, and [web](https://example.com) is kept.

```fortran
program example
end program example
```
""",
        encoding="utf-8",
    )
    monkeypatch.setattr(module, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(module, "DOCS_ROOT", docs_root)
    monkeypatch.setattr(module, "CONTENT_ROOT", content_root)
    monkeypatch.setattr(
        module,
        "DOC_LINKS",
        {"reference/Parameters.md": "parameters", "guide/Other.md": "other"},
    )
    monkeypatch.setattr(module, "source_last_updated", lambda _source: "2026-01-02")

    page = module.Page(
        source="guide/Fixture.md",
        slug="fixture",
        order=42,
        description="fixture description",
    )
    target, content = module.render_page(page)

    assert target == content_root / "fixture.md"
    assert "\n# Fixture title\n" not in content
    assert "[local](/BEACH/parameters.html#run)" in content
    assert "[sibling](/BEACH/other.html)" in content
    assert f"[SPEC]({module.GITHUB_BLOB_ROOT}/SPEC.md#fixture)" in content
    assert f"[source]({module.GITHUB_BLOB_ROOT}/src/example.f90)" in content
    assert f"[ADR]({module.GITHUB_BLOB_ROOT}/docs/adr/0001-example.md)" in content
    assert "![image](/BEACH/images/example.png)" in content
    assert '<img src="/BEACH/media/example.gif">' in content
    assert "[web](https://example.com)" in content
    assert "`docs/OutputGuide.md`" in content
    assert "```fortran-free-form" in content
