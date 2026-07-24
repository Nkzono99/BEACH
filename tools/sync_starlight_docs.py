#!/usr/bin/env python3
"""Prepare Starlight content from the canonical BEACH Markdown docs."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from functools import cache
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
SITE_ROOT = REPO_ROOT / "docs-site"
NAVIGATION_FILE = SITE_ROOT / "navigation.json"
CONTENT_ROOT = SITE_ROOT / "src" / "content" / "docs"
PUBLIC_ROOT = SITE_ROOT / "public"
GITHUB_BLOB_ROOT = "https://github.com/Nkzono99/BEACH/blob/main"
GITHUB_EDIT_ROOT = "https://github.com/Nkzono99/BEACH/edit/main"
CHANGELOG_URL = "https://github.com/Nkzono99/BEACH/blob/main/CHANGELOG.md"
SITE_BASE = "/BEACH"


@dataclass(frozen=True)
class Page:
    source: str
    locale: str
    slug: str
    order: int
    description: str


def load_navigation() -> tuple[tuple[Page, ...], dict[str, str]]:
    """Load the canonical page inventory shared with the Starlight sidebar."""
    navigation = json.loads(NAVIGATION_FILE.read_text(encoding="utf-8"))
    pages: list[Page] = []
    links: dict[str, str] = {}

    for entry in navigation["pages"]:
        slug = entry["slug"]
        order = entry["order"]
        for locale in ("root", "en"):
            source = entry["source"][locale]
            pages.append(
                Page(source, locale, slug, order, entry["description"][locale])
            )
            links[Path(source).stem.removesuffix(".en")] = slug

    page_slugs = {page.slug for page in pages}
    sidebar_page_slugs = {
        entry["slug"] for entry in navigation["pages"] if entry.get("sidebar", True)
    }
    nav_slugs = set(iter_page_slugs(navigation["sections"]))
    if sidebar_page_slugs != nav_slugs or not nav_slugs <= page_slugs:
        missing = sorted(sidebar_page_slugs - nav_slugs)
        unknown = sorted(nav_slugs - page_slugs)
        hidden = sorted(nav_slugs - sidebar_page_slugs)
        raise ValueError(
            "navigation page mismatch: "
            f"missing={missing}, unknown={unknown}, hidden_in_sidebar={hidden}"
        )

    return tuple(pages), links


def iter_page_slugs(groups: list[dict[str, object]]):
    """Yield page slugs from arbitrarily nested sidebar groups."""
    for group in groups:
        for item in group["items"]:
            if item["type"] == "page":
                yield item["slug"]
            elif item["type"] == "group":
                yield from iter_page_slugs([item])


PAGES, DOC_LINKS = load_navigation()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check that generated files are up to date.")
    return parser.parse_args()


def read_doc(path: Path) -> tuple[dict[str, str], str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    metadata: dict[str, str] = {}
    body_start = 0

    for index, line in enumerate(lines):
        if not line.strip():
            body_start = index + 1
            break
        if ":" not in line:
            break
        key, value = line.split(":", 1)
        if not key or any(char.isspace() for char in key):
            break
        metadata[key] = value.strip()
    else:
        body_start = len(lines)

    body = "\n".join(lines[body_start:]).strip() + "\n"
    return metadata, body


def strip_language_line(text: str) -> str:
    lines = [line for line in text.splitlines() if not line.startswith("Lang: ")]
    return "\n".join(lines).strip() + "\n"


def strip_first_h1(text: str) -> str:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if not line.strip():
            continue
        if line.startswith("# "):
            del lines[index]
            if index < len(lines) and not lines[index].strip():
                del lines[index]
        break
    return "\n".join(lines).strip() + "\n"


def replace_doc_links(text: str, locale: str) -> str:
    for source_stem, slug in DOC_LINKS.items():
        current = doc_href(locale, slug)
        other = doc_href("root" if locale == "en" else "en", slug)
        replacements = (
            (f"{source_stem}.en.html", current if locale == "en" else other),
            (f"{source_stem}.html", current if locale != "en" else other),
            (f"{source_stem}.en.md", current if locale == "en" else other),
            (f"{source_stem}.md", current if locale != "en" else other),
        )
        for source, target in replacements:
            # Rewrite Markdown destinations only. Plain-text repository paths such as
            # `docs/OutputGuide.md` must remain intact.
            pattern = rf"(?<=\()\.?/?{re.escape(source)}(?=(?:[?#][^)]*)?\))"
            text = re.sub(pattern, target, text)
    return text


def doc_href(locale: str, slug: str) -> str:
    if locale == "en":
        if slug == "index":
            return f"{SITE_BASE}/en.html"
        return f"{SITE_BASE}/en/{slug}.html"
    if slug == "index":
        return f"{SITE_BASE}/index.html"
    return f"{SITE_BASE}/{slug}.html"


def replace_repo_links(text: str, locale: str) -> str:
    text = text.replace("(../src/", f"({GITHUB_BLOB_ROOT}/src/")
    text = text.replace("(../app/", f"({GITHUB_BLOB_ROOT}/app/")
    text = text.replace("(../tests/", f"({GITHUB_BLOB_ROOT}/tests/")
    text = text.replace("(../schemas/", f"({GITHUB_BLOB_ROOT}/schemas/")
    text = text.replace("(../examples/", f"({GITHUB_BLOB_ROOT}/examples/")
    text = text.replace("(adr/", f"({GITHUB_BLOB_ROOT}/docs/adr/")

    text = text.replace("(../media/", f"({SITE_BASE}/media/")
    text = text.replace("(../images/", f"({SITE_BASE}/images/")
    text = text.replace("(media/", f"({SITE_BASE}/media/")
    text = text.replace("(images/", f"({SITE_BASE}/images/")

    text = text.replace("(../module/", f"({SITE_BASE}/fortran/module/")
    text = text.replace("(../proc/", f"({SITE_BASE}/fortran/proc/")
    text = text.replace("(../type/", f"({SITE_BASE}/fortran/type/")

    return text


def normalize_code_fences(text: str) -> str:
    return text.replace("```fortran\n", "```fortran-free-form\n")


@cache
def source_last_updated(source: str) -> str | None:
    """Return the last committed YYYY-MM-DD date for a canonical source file."""
    result = subprocess.run(
        ["git", "log", "-1", "--format=%cs", "--", f"docs/{source}"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    value = result.stdout.strip()
    return value if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value) else None


def development_banner(locale: str) -> str:
    if locale == "en":
        return (
            "<strong>Development documentation</strong> — This site follows the "
            f"<code>main</code> branch and may differ from an installed release. "
            f'<a href="{CHANGELOG_URL}">View the changelog.</a>'
        )
    return (
        "<strong>開発版ドキュメント</strong> — このサイトは<code>main</code>ブランチ向けで、"
        "インストール済みリリースと異なる場合があります。"
        f'<a href="{CHANGELOG_URL}">変更履歴を確認してください。</a>'
    )


def render_frontmatter(page: Page, title: str) -> str:
    lines = [
        "---",
        f"title: {json.dumps(title, ensure_ascii=False)}",
        f"description: {json.dumps(page.description, ensure_ascii=False)}",
        f"editUrl: {GITHUB_EDIT_ROOT}/docs/{page.source}",
    ]
    last_updated = source_last_updated(page.source)
    if last_updated:
        lines.append(f"lastUpdated: {last_updated}")
    lines.extend(
        [
            "banner:",
            f"  content: {json.dumps(development_banner(page.locale), ensure_ascii=False)}",
            "sidebar:",
            f"  order: {page.order}",
            "---",
            "",
        ]
    )
    return "\n".join(lines)


def render_page(page: Page) -> tuple[Path, str]:
    source = DOCS_ROOT / page.source
    metadata, body = read_doc(source)
    title = metadata.get("title") or infer_title(body) or page.slug
    body = strip_language_line(body)
    body = strip_first_h1(body)
    body = replace_doc_links(body, page.locale)
    body = replace_repo_links(body, page.locale)
    body = normalize_code_fences(body)
    content = render_frontmatter(page, title) + body

    if page.locale == "en":
        target = CONTENT_ROOT / "en" / f"{page.slug}.md"
    else:
        target = CONTENT_ROOT / f"{page.slug}.md"
    return target, content


def infer_title(body: str) -> str | None:
    for line in body.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return None


def sync_assets() -> None:
    for name in ("media", "images"):
        source = DOCS_ROOT / name
        target = PUBLIC_ROOT / name
        if target.exists():
            shutil.rmtree(target)
        if source.exists():
            shutil.copytree(source, target)


def clean_content() -> None:
    if CONTENT_ROOT.exists():
        shutil.rmtree(CONTENT_ROOT)
    CONTENT_ROOT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    args = parse_args()
    expected = [render_page(page) for page in PAGES]

    if args.check:
        missing_or_changed = []
        for target, content in expected:
            if not target.exists() or target.read_text(encoding="utf-8") != content:
                missing_or_changed.append(target)
        if missing_or_changed:
            for target in missing_or_changed:
                print(f"out of date: {target.relative_to(REPO_ROOT)}")
            return 1
        return 0

    clean_content()
    for target, content in expected:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8")
    sync_assets()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
