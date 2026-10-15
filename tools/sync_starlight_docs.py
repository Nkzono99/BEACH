#!/usr/bin/env python3
"""Prepare Starlight content from the canonical BEACH Markdown docs."""

from __future__ import annotations

import argparse
import json
import os
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
# Directories under docs/ that are copied to the site as static assets.
ASSET_DIRS = ("media", "images")
# Relative link destinations in Markdown links and HTML src attributes.
LINK_RE = re.compile(r"(?P<pre>\]\(|src=\")(?P<target>[^)\"\s]+)(?P<post>[)\"])")


@dataclass(frozen=True)
class Page:
    source: str
    slug: str
    order: int
    description: str


def load_navigation() -> tuple[tuple[Page, ...], dict[str, str]]:
    """Load the canonical page inventory shared with the Starlight sidebar.

    Returns the pages and a map from each page source (relative to ``docs/``) to its slug.
    """
    navigation = json.loads(NAVIGATION_FILE.read_text(encoding="utf-8"))
    pages = tuple(
        Page(entry["source"], entry["slug"], entry["order"], entry["description"])
        for entry in navigation["pages"]
    )
    links = {page.source: page.slug for page in pages}

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

    return pages, links


def iter_page_slugs(groups: list[dict[str, object]]):
    """Yield page slugs from arbitrarily nested sidebar groups."""
    for group in groups:
        for item in group.get("items", []):
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


def doc_href(slug: str) -> str:
    if slug == "index":
        return f"{SITE_BASE}/index.html"
    return f"{SITE_BASE}/{slug}.html"


def resolve_link(target: str, source: str) -> str:
    """Map a relative link in ``docs/<source>`` to its site or repository URL."""
    if "://" in target or target.startswith(("#", "/", "mailto:")):
        return target
    path, separator, fragment = target.partition("#")
    if not path:
        return target
    anchor = separator + fragment
    resolved = Path(os.path.normpath(DOCS_ROOT / Path(source).parent / path))
    try:
        in_docs = resolved.relative_to(DOCS_ROOT).as_posix()
    except ValueError:
        in_docs = None

    if in_docs in DOC_LINKS:
        return doc_href(DOC_LINKS[in_docs]) + anchor
    if in_docs is not None and in_docs.split("/", 1)[0] in ASSET_DIRS:
        return f"{SITE_BASE}/{in_docs}{anchor}"
    try:
        in_repo = resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return target
    return f"{GITHUB_BLOB_ROOT}/{in_repo}{anchor}"


def replace_links(text: str, source: str) -> str:
    """Rewrite relative Markdown and image destinations; inline code is never a link destination."""
    return LINK_RE.sub(
        lambda match: match.group("pre") + resolve_link(match.group("target"), source) + match.group("post"),
        text,
    )


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


def development_banner() -> str:
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
            f"  content: {json.dumps(development_banner(), ensure_ascii=False)}",
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
    body = strip_first_h1(body)
    body = replace_links(body, page.source)
    body = normalize_code_fences(body)
    content = render_frontmatter(page, title) + body
    return CONTENT_ROOT / f"{page.slug}.md", content


def infer_title(body: str) -> str | None:
    for line in body.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return None


def sync_assets() -> None:
    for name in ASSET_DIRS:
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
