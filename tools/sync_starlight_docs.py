#!/usr/bin/env python3
"""Prepare Starlight content from the canonical BEACH Markdown docs."""

from __future__ import annotations

import argparse
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
SITE_ROOT = REPO_ROOT / "docs-site"
CONTENT_ROOT = SITE_ROOT / "src" / "content" / "docs"
PUBLIC_ROOT = SITE_ROOT / "public"
GITHUB_BLOB_ROOT = "https://github.com/Nkzono99/BEACH/blob/main"
SITE_BASE = "/BEACH"


@dataclass(frozen=True)
class Page:
    source: str
    locale: str
    slug: str
    order: int
    description: str


PAGES = (
    Page("index.md", "root", "index", 10, "BEACH documentation landing page."),
    Page("index.en.md", "en", "index", 10, "BEACH documentation landing page."),
    Page("OutputGuide.md", "root", "output-guide", 20, "Guide to BEACH output files."),
    Page("OutputGuide.en.md", "en", "output-guide", 20, "Guide to BEACH output files."),
    Page("ConfigurationRecipes.md", "root", "configuration-recipes", 30, "Configuration recipes for common BEACH runs."),
    Page("ConfigurationRecipes.en.md", "en", "configuration-recipes", 30, "Configuration recipes for common BEACH runs."),
    Page("Parameters.md", "root", "parameters", 40, "BEACH input parameter reference."),
    Page("Parameters.en.md", "en", "parameters", 40, "BEACH input parameter reference."),
    Page("Configuration.md", "root", "configuration", 50, "beachx config and high-level notation guide."),
    Page("Configuration.en.md", "en", "configuration", 50, "beachx config and high-level notation guide."),
    Page("PostprocessTutorial.md", "root", "postprocess-tutorial", 60, "First post-processing steps for BEACH outputs."),
    Page("PostprocessTutorial.en.md", "en", "postprocess-tutorial", 60, "First post-processing steps for BEACH outputs."),
    Page("PythonPostprocessAPI.md", "root", "python-postprocess-api", 70, "Python post-processing API reference."),
    Page("PythonPostprocessAPI.en.md", "en", "python-postprocess-api", 70, "Python post-processing API reference."),
    Page("Algorithms.md", "root", "algorithms", 80, "BEACH numerical algorithm overview."),
    Page("Algorithms.en.md", "en", "algorithms", 80, "BEACH numerical algorithm overview."),
    Page("FieldSolvers.md", "root", "field-solvers", 90, "BEACH field solvers and field boundary conditions."),
    Page("FieldSolvers.en.md", "en", "field-solvers", 90, "BEACH field solvers and field boundary conditions."),
    Page("ParticleChargeLoop.md", "root", "particle-charge-loop", 100, "Particle tracking and charge accumulation loop."),
    Page("ParticleChargeLoop.en.md", "en", "particle-charge-loop", 100, "Particle tracking and charge accumulation loop."),
    Page("FMMCore.md", "root", "fmm-core", 110, "Coulomb FMM core details."),
    Page("FMMCore.en.md", "en", "fmm-core", 110, "Coulomb FMM core details."),
    Page("BatchDurationStability.md", "root", "batch-duration-stability", 120, "batch_duration stability and steady value."),
    Page("BatchDurationStability.en.md", "en", "batch-duration-stability", 120, "batch_duration stability and steady value."),
    Page("Workflow.md", "root", "workflow", 130, "BEACH execution and development workflow."),
    Page("Workflow.en.md", "en", "workflow", 130, "BEACH execution and development workflow."),
    Page("FortranDependencyMap.md", "root", "fortran-dependency-map", 140, "Generated Fortran dependency map."),
    Page("FortranDependencyMap.en.md", "en", "fortran-dependency-map", 140, "Generated Fortran dependency map."),
    Page("agent-user-guide.md", "root", "agent-user-guide", 150, "Guide for AI agents operating BEACH simulations."),
    Page("agent-user-guide.en.md", "en", "agent-user-guide", 150, "Guide for AI agents operating BEACH simulations."),
)


DOC_LINKS = {
    "index": "index",
    "Workflow": "workflow",
    "OutputGuide": "output-guide",
    "ConfigurationRecipes": "configuration-recipes",
    "PostprocessTutorial": "postprocess-tutorial",
    "agent-user-guide": "agent-user-guide",
    "Algorithms": "algorithms",
    "FieldSolvers": "field-solvers",
    "ParticleChargeLoop": "particle-charge-loop",
    "FMMCore": "fmm-core",
    "BatchDurationStability": "batch-duration-stability",
    "Configuration": "configuration",
    "Parameters": "parameters",
    "FortranDependencyMap": "fortran-dependency-map",
    "PythonPostprocessAPI": "python-postprocess-api",
}


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
    replacements: list[tuple[str, str]] = []
    for source_stem, slug in DOC_LINKS.items():
        current = doc_href(locale, slug)
        other = doc_href("root" if locale == "en" else "en", slug)

        replacements.extend(
            [
                (f"{source_stem}.en.html", current if locale == "en" else other),
                (f"{source_stem}.html", current if locale != "en" else other),
                (f"{source_stem}.en.md", current if locale == "en" else other),
                (f"{source_stem}.md", current if locale != "en" else other),
            ]
        )

    placeholders: list[tuple[str, str]] = []
    for index, (source, target) in enumerate(replacements):
        placeholder = f"@@BEACH_DOC_LINK_{index}@@"
        text = text.replace(source, placeholder)
        placeholders.append((placeholder, target))
    for placeholder, target in placeholders:
        text = text.replace(placeholder, target)
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


def render_frontmatter(title: str, description: str, order: int) -> str:
    return "\n".join(
        [
            "---",
            f"title: {json.dumps(title, ensure_ascii=False)}",
            f"description: {json.dumps(description, ensure_ascii=False)}",
            "sidebar:",
            f"  order: {order}",
            "---",
            "",
        ]
    )


def render_page(page: Page) -> tuple[Path, str]:
    source = DOCS_ROOT / page.source
    metadata, body = read_doc(source)
    title = metadata.get("title") or infer_title(body) or page.slug
    body = strip_language_line(body)
    body = strip_first_h1(body)
    body = replace_doc_links(body, page.locale)
    body = replace_repo_links(body, page.locale)
    body = normalize_code_fences(body)
    content = render_frontmatter(title, page.description, page.order) + body

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
