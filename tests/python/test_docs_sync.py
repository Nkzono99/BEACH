from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SYNC_PATH = ROOT / "tools" / "sync_starlight_docs.py"


def _load_sync_module():
    spec = importlib.util.spec_from_file_location("beach_docs_sync", SYNC_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_onboarding_pages_are_bilingual_with_localized_descriptions() -> None:
    module = _load_sync_module()
    expected = {
        "installation": {"root", "en"},
        "tutorial": {"root", "en"},
        "validation-guide": {"root", "en"},
        "troubleshooting": {"root", "en"},
        "physics-release-verification": {"root", "en"},
    }

    locales_by_slug: dict[str, set[str]] = {}
    for page in module.PAGES:
        locales_by_slug.setdefault(page.slug, set()).add(page.locale)
        if page.locale == "root":
            assert any("ぁ" <= char <= "龯" for char in page.description)

    for slug, locales in expected.items():
        assert locales_by_slug.get(slug) == locales


def test_sidebar_follows_user_workflow_and_separates_agents() -> None:
    config = (ROOT / "docs-site" / "astro.config.mjs").read_text(encoding="utf-8")
    labels = ["はじめに", "使い方", "リファレンス", "数値アルゴリズム", "開発者向け", "AIエージェント向け"]
    positions = [config.index(f"label: '{label}'") for label in labels]

    assert positions == sorted(positions)
    assert "{ slug: 'installation' }" in config
    assert "{ slug: 'tutorial' }" in config
    assert "{ slug: 'validation-guide' }" in config
    assert "{ slug: 'troubleshooting' }" in config
    assert config.index("{ slug: 'agent-user-guide' }") > config.index("label: 'AIエージェント向け'")
