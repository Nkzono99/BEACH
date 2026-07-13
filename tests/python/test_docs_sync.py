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

    assert "physics-redesign-completion-audit" not in locales_by_slug


def test_sidebar_follows_user_workflow_and_separates_agents() -> None:
    config = (ROOT / "docs-site" / "astro.config.mjs").read_text(encoding="utf-8")
    labels = [
        "はじめに",
        "使い方",
        "リファレンス",
        "数値アルゴリズム",
        "開発者向け",
        "AIエージェント向け",
    ]
    positions = [config.index(f"label: '{label}'") for label in labels]

    assert positions == sorted(positions)
    assert "{ slug: 'installation' }" in config
    assert "{ slug: 'tutorial' }" in config
    assert "{ slug: 'validation-guide' }" in config
    assert "{ slug: 'troubleshooting' }" in config
    assert "{ slug: 'physics-redesign-completion-audit' }" not in config
    assert config.index("{ slug: 'agent-user-guide' }") > config.index(
        "label: 'AIエージェント向け'"
    )


def test_configuration_recipes_cover_production_kinetic_outer_sheath() -> None:
    for name in ("ConfigurationRecipes.md", "ConfigurationRecipes.en.md"):
        text = (ROOT / "docs" / name).read_text(encoding="utf-8")

        assert 'field_periodic_far_correction = "cached_kneq0"' in text
        assert 'nonzero_mode_backend = "cached_kneq0"' in text
        assert 'model = "kinetic_1d"' in text
        assert 'return_model = "kinetic_1d_profile_return"' in text
        assert 'photoelectron_closure = "kinetic_mean"' in text
        assert 'sheath_injection_model = "none"' in text
        assert "examples/periodic2_kinetic_outer.toml" in text


def test_numerics_pages_explain_same_time_boris_and_structure_cached_operator() -> None:
    particle_ja = (ROOT / "docs" / "ParticleChargeLoop.md").read_text(encoding="utf-8")
    particle_en = (ROOT / "docs" / "ParticleChargeLoop.en.md").read_text(
        encoding="utf-8"
    )
    fmm_ja = (ROOT / "docs" / "FMMCore.md").read_text(encoding="utf-8")
    fmm_en = (ROOT / "docs" / "FMMCore.en.md").read_text(encoding="utf-8")
    algorithms_ja = (ROOT / "docs" / "Algorithms.md").read_text(encoding="utf-8")
    algorithms_en = (ROOT / "docs" / "Algorithms.en.md").read_text(encoding="utf-8")

    assert "Boris速度更新を組み込んだ同時刻状態の積分器" in particle_ja
    assert "完全な1ステップを単に「leapfrog」と呼ぶのは不正確" in particle_ja
    assert (
        "same-time state integrator containing a Boris velocity update" in particle_en
    )
    assert (
        "inaccurate to describe the complete BEACH step only as a leapfrog"
        in particle_en
    )
    for section in ("8.1", "8.2", "8.3", "8.4", "8.5"):
        assert f"### {section}" in particle_ja
        assert f"### {section}" in particle_en

    for heading in (
        "演算子の構成",
        "cache lifecycle",
        "MPI/OpenMPによるcold build",
        "SysA測定値",
        "運用指針",
    ):
        assert f"#### {heading}" in fmm_ja
    for heading in (
        "Operator composition",
        "Cache lifecycle",
        "MPI/OpenMP cold build",
        "SysA measurements",
        "Operating guidance",
    ):
        assert f"#### {heading}" in fmm_en
    assert fmm_ja.index("### 10.1 cached periodic nonzero operator") < fmm_ja.index(
        "### 11. 実装との対応"
    )
    assert fmm_en.index("### 10.1 Cached periodic nonzero operator") < fmm_en.index(
        "### 11. Implementation mapping"
    )
    for phrase in (
        "##### point source",
        "##### P0 triangle source",
        "q_i` はすでに要素総電荷",
        "`triangle_p0` では `softening=0` を強制",
    ):
        assert phrase in fmm_ja
    for phrase in (
        "##### Point sources",
        "##### P0 triangle sources",
        "`q_i` is already the total element charge",
        "`triangle_p0` enforces `softening=0`",
    ):
        assert phrase in fmm_en

    expected_ja = [
        "### 2.3 P0 triangle panel field kernel",
        "### 2.4 periodic2 split reference",
        "### 2.5 outer particle interface",
        "### 2.6 periodic2 collision mesh",
        "### 2.7 restart",
    ]
    expected_en = [
        "### 2.3 P0 triangle panel field kernel",
        "### 2.4 periodic2 split reference",
        "### 2.5 Outer particle interface",
        "### 2.6 periodic2 collision mesh",
        "### 2.7 Restart",
    ]
    assert [algorithms_ja.index(heading) for heading in expected_ja] == sorted(
        algorithms_ja.index(heading) for heading in expected_ja
    )
    assert [algorithms_en.index(heading) for heading in expected_en] == sorted(
        algorithms_en.index(heading) for heading in expected_en
    )
