from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path


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
        assert r"\(" not in text, page.source
        assert r"\)" not in text, page.source

        fence_marker: str | None = None
        in_display_math = False
        for line_number, line in enumerate(text.splitlines(), start=1):
            fence_match = FENCE_RE.match(line)
            if fence_match:
                marker = fence_match.group(1)[0]
                if fence_marker is None:
                    fence_marker = marker
                elif fence_marker == marker:
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
            assert UNESCAPED_DOLLAR_RE.search(prose) is None, location
            assert TEX_COMMAND_RE.search(prose) is None, location

        assert fence_marker is None, f"{page.source}: unclosed code fence"
        assert not in_display_math, f"{page.source}: unclosed display math"


def test_japanese_pages_avoid_repetitive_reference_phrasing() -> None:
    module = _load_sync_module()
    japanese_pages = [page for page in module.PAGES if page.locale == "root"]

    for page in japanese_pages:
        assert _read_doc(page.source).count("参照してください") <= 1, page.source

    for name in (
        "BorisPusher.md",
        "FinitePeriodicConfiguration.md",
        "FMM.md",
        "PeriodicElectrostatics.md",
    ):
        assert "[<sup>1</sup>](" in _read_doc(name)


def test_japanese_pages_name_the_batch_fixed_field_concretely() -> None:
    module = _load_sync_module()
    ambiguous_phrases = (
        "field snapshot",
        "electrostatic snapshot",
        "immutable snapshot",
        "固定snapshot",
        "snapshot更新",
    )

    for page in (page for page in module.PAGES if page.locale == "root"):
        text = _read_doc(page.source)
        for phrase in ambiguous_phrases:
            assert phrase not in text, f"{page.source}: {phrase}"


def test_particle_escape_return_starts_from_boundary_ownership() -> None:
    expected_headings = {
        "ParticleEscapeReturn.md": (
            "## 1. `escape`: open面で粒子を除去する",
            "## 2. `potential_barrier`: scalar障壁で反射を判定する",
            "## 3. `linear_debye`: 解析的な1D profileでreturnを写像する",
            "## 4. `kinetic_1d`: 離散sheath profileでreturnを求める",
            "## 5. `unified_linear_response`: 外部3D軌道を積分する",
        ),
        "ParticleEscapeReturn.en.md": (
            "## 1. `escape`: remove a particle at an open face",
            "## 2. `potential_barrier`: decide reflection at a scalar barrier",
            "## 3. `linear_debye`: map return through an analytic 1-D profile",
            "## 4. `kinetic_1d`: obtain return from a discrete sheath profile",
            "## 5. `unified_linear_response`: integrate an external 3-D orbit",
        ),
    }

    for name, headings in expected_headings.items():
        text = _read_doc(name)
        positions = [text.index(heading) for heading in headings]
        assert positions == sorted(positions)
        assert "external_boundary.ordinary_open.model" in text
        assert "external_boundary.field.model" in text
        assert "external_boundary.particles.mode" in text
        assert "sim.open_boundary_model" not in text
        assert "particle_transfer_mode" not in text


def test_config_init_docs_match_official_tutorial_case() -> None:
    for name in ("Configuration.md", "Configuration.en.md"):
        text = _read_doc(name)
        assert "examples/tutorial_insulator.toml" in text
        assert 'field_solver="fmm"' in text
        assert 'field_bc_mode="periodic2"' in text
        assert "field_periodic_image_layers=1" in text
        assert 'field_periodic_far_correction="none"' in text
        assert "run_periodic2" not in text

    for name in ("Workflow.md", "Workflow.en.md"):
        text = _read_doc(name)
        assert "mkdir beach-tutorial" in text
        assert "run_periodic2" not in text


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
    navigation = json.loads(
        (ROOT / "docs-site" / "navigation.json").read_text(encoding="utf-8")
    )
    labels = [section["label"]["root"] for section in navigation["sections"]]

    assert labels == [
        "はじめる",
        "シミュレーションする",
        "リファレンス",
        "モデルと数値手法",
        "開発・検証",
    ]

    slugs = {page["slug"] for page in navigation["pages"]}
    assert {
        "installation",
        "tutorial",
        "validation-guide",
        "troubleshooting",
        "agent-user-guide",
    } <= slugs
    assert "physics-redesign-completion-audit" not in slugs
    assert any(
        item.get("slug") == "agent-user-guide"
        for item in navigation["sections"][-1]["items"]
    )


def test_case_design_and_configuration_pages_have_distinct_tasks() -> None:
    expected_titles = {
        "ConfigurationRecipes.md": (
            "title: シミュレーションケースを設計する",
            "# シミュレーションケースを設計する",
        ),
        "ConfigurationRecipes.en.md": (
            "title: Design a Simulation Case",
            "# Design a Simulation Case",
        ),
        "Configuration.md": (
            "title: beach.tomlを作成・検証する",
            "# `beach.toml`を作成・検証する",
        ),
        "Configuration.en.md": (
            "title: Create and Validate beach.toml",
            "# Create and Validate `beach.toml`",
        ),
    }

    for name, phrases in expected_titles.items():
        text = _read_doc(name)
        for phrase in phrases:
            assert phrase in text

    navigation = json.loads(
        (ROOT / "docs-site" / "navigation.json").read_text(encoding="utf-8")
    )
    simulation_items = navigation["sections"][1]["items"]
    labels = {
        item["slug"]: item["label"]
        for item in simulation_items
        if item.get("slug") in {"configuration-recipes", "configuration"}
    }
    assert labels == {
        "configuration-recipes": {
            "root": "ケースを設計する",
            "en": "Design a case",
        },
        "configuration": {
            "root": "設定ファイルを作成・検証する",
            "en": "Create and validate configuration",
        },
    }


def test_output_lookup_and_validation_pages_have_distinct_tasks() -> None:
    output_pages = {
        "OutputGuide.md": (
            "title: 出力ファイルを調べる",
            "# 出力ファイルを調べる",
            "[計算結果の妥当性確認](ValidationGuide.html)",
        ),
        "OutputGuide.en.md": (
            "title: Inspect Output Files",
            "# Inspect Output Files",
            "[Validating Simulation Results](ValidationGuide.en.html)",
        ),
    }
    for name, phrases in output_pages.items():
        text = _read_doc(name)
        for phrase in phrases:
            assert phrase in text
        assert "## 成功と注意の読み分け" not in text
        assert "## Interpreting Success and Warnings" not in text

    assert "その値を使って計算を受理できるか判定します" in _read_doc(
        "ValidationGuide.md"
    )
    assert "uses those values to decide whether a run" in _read_doc(
        "ValidationGuide.en.md"
    )

    navigation = json.loads(
        (ROOT / "docs-site" / "navigation.json").read_text(encoding="utf-8")
    )
    output_item = next(
        item
        for item in navigation["sections"][1]["items"]
        if item.get("slug") == "output-guide"
    )
    assert output_item["label"] == {
        "root": "出力ファイルを調べる",
        "en": "Inspect output files",
    }


def test_generated_pages_show_development_status_freshness_and_edit_source() -> None:
    module = _load_sync_module()

    for locale in ("root", "en"):
        page = next(
            page
            for page in module.PAGES
            if page.locale == locale and page.slug == "parameters"
        )
        _, content = module.render_page(page)
        assert "lastUpdated:" in content
        assert f"edit/main/docs/{page.source}" in content
        assert "banner:" in content
        if locale == "en":
            assert "Development documentation" in content
        else:
            assert "開発版ドキュメント" in content


def test_configuration_recipes_prioritize_meshes_and_particle_sources() -> None:
    for name in ("ConfigurationRecipes.md", "ConfigurationRecipes.en.md"):
        text = _read_doc(name)

        for kind in (
            "plane",
            "plate_hole",
            "plane_hole",
            "disk",
            "annulus",
            "box",
            "cylinder",
            "sphere",
        ):
            assert f"`{kind}`" in text
        assert 'kind = "sphere"' in text

        for source_mode in ("volume_seed", "reservoir_face", "photo_raycast"):
            assert f'`{source_mode}`' in text
            assert f'source_mode = "{source_mode}"' in text
        assert "target_macro_particles_per_batch" in text
        assert "rays_per_batch" in text
        assert "sim.batch_duration" in text


def test_configuration_recipes_delegate_advanced_outer_coupling() -> None:
    expected_links = {
        "ConfigurationRecipes.md": (
            "InfinitePeriodicOuterConfiguration.html",
            "PhotoelectronEmission.html",
        ),
        "ConfigurationRecipes.en.md": (
            "InfinitePeriodicOuterConfiguration.en.html",
            "PhotoelectronEmission.en.html",
        ),
    }

    for name, links in expected_links.items():
        text = _read_doc(name)

        for link in links:
            assert link in text
        assert "examples/periodic2_kinetic_outer.toml" in text
        assert 'field_periodic_far_correction = "cached_kneq0"' not in text
        assert 'photoelectron_density_model = "kinetic_mean"' not in text


def test_split_detail_pages_cover_migrated_numerics_topics() -> None:
    legacy_pages = (
        "ParticleChargeLoop.md",
        "ParticleChargeLoop.en.md",
        "PeriodicZeroModeOuterPlasma.md",
        "PeriodicZeroModeOuterPlasma.en.md",
        "SheathReservoirBoundary.md",
        "SheathReservoirBoundary.en.md",
    )
    for name in legacy_pages:
        assert not (ROOT / "docs" / name).exists()

    pages = {
        name: _read_doc(name)
        for name in (
            "Algorithms.md",
            "Algorithms.en.md",
            "BorisPusher.md",
            "BorisPusher.en.md",
            "FMMCore.md",
            "FMMCore.en.md",
            "PeriodicElectrostatics.md",
            "PeriodicElectrostatics.en.md",
            "KineticOuterPlasma.md",
            "KineticOuterPlasma.en.md",
            "UnifiedLinearResponse.md",
            "UnifiedLinearResponse.en.md",
            "ReservoirInjection.md",
            "ReservoirInjection.en.md",
            "ParticleEscapeReturn.md",
            "ParticleEscapeReturn.en.md",
            "SheathInjectionClosures.md",
            "SheathInjectionClosures.en.md",
            "PhotoelectronEmission.md",
            "PhotoelectronEmission.en.md",
        )
    }

    for phrase in (
        "## n batchの計算フロー",
        "[Boris粒子更新](BorisPusher.html)",
        "[粒子の衝突・境界イベント](ParticleEvents.html)",
        "[表面電荷更新](SurfaceModels.html)",
    ):
        assert phrase in pages["Algorithms.md"]
    for phrase in (
        "## The n-batch computation flow",
        "[Boris particle update](BorisPusher.en.html)",
        "[Particle collision and boundary events](ParticleEvents.en.html)",
        "[Surface charge update](SurfaceModels.en.html)",
    ):
        assert phrase in pages["Algorithms.en.md"]

    assert "## 予測中点で電場を評価する" in pages["BorisPusher.md"]
    assert "## 台形則で候補位置を作る" in pages["BorisPusher.md"]
    assert "## Sample the field at a predicted midpoint" in pages["BorisPusher.en.md"]
    assert "## Form the candidate position with the trapezoidal rule" in pages[
        "BorisPusher.en.md"
    ]

    for heading in (
        "何を高速化するoperatorか",
        "1回のfield評価で何を足すか",
        "数式との対応",
        "cache lifecycle",
        "MPI/OpenMPによるcold build",
        "cold buildとwarm runの違い",
        "SysA測定値",
        "運用指針",
    ):
        assert f"#### {heading}" in pages["FMMCore.md"]
    for heading in (
        "What the operator accelerates",
        "What one field evaluation adds",
        "Relation to the formula",
        "Cache lifecycle",
        "MPI/OpenMP cold build",
        "Cold versus warm execution",
        "SysA measurements",
        "Operating guidance",
    ):
        assert f"#### {heading}" in pages["FMMCore.en.md"]

    migrated_headings = {
        "PeriodicElectrostatics.md": (
            "## 場を4つの成分に分ける",
            "## Ewald2Pで無限周期の遠方場を分離する",
            "## 物理`k=0`を一度だけ加える",
        ),
        "KineticOuterPlasma.md": (
            "## VDFを電位依存の電荷密度へ写す",
            "## continuation付きNewton法で物理解を追う",
        ),
        "UnifiedLinearResponse.md": (
            "## 表面電荷とplasma応答からzero modeを解く",
            "## 真空nonzero modeをscreened tailへ接続する",
        ),
        "ReservoirInjection.md": (
            "## Maxwell 分布を流入流束で重み付けする",
            "## 1 つの電位差で到達条件と注入面速度を決める",
        ),
        "ParticleEscapeReturn.md": (
            "## linear Debye profileからreturn時間を求める",
            "## outer flightをglobal timeへ加えない近似",
        ),
        "SheathInjectionClosures.md": (
            "## `floating_no_photo`でelectron/ion流入を釣り合わせる",
            "## Zhao closureの無次元量を作る",
        ),
        "PhotoelectronEmission.md": (
            "## 放出から再吸収までを同じbatchで追う",
            "## 放出・再吸収・escapeの電荷収支を確認する",
        ),
    }
    for name, headings in migrated_headings.items():
        for heading in headings:
            assert heading in pages[name]

    migrated_headings_en = {
        "PeriodicElectrostatics.en.md": (
            "## Decompose the field into four components",
            "## Separate the infinite-periodic far field with Ewald2P",
            "## Add the physical `k=0` component exactly once",
        ),
        "KineticOuterPlasma.en.md": (
            "## Map VDFs to potential-dependent charge density",
            "## Follow the physical branch with continued Newton solves",
        ),
        "UnifiedLinearResponse.en.md": (
            "## Solve the zero mode from surface charge and plasma response",
            "## Connect vacuum nonzero modes to a screened tail",
        ),
        "ReservoirInjection.en.md": (
            "## Weight a Maxwell distribution by inflow flux",
            "## Use one potential drop for accessibility and face velocity",
        ),
        "ParticleEscapeReturn.en.md": (
            "## Derive return time from a linear-Debye profile",
            "## Keep outer flight outside global simulation time",
        ),
        "SheathInjectionClosures.en.md": (
            "## Balance electron and ion inflow with `floating_no_photo`",
            "## Form the dimensionless Zhao variables",
        ),
        "PhotoelectronEmission.en.md": (
            "## Track emission through reabsorption in the same batch",
            "## Check charge balance across emission, reabsorption, and escape",
        ),
    }
    for name, headings in migrated_headings_en.items():
        for heading in headings:
            assert heading in pages[name]
