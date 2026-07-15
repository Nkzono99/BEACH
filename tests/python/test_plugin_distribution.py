from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PLUGIN = ROOT / "plugins" / "beach-context"
CLAUDE_COMPAT = ROOT / "plugins" / "beach-context-claude"


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def test_codex_plugin_manifest_and_marketplace_are_consistent() -> None:
    manifest = _load_json(PLUGIN / ".codex-plugin" / "plugin.json")
    marketplace = _load_json(ROOT / ".agents" / "plugins" / "marketplace.json")

    assert manifest["name"] == "beach-context"
    assert re.fullmatch(r"\d+\.\d+\.\d+(?:[-+][0-9A-Za-z.-]+)?", manifest["version"])
    entry = marketplace["plugins"][0]
    assert entry["name"] == manifest["name"]
    assert entry["source"] == {
        "source": "local",
        "path": "./plugins/beach-context",
    }
    assert entry["policy"] == {
        "installation": "AVAILABLE",
        "authentication": "ON_INSTALL",
    }


def test_claude_plugin_manifest_and_marketplace_are_consistent() -> None:
    manifest = _load_json(PLUGIN / ".claude-plugin" / "plugin.json")
    marketplace = _load_json(ROOT / ".claude-plugin" / "marketplace.json")

    assert manifest["name"] == "beach-context"
    assert "version" not in manifest
    assert marketplace["name"] == "beach-claude"
    assert marketplace["owner"]["name"]
    entry = marketplace["plugins"][0]
    assert entry["name"] == manifest["name"]
    assert entry["source"] == "./plugins/beach-context"
    assert (PLUGIN / "skills").is_dir()
    assert (PLUGIN / "agents").is_dir()


def test_claude_compatibility_agents_mirror_shared_plugin_agents() -> None:
    shared_agents = {path.name: path for path in (PLUGIN / "agents").glob("*.md")}
    compat_agents = {
        path.name: path for path in (CLAUDE_COMPAT / "agents").glob("*.md")
    }

    assert shared_agents
    assert compat_agents.keys() == shared_agents.keys()
    for name, shared_path in shared_agents.items():
        assert compat_agents[name].read_bytes() == shared_path.read_bytes()


def test_plugin_docs_use_current_install_and_config_init_contracts() -> None:
    readme = (PLUGIN / "README.md").read_text(encoding="utf-8")
    root_readme = (ROOT / "plugins" / "README.md").read_text(encoding="utf-8")
    usage = (PLUGIN / "docs" / "usage-workflows.md").read_text(encoding="utf-8")

    for text in (readme, root_readme):
        assert "codex plugin add beach-context@beach" in text
        assert "claude plugin install beach-context@beach-claude" in text
    assert 'field_periodic_image_layers=1' in usage
    assert 'field_periodic_far_correction="none"' in usage
