from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from beach.runops_adapter import Adapter, agent_context_path


def _run_dir(tmp_path: Path, config: str = "[sim]\nbatch_count = 10\n") -> Path:
    run = tmp_path / "run"
    (run / "input").mkdir(parents=True)
    (run / "work" / "latest").mkdir(parents=True)
    (run / "input" / "beach.toml").write_text(config, encoding="utf-8")
    return run


def test_build_command() -> None:
    adapter = Adapter()
    assert adapter.name == "beach"
    runtime = SimpleNamespace(executable=Path("/opt/beach/bin/beach"))
    assert adapter.build_command(runtime, Path(".")) == [
        "/opt/beach/bin/beach",
        "input/beach.toml",
    ]
    assert (
        adapter.build_command(SimpleNamespace(executable=Path("")), Path("."))[0]
        == "beach"
    )


def test_detect_status_and_summary(tmp_path: Path) -> None:
    adapter = Adapter()
    run = _run_dir(tmp_path)
    assert adapter.detect_status(run) == "failed"  # work/ has content but no summary
    (run / "manifest.toml").write_text('[job]\njob_id = "42"\n', encoding="utf-8")
    (run / "work" / "stdout.42.log").write_text(
        "batch 3/10\nbatch 10/10\n", encoding="utf-8"
    )
    (run / "work" / "latest" / "summary.txt").write_text(
        "batches = 10\nmean = 1.5\n", encoding="utf-8"
    )
    assert adapter.detect_status(run) == "completed"
    summary = adapter.summarize(run)
    assert summary["batches"] == 10 and summary["mean"] == 1.5
    assert (summary["batch"], summary["batch_count"]) == (10, 10)
    (run / "work" / "stderr.42.log").write_text(
        "forrtl: severe error\n", encoding="utf-8"
    )
    assert adapter.detect_status(run) == "failed"


def test_agent_context_has_plugin_skills() -> None:
    path = agent_context_path()
    assert path is not None
    assert (path / "skills" / "beach-config-review" / "SKILL.md").is_file()
    assert (path / "references").is_dir()
    assert Adapter().agent_context_dir() == path


def test_render_inputs_with_runops(tmp_path: Path) -> None:
    pytest.importorskip("runops.adapters")
    import tomllib

    run = _run_dir(tmp_path)
    assert Adapter().render_inputs({"sim.batch_count": 200}, run) == [
        "input/beach.toml"
    ]
    config = tomllib.loads((run / "input" / "beach.toml").read_text(encoding="utf-8"))
    assert config["sim"]["batch_count"] == 200
    assert config["output"]["dir"] == "work/latest"
