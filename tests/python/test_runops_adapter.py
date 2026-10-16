from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from beach.runops_adapter import Adapter, agent_context_path


GROUPED = """[run]
batch_count = 10

[[mesh.parts]]
kind = "obj"
path = "part.obj"
surface_model = "insulator"
"""


def _run_dir(tmp_path: Path, config: str = GROUPED) -> Path:
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


def test_case_template_is_the_grouped_default() -> None:
    from beach.config import default_config, dump_beach_toml
    from beach.config.schema import load_schema, schema_errors

    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python 3.10
        import tomli as tomllib
    load_toml_text = tomllib.loads

    text = Adapter().case_template()["beach.toml"]
    assert text.startswith("#:schema ")
    assert "runo new CASE --set run.batch_count=200" in text
    schema, _ = load_schema()
    assert schema_errors(load_toml_text(text), schema) == []
    assert load_toml_text(text) == load_toml_text(dump_beach_toml(default_config()))


def test_render_inputs_with_runops(tmp_path: Path) -> None:
    pytest.importorskip("runops.adapters")
    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python 3.10
        import tomli as tomllib

    run = _run_dir(tmp_path)
    (run / "input" / "part.obj").write_text("v 0 0 0\n", encoding="utf-8")
    absolute = tmp_path / "shared" / "plate.obj"
    absolute.parent.mkdir()
    absolute.write_text("v 1 0 0\n", encoding="utf-8")
    config = (run / "input" / "beach.toml").read_text(encoding="utf-8")
    config += f'\n[[mesh.parts]]\nkind = "obj"\npath = "{absolute}"\nsurface_model = "insulator"\n'
    (run / "input" / "beach.toml").write_text(config, encoding="utf-8")

    written = Adapter().render_inputs({"run.batch_count": 200}, run)
    assert written == ["input/beach.toml", "input/part.obj", "input/plate.obj"]
    rendered = tomllib.loads((run / "input" / "beach.toml").read_text(encoding="utf-8"))
    assert rendered["run"]["batch_count"] == 200
    assert rendered["output"]["dir"] == "work/latest"
    assert [p["path"] for p in rendered["mesh"]["parts"]] == [
        "input/part.obj",
        "input/plate.obj",
    ]
    assert (run / "input" / "plate.obj").read_text(encoding="utf-8") == "v 1 0 0\n"


def test_render_inputs_rejects_flat_layout(tmp_path: Path) -> None:
    pytest.importorskip("runops.adapters")
    run = _run_dir(tmp_path, "[sim]\nbatch_count = 10\n")
    with pytest.raises(ValueError, match="pre-3.0 flat layout"):
        Adapter().render_inputs({}, run)


def test_adapter_imports_without_numpy() -> None:
    """runops loads the adapter in its own interpreter, which may lack numpy."""
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[2]
    code = (
        "import sys; sys.modules['numpy'] = None; sys.modules['matplotlib'] = None\n"
        "import beach.runops_adapter, beach.config\n"
        "beach.runops_adapter.Adapter().case_template()\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_public_names_still_resolve() -> None:
    import beach

    assert beach.Beach.__name__ == "Beach"
    assert "load_fortran_result" in dir(beach)
    with pytest.raises(AttributeError):
        beach.no_such_name  # noqa: B018
