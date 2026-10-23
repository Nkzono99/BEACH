from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from beach.tooling.runops_adapter import Adapter, agent_context_path


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


def test_detect_status_needs_the_requested_batches(tmp_path: Path) -> None:
    adapter = Adapter()
    run = _run_dir(tmp_path)  # run.batch_count = 10
    summary = run / "work" / "latest" / "summary.txt"
    summary.write_text("batches=6\n", encoding="utf-8")  # e.g. stopped at the walltime
    assert adapter.detect_status(run) == "failed"
    summary.write_text("batches=10\n", encoding="utf-8")
    assert adapter.detect_status(run) == "completed"


def test_keep_patterns_keep_the_final_checkpoint() -> None:
    import fnmatch

    keep = Adapter().keep_patterns()

    def kept(rel: str) -> bool:
        name = rel.rsplit("/", 1)[-1]
        return any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(name, p) for p in keep)

    for name in (
        "summary.txt",
        "charges.csv",
        "rng_state.txt",
        "rng_state_rank00003.txt",
        "macro_residuals.csv",
        "charge_ledger.csv",
        "checkpoint_complete.txt",
    ):
        assert kept(f"latest/{name}"), name
    assert not kept("latest/checkpoints/slot0/rng_state.txt")
    assert not kept("latest/potential_history.csv")


def test_validate_params_flags_unknown_keys() -> None:
    adapter = Adapter()
    assert adapter.validate_params(
        {
            "run.batch_count": 200,
            "particles.species.0.source.pos_low": [0.0, 0.0, 1.0e-3],
            "particles.species[1].charge_c": -1.6e-19,
            "mesh.parts.0.path": "plate.obj",
            "output.history.stride_batches": 5,
        }
    ) == []
    warnings = adapter.validate_params({"run.batch_cuont": 200, "sim.batch_count": 3})
    assert len(warnings) == 2
    assert "'run.batch_cuont'" in warnings[0]


def _source_run(tmp_path: Path, batches: int = 10) -> Path:
    source = tmp_path / "source"
    (source / "work" / "latest").mkdir(parents=True)
    (source / "work" / "latest" / "summary.txt").write_text(
        f"batches={batches}\n", encoding="utf-8"
    )
    return source


def test_extend_inputs_resumes_from_the_source_checkpoint(tmp_path: Path) -> None:
    pytest.importorskip("runops.adapters")
    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python 3.10
        import tomli as tomllib

    source = _source_run(tmp_path)
    run = _run_dir(tmp_path)
    (run / "input" / "part.obj").write_text("v 0 0 0\n", encoding="utf-8")
    adapter = Adapter()
    adapter.render_inputs({"run.batch_count": 25}, run)
    written = adapter.extend_inputs(source, run, {"run.batch_count": 25})

    assert written == ["input/beach.toml", "input/restart_source.txt", "work/restart"]
    link = run / "work" / "restart"
    assert link.is_symlink() and link.resolve() == (source / "work" / "latest").resolve()
    config = tomllib.loads((run / "input" / "beach.toml").read_text(encoding="utf-8"))
    assert config["run"]["restart"]["from"] == "work/restart"
    assert config["output"]["enabled"] is True
    assert config["output"]["dir"] == "work/latest"
    assert "batches = 10" in (run / "input" / "restart_source.txt").read_text(encoding="utf-8")


def test_extend_inputs_needs_more_batches_and_a_checkpoint(tmp_path: Path) -> None:
    pytest.importorskip("runops.adapters")
    adapter = Adapter()
    run = _run_dir(tmp_path)  # run.batch_count = 10, same as the checkpoint
    with pytest.raises(ValueError, match="more than 10"):
        adapter.extend_inputs(_source_run(tmp_path), run, {})
    empty = tmp_path / "empty"
    (empty / "work").mkdir(parents=True)
    with pytest.raises(ValueError, match="no BEACH checkpoint"):
        adapter.extend_inputs(empty, run, {})


def test_extend_inputs_uses_a_periodic_slot(tmp_path: Path) -> None:
    pytest.importorskip("runops.adapters")
    source = tmp_path / "source"
    slot = source / "work" / "latest" / "checkpoints" / "slot1"
    slot.mkdir(parents=True)
    (slot / "summary.txt").write_text("batches=40\n", encoding="utf-8")
    run = _run_dir(tmp_path, GROUPED.replace("batch_count = 10", "batch_count = 50"))
    Adapter().extend_inputs(source, run, {})
    assert (run / "work" / "restart").resolve() == (source / "work" / "latest").resolve()


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
        "import beach.tooling.runops_adapter, beach.config\n"
        "beach.tooling.runops_adapter.Adapter().case_template()\n"
        "assert beach.tooling.runops_adapter.Adapter().validate_params({'run.batch_cuont': 1})\n"
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
    assert "load_run" in dir(beach)
    with pytest.raises(AttributeError):
        beach.no_such_name  # noqa: B018


def test_progress_counts_batches_from_the_resume_point(tmp_path: Path) -> None:
    adapter = Adapter()
    run = _run_dir(tmp_path)
    assert adapter.progress(run) is None
    (run / "manifest.toml").write_text('[job]\njob_id = "42"\n', encoding="utf-8")
    (run / "work" / "stdout.42.log").write_text(
        "runops: start 2026-10-23T02:43:53+09:00 on xa0146\n"
        "resuming_from_dir=work/restart/checkpoints/slot1\n"
        "resuming_from_batches=104000\n"
        "---------- batch 104001/300000 rel_change=  1.3E-02 ----------\n"
        "---------- batch 104002/300000 rel_change=  1.3E-02 ----------\n",
        encoding="utf-8",
    )
    assert adapter.progress(run) == {
        "done": 104002,
        "total": 300000,
        "unit": "batch",
        "start": 104000,
    }
    (run / "work" / "stdout.41.log").write_text("batch 9/10\n", encoding="utf-8")
    assert adapter.progress(run)["done"] == 104002  # the manifest's job, not the newest file


def test_agent_context_ships_docs() -> None:
    path = agent_context_path()
    assert path is not None and (path / "docs" / "usage-workflows.md").is_file()
