from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

import tools.check_source_text as source_text


def test_source_text_path_filter_includes_fortran_and_project_files() -> None:
    assert source_text.is_source_text_path(Path("src/runtime/bem_simulator.f90"))
    assert source_text.is_source_text_path(Path("src/core/bem_mpi.F90"))
    assert source_text.is_source_text_path(Path("Makefile"))
    assert source_text.is_source_text_path(Path("docs/new-source-format.astro"))
    assert not source_text.is_source_text_path(Path("docs/media/example.PNG"))
    assert not source_text.is_source_text_path(Path("docs/images/example.gif"))
    assert not source_text.is_source_text_path(Path("src/generated_backup.I90"))


def test_main_checks_git_source_set_and_returns_nul_status(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(
        ["git", "init", "--quiet"],
        cwd=repo,
        check=True,
        capture_output=True,
    )

    tracked = repo / "tracked.py"
    untracked = repo / "新規.f90"
    ignored = repo / "ignored.py"
    image = repo / "diagram.PNG"
    backup = repo / "generated.I90"
    (repo / ".gitignore").write_text("ignored.py\n", encoding="utf-8")
    tracked.write_bytes(b"print('tracked')\n\0")
    untracked.write_bytes("! 日本語\n".encode("utf-8") + b"\0")
    ignored.write_bytes(b"\0")
    image.write_bytes(b"\x89PNG\r\n\x1a\n\0")
    backup.write_bytes(b"\0")
    subprocess.run(
        ["git", "add", ".gitignore", tracked.name],
        cwd=repo,
        check=True,
        capture_output=True,
    )
    monkeypatch.setattr(
        source_text,
        "__file__",
        str(repo / "tools" / "check_source_text.py"),
    )

    assert source_text.main() == 1
    failure_output = capsys.readouterr().out
    assert tracked.name in failure_output
    assert untracked.name in failure_output
    assert ignored.name not in failure_output
    assert image.name not in failure_output
    assert backup.name not in failure_output

    tracked.write_text("print('clean')\n", encoding="utf-8")
    untracked.write_text("! 日本語\n", encoding="utf-8")

    assert source_text.main() == 0
