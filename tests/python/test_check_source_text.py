from __future__ import annotations

from pathlib import Path

from tools.check_source_text import find_nul_files, is_source_text_path


def test_source_text_path_filter_includes_fortran_and_project_files() -> None:
    assert is_source_text_path(Path("src/runtime/bem_simulator.f90"))
    assert is_source_text_path(Path("src/core/bem_mpi.F90"))
    assert is_source_text_path(Path("Makefile"))
    assert is_source_text_path(Path("docs/new-source-format.astro"))
    assert not is_source_text_path(Path("docs/media/example.png"))
    assert not is_source_text_path(Path("docs/images/example.gif"))
    assert not is_source_text_path(Path("src/generated_backup.i90"))


def test_find_nul_files_reports_only_contaminated_text(tmp_path: Path) -> None:
    clean = tmp_path / "clean.py"
    contaminated = tmp_path / "contaminated.f90"
    clean.write_bytes(b"print('ok')\n")
    contaminated.write_bytes(b"end module\n\0\0")

    assert find_nul_files([clean, contaminated]) == [contaminated]
