"""Exercise the build wrapper's cache lifecycle with isolated build artifacts."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]


def write_executable(path: Path, body: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\nset -eu\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return path


@pytest.fixture
def project(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    shutil.copy2(ROOT / "build.sh", tmp_path)
    (tmp_path / "fpm.toml").write_text('name = "beach_fortran"\nversion = "1.0.0"\n')
    (tmp_path / "src").mkdir()
    (tmp_path / "src/model.f90").write_text("module model\nend module\n")
    compiler = write_executable(tmp_path / "bin/gfortran", 'echo "compiler ${TEST_FC_VERSION:-1}"\n')
    fpm = write_executable(
        tmp_path / "bin/fpm",
        'printf "%s\\n" "$@" > "$TEST_FPM_LOG"\nexit "${TEST_FPM_STATUS:-0}"\n',
    )
    env = {key: value for key, value in os.environ.items() if not key.startswith(("FPM_", "BEACH_"))}
    env.update(
        FPM=str(fpm), FPM_FC=str(compiler), BEACH_VERSION_MODE="dev",
        TEST_FPM_LOG=str(tmp_path / "fpm.log"),
    )
    return tmp_path, env


def run_build(project: tuple[Path, dict[str, str]], *args: str, **overrides: str) -> subprocess.CompletedProcess[str]:
    root, env = project
    result = subprocess.run(
        [str(root / "build.sh"), *args], cwd=root, env=env | overrides,
        text=True, capture_output=True, timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def cache_artifacts(root: Path, compiler: str = "gfortran") -> list[Path]:
    paths = [
        root / f"build/{compiler}_0123456789ABCDEF/beach_fortran/old.f90.o",
        root / f"build/{compiler}_0123456789ABCDEF/model.mod",
        root / f"build/{compiler}_123456789ABCDEF0/beach_fortran/libbeach_fortran.a",
        root / f"build/{compiler}_23456789ABCDEF01/test/model_test",
        root / f"build/{compiler}_23456789ABCDEF01/app/model",
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    return paths


def test_bootstrap_clears_complete_compiler_cache_and_preserves_products(project: tuple) -> None:
    root, _ = project
    old = cache_artifacts(root)
    kept = cache_artifacts(root, "mpiifx") + [
        root / "build/dependencies/source.f90",
        root / "build/logs/test.log",
        root / "build/docs/index.html",
        root / "build/libbeach_field_kernel.so",
        root / "build/gfortran_notes_0123456789ABCDEF/notes.txt",
    ]
    for path in kept:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    run_build(project)
    assert all(not path.exists() for path in old)
    assert all(path.exists() for path in kept)


@pytest.mark.parametrize("change", ["add", "move", "remove"])
def test_source_layout_changes_invalidate_once(project: tuple, change: str) -> None:
    root, _ = project
    run_build(project)
    old = cache_artifacts(root)
    source = root / "src/model.f90"
    if change == "add":
        (root / "src/new.f90").write_text("! untracked source\n")
    elif change == "move":
        source.rename(root / "src/moved.f90")
    else:
        source.unlink()
    run_build(project)
    assert all(not path.exists() for path in old)
    fresh = cache_artifacts(root)
    run_build(project, "--target", "second", FPM_ACTION="test")
    assert all(path.exists() for path in fresh)


def test_content_flags_profile_and_backup_changes_reuse_cache(project: tuple) -> None:
    root, _ = project
    run_build(project)
    old = cache_artifacts(root)
    (root / "src/model.f90").write_text("module model\ninteger :: value\nend module\n")
    (root / "src/model.i90").write_text("! formatter backup\n")
    run_build(project, "--target", "second", FPM_ACTION="test", FPM_PROFILE="debug", FPM_FFLAGS="-O0")
    assert all(path.exists() for path in old)


@pytest.mark.parametrize("change", ["version", "path"])
def test_compiler_replacement_invalidates_same_name_cache(project: tuple, change: str) -> None:
    root, _ = project
    run_build(project)
    old = cache_artifacts(root)
    if change == "version":
        run_build(project, TEST_FC_VERSION="2")
    else:
        replacement = write_executable(root / "new-bin/gfortran", 'echo "compiler 1"\n')
        run_build(project, FPM_FC=str(replacement))
    assert all(not path.exists() for path in old)


def test_cli_compiler_precedence_and_runtime_arguments(project: tuple) -> None:
    root, env = project
    other = write_executable(root / "bin/mpiifx", 'echo "ifx 1"\n')
    first = cache_artifacts(root)
    second = cache_artifacts(root, "mpiifx")
    run_build(project, "--compiler", str(other), "--", "--compiler", env["FPM_FC"], FPM_ACTION="run")
    assert all(path.exists() for path in first)
    assert all(not path.exists() for path in second)
    forwarded = (root / "fpm.log").read_text().splitlines()
    assert forwarded[-5:] == ["--compiler", str(other), "--", "--compiler", env["FPM_FC"]]


@pytest.mark.parametrize("option", ["--help", "--version", "--list", "--show-model"])
def test_inspection_never_resets_cache(project: tuple, option: str) -> None:
    root, _ = project
    old = cache_artifacts(root)
    run_build(project, option, BEACH_REBUILD="1")
    assert all(path.exists() for path in old)


def test_explicit_rebuild_invalidates_unchanged_cache(project: tuple) -> None:
    root, _ = project
    run_build(project)
    old = cache_artifacts(root)
    run_build(project, BEACH_REBUILD="1")
    assert all(not path.exists() for path in old)


@pytest.mark.parametrize("stream", ["stdout", "stderr"])
def test_kernel_links_fpm_selected_archive_even_with_newer_foreign_build(tmp_path: Path, stream: str) -> None:
    shutil.copy2(ROOT / "Makefile", tmp_path)
    selected = tmp_path / "build/current/beach_fortran/libbeach_fortran.a"
    selected.parent.mkdir(parents=True)
    selected.touch()
    foreign = tmp_path / "build/foreign/beach_fortran/libbeach_fortran.a"
    foreign.parent.mkdir(parents=True)
    foreign.touch()
    os.utime(selected, (1, 1))
    redirect = " >&2" if stream == "stderr" else ""
    fake_build = write_executable(
        tmp_path / "build.sh", 'echo " build/current/beach_fortran/libbeach_fortran.a"' + redirect + "\n",
    )
    fake_compiler = write_executable(tmp_path / "compiler", 'printf "%s\\n" "$@" > kernel-link.log\n')
    result = subprocess.run(
        ["make", "build-kernel", f"BUILD_SH={fake_build}", f"KERNEL_FC={fake_compiler}"],
        cwd=tmp_path, text=True, capture_output=True, timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    args = (tmp_path / "kernel-link.log").read_text().splitlines()
    assert "build/current/beach_fortran/libbeach_fortran.a" in args
    assert all("foreign" not in arg for arg in args)
