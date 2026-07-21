from __future__ import annotations

import tarfile
from pathlib import Path

import pytest

from tools.verify_sdist import SdistContractError, verify_sdist


ROOT = Path(__file__).resolve().parents[2]


def _write_fixture_sdist(tmp_path: Path, *, include_benchmark: bool) -> Path:
    package_root = tmp_path / "beach_bem-1.6.2"
    (package_root / "app").mkdir(parents=True)
    (package_root / "benchmarks" / "fortran").mkdir(parents=True)
    (package_root / "app" / "main.f90").write_text("program main\nend program\n")
    if include_benchmark:
        (package_root / "benchmarks" / "fortran" / "benchmark.f90").write_text(
            "program benchmark\nend program\n"
        )
    (package_root / "fpm.toml").write_text(
        """
[[executable]]
name = "beach"
source-dir = "app"
main = "main.f90"

[[example]]
name = "benchmark"
source-dir = "benchmarks/fortran"
main = "benchmark.f90"
""".lstrip()
    )

    archive_path = tmp_path / "beach_bem-1.6.2.tar.gz"
    with tarfile.open(archive_path, mode="w:gz") as archive:
        archive.add(package_root, arcname=package_root.name)
    return archive_path


def test_verify_sdist_accepts_complete_explicit_targets(tmp_path: Path) -> None:
    verify_sdist(_write_fixture_sdist(tmp_path, include_benchmark=True))


def test_verify_sdist_reports_missing_explicit_target(tmp_path: Path) -> None:
    archive = _write_fixture_sdist(tmp_path, include_benchmark=False)

    with pytest.raises(SdistContractError, match="example benchmark"):
        verify_sdist(archive)


def test_release_metadata_and_sdist_gates_are_synchronized() -> None:
    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python 3.10 compatibility
        import tomli as tomllib

    with (ROOT / "pyproject.toml").open("rb") as stream:
        pyproject = tomllib.load(stream)
        python_version = pyproject["project"]["version"]
    with (ROOT / "fpm.toml").open("rb") as stream:
        fortran_version = tomllib.load(stream)["version"]

    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    publish_workflow = (
        ROOT / ".github" / "workflows" / "publish-pypi.yml"
    ).read_text(encoding="utf-8")

    assert python_version == fortran_version
    assert pyproject["project"]["license"] == "Apache-2.0"
    assert pyproject["project"]["license-files"] == ["LICENSE"]
    assert "setuptools>=77" in pyproject["build-system"]["requires"]
    assert "recursive-include benchmarks *.f90" in manifest
    assert "recursive-include beach/include *.h" in manifest
    assert "include tools/check_source_text.py" in manifest
    assert pyproject["tool"]["setuptools"]["package-data"]["beach"] == [
        "include/*.h"
    ]
    assert "include tools/verify_sdist.py" in manifest
    assert "python tools/verify_sdist.py dist/*.tar.gz" in publish_workflow
    assert "python -m pip wheel --no-deps dist/*.tar.gz" in publish_workflow
