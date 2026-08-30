from __future__ import annotations

import importlib
import re
import tarfile
from pathlib import Path

import pytest

from tools.verify_sdist import SdistContractError, verify_sdist


ROOT = Path(__file__).resolve().parents[2]
EXPLICIT_TARGETS = (
    ("executable", "beach", "app/main.f90"),
    ("test", "contract", "tests/fortran/test_contract.f90"),
    ("example", "benchmark", "benchmarks/fortran/benchmark.f90"),
)


def _load_toml(relative_path: str) -> dict[str, object]:
    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python 3.10 compatibility
        import tomli as tomllib

    with (ROOT / relative_path).open("rb") as stream:
        return tomllib.load(stream)


def _write_fixture_sdist(
    tmp_path: Path,
    *,
    missing_source: str | None = None,
) -> Path:
    package_root = tmp_path / "beach_bem-1.6.2"
    for _, name, source in EXPLICIT_TARGETS:
        if source == missing_source:
            continue
        path = package_root / source
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"program {name}\nend program\n", encoding="utf-8")
    (package_root / "fpm.toml").write_text(
        """
[[executable]]
name = "beach"
source-dir = "app"
main = "main.f90"

[[test]]
name = "contract"
source-dir = "tests/fortran"
main = "test_contract.f90"

[[example]]
name = "benchmark"
source-dir = "benchmarks/fortran"
main = "benchmark.f90"
""".lstrip(),
        encoding="utf-8",
    )

    archive_path = tmp_path / "beach_bem-1.6.2.tar.gz"
    with tarfile.open(archive_path, mode="w:gz") as archive:
        archive.add(package_root, arcname=package_root.name)
    return archive_path


def test_verify_sdist_accepts_complete_explicit_targets(tmp_path: Path) -> None:
    verify_sdist(_write_fixture_sdist(tmp_path))


@pytest.mark.parametrize(
    ("section", "name", "source"),
    EXPLICIT_TARGETS,
    ids=("executable", "test", "example"),
)
def test_verify_sdist_reports_each_missing_explicit_target(
    tmp_path: Path,
    section: str,
    name: str,
    source: str,
) -> None:
    archive = _write_fixture_sdist(tmp_path, missing_source=source)

    expected = rf"{section} {name}: {re.escape(source)}"
    with pytest.raises(SdistContractError, match=expected):
        verify_sdist(archive)


def test_release_metadata_matches_fortran_and_isolated_build_contract() -> None:
    pyproject = _load_toml("pyproject.toml")
    fpm = _load_toml("fpm.toml")
    project = pyproject["project"]
    build_system = pyproject["build-system"]

    assert project["name"] == "beach-bem"
    assert project["version"] == fpm["version"]
    assert project["license"] == "Apache-2.0"
    assert "LICENSE" in project["license-files"]
    assert build_system["build-backend"] == "setuptools.build_meta"
    assert {"setuptools>=77", "wheel", "fpm"} <= set(build_system["requires"])
    assert {"beach", "beach-zhao-response"} <= {
        target["name"] for target in fpm["executable"]
    }


def test_wheel_metadata_exposes_runtime_resources_and_cli_entry_points() -> None:
    pyproject = _load_toml("pyproject.toml")
    project = pyproject["project"]
    setuptools = pyproject["tool"]["setuptools"]
    package_data = setuptools["package-data"]

    assert {
        "beach",
        "beach.cli",
        "beach.config",
        "beach.fortran_results",
    } <= set(setuptools["packages"])
    assert "include/*.h" in package_data["beach"]
    assert "schemas/*.json" in package_data["beach.config"]
    assert (ROOT / "beach/include/beach_field_kernel.h").is_file()
    assert (ROOT / "beach/config/schemas/beach.schema.json").is_file()

    scripts = project["scripts"]
    assert {
        "beachx",
        "beach-inspect",
        "beach-animate-history",
        "beach-plot-coulomb-force-matrix",
        "beach-plot-potential-slices",
        "beach-estimate-workload",
        "beach-plot-performance-profile",
    } <= set(scripts)
    for command, target in scripts.items():
        module_name, separator, attribute_name = target.partition(":")
        assert separator and attribute_name, command
        entry_point = importlib.import_module(module_name)
        for attribute in attribute_name.split("."):
            entry_point = getattr(entry_point, attribute)
        assert callable(entry_point), command
