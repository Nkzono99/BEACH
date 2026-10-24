"""Shared fixtures and gate markers for the BEACH Python tests.

Tests that need a native artifact request it through a fixture, and the fixture name
selects the gate marker, so ``make`` targets choose suites with ``-m`` instead of file
lists: ``field_kernel_library`` marks ``native`` and ``native_checker`` marks
``config_native``. ``field_kernel_cache`` tests are marked explicitly and run only when
``BEACH_RUN_FIELD_KERNEL_CACHE_TESTS=1`` because they generate real periodic caches.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_MARKERS = {
    "field_kernel_library": "native",
    "native_checker": "config_native",
}


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    run_cache_tests = os.environ.get("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS") == "1"
    skip_cache = pytest.mark.skip(
        reason="periodic cache generation is opt-in; run make test-field-kernel-cache"
    )
    for item in items:
        fixtures = set(getattr(item, "fixturenames", ()))
        for fixture, marker in FIXTURE_MARKERS.items():
            if fixture in fixtures:
                item.add_marker(marker)
        if not run_cache_tests and item.get_closest_marker("field_kernel_cache"):
            item.add_marker(skip_cache)


@pytest.fixture(autouse=True)
def close_figures_after_test():
    yield
    pyplot = sys.modules.get("matplotlib.pyplot")
    if pyplot is not None:
        pyplot.close("all")


@pytest.fixture(scope="session")
def field_kernel_library() -> Path:
    """The native field kernel: ``BEACH_FIELD_KERNEL_LIB`` or ``build/``.

    ``make test-python-native`` builds the kernel and names it explicitly; a named
    library that does not exist is a failure. Without the variable, the tests skip
    when ``build/libbeach_field_kernel.so`` is absent. ``make test-python`` does not
    select ``native`` tests, so it never runs them against a stale ``build/`` copy.
    """

    supplied = os.environ.get("BEACH_FIELD_KERNEL_LIB")
    if supplied:
        path = Path(supplied).expanduser().resolve()
        assert path.is_file(), f"missing requested field kernel: {path}"
        return path
    path = REPO_ROOT / "build/libbeach_field_kernel.so"
    if not path.is_file():
        pytest.skip("run make test-python-native for native field-kernel contracts")
    return path


@pytest.fixture(scope="session")
def native_checker() -> Path:
    """The freshly built ``beach`` executable supplied by ``make test-config-contract``.

    Ordinary Python-only runs skip the cases that need it instead of selecting a
    stale binary.
    """

    executable = os.environ.get("BEACH_CONFIG_CHECK_EXE")
    if not executable:
        pytest.skip("run make test-config-contract for native config loading")
    path = Path(executable).resolve()
    assert path.is_file(), f"missing Fortran checker: {path}"
    return path
