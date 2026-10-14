"""Shared fixtures for the BEACH Python tests."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


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
