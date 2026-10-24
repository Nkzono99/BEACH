"""Locate the native kernel used by cross-language contract tests."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


def kernel_library_path() -> Path:
    supplied = os.environ.get("BEACH_FIELD_KERNEL_LIB")
    path = (
        Path(supplied).expanduser().resolve()
        if supplied
        else Path(__file__).resolve().parents[2] / "build/libbeach_field_kernel.so"
    )
    if supplied:
        assert path.is_file(), f"missing requested field kernel: {path}"
    elif not path.is_file():
        pytest.skip("run make test-python-native for native field-kernel contracts")
    return path
