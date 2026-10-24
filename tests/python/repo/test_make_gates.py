"""Makefile gate routing: which pytest markers and Fortran targets each target runs."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def _line_starting_with(text: str, prefix: str) -> str:
    for line in text.splitlines():
        if line.startswith(prefix):
            return line
    raise AssertionError(f"missing line starting with {prefix!r}")


def _make_database(environment: dict[str, str]) -> str:
    result = subprocess.run(
        ["make", "-pn", "--no-print-directory", "-f", "Makefile"],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


def _pytest_marker_expression(environment: dict[str, str], target: str) -> str:
    result = subprocess.run(
        ["make", "-n", "--no-print-directory", "-f", "Makefile", target],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    commands = [line for line in result.stdout.splitlines() if " -m pytest " in line]
    assert len(commands) == 1, (target, result.stdout)
    match = re.search(r" -m '([^']+)'", commands[0])
    assert match is not None, commands[0]
    return match.group(1)


def _selects(expression: str, markers: set[str]) -> bool:
    names = set(re.findall(r"[A-Za-z_]\w*", expression)) - {"and", "or", "not"}
    return bool(eval(expression, {"__builtins__": {}}, {name: name in markers for name in names}))


def test_field_kernel_cache_tests_are_opt_in() -> None:
    environment = os.environ.copy()
    for name in ("MAKEFLAGS", "MAKEOVERRIDES", "MFLAGS"):
        environment.pop(name, None)
    for name in list(environment):
        if name.startswith("FORTRAN_") and name.endswith("_TARGETS"):
            environment.pop(name)
    database = _make_database(environment)

    far_targets = set(
        _line_starting_with(database, "FORTRAN_FAR_CORRECTION_TARGETS = ")
        .split(" = ", maxsplit=1)[1].split()
    )
    assert "test_field_kernel_cache_c" in far_targets
    tier_lines = re.findall(r"^FORTRAN_L\d+_TARGETS = (.*)$", database, re.MULTILINE)
    assert tier_lines
    for line in tier_lines:
        assert "test_field_kernel_cache_c" not in line.split()

    opt_in_assignment = _line_starting_with(
        database,
        "test-field-kernel-cache: BEACH_RUN_FIELD_KERNEL_CACHE_TESTS ",
    )
    assert opt_in_assignment.rsplit("=", maxsplit=1)[1].strip() == "1"

    tiers = sorted(set(re.findall(r"^(test-l\d+):", database, re.MULTILINE)))
    default_routes = subprocess.run(
        ["make", "-n", "--trace", "--no-print-directory", "-f", "Makefile",
         *tiers, "test-physics-release"],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert default_routes.returncode == 0, default_routes.stderr
    assert "target 'test-field-kernel-cache'" not in default_routes.stdout

    cache_test = {"field_kernel_cache", "native"}
    for target in ("test-python", "test-python-native", "test-validation-workflow"):
        expression = _pytest_marker_expression(environment, target)
        assert not _selects(expression, cache_test), (target, expression)
        assert not _selects(expression, cache_test | {"workflow"}), (target, expression)
    assert _selects(
        _pytest_marker_expression(environment, "test-field-kernel-cache"), cache_test
    )
