"""Import-layer contract of the ``beach`` subpackages.

``beach.io`` is the bottom layer, ``beach.config`` and ``beach.tooling`` must run
without numpy (runops loads the adapter in its own interpreter), and the plotting
layer reaches the native kernel only through ``beach.analysis``.
"""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]

_PROBE = """
import importlib, pkgutil, sys
package = importlib.import_module({package!r})
if {walk!r}:
    for info in pkgutil.walk_packages(package.__path__, package.__name__ + "."):
        importlib.import_module(info.name)
print(" ".join(name for name in {forbidden!r} if name in sys.modules))
"""


@pytest.mark.parametrize(
    ("package", "walk", "forbidden"),
    [
        ("beach.io", True, ("matplotlib", "beach.native", "beach.analysis", "beach.plot")),
        ("beach.config", True, ("numpy", "matplotlib")),
        ("beach.tooling", True, ("numpy", "matplotlib")),
        # Kernel-backed figures reach beach.native through beach.analysis, so only
        # the package itself is probed; the static test below covers its modules.
        ("beach.plot", False, ("matplotlib", "beach.native")),
    ],
)
def test_layer_imports_without_heavier_layers(
    package: str, walk: bool, forbidden: tuple[str, ...]
) -> None:
    result = subprocess.run(
        [sys.executable, "-c", _PROBE.format(package=package, walk=walk, forbidden=forbidden)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == []


def _imported_names(path: Path, package: str) -> list[str]:
    names: list[str] = []
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            names += [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                parent = package.split(".")[: len(package.split(".")) - node.level + 1]
                base = ".".join([*parent, base] if base else parent)
            names += [base, *(f"{base}.{alias.name}" for alias in node.names)]
    return names


def test_plot_modules_do_not_import_the_native_kernel() -> None:
    offenders = [
        f"{path.name}: {name}"
        for path in sorted((ROOT / "beach" / "plot").glob("*.py"))
        for name in _imported_names(path, "beach.plot")
        if name == "beach.native" or name.startswith("beach.native.")
    ]

    assert offenders == []
