from __future__ import annotations

import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install

try:
    from setuptools.command.build_scripts import build_scripts as _build_scripts
except Exception:
    from distutils.command.build_scripts import build_scripts as _build_scripts


ROOT_DIR = Path(__file__).resolve().parent
INSTALL_PREFIX = Path(os.environ.get("BEACH_PIP_PREFIX", ROOT_DIR / "build" / "pip-install"))
BIN_BASE = INSTALL_PREFIX / "bin" / "beach"
CANDIDATES = (BIN_BASE, BIN_BASE.with_suffix(".exe"))
_BUILT_BINARY: Path | None = None


def _which(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _built_binary() -> Path | None:
    for path in CANDIDATES:
        if path.exists():
            return path
    return None


def _build_with_make() -> Path:
    if not _which("make"):
        print("\nERROR: 'make' is required to build BEACH.\n", file=sys.stderr)
        sys.exit(1)

    cmd = ["make", f"PREFIX={INSTALL_PREFIX}", "install"]
    install_profile = os.environ.get("INSTALL_PROFILE", "")
    if install_profile:
        cmd.insert(1, f"INSTALL_PROFILE={install_profile}")

    try:
        subprocess.check_call(cmd, cwd=ROOT_DIR)
    except subprocess.CalledProcessError as exc:
        print(
            "\nERROR: failed to build/install Fortran executable via make.\n"
            "       Ensure fpm and a Fortran compiler are available in PATH.\n",
            file=sys.stderr,
        )
        raise SystemExit(exc.returncode) from exc

    binpath = _built_binary()
    if not binpath:
        print(
            f"\nERROR: expected binary not found: {BIN_BASE}(.exe)\n",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        mode = os.stat(binpath).st_mode
        os.chmod(binpath, mode | stat.S_IEXEC)
    except OSError:
        pass
    return binpath


def _ensure_built_binary() -> Path:
    global _BUILT_BINARY
    if _BUILT_BINARY and _BUILT_BINARY.exists():
        return _BUILT_BINARY
    _BUILT_BINARY = _built_binary() or _build_with_make()
    return _BUILT_BINARY


class build_py(_build_py):
    def run(self) -> None:
        self.distribution.scripts = [str(_ensure_built_binary())]
        super().run()


class install(_install):
    def run(self) -> None:
        self.distribution.scripts = [str(_ensure_built_binary())]
        super().run()


class develop(_develop):
    def run(self) -> None:
        self.distribution.scripts = [str(_ensure_built_binary())]
        super().run()


class build_scripts(_build_scripts):
    def copy_scripts(self) -> None:
        if not hasattr(self, "outfiles") or self.outfiles is None:
            self.outfiles = []

        self.mkpath(self.build_dir)
        for script in self.scripts:
            src = Path(script)
            dst = Path(self.build_dir) / src.name
            with open(src, "rb") as fsrc, open(dst, "wb") as fdst:
                fdst.write(fsrc.read())
            try:
                mode = os.stat(dst).st_mode
                os.chmod(dst, mode | stat.S_IEXEC)
            except OSError:
                pass
            self.outfiles.append(str(dst))


setup(
    cmdclass={
        "build_py": build_py,
        "install": install,
        "develop": develop,
        "build_scripts": build_scripts,
    }
)
