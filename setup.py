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
BIN_NAMES = ("beach", "beach-zhao-response")
BIN_BASES = tuple(INSTALL_PREFIX / "bin" / name for name in BIN_NAMES)
_BUILT_BINARIES: tuple[Path, ...] | None = None


def _which(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _built_binaries() -> tuple[Path, ...] | None:
    resolved: list[Path] = []
    for base in BIN_BASES:
        candidates = (base, base.with_suffix(".exe"))
        path = next((candidate for candidate in candidates if candidate.exists()), None)
        if path is None:
            return None
        resolved.append(path)
    return tuple(resolved)


def _build_with_make() -> tuple[Path, ...]:
    if not _which("make"):
        print("\nERROR: 'make' is required to build BEACH.\n", file=sys.stderr)
        sys.exit(1)

    # Default to auto profile for pip builds.
    install_profile = os.environ.get("INSTALL_PROFILE", "auto")
    cmd = ["make", f"INSTALL_PROFILE={install_profile}", f"PREFIX={INSTALL_PREFIX}", "install"]

    try:
        subprocess.check_call(cmd, cwd=ROOT_DIR)
    except subprocess.CalledProcessError as exc:
        # Optional fallback for constrained build environments.
        use_fallback = os.environ.get("BEACH_PIP_FALLBACK_GENERIC", "1") == "1"
        if use_fallback and install_profile == "auto" and "INSTALL_PROFILE" not in os.environ:
            retry_cmd = ["make", "INSTALL_PROFILE=generic", f"PREFIX={INSTALL_PREFIX}", "install"]
            print(
                "\nWARN: auto profile build failed; retrying with INSTALL_PROFILE=generic.\n",
                file=sys.stderr,
            )
            try:
                subprocess.check_call(retry_cmd, cwd=ROOT_DIR)
            except subprocess.CalledProcessError as retry_exc:
                print(
                    "\nERROR: failed to build/install Fortran executable via make.\n"
                    "       Ensure fpm and a Fortran compiler are available in PATH.\n",
                    file=sys.stderr,
                )
                raise SystemExit(retry_exc.returncode) from retry_exc
        else:
            print(
                "\nERROR: failed to build/install Fortran executable via make.\n"
                "       Ensure fpm and a Fortran compiler are available in PATH.\n",
                file=sys.stderr,
            )
            raise SystemExit(exc.returncode) from exc

    binpaths = _built_binaries()
    if not binpaths:
        print(
            "\nERROR: expected binaries not found: "
            + ", ".join(f"{base}(.exe)" for base in BIN_BASES)
            + "\n",
            file=sys.stderr,
        )
        sys.exit(1)

    for binpath in binpaths:
        try:
            mode = os.stat(binpath).st_mode
            os.chmod(binpath, mode | stat.S_IEXEC)
        except OSError:
            pass
    return binpaths


def _ensure_built_binaries() -> tuple[Path, ...]:
    global _BUILT_BINARIES
    if _BUILT_BINARIES and all(path.exists() for path in _BUILT_BINARIES):
        return _BUILT_BINARIES
    # Always build once per invocation to avoid reusing stale binaries from prior profile builds.
    _BUILT_BINARIES = _build_with_make()
    return _BUILT_BINARIES


class build_py(_build_py):
    def run(self) -> None:
        self.distribution.scripts = [str(path) for path in _ensure_built_binaries()]
        super().run()


class install(_install):
    def run(self) -> None:
        self.distribution.scripts = [str(path) for path in _ensure_built_binaries()]
        super().run()


class develop(_develop):
    def run(self) -> None:
        self.distribution.scripts = [str(path) for path in _ensure_built_binaries()]
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
