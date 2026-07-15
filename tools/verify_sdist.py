#!/usr/bin/env python3
"""Verify that a BEACH source distribution contains every explicit fpm target."""

from __future__ import annotations

import argparse
import tarfile
from pathlib import Path, PurePosixPath
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 compatibility
    import tomli as tomllib


TARGET_SECTIONS = ("executable", "test", "example")


class SdistContractError(RuntimeError):
    """Raised when an sdist cannot satisfy its fpm manifest."""


def _read_fpm_manifest(
    archive: tarfile.TarFile,
) -> tuple[PurePosixPath, dict[str, Any], set[PurePosixPath]]:
    files = {
        PurePosixPath(member.name): member
        for member in archive.getmembers()
        if member.isfile()
    }
    manifests = [
        path for path in files if path.name == "fpm.toml" and len(path.parts) == 2
    ]
    if len(manifests) != 1:
        raise SdistContractError(
            f"expected one top-level fpm.toml in sdist, found {len(manifests)}"
        )

    manifest_path = manifests[0]
    stream = archive.extractfile(files[manifest_path])
    if stream is None:
        raise SdistContractError(f"could not read {manifest_path}")
    manifest = tomllib.loads(stream.read().decode("utf-8"))
    return manifest_path.parent, manifest, set(files)


def verify_sdist(path: Path) -> None:
    """Raise ``SdistContractError`` if an explicit fpm target is missing."""

    with tarfile.open(path, mode="r:*") as archive:
        root, manifest, files = _read_fpm_manifest(archive)

    missing: list[str] = []
    for section in TARGET_SECTIONS:
        for target in manifest.get(section, []):
            main = target.get("main")
            if not main:
                continue
            source_dir = target.get("source-dir", ".")
            expected = root / source_dir / main
            if expected not in files:
                name = target.get("name", "<unnamed>")
                missing.append(f"{section} {name}: {expected.relative_to(root)}")

    if missing:
        details = "\n".join(f"- {item}" for item in missing)
        raise SdistContractError(f"sdist is missing fpm target sources:\n{details}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", nargs="+", type=Path)
    args = parser.parse_args()

    for archive in args.archives:
        verify_sdist(archive)
        print(f"sdist contract ok: {archive}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
