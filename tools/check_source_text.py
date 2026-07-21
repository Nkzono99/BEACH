#!/usr/bin/env python3
"""Reject NUL bytes in tracked or untracked BEACH source-text files."""

from __future__ import annotations

import os
import subprocess
from collections.abc import Iterable
from pathlib import Path


BINARY_SUFFIXES = frozenset({".gif", ".png"})
IGNORED_SUFFIXES = frozenset({".i90"})


def is_source_text_path(path: Path) -> bool:
    """Return whether ``path`` is expected to contain source text."""

    suffix = path.suffix.lower()
    return suffix not in BINARY_SUFFIXES and suffix not in IGNORED_SUFFIXES


def repository_source_text_paths(root: Path) -> list[Path]:
    """List tracked and non-ignored untracked source-text paths under ``root``."""

    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
    )
    paths: list[Path] = []
    for raw_path in result.stdout.split(b"\0"):
        if not raw_path:
            continue
        relative_path = Path(os.fsdecode(raw_path))
        path = root / relative_path
        if is_source_text_path(relative_path) and path.is_file():
            paths.append(path)
    return paths


def find_nul_files(paths: Iterable[Path]) -> list[Path]:
    """Return every path whose contents contain at least one NUL byte."""

    return [path for path in paths if b"\0" in path.read_bytes()]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    nul_files = find_nul_files(repository_source_text_paths(root))
    if not nul_files:
        print("source text check passed: no NUL bytes found")
        return 0

    print("source text check failed: NUL bytes found in:")
    for path in nul_files:
        print(f"- {path.relative_to(root)}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
