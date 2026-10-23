"""Shared user cache locations for native and Python entry points."""

import os
from pathlib import Path


def default_periodic_cache_dir() -> str:
    """Resolve the cache root at use time, following the native runtime precedence."""
    if root := os.environ.get("BEACH_CACHE_DIR"):
        return str(Path(root) / "periodic2")
    if root := os.environ.get("XDG_CACHE_HOME"):
        return str(Path(root) / "beach" / "periodic2")
    return str(Path.home() / ".cache" / "beach" / "periodic2")
