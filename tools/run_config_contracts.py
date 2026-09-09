"""fpm runner: test the configuration contract with this exact BEACH binary."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_config_contracts.py BEACH_EXECUTABLE")
    environment = os.environ.copy()
    environment["BEACH_CONFIG_CHECK_EXE"] = str(Path(sys.argv[1]).resolve())
    return subprocess.call(
        [sys.executable, "-m", "pytest", "-q", "tests/python/test_config_contract.py"],
        env=environment,
    )


if __name__ == "__main__":
    raise SystemExit(main())
