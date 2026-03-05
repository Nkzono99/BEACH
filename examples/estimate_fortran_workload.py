"""Estimate per-batch particle workload and per-thread load from a Fortran TOML config."""

from __future__ import annotations

from pathlib import Path
import sys


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    from beach.cli_estimate_fortran_workload import main as cli_main

    cli_main()


if __name__ == "__main__":
    main()
