"""Plot XY/YZ/XZ electric-potential slices in the simulation box."""

from __future__ import annotations

from pathlib import Path
import sys


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    from beach.cli.main import main as cli_main

    cli_main(["slices", *sys.argv[1:]])


if __name__ == "__main__":
    main()
