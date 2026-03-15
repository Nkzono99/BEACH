"""Plot a BEACH performance profile from ``performance_profile.csv``."""

from __future__ import annotations

from pathlib import Path
import sys


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    from beach.cli_plot_performance_profile import main as cli_main

    cli_main()


if __name__ == "__main__":
    main()
