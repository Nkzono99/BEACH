"""``beachx config init``: write a runnable ``beach.toml``."""

from __future__ import annotations

import argparse
from pathlib import Path

from beach.config import CONFIG_FILENAME, default_config, dump_beach_toml


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``config init`` options."""

    parser.add_argument(
        "output",
        nargs="?",
        default=Path(CONFIG_FILENAME),
        type=Path,
        help=f"destination config path (default: ./{CONFIG_FILENAME})",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite the destination if it already exists",
    )


def run(args: argparse.Namespace) -> None:
    """Create one new ``beach.toml`` file."""

    destination = args.output
    if destination.exists() and not args.force:
        raise SystemExit(f"config file already exists: {destination}")

    text = dump_beach_toml(default_config())
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    print(f"saved={destination}")
