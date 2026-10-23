"""``beachx config diff``: semantic differences between two configs."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from beach.config import ConfigError, load_config_file, semantic_diff
from beach.config._toml import load_toml_file


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``config diff`` options."""

    parser.add_argument("left", type=Path, help="left-hand beach.toml")
    parser.add_argument("right", type=Path, help="right-hand beach.toml")
    parser.add_argument(
        "--raw",
        action="store_true",
        help="diff raw TOML payloads without high-level normalization",
    )


def run(args: argparse.Namespace) -> None:
    """Show semantic differences between two configs."""

    try:
        if args.raw:
            left_document: Any = load_toml_file(args.left)
            right_document: Any = load_toml_file(args.right)
        else:
            left_document = load_config_file(args.left)
            right_document = load_config_file(args.right)
    except FileNotFoundError as exc:
        raise SystemExit(f"config file not found: {exc.filename}") from exc
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc
    except ValueError as exc:
        raise SystemExit(f"TOML parse error: {exc}") from exc

    lines = semantic_diff(left_document, right_document)
    if not lines:
        print("no differences.")
        return
    for line in lines:
        print(line)
