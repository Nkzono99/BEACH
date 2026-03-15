"""Shared helpers for BEACH CLI entry points."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable


def configure_entry_parser(
    parser: argparse.ArgumentParser,
    handler: Callable[[argparse.Namespace], None],
) -> argparse.ArgumentParser:
    """Attach the parsed-command handler to one parser instance."""

    parser.set_defaults(func=handler, _parser=parser)
    return parser


def print_legacy_warning(legacy_name: str, subcommand: str) -> None:
    """Print a one-line deprecation warning for a legacy CLI alias."""

    print(
        f"WARNING: `{legacy_name}` is deprecated; use `beachx {subcommand}` instead.",
        file=sys.stderr,
    )
