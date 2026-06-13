"""Shared helpers for BEACH CLI entry points."""

from __future__ import annotations

import argparse
import math
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


def require_finite(
    parser: argparse.ArgumentParser,
    value: float | None,
    option: str,
) -> None:
    """Reject non-finite numeric CLI options."""

    if value is None:
        return
    if not math.isfinite(float(value)):
        parser.error(f"{option} must be finite.")


def require_nonnegative_finite(
    parser: argparse.ArgumentParser,
    value: float | None,
    option: str,
) -> None:
    """Reject negative or non-finite numeric CLI options."""

    require_finite(parser, value, option)
    if value is not None and float(value) < 0.0:
        parser.error(f"{option} must be >= 0.")


def require_positive_finite(
    parser: argparse.ArgumentParser,
    value: float | None,
    option: str,
) -> None:
    """Reject non-positive or non-finite numeric CLI options."""

    require_finite(parser, value, option)
    if value is not None and float(value) <= 0.0:
        parser.error(f"{option} must be > 0.")
