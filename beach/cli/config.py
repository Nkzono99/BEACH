"""Direct config management commands for ``beachx``."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Sequence

from beach.config import (
    CONFIG_FILENAME,
    ConfigError,
    default_config,
    load_config_file,
    dump_beach_toml,
    semantic_diff,
)
from beach.config._toml import load_toml_file

from ._shared import configure_entry_parser

COMMAND_NAME = "config"


def build_parser(*, prog: str | None = None) -> argparse.ArgumentParser:
    """Build the parser for the ``beachx config`` group."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_group_parser(parser)
    return parser


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register the ``config`` command group under the root CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="create and validate beach.toml",
        description="Create, validate, and diff direct BEACH config files.",
    )
    _configure_group_parser(parser)
    return parser


def _configure_group_parser(parser: argparse.ArgumentParser) -> None:
    config_subparsers = parser.add_subparsers(
        dest="config_command",
        metavar="subcommand",
        required=True,
    )
    _add_init_subparser(config_subparsers)
    _add_validate_subparser(config_subparsers)
    _add_diff_subparser(config_subparsers)


def _add_init_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "init",
        help="create a runnable beach.toml",
    )
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
    configure_entry_parser(parser, run_init)


def _add_validate_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "validate",
        help="validate a beach.toml config",
    )
    parser.add_argument(
        "config_path",
        nargs="?",
        default=Path(CONFIG_FILENAME),
        type=Path,
        help=f"input config file (default: ./{CONFIG_FILENAME})",
    )
    configure_entry_parser(parser, run_validate)


def _add_diff_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "diff",
        help="show semantic differences between two configs",
    )
    parser.add_argument("left", type=Path, help="left-hand beach.toml")
    parser.add_argument("right", type=Path, help="right-hand beach.toml")
    parser.add_argument(
        "--raw",
        action="store_true",
        help="diff raw TOML payloads without high-level normalization",
    )
    configure_entry_parser(parser, run_diff)


def run_init(args: argparse.Namespace) -> None:
    """Create one new ``beach.toml`` file."""

    destination = args.output
    if destination.exists() and not args.force:
        raise SystemExit(f"config file already exists: {destination}")

    text = dump_beach_toml(default_config())
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    print(f"saved={destination}")


def run_validate(args: argparse.Namespace) -> None:
    """Validate one direct config file."""

    try:
        load_config_file(args.config_path)
    except FileNotFoundError as exc:
        raise SystemExit(f"config file not found: {args.config_path}") from exc
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc

    print(f"config={args.config_path}")
    print("status=ok")


def run_diff(args: argparse.Namespace) -> None:
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

    lines = semantic_diff(left_document, right_document)
    if not lines:
        print("no differences.")
        return
    for line in lines:
        print(line)


def main(argv: Sequence[str] | None = None) -> None:
    """Run the ``beachx config`` group as a standalone entry."""

    args = build_parser(prog="beachx config").parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
