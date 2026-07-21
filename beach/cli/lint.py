"""Config lint command for ``beachx``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from beach.config import (
    CONFIG_FILENAME,
    ConfigError,
    normalize_config_document,
    normalize_high_level_config,
)
from beach.config._toml import load_toml_file
from beach.config.schema import load_schema as _load_schema
from beach.config.schema import schema_errors as _schema_errors

from ._shared import configure_entry_parser

COMMAND_NAME = "lint"


def build_parser(*, prog: str | None = None) -> argparse.ArgumentParser:
    """Build the parser for ``beachx lint``."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return parser


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register the ``lint`` command under the root CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="lint beach.toml with schema and BEACH semantic checks",
        description="Validate a BEACH config with TOML parsing, JSON Schema, and semantic checks.",
    )
    _configure_parser(parser)
    return parser


def _configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "config_path",
        nargs="?",
        default=Path(CONFIG_FILENAME),
        type=Path,
        help=f"input config file (default: ./{CONFIG_FILENAME})",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        help="JSON Schema path (default: packaged beach.schema.json)",
    )
    parser.add_argument(
        "--max-errors",
        type=int,
        default=20,
        help="maximum schema errors to print (default: 20)",
    )
    configure_entry_parser(parser, run_lint)


def run_lint(args: argparse.Namespace) -> None:
    """Lint one BEACH config file."""

    if args.max_errors < 1:
        raise SystemExit("--max-errors must be >= 1")

    try:
        raw_config = load_toml_file(args.config_path)
    except FileNotFoundError as exc:
        raise SystemExit(f"config file not found: {exc.filename}") from exc
    except ValueError as exc:
        raise SystemExit(f"TOML parse error: {exc}") from exc

    try:
        schema, schema_label = _load_schema(args.schema)
    except FileNotFoundError as exc:
        raise SystemExit(f"schema file not found: {exc.filename}") from exc
    except json.JSONDecodeError as exc:
        raise SystemExit(f"schema JSON parse error: {exc}") from exc
    except ValueError as exc:
        raise SystemExit(f"schema error: {exc}") from exc

    raw_schema_errors = _schema_errors(raw_config, schema)
    if raw_schema_errors:
        _raise_schema_errors(
            path=args.config_path,
            phase="authoring",
            errors=raw_schema_errors,
            max_errors=args.max_errors,
        )

    try:
        normalized_config = normalize_high_level_config(raw_config)
    except (ConfigError, TypeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    schema_errors = _schema_errors(normalized_config, schema)
    if schema_errors:
        _raise_schema_errors(
            path=args.config_path,
            phase="normalized",
            errors=schema_errors,
            max_errors=args.max_errors,
        )

    try:
        normalize_config_document(raw_config)
    except (ConfigError, TypeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    print(f"config={args.config_path}")
    print(f"schema={schema_label}")
    print("checks=toml,schema,semantic")
    print("status=ok")


def _raise_schema_errors(
    *,
    path: Path,
    phase: str,
    errors: Sequence[str],
    max_errors: int,
) -> None:
    shown = list(errors[:max_errors])
    lines = [
        f"schema validation failed: {path}",
        f"schema phase={phase}",
        *shown,
    ]
    hidden_count = len(errors) - len(shown)
    if hidden_count > 0:
        lines.append(f"... {hidden_count} more schema error(s) omitted.")
    raise SystemExit("\n".join(lines))


def main(argv: Sequence[str] | None = None) -> None:
    """Run ``beachx lint`` as a standalone entry."""

    args = build_parser(prog="beachx lint").parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
