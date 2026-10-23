"""Config lint command for ``beachx``."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Sequence

from beach.config import (
    CONFIG_FILENAME,
    ConfigError,
)
from beach.config.core import _load_config_file
from beach.config.schema import ConfigSchemaError
from beach.config.schema import load_schema as _load_schema

NATIVE_CHECKER_ENV = "BEACH_CONFIG_CHECK_EXE"
NATIVE_CHECKER_NAME = "beach"


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``lint`` options."""

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
        help="additional JSON Schema constraints (the packaged BEACH contract always applies)",
    )
    parser.add_argument(
        "--max-errors",
        type=int,
        default=20,
        help="maximum schema errors to print (default: 20)",
    )
    parser.add_argument(
        "--beach",
        type=Path,
        help=(
            "beach executable for the preflight check "
            f"(default: ${NATIVE_CHECKER_ENV}, then `{NATIVE_CHECKER_NAME}` on PATH)"
        ),
    )


def resolve_native_checker(explicit: Path | None = None) -> Path:
    """Locate the ``beach`` executable that owns the semantic configuration checks."""

    if explicit is not None:
        path = explicit
    elif os.environ.get(NATIVE_CHECKER_ENV):
        path = Path(os.environ[NATIVE_CHECKER_ENV])
    else:
        found = shutil.which(NATIVE_CHECKER_NAME)
        if found is None:
            raise SystemExit(
                f"{NATIVE_CHECKER_NAME} executable not found: install BEACH, "
                f"set {NATIVE_CHECKER_ENV}, or pass --beach"
            )
        path = Path(found)
    if not path.is_file():
        raise SystemExit(f"{NATIVE_CHECKER_NAME} executable not found: {path}")
    return path


def run_native_check(checker: Path, config_path: Path) -> None:
    """Run ``beach --check-config`` and relay its diagnostics on failure."""

    result = subprocess.run(
        [str(checker), "--check-config", str(config_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode == 0 and "status=ok" in result.stdout:
        return
    diagnostics = [
        line
        for line in result.stderr.splitlines()
        if line.strip()
        and not line.startswith("#")
        and not line.startswith("Error termination")
        and not line.startswith("ERROR STOP 1")
    ]
    if not diagnostics:
        diagnostics = [f"exit code {result.returncode}"]
    raise SystemExit(
        "\n".join([f"preflight check failed: {config_path}", f"beach={checker}", *diagnostics])
    )


def run(args: argparse.Namespace) -> None:
    """Lint one BEACH config file."""

    if args.max_errors < 1:
        raise SystemExit("--max-errors must be >= 1")

    try:
        schema, schema_label = _load_schema()
        additional_schema = None
        if args.schema is not None:
            additional_schema, schema_label = _load_schema(args.schema)
    except FileNotFoundError as exc:
        raise SystemExit(f"schema file not found: {exc.filename}") from exc
    except json.JSONDecodeError as exc:
        raise SystemExit(f"schema JSON parse error: {exc}") from exc
    except ValueError as exc:
        raise SystemExit(f"schema error: {exc}") from exc

    try:
        _load_config_file(
            args.config_path, schema=schema, additional_schema=additional_schema,
        )
    except ConfigSchemaError as exc:
        _raise_schema_errors(
            path=args.config_path,
            phase=exc.phase,
            errors=exc.errors,
            max_errors=args.max_errors,
        )
    except FileNotFoundError as exc:
        raise SystemExit(f"config file not found: {exc.filename}") from exc
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc
    except ValueError as exc:
        raise SystemExit(f"TOML parse error: {exc}") from exc

    checker = resolve_native_checker(getattr(args, "beach", None))
    run_native_check(checker, args.config_path)

    print(f"config={args.config_path}")
    print(f"schema={schema_label}")
    print(f"beach={checker}")
    print("checks=toml,schema,preflight")
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
