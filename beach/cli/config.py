"""Preset-based config management commands for ``beachx``."""

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Sequence

from beach.config import (
    CASE_FILENAME,
    DEFAULT_PRESET_NAMES,
    RENDERED_FILENAME,
    CaseSpecError,
    ConfigError,
    default_case_document,
    list_saved_cases,
    load_case_document,
    load_saved_case_path,
    render_beach_toml,
    render_case_file,
    render_case_toml,
    save_case_file,
    semantic_diff,
    validate_rendered_config,
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
        help="compose beach.toml from presets and case.toml overrides",
        description="Compose and validate BEACH case files from reusable presets.",
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
    _add_render_subparser(config_subparsers)
    _add_validate_subparser(config_subparsers)
    _add_diff_subparser(config_subparsers)
    _add_save_subparser(config_subparsers)
    _add_list_saved_subparser(config_subparsers)


def _add_init_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "init",
        help="create a new editable case.toml",
    )
    parser.add_argument(
        "output",
        nargs="?",
        default=Path(CASE_FILENAME),
        type=Path,
        help=f"destination case path (default: ./{CASE_FILENAME})",
    )
    parser.add_argument(
        "--preset",
        action="append",
        dest="presets",
        help="preset name to include; repeat to combine multiple presets",
    )
    parser.add_argument(
        "--title",
        help="optional title saved into case.toml",
    )
    parser.add_argument(
        "--from",
        dest="from_saved",
        help="copy a saved case template from ~/.config/beachx/cases",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite the destination if it already exists",
    )
    configure_entry_parser(parser, run_init)


def _add_render_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "render",
        help="render case.toml into the final beach.toml",
    )
    parser.add_argument(
        "case_path",
        nargs="?",
        default=Path(CASE_FILENAME),
        type=Path,
        help=f"input case file (default: ./{CASE_FILENAME})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="destination beach.toml path (default: alongside case.toml)",
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="print rendered TOML instead of writing a file",
    )
    configure_entry_parser(parser, run_render)


def _add_validate_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "validate",
        help="validate case.toml and the rendered beach.toml structure",
    )
    parser.add_argument(
        "case_path",
        nargs="?",
        default=Path(CASE_FILENAME),
        type=Path,
        help=f"input case file (default: ./{CASE_FILENAME})",
    )
    configure_entry_parser(parser, run_validate)


def _add_diff_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "diff",
        help="show semantic differences between two case/rendered configs",
    )
    parser.add_argument("left", type=Path, help="left-hand case.toml or beach.toml")
    parser.add_argument("right", type=Path, help="right-hand case.toml or beach.toml")
    parser.add_argument(
        "--rendered",
        action="store_true",
        help="render both inputs first, then diff the final beach.toml payloads",
    )
    configure_entry_parser(parser, run_diff)


def _add_save_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "save",
        help="save the current case.toml into the local BEACHX case store",
    )
    parser.add_argument("name", help="saved case name, for example cavity-base")
    parser.add_argument(
        "case_path",
        nargs="?",
        default=Path(CASE_FILENAME),
        type=Path,
        help=f"source case file (default: ./{CASE_FILENAME})",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing saved case",
    )
    configure_entry_parser(parser, run_save)


def _add_list_saved_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "list-saved",
        help="list saved case templates under ~/.config/beachx/cases",
    )
    configure_entry_parser(parser, run_list_saved)


def run_init(args: argparse.Namespace) -> None:
    """Create one new ``case.toml`` file."""

    destination = args.output
    if destination.exists() and not args.force:
        raise SystemExit(f"case file already exists: {destination}")

    if args.from_saved is not None and args.presets:
        raise SystemExit("cannot combine --from with explicit --preset values")

    if args.from_saved is not None:
        try:
            saved_path = load_saved_case_path(args.from_saved)
        except (CaseSpecError, FileNotFoundError) as exc:
            raise SystemExit(str(exc)) from exc
        if args.title is None:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(saved_path, destination)
            print(f"saved={destination}")
            print(f"source_saved_case={saved_path}")
            return

        try:
            base_case = load_case_document(saved_path)
        except ConfigError as exc:
            raise SystemExit(str(exc)) from exc
        case_document = replace(base_case, title=args.title, source_path=None)
    else:
        presets = tuple(args.presets) if args.presets else DEFAULT_PRESET_NAMES
        try:
            case_document = default_case_document(use_presets=presets, title=args.title)
        except CaseSpecError as exc:
            raise SystemExit(str(exc)) from exc

    text = render_case_toml(case_document)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    print(f"saved={destination}")
    print(f"use_presets={list(case_document.use_presets)}")


def run_render(args: argparse.Namespace) -> None:
    """Render one case file into the final BEACH config."""

    try:
        result = render_case_file(args.case_path)
    except FileNotFoundError as exc:
        raise SystemExit(f"case file not found: {args.case_path}") from exc
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc

    for warning in result.warnings:
        print(warning, file=sys.stderr)

    text = render_beach_toml(result.config, source_case=args.case_path)
    if args.stdout:
        print(text, end="")
        return

    output_path = args.output
    if output_path is None:
        output_path = args.case_path.parent / RENDERED_FILENAME
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")
    print(f"saved={output_path}")
    print(f"preset_count={len(result.presets)}")


def run_validate(args: argparse.Namespace) -> None:
    """Validate one case file and the rendered final config."""

    try:
        result = render_case_file(args.case_path)
    except FileNotFoundError as exc:
        raise SystemExit(f"case file not found: {args.case_path}") from exc
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc

    for warning in result.warnings:
        print(warning, file=sys.stderr)

    print(f"case={args.case_path}")
    print("status=ok")
    print(f"preset_count={len(result.presets)}")


def run_diff(args: argparse.Namespace) -> None:
    """Show semantic differences between two configs."""

    try:
        left_kind, left_document = _load_diff_document(args.left, rendered=args.rendered)
        right_kind, right_document = _load_diff_document(args.right, rendered=args.rendered)
    except FileNotFoundError as exc:
        raise SystemExit(f"config file not found: {exc.filename}") from exc
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc

    if left_kind != right_kind:
        raise SystemExit(
            "diff error: inputs must both be case documents or both be rendered configs. "
            "Use --rendered to compare the final merged beach.toml payloads."
        )

    lines = semantic_diff(left_document, right_document)
    if not lines:
        print("no differences.")
        return
    for line in lines:
        print(line)


def run_save(args: argparse.Namespace) -> None:
    """Save the current case file to the BEACHX case store."""

    try:
        destination = save_case_file(args.case_path, args.name, force=args.force)
    except FileNotFoundError as exc:
        raise SystemExit(f"case file not found: {args.case_path}") from exc
    except (CaseSpecError, FileExistsError) as exc:
        raise SystemExit(str(exc)) from exc

    print(f"saved={destination}")


def run_list_saved(args: argparse.Namespace) -> None:
    """List saved case templates."""

    del args
    names = list_saved_cases()
    if not names:
        print("no saved cases.")
        return
    for name in names:
        print(name)


def main(argv: Sequence[str] | None = None) -> None:
    """Run the ``beachx config`` group as a standalone entry."""

    args = build_parser(prog="beachx config").parse_args(argv)
    args.func(args)


def _load_diff_document(path: Path, *, rendered: bool) -> tuple[str, Any]:
    if rendered:
        result = render_case_file(path)
        return "rendered", result.config

    raw = load_toml_file(path)
    if _looks_like_case_document(raw):
        case_document = load_case_document(path)
        return "case", case_document.to_dict()

    validate_rendered_config(raw)
    return "rendered", raw


def _looks_like_case_document(document: dict[str, Any]) -> bool:
    return any(
        key in document
        for key in ("schema_version", "title", "use_presets", "override")
    )


if __name__ == "__main__":
    main()
