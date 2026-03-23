"""User preset management commands for ``beachx``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from beach.config import (
    CaseSpecError,
    ConfigError,
    default_preset_fragment,
    list_presets,
    render_preset_toml,
    resolve_preset,
    save_preset_fragment,
)

from ._shared import configure_entry_parser

COMMAND_NAME = "preset"


def build_parser(*, prog: str | None = None) -> argparse.ArgumentParser:
    """Build the parser for the ``beachx preset`` group."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_group_parser(parser)
    return parser


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register the ``preset`` group under the root CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="create and inspect reusable config presets",
        description="Manage project-local and user-level BEACH preset fragments.",
    )
    _configure_group_parser(parser)
    return parser


def _configure_group_parser(parser: argparse.ArgumentParser) -> None:
    preset_subparsers = parser.add_subparsers(
        dest="preset_command",
        metavar="subcommand",
        required=True,
    )
    _add_new_subparser(preset_subparsers)
    _add_list_subparser(preset_subparsers)
    _add_show_subparser(preset_subparsers)
    _add_path_subparser(preset_subparsers)
    _add_validate_subparser(preset_subparsers)


def _add_new_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "new",
        help="create a new user or project-local preset",
    )
    parser.add_argument("name", help="preset name such as sim/lab/periodic2_fast")
    parser.add_argument(
        "--from",
        dest="from_preset",
        help="clone an existing preset into the new location",
    )
    parser.add_argument(
        "--local",
        action="store_true",
        help="save under ./.beachx/presets instead of ~/.config/beachx/presets",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite the destination if it already exists",
    )
    configure_entry_parser(parser, run_new)


def _add_list_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "list",
        help="list visible presets after precedence resolution",
    )
    configure_entry_parser(parser, run_list)


def _add_show_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "show",
        help="print one resolved preset fragment",
    )
    parser.add_argument("name", help="preset name")
    configure_entry_parser(parser, run_show)


def _add_path_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "path",
        help="print the resolved file path of one preset",
    )
    parser.add_argument("name", help="preset name")
    configure_entry_parser(parser, run_path)


def _add_validate_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "validate",
        help="validate one preset fragment",
    )
    parser.add_argument("name", help="preset name")
    configure_entry_parser(parser, run_validate)


def run_new(args: argparse.Namespace) -> None:
    """Create one new preset fragment."""

    scope = "project-local" if args.local else "user"
    try:
        if args.from_preset is not None:
            resolved = resolve_preset(args.from_preset, case_dir=Path.cwd())
            fragment = resolved.data
        else:
            resolved = None
            fragment = default_preset_fragment(args.name)
        destination = save_preset_fragment(
            args.name,
            fragment,
            scope=scope,
            case_dir=Path.cwd(),
            force=args.force,
        )
    except (CaseSpecError, ConfigError, FileExistsError) as exc:
        raise SystemExit(str(exc)) from exc

    print(f"saved={destination}")
    print(f"scope={scope}")
    if resolved is not None:
        print(f"from={resolved.name}")


def run_list(args: argparse.Namespace) -> None:
    """List visible presets."""

    del args
    presets = list_presets(case_dir=Path.cwd())
    if not presets:
        print("no presets.")
        return
    for preset in presets:
        print(f"{preset.name}\t{preset.scope}\t{preset.source_path}")


def run_show(args: argparse.Namespace) -> None:
    """Print the resolved preset fragment."""

    try:
        preset = resolve_preset(args.name, case_dir=Path.cwd())
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc

    for source in preset.shadowed_sources:
        print(
            f'warning: preset "{preset.name}" resolved to {preset.source_label} and shadows {source}',
            file=sys.stderr,
        )
    print(
        render_preset_toml(
            preset.data,
            preset_name=preset.name,
            source_label=preset.source_label,
        ),
        end="",
    )


def run_path(args: argparse.Namespace) -> None:
    """Print the resolved preset path."""

    try:
        preset = resolve_preset(args.name, case_dir=Path.cwd())
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc
    print(preset.source_label)


def run_validate(args: argparse.Namespace) -> None:
    """Validate the resolved preset fragment."""

    try:
        preset = resolve_preset(args.name, case_dir=Path.cwd())
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"preset={preset.name}")
    print(f"scope={preset.scope}")
    print("status=ok")


def main(argv: Sequence[str] | None = None) -> None:
    """Run the ``beachx preset`` group as a standalone entry."""

    args = build_parser(prog="beachx preset").parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
