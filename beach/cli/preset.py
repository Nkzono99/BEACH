"""User preset management commands for ``beachx``."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence

from beach.config import (
    CASE_FILENAME,
    RENDERED_FILENAME,
    CaseSpecError,
    ConfigError,
    default_preset_fragment,
    extract_preset_fragment,
    list_presets,
    load_case_document,
    normalize_preset_section,
    preset_category,
    render_preset_toml,
    render_case_document,
    resolve_preset,
    save_preset_fragment,
    validate_rendered_config,
)
from beach.config._toml import load_toml_file

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
    _add_save_subparser(preset_subparsers)
    _add_list_subparser(preset_subparsers)
    _add_edit_subparser(preset_subparsers)
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


def _add_save_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "save",
        help="extract one section from case.toml or beach.toml into a preset",
    )
    parser.add_argument("name", help="preset name such as sim/lab/my_run_baseline")
    parser.add_argument(
        "--section",
        required=True,
        help="section to extract: sim, mesh, output, species, or mesh.templates",
    )
    parser.add_argument(
        "--index",
        type=int,
        help="1-based item index for species or mesh.templates extraction",
    )
    parser.add_argument(
        "source_path",
        nargs="?",
        type=Path,
        help="source case.toml or rendered beach.toml (default: auto-detect in current directory)",
    )
    parser.add_argument(
        "--rendered",
        action="store_true",
        help="treat the source as a rendered beach.toml instead of auto-detecting",
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
    configure_entry_parser(parser, run_save)


def _add_list_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "list",
        help="list visible presets after precedence resolution",
    )
    configure_entry_parser(parser, run_list)


def _add_edit_subparser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "edit",
        help="open one editable preset in $VISUAL or $EDITOR",
    )
    parser.add_argument("name", help="preset name")
    configure_entry_parser(parser, run_edit)


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
    except (CaseSpecError, ConfigError, FileExistsError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    print(f"saved={destination}")
    print(f"scope={scope}")
    if resolved is not None:
        print(f"from={resolved.name}")


def run_save(args: argparse.Namespace) -> None:
    """Extract and save one preset fragment from case.toml or beach.toml."""

    scope = "project-local" if args.local else "user"
    try:
        section = normalize_preset_section(args.section)
        _validate_name_matches_section(args.name, section=section)
        source_path, config = _load_preset_save_source(
            args.source_path,
            rendered=args.rendered,
        )
        fragment = extract_preset_fragment(
            config,
            section=section,
            index=args.index,
        )
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
    print(f"source={source_path}")
    print(f"section={section}")


def run_list(args: argparse.Namespace) -> None:
    """List visible presets."""

    del args
    presets = list_presets(case_dir=Path.cwd())
    if not presets:
        print("no presets.")
        return
    for preset in presets:
        print(f"{preset.name}\t{preset.scope}\t{preset.source_path}")


def run_edit(args: argparse.Namespace) -> None:
    """Open one editable preset in the configured editor."""

    try:
        preset = resolve_preset(args.name, case_dir=Path.cwd())
    except ConfigError as exc:
        raise SystemExit(str(exc)) from exc

    if preset.scope == "built-in":
        raise SystemExit(
            "built-in presets are read-only. Copy one first with "
            f"`beachx preset new {preset.name} --from {preset.name}`."
        )

    editor_command = _resolve_editor_command()
    completed = subprocess.run([*editor_command, preset.source_label], check=False)
    if completed.returncode != 0:
        raise SystemExit(
            f"editor exited with status {completed.returncode}: {' '.join(editor_command)}"
        )
    print(f"edited={preset.source_label}")


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


def _load_preset_save_source(
    source_path: Path | None,
    *,
    rendered: bool,
) -> tuple[Path, dict[str, Any]]:
    resolved_path = _resolve_source_path(source_path, rendered=rendered)
    if rendered:
        raw = load_toml_file(resolved_path)
        validate_rendered_config(raw)
        return resolved_path, raw

    raw = load_toml_file(resolved_path)
    if _looks_like_case_document(raw):
        return resolved_path, render_case_document(load_case_document(resolved_path)).config
    validate_rendered_config(raw)
    return resolved_path, raw


def _resolve_source_path(source_path: Path | None, *, rendered: bool) -> Path:
    if source_path is not None:
        if not source_path.exists():
            raise ConfigError(f"preset save error: source file not found: {source_path}")
        return source_path

    preferred = Path(RENDERED_FILENAME) if rendered else Path(CASE_FILENAME)
    if preferred.exists():
        return preferred

    fallback = Path(CASE_FILENAME if rendered else RENDERED_FILENAME)
    if fallback.exists():
        return fallback

    raise ConfigError(
        f"preset save error: neither {CASE_FILENAME} nor {RENDERED_FILENAME} exists in the current directory."
    )


def _looks_like_case_document(document: dict[str, Any]) -> bool:
    return any(
        key in document
        for key in ("schema_version", "title", "use_presets", "override")
    )


def _validate_name_matches_section(name: str, *, section: str) -> None:
    category = preset_category(name)
    allowed = {
        "sim": {"sim"},
        "mesh": {"mesh"},
        "output": {"output"},
        "species": {"species"},
        "mesh.templates": {"mesh"},
    }[section]
    if category not in allowed:
        raise CaseSpecError(
            f"case spec error: preset name {name!r} does not match --section {section!r}."
        )


def _resolve_editor_command() -> list[str]:
    for env_key in ("VISUAL", "EDITOR"):
        value = os.environ.get(env_key, "").strip()
        if value:
            return shlex.split(value)
    for command in ("nano", "vim", "vi"):
        path = shutil.which(command)
        if path is not None:
            return [path]
    raise SystemExit(
        "no editor found. Set $VISUAL or $EDITOR, or install vi/nano/vim."
    )


if __name__ == "__main__":
    main()
