"""Command table, parser construction, and dispatch for ``beachx``.

Each leaf command is one module under :mod:`beach.cli` that defines
``configure(parser)`` and ``run(args)``. A command whose ``configure`` adds the
shared ``output_dir`` argument (:func:`beach.cli._options.add_output_dir`) instead
defines ``run(args, beach)``: the registry
opens the run directory once as a :class:`~beach.Beach` session (with the
command's ``--config`` and ``--library``), reports missing output files with one
message, and closes the session's kernels after ``run`` returns.
The order of :data:`COMMANDS` is the order in ``beachx --help``.
"""

from __future__ import annotations

import argparse
import importlib
from collections.abc import Sequence
from dataclasses import dataclass

_MISSING_OUTPUT = (
    'Fortran output files are missing under "{}". '
    "Expected at least summary.txt and charges.csv."
)
_INSTALL_HINT = "Install dependencies with `python -m pip install -e . --no-build-isolation`."


@dataclass(frozen=True)
class Command:
    """One ``beachx`` subcommand: a leaf module or a group of subcommands."""

    name: str
    help: str
    target: str | tuple[Command, ...]
    description: str | None = None


COMMANDS: tuple[Command, ...] = (
    Command("inspect", "inspect Fortran output files", "inspect_fortran_output"),
    Command("animate", "animate charge or potential history", "animate_fortran_history"),
    Command("coulomb", "plot object-wise Coulomb force matrix", "plot_coulomb_force_matrix"),
    Command(
        "mobility",
        "analyze object-wise Coulomb mobility indicators",
        "analyze_coulomb_mobility",
    ),
    Command(
        "kernel-forces",
        "write object-wise net forces from the Fortran field kernel",
        "kernel_forces",
    ),
    Command(
        "object-detachment",
        "evaluate periodic object force, work, and release speed",
        "object_detachment",
    ),
    Command(
        "slices",
        "plot potential slices inside the simulation box",
        "plot_fortran_potential_slices",
    ),
    Command(
        "workload",
        "estimate per-batch particle workload",
        "estimate_fortran_workload",
        description=(
            "Estimate particle workload from Fortran TOML config: "
            "per-batch local(rank) particles and per-thread particle counts."
        ),
    ),
    Command("profile", "plot one performance_profile.csv summary", "plot_performance_profile"),
    Command(
        "lint",
        "lint beach.toml with the JSON Schema and the Fortran preflight",
        "lint",
        description=(
            "Validate a BEACH config with TOML parsing, the packaged JSON Schema, "
            "and the Fortran preflight (beach --check-config)."
        ),
    ),
    Command(
        "config",
        "create and diff beach.toml",
        (
            Command("init", "create a runnable beach.toml", "config_init"),
            Command("diff", "show semantic differences between two configs", "config_diff"),
        ),
        description="Create and diff BEACH config files.",
    ),
    Command(
        "model",
        "generate parameterized model configs",
        (
            Command(
                "close-pack",
                "generate a close-packed sphere-stack beach.toml",
                "model_close_pack",
            ),
        ),
        description="Generate parameterized BEACH configuration files.",
    ),
)


def build_parser() -> argparse.ArgumentParser:
    """Build the root parser for ``beachx``."""

    parser = argparse.ArgumentParser(
        prog="beachx",
        description="BEACH tools for post-processing, visualization, and model generation.",
    )
    _add_commands(parser, COMMANDS, dest="command")
    return parser


def _add_commands(
    parser: argparse.ArgumentParser,
    commands: Sequence[Command],
    *,
    dest: str,
) -> None:
    subparsers = parser.add_subparsers(dest=dest, metavar="command", required=True)
    for command in commands:
        subparser = subparsers.add_parser(
            command.name,
            help=command.help,
            description=command.description,
        )
        if isinstance(command.target, tuple):
            _add_commands(subparser, command.target, dest=f"{dest}_{command.name}")
            continue
        module = importlib.import_module(f"{__package__}.{command.target}")
        module.configure(subparser)
        subparser.set_defaults(_run=module.run, _parser=subparser)


def main(argv: Sequence[str] | None = None) -> None:
    """Parse ``argv`` and run the selected ``beachx`` subcommand."""

    from beach.native.kernel import FieldKernelError

    args = build_parser().parse_args(argv)
    try:
        if getattr(args, "_opens_run", False):
            with _open_run(args) as beach:
                args._run(args, beach)
        else:
            args._run(args)
    except ModuleNotFoundError as exc:
        package = (exc.name or "").partition(".")[0]
        if package in {"matplotlib", "PIL"}:
            name = "pillow" if package == "PIL" else package
            raise SystemExit(f"{name} is required for this command. {_INSTALL_HINT}") from exc
        raise
    except (ValueError, FieldKernelError) as exc:
        raise SystemExit(str(exc)) from exc


def _open_run(args: argparse.Namespace):
    from beach import Beach

    beach = Beach(
        args.output_dir,
        config_path=getattr(args, "config", None),
        library_path=getattr(args, "library", None),
    )
    try:
        beach.result
    except FileNotFoundError as exc:
        raise SystemExit(_MISSING_OUTPUT.format(args.output_dir)) from exc
    return beach
