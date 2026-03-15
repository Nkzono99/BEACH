from __future__ import annotations

import argparse

import pytest

from beach.cli import (
    animate_fortran_history,
    estimate_fortran_workload,
    inspect_fortran_output,
    legacy,
    plot_coulomb_force_matrix,
    plot_fortran_potential_slices,
    plot_performance_profile,
)
from beach.cli.main import build_parser as build_beachx_parser
from beach.cli.main import main as beachx_main


CLI_CASES = [
    (
        "inspect",
        inspect_fortran_output.build_parser,
        legacy.inspect_main,
        ["no_such_dir"],
        'Fortran output files are missing under "no_such_dir". '
        "Expected at least summary.txt and charges.csv.",
        "WARNING: `beach-inspect` is deprecated; use `beachx inspect` instead.",
    ),
    (
        "animate",
        animate_fortran_history.build_parser,
        legacy.animate_main,
        ["no_such_dir"],
        'Fortran output files are missing under "no_such_dir". '
        "Expected at least summary.txt and charges.csv.",
        "WARNING: `beach-animate-history` is deprecated; use `beachx animate` instead.",
    ),
    (
        "coulomb",
        plot_coulomb_force_matrix.build_parser,
        legacy.coulomb_main,
        ["no_such_dir"],
        'Fortran output files are missing under "no_such_dir". '
        "Expected at least summary.txt and charges.csv.",
        "WARNING: `beach-plot-coulomb-force-matrix` is deprecated; use `beachx coulomb` instead.",
    ),
    (
        "slices",
        plot_fortran_potential_slices.build_parser,
        legacy.slices_main,
        ["no_such_dir"],
        'Fortran output files are missing under "no_such_dir". '
        "Expected at least summary.txt and charges.csv.",
        "WARNING: `beach-plot-potential-slices` is deprecated; use `beachx slices` instead.",
    ),
    (
        "workload",
        estimate_fortran_workload.build_parser,
        legacy.workload_main,
        ["no_such_file.toml"],
        "config file not found: no_such_file.toml",
        "WARNING: `beach-estimate-workload` is deprecated; use `beachx workload` instead.",
    ),
    (
        "profile",
        plot_performance_profile.build_parser,
        legacy.profile_main,
        ["no_such_dir"],
        'Performance profile file is missing under "no_such_dir/performance_profile.csv". '
        "Expected performance_profile.csv.",
        "WARNING: `beach-plot-performance-profile` is deprecated; use `beachx profile` instead.",
    ),
]


def _get_subparser(command: str) -> argparse.ArgumentParser:
    parser = build_beachx_parser()
    subparsers = next(
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices[command]


def _parser_signature(parser: argparse.ArgumentParser) -> list[tuple[object, ...]]:
    signature: list[tuple[object, ...]] = []
    for action in parser._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        choices = action.choices
        if choices is None:
            normalized_choices: object = None
        elif isinstance(choices, dict):
            normalized_choices = tuple(sorted(choices))
        else:
            normalized_choices = tuple(choices)
        signature.append(
            (
                type(action).__name__,
                tuple(action.option_strings),
                action.dest,
                action.nargs,
                action.required,
                normalized_choices,
            )
        )
    return signature


def test_beachx_help_lists_all_subcommands(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["--help"])

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "inspect" in captured.out
    assert "animate" in captured.out
    assert "coulomb" in captured.out
    assert "slices" in captured.out
    assert "workload" in captured.out
    assert "profile" in captured.out


@pytest.mark.parametrize(
    ("command", "build_legacy_parser", "legacy_main", "_argv", "_message", "_warning"),
    CLI_CASES,
)
def test_beachx_subparser_matches_legacy_parser_shape(
    command: str,
    build_legacy_parser,
    legacy_main,
    _argv: list[str],
    _message: str,
    _warning: str,
) -> None:
    del legacy_main, _argv, _message, _warning
    assert _parser_signature(_get_subparser(command)) == _parser_signature(
        build_legacy_parser()
    )


@pytest.mark.parametrize(
    ("command", "_build_legacy_parser", "legacy_main", "argv", "message", "warning"),
    CLI_CASES,
)
def test_legacy_alias_warns_and_matches_beachx_errors(
    capsys: pytest.CaptureFixture[str],
    command: str,
    _build_legacy_parser,
    legacy_main,
    argv: list[str],
    message: str,
    warning: str,
) -> None:
    del _build_legacy_parser

    with pytest.raises(SystemExit) as legacy_exc:
        legacy_main(argv)
    legacy_streams = capsys.readouterr()

    with pytest.raises(SystemExit) as beachx_exc:
        beachx_main([command, *argv])
    beachx_streams = capsys.readouterr()

    assert str(legacy_exc.value) == message
    assert str(beachx_exc.value) == message
    assert legacy_streams.err.strip() == warning
    assert beachx_streams.err == ""
