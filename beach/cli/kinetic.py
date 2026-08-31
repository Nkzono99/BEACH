"""Offline matching-plane kinetic-oracle commands."""

from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path
from typing import Sequence

from beach.outer_kinetic import (
    KineticConfigError,
    QUERY_HEADER,
    convert_kinetic_table,
    load_outer_kinetic_config,
    read_kinetic_queries,
    write_kinetic_atlas,
)

from ._shared import configure_entry_parser

RESPONSE_COMMAND = "kinetic-response"
TABLE_COMMAND = "kinetic-table"


def add_response_subparser(
    subparsers: argparse._SubParsersAction,
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        RESPONSE_COMMAND,
        help="run the offline 1D1V matching-plane kinetic oracle",
        description=(
            "Run fixed matching-plane queries with the offline 1D1V Vlasov "
            "oracle and write certification diagnostics."
        ),
    )
    _configure_response_parser(parser)
    return parser


def add_table_subparser(
    subparsers: argparse._SubParsersAction,
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        TABLE_COMMAND,
        help="convert a certified kinetic atlas to a response table",
        description=(
            "Convert a complete Cartesian subset of a certified kinetic atlas "
            "to matching-plane response CSV v1."
        ),
    )
    _configure_table_parser(parser)
    return parser


def _configure_response_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("config", type=Path, help="outer-kinetic TOML config")
    parser.add_argument("queries", type=Path, help="matching query CSV")
    parser.add_argument("output_directory", type=Path, help="new output directory")
    configure_entry_parser(parser, run_response)


def _configure_table_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("raw_csv", type=Path, help="kinetic_response_raw.csv")
    parser.add_argument("output_csv", type=Path, help="response table CSV to create")
    parser.add_argument(
        "--manifest",
        type=Path,
        help=(
            "kinetic manifest (default: kinetic_response_manifest.json beside raw CSV)"
        ),
    )
    parser.add_argument(
        "--range",
        dest="ranges",
        action="append",
        default=[],
        metavar="AXIS=MIN:MAX",
        help="retain an inclusive range on one query axis; repeatable",
    )
    parser.add_argument(
        "--allow-stationary-average",
        action="store_true",
        help="allow stationary_average rows in addition to steady rows",
    )
    configure_entry_parser(parser, run_table)


def _parse_ranges(values: list[str]) -> dict[str, tuple[float, float]]:
    ranges: dict[str, tuple[float, float]] = {}
    for value in values:
        try:
            axis, encoded = value.split("=", 1)
            lower_text, upper_text = encoded.split(":", 1)
            bounds = (float(lower_text), float(upper_text))
        except ValueError as exc:
            raise ValueError(
                f"invalid --range {value!r}; expected AXIS=MIN:MAX"
            ) from exc
        if axis not in QUERY_HEADER:
            raise ValueError(f"unknown kinetic table range axis: {axis}")
        if not all(math.isfinite(bound) for bound in bounds):
            raise ValueError(f"invalid --range {value!r}; bounds must be finite")
        if bounds[0] > bounds[1]:
            raise ValueError(f"invalid --range {value!r}; MIN must be <= MAX")
        if axis in ranges:
            raise ValueError(f"duplicate --range axis: {axis}")
        ranges[axis] = bounds
    return ranges


def run_response(args: argparse.Namespace) -> None:
    try:
        config = load_outer_kinetic_config(args.config)
        queries = read_kinetic_queries(args.queries)
        results = write_kinetic_atlas(config, queries, args.output_directory)
    except (FileNotFoundError, KineticConfigError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    counts = Counter(result.classification for result in results)
    print(f"output_directory={args.output_directory}")
    print(f"query_count={len(results)}")
    print(
        "classifications="
        + ",".join(f"{name}:{counts[name]}" for name in sorted(counts))
    )


def run_table(args: argparse.Namespace) -> None:
    manifest = args.manifest or args.raw_csv.with_name(
        "kinetic_response_manifest.json"
    )
    try:
        ranges = _parse_ranges(args.ranges)
        row_count = convert_kinetic_table(
            args.raw_csv,
            manifest,
            args.output_csv,
            ranges=ranges,
            allow_stationary_average=args.allow_stationary_average,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    print(f"output_csv={args.output_csv}")
    print(f"rows={row_count}")


def response_main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="beach-kinetic-response")
    _configure_response_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


def table_main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="beach-kinetic-table")
    _configure_table_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)
