"""Grouped model/config generators for ``beachx``."""

from __future__ import annotations

import argparse

from . import model_close_pack

COMMAND_NAME = "model"


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register the ``model`` command group under ``beachx``."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="generate parameterized model configs",
        description="Generate parameterized BEACH configuration files.",
    )
    model_subparsers = parser.add_subparsers(
        dest="model_name",
        metavar="model",
        required=True,
    )
    model_close_pack.add_subparser(model_subparsers)
    return parser
