"""Backward-compatible wrapper for :mod:`beach.cli.inspect_fortran_output`."""

from .cli.inspect_fortran_output import build_parser, main

__all__ = ["build_parser", "main"]
