"""Backward-compatible wrapper for :mod:`beach.cli.animate_fortran_history`."""

from .cli.animate_fortran_history import build_parser, main

__all__ = ["build_parser", "main"]
