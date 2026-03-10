"""Backward-compatible wrapper for :mod:`beach.cli.plot_fortran_potential_slices`."""

from .cli.plot_fortran_potential_slices import build_parser, main

__all__ = ["build_parser", "main"]
