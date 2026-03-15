"""Backward-compatible wrapper for :mod:`beach.cli.plot_coulomb_force_matrix`."""

from .cli.plot_coulomb_force_matrix import build_parser, main

__all__ = ["build_parser", "main"]
