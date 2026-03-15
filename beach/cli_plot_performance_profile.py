"""Backward-compatible wrapper for :mod:`beach.cli.plot_performance_profile`."""

from .cli.plot_performance_profile import build_parser as build_parser
from .cli.plot_performance_profile import main as main

__all__ = ["build_parser", "main"]
