"""Backward-compatible wrapper for :mod:`beach.cli.estimate_fortran_workload`."""

from .cli.estimate_fortran_workload import (
    ALLOWED_SIM_KEYS,
    ALLOWED_SPECIES_KEYS,
    DEFAULT_SIM,
    DEFAULT_SPECIES,
    build_parser,
    estimate_workload,
    load_toml,
    main,
    read_macro_residuals,
)

__all__ = [
    "ALLOWED_SIM_KEYS",
    "ALLOWED_SPECIES_KEYS",
    "DEFAULT_SIM",
    "DEFAULT_SPECIES",
    "load_toml",
    "read_macro_residuals",
    "estimate_workload",
    "build_parser",
    "main",
]
