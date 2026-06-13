"""Backward-compatible wrapper for :mod:`beach.cli.estimate_fortran_workload`."""

from .cli.estimate_fortran_workload import (
    ALLOWED_SIM_KEYS,
    ALLOWED_SPECIES_KEYS,
    DEFAULT_SIM,
    DEFAULT_SPECIES,
    build_parser,
    completed_batches_from_resume_config,
    estimate_workload,
    load_toml,
    main,
    read_macro_residuals,
    read_summary_batches,
)

__all__ = [
    "ALLOWED_SIM_KEYS",
    "ALLOWED_SPECIES_KEYS",
    "DEFAULT_SIM",
    "DEFAULT_SPECIES",
    "load_toml",
    "read_macro_residuals",
    "read_summary_batches",
    "completed_batches_from_resume_config",
    "estimate_workload",
    "build_parser",
    "main",
]
