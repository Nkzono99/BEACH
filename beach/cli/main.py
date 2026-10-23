"""``beachx`` entry point; the command table lives in :mod:`beach.cli._registry`."""

from __future__ import annotations

from ._registry import build_parser, main

__all__ = ["build_parser", "main"]


if __name__ == "__main__":
    main()
