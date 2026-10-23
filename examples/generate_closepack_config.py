#!/usr/bin/env python3
"""Backward-compatible wrapper for ``beachx model close-pack``."""

import sys

from beach.cli.main import main


if __name__ == "__main__":
    main(["model", "close-pack", *sys.argv[1:]])
