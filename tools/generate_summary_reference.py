#!/usr/bin/env python3
"""Update the generated summary scalar table in OutputReference.md."""

from __future__ import annotations

import argparse
from pathlib import Path

from beach.io.summary import summary_fields_markdown

_START = "<!-- BEGIN GENERATED SUMMARY FIELDS -->"
_END = "<!-- END GENERATED SUMMARY FIELDS -->"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="fail if the generated table is stale"
    )
    args = parser.parse_args()
    path = Path(__file__).resolve().parents[1] / "docs/reference/OutputReference.md"
    text = path.read_text(encoding="utf-8")
    begin = text.index(_START) + len(_START)
    end = text.index(_END, begin)
    updated = text[:begin] + "\n\n" + summary_fields_markdown() + "\n\n" + text[end:]
    if args.check:
        if updated != text:
            raise SystemExit(
                "summary reference is stale; run tools/generate_summary_reference.py"
            )
    else:
        path.write_text(updated, encoding="utf-8")


if __name__ == "__main__":
    main()
