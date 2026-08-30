from __future__ import annotations

from pathlib import Path

import pytest

from beach.summary import (
    CORE_SUMMARY_REQUIRED_KEYS,
    load_summary_file,
    parse_summary_text,
)


CORE_SUMMARY_TEXT = """\
mesh_nelem=1
processed_particles=4
absorbed=2
escaped=2
batches=1
last_rel_change=0.5
"""


def test_summary_parser_normalizes_whitespace_and_preserves_value_delimiters() -> None:
    assert parse_summary_text(
        "\n mesh_nelem = 2 \nperiodic2_cache_path = cache=generated\n"
    ) == {
        "mesh_nelem": "2",
        "periodic2_cache_path": "cache=generated",
    }


@pytest.mark.parametrize(
    ("text", "match"),
    [
        ("batches=1\n batches =2\n", "duplicate key 'batches'"),
        (" = value\n", "empty key"),
    ],
)
def test_summary_parser_rejects_invalid_keys(text: str, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        parse_summary_text(text)


def test_summary_loader_accepts_core_and_unknown_optional_keys(tmp_path: Path) -> None:
    path = tmp_path / "summary.txt"
    path.write_text(
        CORE_SUMMARY_TEXT + "future_receipt=enabled\n",
        encoding="utf-8",
    )

    values = load_summary_file(path, required_keys=CORE_SUMMARY_REQUIRED_KEYS)

    assert values["future_receipt"] == "enabled"


def test_summary_loader_enforces_required_keys(tmp_path: Path) -> None:
    path = tmp_path / "summary.txt"
    path.write_text(CORE_SUMMARY_TEXT.replace("escaped=2\n", ""), encoding="utf-8")

    with pytest.raises(ValueError, match="escaped"):
        load_summary_file(path, required_keys=CORE_SUMMARY_REQUIRED_KEYS)
