"""Small TOML helpers for BEACH config tooling."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import tomli_w


def load_toml_file(path: str | Path) -> dict[str, Any]:
    """Load one TOML file into a plain Python mapping."""

    resolved = Path(path)
    try:
        import tomllib  # py311+

        with resolved.open("rb") as stream:
            data = tomllib.load(stream)
    except ModuleNotFoundError:
        try:
            import tomli  # type: ignore

            with resolved.open("rb") as stream:
                data = tomli.load(stream)
        except ModuleNotFoundError as exc:
            raise SystemExit(
                "TOML parser is missing. Use Python 3.11+ or install tomli: "
                "`python -m pip install tomli`."
            ) from exc
    if not isinstance(data, dict):
        raise ValueError(f"TOML document must decode to a table: {resolved}")
    return data


def render_toml_document(
    document: Mapping[str, Any],
    *,
    header_comments: Sequence[str] | None = None,
    top_level_order: Sequence[str] | None = None,
) -> str:
    """Render one TOML document from nested Python mappings/lists."""

    lines: list[str] = []
    if header_comments:
        lines.extend(header_comments)
        lines.append("")

    ordered = {key: document[key] for key in _ordered_keys(document, top_level_order)}
    lines.append(tomli_w.dumps(ordered).rstrip())
    return "\n".join(lines).rstrip() + "\n"


def _ordered_keys(
    document: Mapping[str, Any],
    preferred_order: Sequence[str] | None,
) -> list[str]:
    if preferred_order is None:
        return list(document.keys())
    remaining = [key for key in document if key not in preferred_order]
    present_preferred = [key for key in preferred_order if key in document]
    return [*present_preferred, *remaining]
