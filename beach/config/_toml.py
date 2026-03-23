"""Small TOML helpers for BEACH config tooling."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


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

    ordered_keys = _ordered_keys(document, top_level_order)
    scalars: list[tuple[str, Any]] = []
    nested: list[tuple[str, Any]] = []
    for key in ordered_keys:
        value = document[key]
        if isinstance(value, Mapping) or _is_array_of_tables(value):
            nested.append((key, value))
        else:
            scalars.append((key, value))

    for key, value in scalars:
        lines.append(f"{key} = {_format_toml_value(value)}")

    for index, (key, value) in enumerate(nested):
        if lines and (index > 0 or scalars):
            lines.append("")
        if isinstance(value, Mapping):
            _write_table(lines, (key,), value)
        else:
            _write_array_tables(lines, (key,), value)

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


def _write_table(lines: list[str], path: tuple[str, ...], table: Mapping[str, Any]) -> None:
    lines.append(f"[{'.'.join(path)}]")
    nested_entries: list[tuple[str, str, Any]] = []
    for key, value in table.items():
        if isinstance(value, Mapping) or _is_array_of_tables(value):
            nested_entries.append((key, "nested", value))
            continue
        lines.append(f"{key} = {_format_toml_value(value)}")
    for key, _, value in nested_entries:
        lines.append("")
        if isinstance(value, Mapping):
            _write_table(lines, (*path, key), value)
        else:
            _write_array_tables(lines, (*path, key), value)


def _write_array_tables(
    lines: list[str],
    path: tuple[str, ...],
    items: Sequence[Mapping[str, Any]],
) -> None:
    for index, item in enumerate(items):
        if index > 0:
            lines.append("")
        lines.append(f"[[{'.'.join(path)}]]")
        nested_entries: list[tuple[str, Any]] = []
        for key, value in item.items():
            if isinstance(value, Mapping) or _is_array_of_tables(value):
                nested_entries.append((key, value))
                continue
            lines.append(f"{key} = {_format_toml_value(value)}")
        for key, value in nested_entries:
            lines.append("")
            if isinstance(value, Mapping):
                _write_table(lines, (*path, key), value)
            else:
                _write_array_tables(lines, (*path, key), value)


def _is_array_of_tables(value: object) -> bool:
    return isinstance(value, list) and len(value) > 0 and all(
        isinstance(item, Mapping) for item in value
    )


def _format_toml_value(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped}"'
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return "[" + ", ".join(_format_toml_value(item) for item in value) + "]"
    raise TypeError(f"Unsupported TOML value: {value!r}")
