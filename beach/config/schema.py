"""Shared JSON Schema access and validation helpers for BEACH configs."""

from __future__ import annotations

import json
from collections.abc import Mapping
from importlib import resources
from pathlib import Path
from typing import Any


DEFAULT_SCHEMA_RESOURCE = "schemas/beach.schema.json"


def load_schema(path: Path | None = None) -> tuple[dict[str, Any], str]:
    """Load an explicit schema or the schema packaged with ``beach.config``."""

    if path is not None:
        return load_json_schema(path), str(path)

    schema_resource = resources.files("beach.config").joinpath(DEFAULT_SCHEMA_RESOURCE)
    with schema_resource.open("r", encoding="utf-8") as stream:
        schema = json.load(stream)
    if not isinstance(schema, dict):
        raise ValueError("packaged BEACH schema must decode to a JSON object")
    return schema, f"package:beach.config/{DEFAULT_SCHEMA_RESOURCE}"


def load_json_schema(path: Path) -> dict[str, Any]:
    """Load a JSON Schema document from ``path``."""

    with path.open("r", encoding="utf-8") as stream:
        schema = json.load(stream)
    if not isinstance(schema, dict):
        raise ValueError(f"schema must decode to a JSON object: {path}")
    return schema


def schema_definition_property_names(definition: str) -> frozenset[str]:
    """Return property names declared by one packaged-schema definition."""

    schema, _ = load_schema()
    definitions = schema.get("$defs")
    if not isinstance(definitions, Mapping):
        raise ValueError("BEACH schema is missing the $defs table")
    definition_schema = definitions.get(definition)
    if not isinstance(definition_schema, Mapping):
        raise ValueError(f"BEACH schema is missing $defs.{definition}")
    properties = definition_schema.get("properties")
    if not isinstance(properties, Mapping):
        raise ValueError(f"BEACH schema is missing $defs.{definition}.properties")
    return frozenset(str(key) for key in properties)


def schema_errors(config: Mapping[str, Any], schema: Mapping[str, Any]) -> list[str]:
    """Return stable, path-qualified JSON Schema validation errors."""

    try:
        from jsonschema import Draft7Validator
        from jsonschema.exceptions import SchemaError
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "jsonschema is required for `beachx lint`. "
            "Install BEACH dependencies or run `python -m pip install jsonschema`."
        ) from exc

    try:
        Draft7Validator.check_schema(schema)
    except SchemaError as exc:
        raise SystemExit(f"schema file is invalid: {exc.message}") from exc
    validator = Draft7Validator(schema)
    errors = sorted(
        validator.iter_errors(config),
        key=lambda error: (
            tuple(error.absolute_path),
            tuple(error.absolute_schema_path),
        ),
    )
    return [_format_schema_error(error) for error in errors]


def _format_schema_error(error: Any) -> str:
    path = _format_json_path(tuple(error.absolute_path))
    return f"schema error at {path}: {error.message}"


def _format_json_path(path: tuple[Any, ...]) -> str:
    if not path:
        return "<root>"
    parts: list[str] = []
    for item in path:
        if isinstance(item, int):
            parts.append(f"[{item}]")
        else:
            if parts:
                parts.append(".")
            parts.append(str(item))
    return "".join(parts)
