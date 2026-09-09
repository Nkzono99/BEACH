"""Shared JSON Schema access and validation helpers for BEACH configs."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from importlib import resources
from pathlib import Path
from typing import Any

from ._shared import ConfigValidationError


DEFAULT_SCHEMA_RESOURCE = "schemas/beach.schema.json"


class ConfigSchemaError(ConfigValidationError):
    """A configuration rejected by a named phase of the shared schema check."""

    def __init__(self, errors: list[Any], *, phase: str) -> None:
        self.phase = phase
        self.errors = [_format_schema_error(error) for error in errors]
        super().__init__("\n".join(self.errors))


def prepare_schema_document(
    config: Mapping[str, Any], schema: Mapping[str, Any]
) -> dict[str, Any]:
    """Copy input, normalize declared identifiers, and reject nonfinite numbers."""

    def copy_value(value: Any, rule: Mapping[str, Any], path: tuple[Any, ...]) -> Any:
        while "$ref" in rule:
            reference = rule["$ref"]
            if not reference.startswith("#/"):
                break
            rule = schema
            for part in reference[2:].split("/"):
                rule = rule[part.replace("~1", "/").replace("~0", "~")]
        if isinstance(value, Mapping):
            properties = rule.get("properties", {})
            additional = rule.get("additionalProperties", {})
            if not isinstance(additional, Mapping):
                additional = {}
            result = {}
            for key, item in value.items():
                normalized_key = key.lower() if isinstance(key, str) else key
                if normalized_key not in properties:
                    normalized_key = key
                if normalized_key in result:
                    raise ConfigValidationError(
                        f"config error at {_format_json_path(path)}: duplicate key {normalized_key!r}."
                    )
                result[normalized_key] = copy_value(
                    item, properties.get(normalized_key, additional), (*path, normalized_key)
                )
            return result
        if isinstance(value, (list, tuple)):
            return [copy_value(item, rule.get("items", {}), (*path, index))
                    for index, item in enumerate(value)]
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            try:
                finite = math.isfinite(value)
            except OverflowError:
                finite = False
            if not finite:
                raise ConfigValidationError(
                    f"config error at {_format_json_path(path)}: number must be finite."
                )
        if isinstance(value, str):
            value = value.rstrip(" ")
            choices = rule.get("enum", [rule.get("const")])
            normalized = value.lower()
            if normalized in choices:
                return normalized
        return value

    return copy_value(config, schema, ())


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

    return [_format_schema_error(error) for error in validation_errors(config, schema)]


def validation_errors(config: Mapping[str, Any], schema: Mapping[str, Any]) -> list[Any]:
    """Return schema errors for the shared authoring/runtime validation pipeline."""

    try:
        from jsonschema import Draft7Validator, validators
        from jsonschema.exceptions import SchemaError
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "jsonschema is required for BEACH configuration validation. "
            "Install BEACH dependencies or run `python -m pip install jsonschema`."
        ) from exc

    try:
        Draft7Validator.check_schema(schema)
    except SchemaError as exc:
        raise SystemExit(f"schema file is invalid: {exc.message}") from exc
    # TOML integers are distinct from reals, including mathematically integral 1.0.
    integer_types = Draft7Validator.TYPE_CHECKER.redefine(
        "integer", lambda checker, value: isinstance(value, int)
        and not isinstance(value, bool) and -(2**31) <= value < 2**31
    )
    validator = validators.extend(Draft7Validator, type_checker=integer_types)(schema)
    errors = sorted(
        validator.iter_errors(config),
        key=lambda error: (
            tuple(str(part) for part in error.absolute_path),
            tuple(str(part) for part in error.absolute_schema_path),
        ),
    )
    return errors


def reject_basic_schema_errors(errors: list[Any], *, phase: str) -> None:
    """Check values before interpretation; semantic validators explain combinations."""

    composition = {"allOf", "anyOf", "oneOf", "not", "if", "then", "else"}
    basic_errors = [error for error in errors
                    if not composition.intersection(error.absolute_schema_path)]
    if basic_errors:
        raise ConfigSchemaError(basic_errors, phase=phase)


def _format_schema_error(error: Any) -> str:
    path = _format_json_path(tuple(error.absolute_path))
    if (error.validator == "type" and error.validator_value == "integer"
            and isinstance(error.instance, int) and not isinstance(error.instance, bool)):
        return f"schema error at {path}: integer must fit the signed 32-bit range."
    if error.validator == "maxLength":
        return f"schema error at {path}: string must contain at most {error.validator_value} characters."
    if error.validator == "minLength" and error.validator_value == 1:
        return f"schema error at {path}: expected a non-empty string."
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
