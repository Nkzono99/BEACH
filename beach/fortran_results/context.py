"""Shared run-result and nearby-config resolution."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Mapping

import numpy as np

from .history import FortranChargeHistory
from .types import FortranRunResult


def resolve_result(value: FortranRunResult | object) -> FortranRunResult:
    """Resolve a result object while preserving the accepted Beach-like protocol."""

    if isinstance(value, FortranRunResult):
        return value

    candidate = getattr(value, "result", None)
    if isinstance(candidate, FortranRunResult):
        return candidate

    raise TypeError("result must be FortranRunResult or Beach.")


def find_config_path_near_output(output_dir: Path) -> Path | None:
    """Find the nearest ``beach.toml`` using the historical search order."""

    candidates = (
        output_dir / "beach.toml",
        output_dir.parent / "beach.toml",
        output_dir.parent.parent / "beach.toml",
    )
    return next((candidate for candidate in candidates if candidate.exists()), None)


def resolve_config_path(
    output_dir: Path,
    *,
    config_path: str | Path | None,
) -> Path | None:
    """Resolve an explicit or nearby config path with stable error behavior."""

    if config_path is None:
        return find_config_path_near_output(output_dir)

    path = Path(config_path)
    if not path.exists():
        raise ValueError(f'config file is not found: "{path}".')
    return path


def load_toml(path: Path) -> dict[str, object]:
    """Load TOML for post-processing helpers using their historical errors."""

    try:
        import tomllib  # py311+

        with path.open("rb") as stream:
            return tomllib.load(stream)
    except ModuleNotFoundError:
        try:
            import tomli  # type: ignore

            with path.open("rb") as stream:
                return tomli.load(stream)
        except ModuleNotFoundError as exc:
            raise ValueError(
                "TOML parser is missing. Use Python 3.11+ or install tomli: "
                "`python -m pip install tomli`."
            ) from exc


def load_config_for_output(
    output_dir: Path,
    *,
    config_path: str | Path | None,
) -> Mapping[str, object] | None:
    """Load an explicit or auto-discovered config for one output directory."""

    path = resolve_config_path(output_dir, config_path=config_path)
    return None if path is None else load_toml(path)


@dataclass(frozen=True)
class RunContext:
    """Resolved run plus lazily loaded, cached nearby configuration."""

    result: FortranRunResult
    requested_config_path: Path | None = None

    @property
    def output_dir(self) -> Path:
        return self.result.directory

    @classmethod
    def from_value(
        cls,
        value: FortranRunResult | object,
        *,
        config_path: str | Path | None = None,
    ) -> "RunContext":
        if isinstance(value, cls):
            if config_path is None:
                return value
            return cls(result=value.result, requested_config_path=Path(config_path))
        inherited_config_path = (
            getattr(value, "config_path", None) if config_path is None else None
        )
        requested_config_path = (
            config_path if config_path is not None else inherited_config_path
        )
        return cls(
            result=resolve_result(value),
            requested_config_path=(
                None if requested_config_path is None else Path(requested_config_path)
            ),
        )

    @cached_property
    def config_path(self) -> Path | None:
        return resolve_config_path(
            self.result.directory,
            config_path=self.requested_config_path,
        )

    @cached_property
    def config(self) -> Mapping[str, object] | None:
        return None if self.config_path is None else load_toml(self.config_path)

    @cached_property
    def sim(self) -> Mapping[str, object] | None:
        if self.config is None:
            return None
        sim = self.config.get("sim")
        return sim if isinstance(sim, Mapping) else None

    def charges_at(self, step: int | None) -> np.ndarray:
        if step is None:
            return self.result.charges

        history = self.result.history
        if history is None or not history.has_data:
            if step == -1:
                return self.result.charges
            raise ValueError(
                "charge_history.csv is required when step is specified and must not be empty."
            )
        if step == -1:
            return history.get_step(-1)
        return history.get_step(step)

    def mesh_ids_or_default(self) -> np.ndarray:
        result = self.result
        if result.mesh_ids is None or result.mesh_ids.size != result.mesh_nelem:
            return np.ones(result.mesh_nelem, dtype=np.int64)
        return result.mesh_ids.astype(np.int64, copy=False)

    def require_triangles(self) -> np.ndarray:
        if self.result.triangles is None:
            raise ValueError(
                "mesh_triangles.csv is not found. Re-run Fortran with latest output format."
            )
        return self.result.triangles

    def require_history(self) -> FortranChargeHistory:
        history = self.result.history
        if history is None or not history.has_data:
            raise ValueError(
                "charge_history.csv is not found or empty. Enable history output and rerun."
            )
        return history
