from __future__ import annotations

import re
from pathlib import Path

import pytest

from beach.summary import (
    CORE_SUMMARY_REQUIRED_KEYS,
    CURRENT_CHECKPOINT_SCHEMA_VERSION,
    load_summary_file,
    parse_summary_text,
)


ROOT = Path(__file__).resolve().parents[2]
RESOLVED_EXTERNAL_BOUNDARY_SUMMARY_KEYS = (
    "external_inflow_map",
    "external_ordinary_open_model",
    "external_interface_transport",
    "outer_particle_mode_resolved",
)


def test_summary_parser_rejects_duplicate_normalized_keys() -> None:
    with pytest.raises(ValueError, match="duplicate key 'batches'"):
        parse_summary_text("batches=1\n batches =2\n")


def test_summary_parser_enforces_required_keys(tmp_path: Path) -> None:
    path = tmp_path / "summary.txt"
    path.write_text("mesh_nelem=1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="processed_particles"):
        load_summary_file(path, required_keys=CORE_SUMMARY_REQUIRED_KEYS)


def test_fortran_writer_matches_python_core_summary_contract() -> None:
    writer = (ROOT / "src/runtime/bem_output_writer.f90").read_text(encoding="utf-8")
    checkpoint_contract = (ROOT / "src/runtime/bem_checkpoint_contract.f90").read_text(
        encoding="utf-8"
    )

    for key in CORE_SUMMARY_REQUIRED_KEYS:
        assert f"'{key}='" in writer
    match = re.search(
        r"checkpoint_schema_version_current\s*=\s*(\d+)_i32",
        checkpoint_contract,
    )
    assert match is not None
    assert int(match.group(1)) == CURRENT_CHECKPOINT_SCHEMA_VERSION


def test_fortran_writer_emits_resolved_external_boundary_receipt() -> None:
    writer = (ROOT / "src/runtime/bem_output_writer.f90").read_text(encoding="utf-8")

    for key in RESOLVED_EXTERNAL_BOUNDARY_SUMMARY_KEYS:
        assert f"'{key}='" in writer
