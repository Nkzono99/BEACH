"""Clean-source and build-origin provenance of the periodic validation tool."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest
from validation_support import (
    ROOT,
    _canonical_build_info,
    _load_tool,
    _stage_clean_validation,
)

pytestmark = pytest.mark.workflow


def test_stage_can_require_clean_source(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    monkeypatch.setattr(
        tool,
        "_git_provenance",
        lambda: {
            "commit": "deadbeef",
            "describe": "deadbeef-dirty",
            "dirty": True,
            "status_porcelain": [" M src/example.f90"],
        },
    )

    with pytest.raises(tool.ValidationError, match="clean source"):
        tool.stage_validation(
            archive_run,
            tmp_path / "validation",
            binary,
            require_clean_source=True,
        )


def test_clean_production_stage_requires_analysis_library(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    monkeypatch.setattr(
        tool,
        "_git_provenance",
        lambda: {
            "commit": "deadbeef",
            "describe": "deadbeef",
            "dirty": False,
            "status_porcelain": [],
        },
    )

    with pytest.raises(tool.ValidationError, match="analysis library"):
        tool.stage_validation(
            archive_run,
            tmp_path / "validation",
            binary,
            require_clean_source=True,
        )


@pytest.mark.parametrize(
    "relative",
    [
        "input/release_kernel_base.toml",
        "analysis/local_release/release_model_summary.json",
    ],
)
def test_clean_production_stage_requires_explicit_release_mechanics(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    relative: str,
) -> None:
    tool = _load_tool()
    (archive_run / relative).unlink()

    with pytest.raises(tool.ValidationError, match="production mechanics.*missing"):
        _stage_clean_validation(
            tool,
            archive_run,
            tmp_path / "validation",
            binary,
            monkeypatch,
        )


@pytest.mark.parametrize("mutation", ["summary_nonfinite", "adhesion_schema"])
def test_strict_analysis_revalidates_release_mechanics_schema_and_finiteness(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool,
        archive_run,
        validation_root,
        binary,
        monkeypatch,
    )
    if mutation == "summary_nonfinite":
        relative = "analysis/local_release/release_model_summary.json"
        path = archive_run / relative
        value = json.loads(path.read_text(encoding="utf-8"))
        value["assumptions"]["dust_density_kg_m3"] = float("nan")
        path.write_text(json.dumps(value), encoding="utf-8")
    else:
        relative = "input/release_kernel_base.toml"
        path = archive_run / relative
        path.write_text(
            '[adhesion]\nmodel = "vdw_work"\nhamaker_constant = 1.0e-19\n',
            encoding="utf-8",
        )
    manifest["archive"]["analysis_inputs"][relative]["sha256"] = tool._sha256(path)
    (validation_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="production mechanics"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=Path(manifest["analysis_library"]["staged_path"]),
            require_complete=True,
        )


@pytest.mark.parametrize(
    "mutation",
    ["binary_commit", "library_dirty", "library_version"],
)
def test_clean_production_stage_rejects_artifact_build_origin_mismatch(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    tool = _load_tool()
    commit = "a" * 40
    binary_info = _canonical_build_info(commit)
    library_info = _canonical_build_info(commit)
    if mutation == "binary_commit":
        binary_info["build_source_commit"] = "b" * 40
        binary_info["build_id"] = f"{'b' * 40}:clean"
    elif mutation == "library_dirty":
        library_info["build_id"] = f"{commit}:dirty"
    else:
        library_info["build_version"] = "1.5.0-stale"
    monkeypatch.setattr(
        tool,
        "_git_provenance",
        lambda: {
            "commit": commit,
            "describe": commit,
            "dirty": False,
            "status_porcelain": [],
        },
    )
    monkeypatch.setattr(
        tool,
        "_binary_build_info",
        lambda _path: binary_info,
        raising=False,
    )
    monkeypatch.setattr(
        tool,
        "_library_build_info",
        lambda _path: library_info,
        raising=False,
    )

    with pytest.raises(tool.ValidationError, match="build origin"):
        tool.stage_validation(
            archive_run,
            tmp_path / "validation",
            binary,
            library=binary,
            require_clean_source=True,
        )


def test_clean_production_stage_and_verify_inputs_attest_artifact_build_origin(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    commit = "a" * 40
    current = _canonical_build_info(commit)
    monkeypatch.setattr(
        tool,
        "_git_provenance",
        lambda: {
            "commit": commit,
            "describe": commit,
            "dirty": False,
            "status_porcelain": [],
        },
    )
    monkeypatch.setattr(
        tool,
        "_binary_build_info",
        lambda _path: dict(current),
        raising=False,
    )
    monkeypatch.setattr(
        tool,
        "_library_build_info",
        lambda _path: dict(current),
        raising=False,
    )
    validation_root = tmp_path / "validation"

    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
        require_clean_source=True,
    )

    assert manifest["build_origin"] == current
    assert tool.verify_inputs(validation_root)["status"] == "ok"
    current["build_id"] = f"{commit}:dirty"
    with pytest.raises(tool.ValidationError, match="build origin"):
        tool.verify_inputs(validation_root)


def test_strict_analysis_requires_attested_build_origin(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )

    with pytest.raises(tool.ValidationError, match="clean build origin"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=Path(manifest["analysis_library"]["staged_path"]),
            require_complete=True,
        )


def test_strict_analysis_rejects_build_origin_policy_downgrade(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool,
        archive_run,
        validation_root,
        binary,
        monkeypatch,
    )
    manifest["execution"]["require_clean_source"] = False
    manifest["build_origin"] = None
    (validation_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="clean build origin"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=Path(manifest["analysis_library"]["staged_path"]),
            require_complete=True,
        )


@pytest.mark.skipif(
    not os.environ.get("BEACH_VALIDATION_BINARY"),
    reason="native build-origin integration requires BEACH_VALIDATION_BINARY",
)
def test_native_binary_and_library_match_clean_current_source(
    field_kernel_library: Path,
) -> None:
    tool = _load_tool()
    binary_value = os.environ.get("BEACH_VALIDATION_BINARY")
    assert binary_value is not None
    binary_path = Path(binary_value).expanduser().resolve()
    library_path = field_kernel_library

    assert binary_path.is_file()
    source_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()

    build_origin = tool._validate_production_build_origin(
        binary_path,
        library_path,
        source_commit,
    )
    assert build_origin["build_source_commit"] == source_commit
    assert build_origin["build_id"] == f"{source_commit}:clean"


def _probe_kernel(electric_field, potential, calls: list[Path] | None = None):
    """A kernel double that returns the given E and phi for the probe's one target."""

    class FakeKernel:
        def __init__(self, triangles, charges, *, library_path):
            assert np.asarray(triangles).shape == (1, 3, 3)
            assert np.asarray(charges).shape == (1,)
            if calls is not None:
                calls.append(Path(library_path))

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def eval_e(self, points):
            assert np.asarray(points).shape == (1, 3)
            return np.asarray(electric_field, dtype=float)

        def eval_phi(self, points):
            assert np.asarray(points).shape == (1, 3)
            return np.asarray(potential, dtype=float)

    return FakeKernel


def test_probe_library_executes_basic_field_and_potential(
    binary: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    calls: list[Path] = []
    kernel = _probe_kernel([[1.5, -2.5, 3.5]], [7.25], calls)
    monkeypatch.setattr(tool, "FieldKernel", kernel, raising=False)

    report = tool.probe_library(binary)

    assert report["status"] == "ok"
    assert report["library_sha256"] == tool._sha256(binary)
    assert calls == [binary.resolve()]
    assert report["electric_field_V_m"] == [1.5, -2.5, 3.5]
    assert report["potential_V"] == 7.25


@pytest.mark.parametrize(
    ("electric_field", "potential"),
    [
        pytest.param([[0.0, 0.0, 0.0]], [7.25], id="zero-field"),
        pytest.param([[1.5, -2.5, 3.5]], [0.0], id="zero-potential"),
        pytest.param([[np.nan, 0.0, 1.0]], [7.25], id="nan-field"),
        pytest.param([[1.5, -2.5, 3.5]], [np.inf], id="inf-potential"),
        pytest.param([1.5, -2.5, 3.5], [7.25], id="field-shape"),
        pytest.param([[1.5, -2.5, 3.5]], [7.25, 7.25], id="potential-shape"),
    ],
)
def test_probe_library_rejects_invalid_field_or_potential(
    binary: Path,
    monkeypatch: pytest.MonkeyPatch,
    electric_field,
    potential,
) -> None:
    tool = _load_tool()
    monkeypatch.setattr(
        tool, "FieldKernel", _probe_kernel(electric_field, potential), raising=False
    )

    with pytest.raises(tool.ValidationError, match="returned invalid values"):
        tool.probe_library(binary)


def test_verify_inputs_hashes_archive_analysis_assumptions(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    release = archive_run / "input/release_kernel_base.toml"
    release.write_text('[adhesion]\nmodel = "none"\n', encoding="utf-8")
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    release.write_text('[adhesion]\nmodel = "vdw_work"\n', encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="analysis input hash mismatch"):
        tool.verify_inputs(validation_root)


def test_verify_inputs_hashes_archive_analysis_outputs(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    output = archive_run / "work/latest"
    (output / "summary.txt").write_text(
        "mesh_nelem=1\nmesh_count=1\n", encoding="utf-8"
    )
    (output / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-15\n", encoding="utf-8"
    )
    (output / "mesh_triangles.csv").write_text(
        "elem_idx,v0_x_m\n1,0.0\n", encoding="utf-8"
    )
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(archive_run, validation_root, binary)
    assert "work/latest/charges.csv" in manifest["archive"]["analysis_outputs"]
    (output / "charges.csv").write_text(
        "elem_idx,charge_C\n1,2.0e-15\n", encoding="utf-8"
    )

    with pytest.raises(tool.ValidationError, match="archive analysis output hash"):
        tool.verify_inputs(validation_root)
