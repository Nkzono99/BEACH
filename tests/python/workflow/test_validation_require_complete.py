"""analyze --require-complete over every staged case."""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from validation_support import (
    _TEST_CASE_PRODUCER_ROLES,
    _load_tool,
    _stage_clean_validation,
    _write_periodic_oracle_receipt,
    _write_run_output,
    _write_submission_provenance,
)

from beach import RunDiagnostics, RunField, RunHistory, RunMesh, RunReceipts, RunSummary

pytestmark = pytest.mark.workflow


def test_analyze_requires_existing_library(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)

    with pytest.raises(tool.ValidationError, match="library does not exist"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=tmp_path / "missing.so",
        )


def test_analyze_require_complete_fails_closed_for_missing_full_runs(
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

    with pytest.raises(tool.ValidationError, match="complete analysis requires"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=Path(manifest["analysis_library"]["staged_path"]),
            require_complete=True,
        )


def test_analyze_require_complete_requires_periodic_oracle_receipt(
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
    monkeypatch.setattr(
        tool,
        "_verify_submission_provenance",
        lambda *_args, **_kwargs: {
            "status": "verified",
            "job_ids": {"analysis": "9006"},
        },
    )
    monkeypatch.setattr(
        tool,
        "_verify_complete_runs",
        lambda *_args, **_kwargs: {},
    )

    with pytest.raises(
        tool.ValidationError,
        match="periodic plane-oracle receipt is missing",
    ):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=Path(manifest["analysis_library"]["staged_path"]),
            require_complete=True,
        )


def test_analyze_require_complete_rechecks_staged_archive_input(
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
    archive_input = archive_run / "input/beach.toml"
    archive_input.write_text(
        archive_input.read_text(encoding="utf-8") + "\n# changed after staging\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="archived input hash mismatch"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=Path(manifest["analysis_library"]["staged_path"]),
            require_complete=True,
        )


def test_analyze_require_complete_rejects_unstaged_library(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    _stage_clean_validation(
        tool,
        archive_run,
        validation_root,
        binary,
        monkeypatch,
    )
    other = tmp_path / "other.so"
    other.write_bytes(binary.read_bytes())

    with pytest.raises(tool.ValidationError, match="staged analysis library"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=other,
            require_complete=True,
        )


def test_analyze_require_complete_verifies_all_cases_and_physics(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    archive_output = archive_run / "work/latest"
    (archive_output / "summary.txt").write_text(
        "mesh_nelem=1\nmesh_count=1\nprocessed_particles=10\nabsorbed=6\n"
        "escaped=4\nbatches=280000\nlast_rel_change=1.0e-3\n",
        encoding="utf-8",
    )
    (archive_output / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-15\n", encoding="utf-8"
    )
    (archive_output / "mesh_triangles.csv").write_text(
        "elem_idx,v0_x_m,v0_y_m,v0_z_m,v1_x_m,v1_y_m,v1_z_m,"
        "v2_x_m,v2_y_m,v2_z_m,charge_C,mesh_id\n"
        "1,0,0,0,1,0,0,0,1,0,1.0e-15,1\n",
        encoding="utf-8",
    )
    (archive_output / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,sphere,insulator,1.0,1\n",
        encoding="utf-8",
    )
    (archive_output / "mesh_potential.csv").write_text(
        "elem_idx,potential_V\n1,2.0\n", encoding="utf-8"
    )
    (archive_output / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        "149001,10,1.0e-3,1,1.0e-15\n"
        "180001,10,1.0e-3,1,1.0e-15\n"
        "279001,10,1.0e-3,1,1.0e-15\n",
        encoding="utf-8",
    )
    manifest = _stage_clean_validation(
        tool,
        archive_run,
        validation_root,
        binary,
        monkeypatch,
    )
    job_ids = _write_submission_provenance(validation_root, manifest)
    definitions = (
        ("cache_prime", 1, False, 1),
        ("smoke_finite_configured", 100, None, None),
        ("smoke_infinite_physical", 100, True, 0),
        ("full_finite_configured_140000", 140000, None, None),
        ("full_finite_configured_280000", 280000, None, None),
        ("full_infinite_physical_140000", 140000, True, 0),
        ("full_infinite_physical_280000", 280000, True, 0),
    )
    for name, batches, hit, count in definitions:
        output = _write_run_output(
            validation_root,
            name,
            batches=batches,
            cache_hit=hit,
            build_count=count,
        )
        producer_role = _TEST_CASE_PRODUCER_ROLES[name]
        monkeypatch.setenv("SLURM_JOB_ID", job_ids[producer_role])
        tool.verify_run(
            output,
            batches,
            producer_job_role=producer_role,
        )
    staged_library = Path(manifest["analysis_library"]["staged_path"])
    _write_periodic_oracle_receipt(tool, validation_root, staged_library)
    monkeypatch.setenv("SLURM_JOB_ID", job_ids["analysis"])

    class FakeBeach:
        def __init__(self, _output_dir, *, config_path=None):
            assert config_path is not None
            self.result = SimpleNamespace(
                charges=np.array([1.0e-15]),
                ledger=None,
                summary=RunSummary(
                    mesh_nelem=1,
                    processed_particles=10,
                    absorbed=6,
                    escaped=4,
                    batches=280000,
                    escaped_boundary=3,
                    survived_max_step=1,
                    last_rel_change=1.0e-3,
                ),
                mesh=RunMesh(
                    mesh_ids=np.array([1]),
                    triangles=np.array(
                        [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
                    ),
                    mesh_potential_v=None,
                ),
                diagnostics=RunDiagnostics(charge_ledger_residual_c=0.0),
                receipts=RunReceipts(
                    model_fingerprint="0123456789ABCDEF",
                    mesh_fingerprint="1111111111111111",
                    species_fingerprint="2222222222222222",
                ),
                history=RunHistory(),
                field=RunField(),
            )

    monkeypatch.setattr(tool, "Beach", FakeBeach)
    monkeypatch.setattr(
        tool,
        "_evaluate_object_physics",
        lambda **_kwargs: (
            [],
            [],
            [],
            [],
            {
                "status": "available",
                "successful_target_models": 1,
                "failures": [],
            },
        ),
    )
    monkeypatch.setattr(
        tool,
        "_local_model_numerical_qualification",
        lambda *_args, **_kwargs: {
            "status": "qualified",
            "claim_scope": "test",
        },
    )
    monkeypatch.setattr(
        tool,
        "_comparison_artifact_contract",
        lambda *_args, **_kwargs: {"status": "complete"},
    )
    monkeypatch.setattr(
        tool,
        "_legacy_estimator_audit",
        lambda *_args, **_kwargs: (
            [{"status": "computed", "interpretation": "test"}],
            {
                "status": "complete",
                "comparison_row_count": 1,
                "covered_native_keys": [],
                "closeness_is_a_gate": False,
            },
        ),
    )
    report = tool.analyze_validation(
        archive_run,
        validation_root,
        library=staged_library,
        require_complete=True,
    )

    assert report["strict_validation"]["required"] is True
    assert set(report["strict_validation"]["verified_cases"]) == {
        name for name, *_rest in definitions
    }
    assert report["strict_validation"]["periodic_plane_oracles"]["status"] == (
        "qualified"
    )
    assert report["physics_evaluation"]["status"] == "available"
    assert report["legacy_estimator_audit"]["status"] == "complete"
    assert (validation_root / "analysis/legacy_estimator_comparison.csv").is_file()
    with (validation_root / "analysis/comparison_matrix.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        comparison_rows = list(csv.DictReader(stream))
    assert any(row["metric"] == "processed_particles" for row in comparison_rows)
    assert any(row["metric"] == "rel_change" for row in comparison_rows)

    monkeypatch.setattr(
        tool,
        "_evaluate_object_physics",
        lambda **_kwargs: (
            [],
            [],
            [],
            [],
            {"status": "partial", "successful_target_models": 0, "failures": [{}]},
        ),
    )
    with pytest.raises(tool.ValidationError, match="already published"):
        tool.analyze_validation(
            archive_run,
            validation_root,
            library=staged_library,
            require_complete=True,
        )


@pytest.mark.parametrize(
    ("require_complete", "report"),
    [
        pytest.param(
            False,
            {"physics_evaluation": {"status": "partial"}},
            id="partial-physics",
        ),
        pytest.param(
            True,
            {"physics_evaluation": {"status": "not_evaluated"}},
            id="complete-physics-not-evaluated",
        ),
        pytest.param(
            True,
            {
                "physics_evaluation": {"status": "available"},
                "numerical_qualification_for_local_frozen_model": {
                    "status": "not_qualified"
                },
            },
            id="complete-local-model-not-qualified",
        ),
    ],
)
def test_analyze_cli_returns_nonzero_for_incomplete_result(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    require_complete: bool,
    report: dict[str, object],
) -> None:
    tool = _load_tool()
    monkeypatch.setattr(
        tool,
        "analyze_validation",
        lambda *_args, **_kwargs: report,
    )

    arguments = [
        "analyze",
        "--archive-run",
        str(tmp_path / "archive"),
        "--validation-root",
        str(tmp_path / "validation"),
        "--library",
        str(tmp_path / "kernel.so"),
    ]
    if require_complete:
        arguments.append("--require-complete")

    status = tool.main(arguments)

    assert status == 2
