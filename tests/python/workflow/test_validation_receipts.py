"""Execution receipts and submission provenance of the periodic validation tool."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from validation_support import (
    _TEST_CASE_PRODUCER_ROLES,
    _TEST_SUBMITTED_JOBS,
    _load_tool,
    _stage_clean_validation,
    _write_run_output,
    _write_submission_provenance,
)

pytestmark = pytest.mark.workflow


def _write_partial_producer_provenance(
    validation_root: Path,
    manifest: dict[str, object],
    case_name: str,
    *,
    job_id: str | None = None,
    config_line: str | None = None,
) -> tuple[str, str, Path]:
    role = _TEST_CASE_PRODUCER_ROLES[case_name]
    expected_job_id, job_name, _resources = _TEST_SUBMITTED_JOBS[role]
    selected_job_id = expected_job_id if job_id is None else job_id
    journal = {
        "submission_complete": False,
        **{name: "" for name in _TEST_SUBMITTED_JOBS},
    }
    journal[role] = selected_job_id
    (validation_root / "submit/job_ids.json").write_text(
        json.dumps(journal), encoding="utf-8"
    )
    cases = manifest["cases"]
    assert isinstance(cases, dict)
    config_path = str(cases[case_name]["config_path"])
    line = f"{config_path}: OK" if config_line is None else config_line
    hash_path = (
        validation_root / "provenance/hashes" / f"{selected_job_id}.{job_name}.sha256"
    )
    hash_path.parent.mkdir(parents=True, exist_ok=True)
    hash_path.write_text(line + "\n", encoding="utf-8")
    return role, selected_job_id, hash_path


def test_clean_first_receipt_accepts_bound_partial_producer_journal(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool, archive_run, validation_root, binary, monkeypatch
    )
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    role, job_id, _hash_path = _write_partial_producer_provenance(
        validation_root, manifest, "smoke_finite_configured"
    )
    monkeypatch.setenv("SLURM_JOB_ID", job_id)

    report = tool.verify_run(
        output,
        100,
        producer_job_role=role,
    )

    assert report["execution_producer"] == {
        "job_role": "smoke",
        "job_id": "9001",
        "config_sha256": manifest["cases"]["smoke_finite_configured"]["config_sha256"],
    }


@pytest.mark.parametrize(
    ("producer_role", "current_job", "hash_suffix", "message"),
    [
        (None, "9001", "", "producer job role"),
        ("finite_140000", "9001", "", "producer job role"),
        ("smoke", None, "", "SLURM_JOB_ID"),
        ("smoke", "9999", "", "SLURM_JOB_ID"),
        ("smoke", "9001", ".backup", "config hash log"),
    ],
)
def test_clean_first_receipt_rejects_unbound_producer(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    producer_role: str | None,
    current_job: str | None,
    hash_suffix: str,
    message: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool, archive_run, validation_root, binary, monkeypatch
    )
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    config_path = manifest["cases"]["smoke_finite_configured"]["config_path"]
    _write_partial_producer_provenance(
        validation_root,
        manifest,
        "smoke_finite_configured",
        config_line=f"{config_path}{hash_suffix}: OK",
    )
    if current_job is None:
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    else:
        monkeypatch.setenv("SLURM_JOB_ID", current_job)

    with pytest.raises(tool.ValidationError, match=message):
        tool.verify_run(
            output,
            100,
            producer_job_role=producer_role,
        )


@pytest.mark.parametrize(
    "mutation",
    [
        "nonboolean_complete",
        "extra_key",
        "missing_role_id",
        "numeric_job_id",
        "duplicate_job_id",
    ],
)
def test_clean_first_receipt_rejects_invalid_partial_journal_schema(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool, archive_run, validation_root, binary, monkeypatch
    )
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    role, job_id, _hash_path = _write_partial_producer_provenance(
        validation_root, manifest, "smoke_finite_configured"
    )
    journal_path = validation_root / "submit/job_ids.json"
    journal = json.loads(journal_path.read_text(encoding="utf-8"))
    if mutation == "nonboolean_complete":
        journal["submission_complete"] = 0
    elif mutation == "extra_key":
        journal["unexpected"] = ""
    elif mutation == "missing_role_id":
        journal[role] = ""
    elif mutation == "numeric_job_id":
        journal[role] = int(job_id)
    else:
        journal["finite_140000"] = job_id
    journal_path.write_text(json.dumps(journal), encoding="utf-8")
    monkeypatch.setenv("SLURM_JOB_ID", job_id)

    with pytest.raises(tool.ValidationError, match="journal|job ID"):
        tool.verify_run(output, 100, producer_job_role=role)


def test_clean_existing_receipt_requires_complete_journal_but_ignores_child_job_id(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool, archive_run, validation_root, binary, monkeypatch
    )
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    role, job_id, _hash_path = _write_partial_producer_provenance(
        validation_root, manifest, "smoke_finite_configured"
    )
    monkeypatch.setenv("SLURM_JOB_ID", job_id)
    tool.verify_run(output, 100, producer_job_role=role)
    receipt = validation_root / "provenance/verified/smoke_finite_configured.json"
    original = receipt.read_bytes()

    monkeypatch.setenv("SLURM_JOB_ID", "9006")
    with pytest.raises(tool.ValidationError, match="journal is incomplete"):
        tool.verify_run(output, 100, require_existing_receipt=True)

    _write_submission_provenance(validation_root, manifest)
    report = tool.verify_run(output, 100, require_existing_receipt=True)
    assert report["execution_producer"]["job_id"] == job_id
    assert receipt.read_bytes() == original


@pytest.mark.parametrize(
    "mutation",
    ["journal_job_id", "config_hash_line", "receipt_binding"],
)
def test_clean_existing_receipt_rechecks_immutable_producer_binding(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool, archive_run, validation_root, binary, monkeypatch
    )
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    job_ids = _write_submission_provenance(validation_root, manifest)
    monkeypatch.setenv("SLURM_JOB_ID", job_ids["smoke"])
    tool.verify_run(output, 100, producer_job_role="smoke")
    monkeypatch.setenv("SLURM_JOB_ID", job_ids["analysis"])

    if mutation == "journal_job_id":
        journal_path = validation_root / "submit/job_ids.json"
        journal = json.loads(journal_path.read_text(encoding="utf-8"))
        journal["smoke"] = "9011"
        journal_path.write_text(json.dumps(journal), encoding="utf-8")
    elif mutation == "config_hash_line":
        hash_path = (
            validation_root / "provenance/hashes/9001.beach-periodic-smoke.sha256"
        )
        config_path = manifest["cases"]["smoke_finite_configured"]["config_path"]
        hash_path.write_text(
            hash_path.read_text(encoding="utf-8").replace(
                f"{config_path}: OK", f"{config_path}.backup: OK"
            ),
            encoding="utf-8",
        )
    else:
        receipt_path = (
            validation_root / "provenance/verified/smoke_finite_configured.json"
        )
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        del receipt["execution_producer"]
        receipt["receipt_payload_sha256"] = tool._receipt_payload_sha256(receipt)
        receipt_path.chmod(0o644)
        receipt_path.write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        receipt_path.chmod(0o444)

    with pytest.raises(tool.ValidationError, match="producer|config hash log"):
        tool.verify_run(output, 100, require_existing_receipt=True)


def test_clean_canonical_case_cannot_publish_production_receipt(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = _stage_clean_validation(
        tool, archive_run, validation_root, binary, monkeypatch
    )
    output = _write_run_output(
        validation_root,
        "full_finite_configured",
        batches=280000,
        cache_hit=None,
        build_count=None,
    )
    job_ids = _write_submission_provenance(validation_root, manifest)
    monkeypatch.setenv("SLURM_JOB_ID", job_ids["finite_280000"])

    with pytest.raises(tool.ValidationError, match="production job graph"):
        tool.verify_run(
            output,
            280000,
            producer_job_role="finite_280000",
        )


def test_nonstrict_receipt_remains_backward_compatible_with_generated_role(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )

    report = tool.verify_run(output, 100, producer_job_role="smoke")

    assert "execution_producer" not in report


@pytest.mark.parametrize(
    "foreign_module",
    ["SysB/2022", "SysC/2022", "SysCL/2022", "SysG/2022"],
)
def test_submission_provenance_rejects_coexisting_system_module(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    foreign_module: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    _write_submission_provenance(validation_root, manifest)
    module_path = (
        validation_root / "provenance/modules/9006.beach-periodic-analysis.txt"
    )
    module_path.write_text(
        module_path.read_text(encoding="utf-8") + foreign_module + "\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("SLURM_JOB_ID", "9006")

    with pytest.raises(tool.ValidationError, match="SysA/2022 only"):
        tool._verify_submission_provenance(validation_root, manifest)


@pytest.mark.parametrize(
    "missing_module",
    ["intel/2023.2", "intelmpi/2023.2"],
)
def test_submission_provenance_requires_exact_intel_modules(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    missing_module: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    _write_submission_provenance(validation_root, manifest)
    module_path = (
        validation_root / "provenance/modules/9006.beach-periodic-analysis.txt"
    )
    module_text = module_path.read_text(encoding="utf-8")
    module_path.write_text(
        module_text.replace(f"{missing_module}\n", ""),
        encoding="utf-8",
    )
    monkeypatch.setenv("SLURM_JOB_ID", "9006")

    with pytest.raises(
        tool.ValidationError,
        match="intel/2023.2 and intelmpi/2023.2",
    ):
        tool._verify_submission_provenance(validation_root, manifest)


def test_submission_provenance_requires_exact_role_config_hash_lines(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    job_ids = _write_submission_provenance(validation_root, manifest)
    hash_path = validation_root / "provenance/hashes/9001.beach-periodic-smoke.sha256"
    config_path = manifest["cases"]["smoke_finite_configured"]["config_path"]
    hash_path.write_text(
        hash_path.read_text(encoding="utf-8").replace(
            f"{config_path}: OK", f"{config_path}.backup: OK"
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("SLURM_JOB_ID", job_ids["analysis"])

    with pytest.raises(tool.ValidationError, match="hash log is incomplete"):
        tool._verify_submission_provenance(validation_root, manifest)


def test_submission_provenance_rejects_numeric_job_id_value(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    job_ids = _write_submission_provenance(validation_root, manifest)
    journal_path = validation_root / "submit/job_ids.json"
    journal = json.loads(journal_path.read_text(encoding="utf-8"))
    journal["smoke"] = 9001
    journal_path.write_text(json.dumps(journal), encoding="utf-8")
    monkeypatch.setenv("SLURM_JOB_ID", job_ids["analysis"])

    with pytest.raises(tool.ValidationError, match="job IDs"):
        tool._verify_submission_provenance(validation_root, manifest)
