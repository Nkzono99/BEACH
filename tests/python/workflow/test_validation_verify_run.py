"""verify-run checks of producer outputs, caches, restarts, and checkpoints."""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import pytest
from validation_support import (
    ROOT,
    _canonical_build_info,
    _load_tool,
    _write_run_output,
)

pytestmark = pytest.mark.workflow


def _replace_summary_value(output: Path, key: str, value: str | None) -> None:
    summary = output / "summary.txt"
    prefix = f"{key}="
    lines = summary.read_text(encoding="utf-8").splitlines()
    assert sum(line.startswith(prefix) for line in lines) == 1
    rewritten = [
        *(line for line in lines if not line.startswith(prefix)),
        *([] if value is None else [f"{key}={value}"]),
    ]
    summary.write_text("\n".join(rewritten) + "\n", encoding="utf-8")


@pytest.mark.parametrize(
    ("key", "tampered_value"),
    [
        ("field_backend", None),
        ("field_backend", "direct"),
        ("field_normalization", None),
        ("field_normalization", "box"),
        ("field_source_model", None),
        ("field_source_model", "point"),
        ("field_kernel_id", None),
        ("field_kernel_id", "softened_point"),
    ],
)
def test_verify_run_rejects_nonproduction_field_execution_contract(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    key: str,
    tampered_value: str | None,
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
    _replace_summary_value(output, key, tampered_value)

    with pytest.raises(tool.ValidationError, match=key):
        tool.verify_run(output, 100)


def test_verify_run_rejects_duplicate_summary_keys(
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
    summary = output / "summary.txt"
    summary.write_text(
        summary.read_text(encoding="utf-8") + "mesh_nelem=1\n",
        encoding="utf-8",
    )

    with pytest.raises(
        tool.ValidationError,
        match=r"duplicate.*summary|summary.*duplicate",
    ):
        tool.verify_run(output, 100)


def test_verify_run_rejects_summary_build_origin_mismatch(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    commit = "a" * 40
    build_origin = _canonical_build_info(commit)
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["build_origin"] = build_origin
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        tool,
        "_binary_build_info",
        lambda _path: dict(build_origin),
        raising=False,
    )
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    _replace_summary_value(output, "build_id", f"{commit}:dirty")

    with pytest.raises(tool.ValidationError, match="build origin"):
        tool.verify_run(output, 100)


def test_verify_run_rejects_per_mesh_element_count_mismatch(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    (archive_run / "work/latest/mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,plane,insulator,1.0,2\n"
        "2,template,sphere,insulator,1.0,1\n",
        encoding="utf-8",
    )
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
        mesh_sources=[
            {
                "mesh_id": 1,
                "source_kind": "template",
                "template_kind": "plane",
                "surface_model": "insulator",
                "epsilon_r": 1.0,
                "elem_count": 2,
            },
            {
                "mesh_id": 2,
                "source_kind": "template",
                "template_kind": "sphere",
                "surface_model": "insulator",
                "epsilon_r": 1.0,
                "elem_count": 1,
            },
        ],
    )
    triangles = output / "mesh_triangles.csv"
    rows = triangles.read_text(encoding="utf-8").splitlines()
    prefix, mesh_id = rows[2].rsplit(",", maxsplit=1)
    assert mesh_id == "1"
    rows[2] = f"{prefix},2"
    triangles.write_text("\n".join(rows) + "\n", encoding="utf-8")

    with pytest.raises(
        tool.ValidationError,
        match=r"mesh_sources\.csv.*mesh_triangles\.csv|per-mesh.*count",
    ):
        tool.verify_run(output, 100)


@pytest.mark.parametrize(
    "mutation",
    [
        "row_order",
        "per_mesh_element_count",
        "source_kind",
        "template_kind",
        "surface_model",
        "epsilon_r",
    ],
)
def test_verify_run_rejects_archive_mesh_source_semantic_mismatch(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    archive_sources = [
        {
            "mesh_id": 1,
            "source_kind": "template",
            "template_kind": "plane",
            "surface_model": "insulator",
            "epsilon_r": 1.0,
            "elem_count": 2,
        },
        {
            "mesh_id": 2,
            "source_kind": "template",
            "template_kind": "sphere",
            "surface_model": "insulator",
            "epsilon_r": 1.0,
            "elem_count": 1,
        },
    ]
    archive_source_path = archive_run / "work/latest/mesh_sources.csv"
    archive_source_path.write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,plane,insulator,1.0,2\n"
        "2,template,sphere,insulator,1.0,1\n",
        encoding="utf-8",
    )
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    actual_sources = [dict(source) for source in archive_sources]
    if mutation == "row_order":
        actual_sources.reverse()
    elif mutation == "per_mesh_element_count":
        actual_sources[0]["elem_count"] = 1
        actual_sources[1]["elem_count"] = 2
    elif mutation == "source_kind":
        actual_sources[1]["source_kind"] = "file"
    elif mutation == "template_kind":
        actual_sources[1]["template_kind"] = "plane"
    elif mutation == "surface_model":
        actual_sources[1]["surface_model"] = "conductor"
    elif mutation == "epsilon_r":
        actual_sources[1]["epsilon_r"] = 2.5
    else:  # pragma: no cover - parameter list is exhaustive
        raise AssertionError(mutation)
    output = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
        mesh_sources=actual_sources,
    )

    with pytest.raises(tool.ValidationError, match="archived mesh source contract"):
        tool.verify_run(output, 100)


def test_verify_run_gates_cache_prime_warm_reuse_and_finite_case(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    prime = _write_run_output(
        validation_root,
        "cache_prime",
        batches=1,
        cache_hit=False,
        build_count=1,
    )
    warm = _write_run_output(
        validation_root,
        "smoke_infinite_physical",
        batches=100,
        cache_hit=True,
        build_count=0,
    )
    finite = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )

    prime_report = tool.verify_run(prime, 1)
    warm_report = tool.verify_run(warm, 100)
    assert prime_report["status"] == "ok"
    assert warm_report["cache"]["hit"] is True
    assert warm_report["cache"]["path_sha256"] == prime_report["cache"]["path_sha256"]
    receipt = validation_root / "provenance/verified/smoke_infinite_physical.json"
    assert json.loads(receipt.read_text(encoding="utf-8"))["case"] == (
        "smoke_infinite_physical"
    )
    assert tool.verify_run(finite, 100)["cache"]["expectation"] == "none"


def test_verify_run_keeps_execution_receipt_immutable_and_detects_output_change(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    out = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    timestamps = iter(("execution-time", "revalidation-time", "changed-time"))
    monkeypatch.setattr(tool, "_utc_now", lambda: next(timestamps))

    tool.verify_run(out, 100)
    receipt = validation_root / "provenance/verified/smoke_finite_configured.json"
    original = receipt.read_bytes()
    tool.verify_run(out, 100)
    assert receipt.read_bytes() == original

    charges = out / "charges.csv"
    charges.write_text(
        charges.read_text(encoding="utf-8").replace("1.0e-15", "2.0e-15"),
        encoding="utf-8",
    )
    with pytest.raises(tool.ValidationError, match="immutable execution receipt"):
        tool.verify_run(out, 100)
    assert receipt.read_bytes() == original


def test_verify_run_requires_prime_summary_and_cached_operator(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    warm = _write_run_output(
        validation_root,
        "smoke_infinite_physical",
        batches=100,
        cache_hit=True,
        build_count=0,
    )

    with pytest.raises(tool.ValidationError, match="cache prime summary"):
        tool.verify_run(warm, 100)

    prime = _write_run_output(
        validation_root,
        "cache_prime",
        batches=1,
        cache_hit=False,
        build_count=1,
    )
    (validation_root / "cache/periodic2/operator.bin").unlink()
    with pytest.raises(tool.ValidationError, match="cached operator"):
        tool.verify_run(prime, 1)


def test_restart_segments_verify_and_history_combines_without_duplicates(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    first = _write_run_output(
        validation_root,
        "full_finite_configured_140000",
        batches=140000,
        cache_hit=None,
        build_count=None,
    )
    second = _write_run_output(
        validation_root,
        "full_finite_configured_280000",
        batches=280000,
        cache_hit=None,
        build_count=None,
    )

    assert tool.verify_run(first, 140000)["status"] == "ok"
    assert tool.verify_run(second, 280000)["status"] == "ok"
    rows = tool._history_rows(
        "new_finite_configured",
        [first / "charge_history.csv", second / "charge_history.csv"],
    )
    assert [int(row["batch"]) for row in rows] == list(range(1, 280001, 1000))

    summary = (first / "summary.txt").read_text(encoding="utf-8")
    (first / "summary.txt").write_text(
        summary.replace(
            "mesh_fingerprint=1111111111111111",
            "mesh_fingerprint=FFFFFFFFFFFFFFFF",
        ),
        encoding="utf-8",
    )
    with pytest.raises(tool.ValidationError, match="restart parent"):
        tool.verify_run(second, 280000)


def test_history_total_charge_uses_compensated_sum(tmp_path: Path) -> None:
    tool = _load_tool()
    history = tmp_path / "segment" / "charge_history.csv"
    history.parent.mkdir()
    history.write_text(
        "batch,processed_particles,rel_change,charge_C\n"
        "1,10,0.0,1e-15\n"
        "1,10,0.0,1e-32\n"
        "1,10,0.0,-1e-15\n",
        encoding="utf-8",
    )

    rows = tool._history_rows("case", [history])

    # Naive left-to-right addition loses 1e-32 against 1e-15 and returns 0.0.
    assert [row["total_charge_C"] for row in rows] == [1e-32]


def test_restart_child_revalidation_detects_parent_output_change(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    first = _write_run_output(
        validation_root,
        "full_finite_configured_140000",
        batches=140000,
        cache_hit=None,
        build_count=None,
    )
    second = _write_run_output(
        validation_root,
        "full_finite_configured_280000",
        batches=280000,
        cache_hit=None,
        build_count=None,
    )
    tool.verify_run(first, 140000)
    tool.verify_run(second, 280000)
    charges = first / "charges.csv"
    charges.write_text(
        charges.read_text(encoding="utf-8").replace("1.0e-15", "2.0e-15"),
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="restart parent.*receipt"):
        tool.verify_run(second, 280000)


@pytest.mark.parametrize(
    "mutation",
    ["cache_expectation", "processed_particles", "outcome_counts"],
)
def test_verify_run_rejects_cache_and_particle_ledger_mismatch(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    if mutation == "cache_expectation":
        output = _write_run_output(
            validation_root,
            "cache_prime",
            batches=1,
            cache_hit=True,
            build_count=0,
        )
        expected_batches = 1
        message = "cache miss"
    else:
        output = _write_run_output(
            validation_root,
            "smoke_finite_configured",
            batches=100,
            cache_hit=None,
            build_count=None,
        )
        expected_batches = 100
        if mutation == "processed_particles":
            summary_path = output / "summary.txt"
            summary = summary_path.read_text(encoding="utf-8")
            summary_path.write_text(
                summary.replace("processed_particles=10", "processed_particles=11"),
                encoding="utf-8",
            )
            message = "processed_particles"
        else:
            ledger_path = output / "charge_ledger.csv"
            ledger = ledger_path.read_text(encoding="utf-8")
            ledger_path.write_text(
                ledger.replace(",6,3,1\n", ",6,4,0\n"),
                encoding="utf-8",
            )
            message = "outcome counts"

    with pytest.raises(tool.ValidationError, match=message):
        tool.verify_run(output, expected_batches)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("nonfinite_mesh", "mesh_triangles.csv"),
        ("history_gap", "history batch sequence"),
    ],
)
def test_verify_run_rejects_nonfinite_mesh_and_history_gap(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
    message: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    out = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    if mutation == "nonfinite_mesh":
        triangles_path = out / "mesh_triangles.csv"
        triangles = triangles_path.read_text(encoding="utf-8")
        triangles_path.write_text(
            triangles.replace("1,0,0,0", "1,nan,0,0"), encoding="utf-8"
        )
    else:
        history_path = out / "charge_history.csv"
        history = history_path.read_text(encoding="utf-8").splitlines()
        history_path.write_text(
            "\n".join([history[0], *history[2:]]) + "\n", encoding="utf-8"
        )

    with pytest.raises(tool.ValidationError, match=message):
        tool.verify_run(out, 100)


@pytest.mark.parametrize("missing_name", ["mesh_potential.csv", "mesh_sources.csv"])
def test_verify_run_requires_configured_mesh_outputs(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    missing_name: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    out = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    (out / missing_name).unlink()

    with pytest.raises(tool.ValidationError, match=re.escape(missing_name)):
        tool.verify_run(out, 100)


def _legacy_schema3_checkpoint(out: Path) -> None:
    _replace_summary_value(out, "checkpoint_schema_version", "3")
    for rank in range(6):
        (out / f"rng_state_rank{rank:05d}.txt").write_text("1\n12345\n")


def _rewrite_manifest(out: Path, key: str, value: str) -> None:
    manifest = out / "checkpoint_complete.txt"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    manifest.write_text(
        "\n".join(
            f"{key}={value}" if line.startswith(f"{key}=") else line for line in lines
        )
        + "\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda out: (out / "checkpoint_complete.txt").unlink(), "checkpoint_complete"),
        (lambda out: (out / "macro_residuals.csv").unlink(), "macro_residuals"),
        (lambda out: _rewrite_manifest(out, "state", "in_progress"), "state"),
        (lambda out: _rewrite_manifest(out, "batches", "0"), "batches"),
        (lambda out: _rewrite_manifest(out, "mpi_world_size", "1"), "mpi_world_size"),
        (
            lambda out: _rewrite_manifest(out, "macro_residuals_present", "F"),
            "macro_residuals_present",
        ),
        (_legacy_schema3_checkpoint, "current runtime schema 11"),
    ],
    ids=[
        "no_manifest",
        "no_residuals",
        "in_progress",
        "batch_mismatch",
        "rank_mismatch",
        "residuals_undeclared",
        "legacy_schema3",
    ],
)
def test_verify_run_requires_complete_schema11_checkpoint(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation,
    message: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    out = _write_run_output(
        validation_root,
        "cache_prime",
        batches=1,
        cache_hit=False,
        build_count=1,
    )
    mutation(out)

    with pytest.raises(tool.ValidationError, match=message):
        tool.verify_run(out, 1)


def test_native_writer_checkpoint_satisfies_validator(
    native_checker: Path, tmp_path: Path
) -> None:
    """Producer/consumer contract: the checkpoint the real Fortran writer publishes
    passes the validator's checkpoint check, and losing its manifest does not."""
    import tomllib

    import tomli_w

    config = tomllib.loads((ROOT / "examples/panel_direct.toml").read_text())
    out = tmp_path / "out"
    config["output"] = {"enabled": True, "dir": str(out)}
    path = tmp_path / "beach.toml"
    path.write_text(tomli_w.dumps(config))
    result = subprocess.run(
        [str(native_checker), str(path)],
        cwd=tmp_path,
        env=os.environ | {"OMP_NUM_THREADS": "1"},
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    tool = _load_tool()
    summary = tool._summary(out / "summary.txt")
    tool._validate_checkpoint(out, summary)
    assert not list(out.glob("rng_state*.txt"))
    (out / "checkpoint_complete.txt").unlink()
    with pytest.raises(tool.ValidationError, match="checkpoint_complete"):
        tool._validate_checkpoint(out, summary)


def test_verify_run_requires_complete_cached_electrostatic_state(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    out = _write_run_output(
        validation_root,
        "cache_prime",
        batches=1,
        cache_hit=False,
        build_count=1,
    )
    summary = (out / "summary.txt").read_text(encoding="utf-8")
    (out / "summary.txt").write_text(
        "\n".join(
            line
            for line in summary.splitlines()
            if not line.startswith("gauss_residual_C=")
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="cached electrostatic state"):
        tool.verify_run(out, 1)
