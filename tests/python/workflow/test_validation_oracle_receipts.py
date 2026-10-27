"""Periodic-oracle probe payloads and receipt binding."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from validation_support import (
    _canonical_build_info,
    _fake_periodic_kernel_oracle,
    _fake_periodic_kernel_oracles,
    _load_tool,
    _write_periodic_oracle_receipt,
)

pytestmark = pytest.mark.workflow


def test_probe_periodic_oracles_rejects_stale_point_config(
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
    library = Path(manifest["analysis_library"]["staged_path"])
    stale = validation_root / "provenance/oracles/periodic_plane_point.toml"
    stale.write_text('[sim]\nfield_source_model = "point"\n', encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="must be empty"):
        tool.probe_periodic_oracles(validation_root, library)


@pytest.mark.parametrize("mutation", ["quadrature", "neutrality", "component"])
def test_probe_periodic_oracles_deeply_validates_payload_before_publish(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    oracle = _fake_periodic_kernel_oracle()
    if mutation == "quadrature":
        oracle["uniform_plane"]["quadrature_relative_difference"] = 0.02
    elif mutation == "neutrality":
        oracle["neutral_cosine_plane"]["errors"][1]["charge_neutrality_ratio"] = 1.0e-8
    else:
        row = oracle["uniform_plane"]["wrench_refinement"][0]
        row["total_minus_images_normalized_absolute"] = 1.0e-5
        row["component_consistency_relative_error"] = 1.0e-5

    def fake_run(*, cache_dir: Path, **_kwargs):
        result = copy.deepcopy(oracle)
        identities = result["cache_identities"]
        for group, identity in identities.items():
            cache_path = cache_dir / f"operator-triangle_p0-{group}.bin"
            cache_path.write_bytes(f"cache:triangle_p0:{group}\n".encode())
            identity["path"] = str(cache_path.resolve())
            identity["sha256"] = tool._sha256(cache_path)
        result["cache_diagnostics"] = copy.deepcopy(identities["4"])
        group_for_evaluation = {
            evaluation_label: group
            for group, evaluation_labels in tool.ORACLE_CACHE_EVALUATION_GROUPS.items()
            for evaluation_label in evaluation_labels
        }
        for evaluation in result["cache_evaluations"]:
            identity = identities[group_for_evaluation[evaluation["label"]]]
            evaluation["fingerprint"] = identity["fingerprint"]
            evaluation["path"] = identity["path"]
            evaluation["sha256"] = identity["sha256"]
        return result

    monkeypatch.setattr(tool, "_run_periodic_plane_kernel_oracle", fake_run)
    receipt = validation_root / "provenance/oracles/periodic_plane.json"

    with pytest.raises(
        tool.ValidationError,
        match="threshold|integration|neutrality|triangle_p0",
    ):
        tool.probe_periodic_oracles(validation_root, library)

    assert not receipt.exists()


@pytest.mark.parametrize("shared_operator", [False, True])
def test_periodic_oracle_receipt_is_bound_to_library_and_cache_inventory(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    shared_operator: bool,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    receipt = _write_periodic_oracle_receipt(
        tool, validation_root, library, shared_operator=shared_operator,
    )

    report = tool._verify_periodic_oracle_receipt(validation_root, library)

    assert report["status"] == "qualified"
    assert report["receipt_sha256"] == tool._sha256(receipt)
    assert set(report["kernel_oracles"]) == {"triangle_p0"}
    oracle = report["kernel_oracles"]["triangle_p0"]
    assert oracle["field_source_model"] == "triangle_p0"
    assert oracle["field_kernel_id"] == "triangle_p0_exact_p2m_near"
    assert [row["label"] for row in oracle["cache_evaluations"]] == [
        "uniform_4",
        "cosine_4",
        "cosine_8",
    ]
    identities = oracle["cache_identities"]
    assert oracle["cache_diagnostics"] == identities["4"]
    assert (identities["4"]["fingerprint"] == identities["8"]["fingerprint"]) is shared_operator
    assert (identities["4"]["path"] == identities["8"]["path"]) is shared_operator


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_sample_abs_z",
        "tampered_sample_abs_z",
        "missing_ratio_error",
        "tampered_ratio",
        "ratio_threshold",
        "missing_ratio_tolerance",
        "tampered_ratio_tolerance",
        "charge_neutrality",
    ],
)
def test_periodic_oracle_receipt_rejects_invalid_multiheight_decay_contract(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    cosine = kernel_oracles["triangle_p0"]["neutral_cosine_plane"]
    fine = cosine["errors"][1]
    if mutation == "missing_sample_abs_z":
        del cosine["sample_abs_z_m"]
    elif mutation == "tampered_sample_abs_z":
        cosine["sample_abs_z_m"] = [0.25, 0.55]
    elif mutation == "missing_ratio_error":
        del fine["field_decay_ratio_relative_error"]
    elif mutation == "tampered_ratio":
        fine["field_decay_ratio"] *= 1.01
    elif mutation == "ratio_threshold":
        fine["field_decay_ratio"] = cosine["expected_decay_ratio"] * 1.19
        fine["field_decay_ratio_relative_error"] = 0.19
    elif mutation == "missing_ratio_tolerance":
        del cosine["decay_ratio_relative_tolerance"]
    elif mutation == "tampered_ratio_tolerance":
        cosine["decay_ratio_relative_tolerance"] = 0.19
    else:
        fine["charge_neutrality_ratio"] = 1.0e-8
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(
        tool.ValidationError,
        match="cosine|decay|sample|threshold|contract|neutrality",
    ):
        tool._verify_periodic_oracle_receipt(validation_root, library)


def test_periodic_oracle_receipt_rechecks_library_build_origin(
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
    build_origin = _canonical_build_info("a" * 40)
    manifest["build_origin"] = build_origin
    (validation_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    current = dict(build_origin)
    monkeypatch.setattr(
        tool,
        "_library_build_info",
        lambda _path: dict(current),
        raising=False,
    )
    _write_periodic_oracle_receipt(tool, validation_root, library)

    assert (
        tool._verify_periodic_oracle_receipt(validation_root, library)["status"]
        == "qualified"
    )
    current["build_id"] = f"{'a' * 40}:dirty"
    with pytest.raises(tool.ValidationError, match="build origin"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


def test_periodic_oracle_receipt_rejects_legacy_point_schema(
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
    library = Path(manifest["analysis_library"]["staged_path"])
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        overrides={"oracle_schema_version": 2},
    )

    with pytest.raises(tool.ValidationError, match="no longer matches staging"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    ("field_name", "tampered_value"),
    [
        ("requested_periodic_model", "finite_images"),
        ("configured_far_correction", "none"),
        ("effective_far_correction", "none"),
        ("nonzero_mode_backend", "legacy_finite_images"),
        ("zero_mode_policy", "legacy_not_decomposed"),
        ("lower_boundary_model", "legacy_implicit"),
        ("cache_hit", False),
        ("cache_build_count", 1),
        ("cache_fingerprint", ""),
        ("cache_path", "/tmp/outside-oracle-cache.bin"),
    ],
)
def test_periodic_oracle_receipt_gates_effective_backend_and_cache_diagnostics(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    field_name: str,
    tampered_value: object,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    oracle = kernel_oracles["triangle_p0"]
    if field_name in oracle["effective_field_contract"]:
        oracle["effective_field_contract"][field_name] = tampered_value
    else:
        key = field_name.removeprefix("cache_")
        oracle["cache_identities"]["4"][key] = tampered_value
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(
        tool.ValidationError,
        match="effective.*contract|cache diagnostics|oracle cache|cache group",
    ):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize("mutation", ["missing", "extra_legacy"])
def test_periodic_oracle_receipt_requires_only_triangle_kernel_model(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    if mutation == "missing":
        del kernel_oracles["triangle_p0"]
    else:
        kernel_oracles["legacy_point"] = copy.deepcopy(kernel_oracles["triangle_p0"])
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="kernel model set mismatch"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


def test_periodic_oracle_receipt_rejects_triangle_config_tampering(
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
    library = Path(manifest["analysis_library"]["staged_path"])
    _write_periodic_oracle_receipt(tool, validation_root, library)
    config = validation_root / "provenance/oracles/periodic_plane.toml"
    config.write_text(
        config.read_text(encoding="utf-8") + "\n# tampered\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="config"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    ("field_name", "tampered_value"),
    [
        ("field_source_model", "point"),
        ("field_kernel_id", "softened_point"),
    ],
)
def test_periodic_oracle_receipt_rejects_legacy_kernel_identity(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    field_name: str,
    tampered_value: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    kernel_oracles["triangle_p0"][field_name] = tampered_value
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match=field_name):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize("mutation", ["integration", "orders", "wrench_count"])
def test_periodic_oracle_receipt_gates_triangle_integration_structure(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    triangle = kernel_oracles["triangle_p0"]["uniform_plane"]
    if mutation == "integration":
        triangle["target_integration"] = "centroid"
    elif mutation == "orders":
        triangle["quadrature_order"] = [7, 7]
    else:
        triangle["wrench_refinement"].pop()
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="triangle_p0.*integration"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize("mutation", ["missing_component", "nonzero_identity"])
def test_periodic_oracle_receipt_gates_triangle_component_identities(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    row = kernel_oracles["triangle_p0"]["uniform_plane"]["wrench_refinement"][1]
    if mutation == "missing_component":
        del row["other_objects_normalized_absolute"]
    else:
        row["total_minus_images_normalized_absolute"] = 1.0e-5
        row["component_consistency_relative_error"] = 1.0e-5
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="triangle_p0"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize("oracle_kind", ["uniform", "cosine"])
def test_periodic_oracle_receipt_rejects_triangle_threshold_excess(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    oracle_kind: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    oracle = kernel_oracles["triangle_p0"]
    if oracle_kind == "uniform":
        oracle["uniform_plane"]["nonzero_relative_error"] = 0.13
    else:
        oracle["neutral_cosine_plane"]["errors"][1]["field_relative_error"] = 0.09
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="threshold"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_label",
        "duplicate_label",
        "fingerprint_drift",
        "path_drift",
        "sha256_drift",
        "cache_miss",
        "cache_build",
    ],
)
def test_periodic_oracle_receipt_binds_every_triangle_cache_evaluation(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
    )
    library = Path(manifest["analysis_library"]["staged_path"])
    kernel_oracles = _fake_periodic_kernel_oracles()
    oracle = kernel_oracles["triangle_p0"]
    evaluations = oracle["cache_evaluations"]
    if mutation == "missing_label":
        evaluations.pop()
    elif mutation == "duplicate_label":
        evaluations[-1]["label"] = evaluations[0]["label"]
    elif mutation == "fingerprint_drift":
        evaluations[-1]["fingerprint"] = "drifted-triangle-cache"
    elif mutation == "path_drift":
        evaluations[-1]["path"] = str(
            (validation_root / "cache/oracles/operator-unbound.bin").resolve()
        )
    elif mutation == "sha256_drift":
        evaluations[-1]["sha256"] = "0" * 64
    elif mutation == "cache_miss":
        evaluations[-1]["hit"] = False
    elif mutation == "cache_build":
        evaluations[-1]["build_count"] = 1
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="cache evaluation"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


def test_periodic_oracle_receipt_rejects_self_hash_tampering(
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
    library = Path(manifest["analysis_library"]["staged_path"])
    receipt = _write_periodic_oracle_receipt(tool, validation_root, library)
    value = json.loads(receipt.read_text(encoding="utf-8"))
    value["status"] = "tampered"
    receipt.chmod(0o644)
    receipt.write_text(
        json.dumps(value, sort_keys=True, separators=(",", ":")),
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="self-hash mismatch"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


def test_periodic_oracle_receipt_rejects_library_hash_mismatch(
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
    library = Path(manifest["analysis_library"]["staged_path"])
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        overrides={"library_sha256": "0" * 64},
    )

    with pytest.raises(tool.ValidationError, match="no longer matches staging"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


def test_periodic_oracle_receipt_rejects_cache_inventory_change(
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
    library = Path(manifest["analysis_library"]["staged_path"])
    _write_periodic_oracle_receipt(tool, validation_root, library)
    (validation_root / "cache/oracles/unrecorded.bin").write_bytes(b"changed\n")

    with pytest.raises(tool.ValidationError, match="no longer matches staging"):
        tool._verify_periodic_oracle_receipt(validation_root, library)
