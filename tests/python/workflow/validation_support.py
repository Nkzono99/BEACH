"""Shared fixture writers and loaders for the periodic validation tool tests."""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]


TOOL_PATH = ROOT / "tools" / "periodic_object_validation.py"


LEGACY_NATIVE_KEYS = ((149001, 7), (180001, 6), (279001, 6), (279001, 7))


def _write_release_mechanics_fixture(run: Path) -> None:
    (run / "input/release_kernel_base.toml").write_text(
        """
[adhesion]
model = "vdw_work"
hamaker_constant = 1.0e-19
contact_distance = 0.4e-9
cutoff_distance = 10.0e-9
roughness_factor = 0.1
contact_geometry = "sphere_sphere"
contact_count = 3.0
peel_factor = 0.5
""".strip()
        + "\n",
        encoding="utf-8",
    )
    analysis = run / "analysis/local_release"
    analysis.mkdir(parents=True, exist_ok=True)
    (analysis / "release_model_summary.json").write_text(
        json.dumps(
            {
                "assumptions": {
                    "radius_m": 3.5e-5,
                    "dust_density_kg_m3": 3000.0,
                    "moon_gravity_m_s2": 1.62,
                    "energy_partition": 0.5,
                }
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def _write_legacy_estimator_fixture(run: Path) -> None:
    analysis = run / "analysis/local_release"
    analysis.mkdir(parents=True, exist_ok=True)
    radius = 3.5e-5
    z_max = 2.0 * radius
    (analysis / "moving_sphere_model_summary.json").write_text(
        json.dumps(
            {
                "status": "ok",
                "model": "moving_top_mesh_pairwise_coulomb",
                "samples": 280,
                "z_samples": 2,
                "z_max_m": z_max,
                "pair_softening_m": 0.0,
                "radius_m": radius,
                "force_curves_csv": "local_release/moving_sphere_force_curves.csv",
                "release_summary_csv": "local_release/moving_sphere_release_summary.csv",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    curve_header = (
        "target_mesh_id,batch,time_s,target_charge_multiplier,"
        "displacement_z_m,f_coulomb_z_n,f_resist_n,net_force_z_n\n"
    )
    curve_lines = [curve_header]
    release_header = (
        "target_mesh_id,batch,time_s,processed_particles,rel_change,"
        "target_charge_multiplier,radius_m,z_max_m,z_samples,q_target_base_c,"
        "q_target_effective_c,q_source_c,f_coulomb_initial_z_n,"
        "f_coulomb_final_z_n,f_coulomb_max_z_n,f_coulomb_min_z_n,f_adh_n,"
        "f_gravity_n,f_resist_n,w_coulomb_signed_j,w_resist_j,w_net_signed_j,"
        "w_net_positive_part_j,v_release_signed_net_m_per_s,"
        "v_release_positive_part_m_per_s,crossed_force_threshold,"
        "first_crossing_z_m\n"
    )
    release_lines = [release_header]
    timeseries_header = (
        "top_mesh_id,radius_m,local_charge_multiplier,batch,time_s,"
        "processed_particles,rel_change,q_base_signed_c,q_base_abs_c,"
        "q_effective_abs_c,f_direct_object_z_n,f_local_pair_proxy_n,f_adh_n,"
        "f_gravity_n,f_resist_n,net_direct_z_n,net_local_pair_n,w_elec_proxy_j,"
        "w_barrier_j,w_release_proxy_j,v_release_proxy_m_per_s\n"
    )
    timeseries_lines = [timeseries_header]
    for batch, mesh_id in LEGACY_NATIVE_KEYS:
        charge = (mesh_id - 5) * 1.0e-15
        force0 = mesh_id * 1.0e-12 + batch * 1.0e-18
        force1 = force0 + 1.0e-12
        work = 0.5 * (force0 + force1) * z_max
        for displacement, force in ((0.0, force0), (z_max, force1)):
            curve_lines.append(
                f"{mesh_id},{batch},1.0,1.0,{displacement:.17g},"
                f"{force:.17g},2e-12,{force - 2.0e-12:.17g}\n"
            )
        release_lines.append(
            f"{mesh_id},{batch},1.0,1,0.0,1.0,{radius:.17g},{z_max:.17g},2,"
            f"{charge:.17g},{charge:.17g},3e-15,{force0:.17g},{force1:.17g},"
            f"{force1:.17g},{force0:.17g},1e-12,1e-13,1.1e-12,{work:.17g},"
            "1e-18,1e-18,1e-18,1.0,1.0,1,0.0\n"
        )
        direct = 0.9 * force0
        proxy = 1.1 * abs(force0)
        timeseries_lines.append(
            f"{mesh_id},{radius:.17g},1.0,{batch},1.0,1,0.0,{charge:.17g},"
            f"{abs(charge):.17g},{abs(charge):.17g},{direct:.17g},"
            f"{proxy:.17g},1e-12,1e-13,1.1e-12,{direct - 1.1e-12:.17g},"
            f"{proxy - 1.1e-12:.17g},1e-18,1e-18,0.0,0.0\n"
        )
    (analysis / "moving_sphere_force_curves.csv").write_text(
        "".join(curve_lines), encoding="utf-8"
    )
    (analysis / "moving_sphere_release_summary.csv").write_text(
        "".join(release_lines), encoding="utf-8"
    )
    (analysis / "force_timeseries.csv").write_text(
        "".join(timeseries_lines), encoding="utf-8"
    )


def _load_tool():
    spec = importlib.util.spec_from_file_location(
        "beach_periodic_object_validation", TOOL_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _canonical_build_info(
    commit: str,
    *,
    version: str = "1.5.0-test",
    state: str = "clean",
) -> dict[str, object]:
    return {
        "build_info_schema_version": 1,
        "build_version": version,
        "build_version_mode": "git",
        "build_source_commit": commit,
        "build_id": f"{commit}:{state}",
    }


def _stage_clean_validation(
    tool,
    archive_run: Path,
    validation_root: Path,
    binary: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> dict[str, object]:
    commit = "a" * 40
    build_origin = _canonical_build_info(commit)
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
        lambda _path: dict(build_origin),
        raising=False,
    )
    monkeypatch.setattr(
        tool,
        "_library_build_info",
        lambda _path: dict(build_origin),
        raising=False,
    )
    return tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=binary,
        require_clean_source=True,
    )


def _load_toml(path: Path) -> dict[str, object]:
    import tomllib

    with path.open("rb") as stream:
        document = tomllib.load(stream)
    from beach.config import to_runtime_layout
    from beach.config._layout import is_grouped
    if is_grouped(document):
        document = to_runtime_layout(document)
        document.setdefault("output", {}).setdefault("resume", False)
    return document


_TEST_CASE_PRODUCER_ROLES = {
    "cache_prime": "smoke",
    "smoke_finite_configured": "smoke",
    "smoke_infinite_physical": "smoke",
    "full_finite_configured_140000": "finite_140000",
    "full_finite_configured_280000": "finite_280000",
    "full_infinite_physical_140000": "infinite_140000",
    "full_infinite_physical_280000": "infinite_280000",
}


_TEST_SUBMITTED_JOBS = {
    "smoke": ("9001", "beach-periodic-smoke", "p=6:t=112:c=112"),
    "finite_140000": ("9002", "beach-finite-140k", "p=6:t=112:c=112"),
    "finite_280000": ("9003", "beach-finite-280k", "p=6:t=112:c=112"),
    "infinite_140000": ("9004", "beach-infinite-140k", "p=6:t=112:c=112"),
    "infinite_280000": ("9005", "beach-infinite-280k", "p=6:t=112:c=112"),
    "analysis": ("9006", "beach-periodic-analysis", "p=1:t=28:c=28"),
}


def _write_run_output(
    validation_root: Path,
    case_name: str,
    *,
    batches: int,
    cache_hit: bool | None,
    build_count: int | None,
    fingerprint: str = "cache-fingerprint-a",
    mesh_sources: list[dict[str, object]] | None = None,
) -> Path:
    manifest = json.loads(
        (validation_root / "manifest.json").read_text(encoding="utf-8")
    )
    case = manifest["cases"][case_name]
    out = Path(case["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    if mesh_sources is None:
        mesh_sources = [
            {
                "mesh_id": 1,
                "source_kind": "template",
                "template_kind": "sphere",
                "surface_model": "insulator",
                "epsilon_r": 1.0,
                "elem_count": 1,
            }
        ]
    element_mesh_ids = [
        int(source["mesh_id"])
        for source in mesh_sources
        for _ in range(int(source["elem_count"]))
    ]
    mesh_nelem = len(element_mesh_ids)
    mesh_count = len(mesh_sources)
    backend = (
        "cached_kneq0"
        if case["periodic_model"] == "infinite_physical"
        else "legacy_finite_images"
    )
    zero_mode = (
        "exclude_k0"
        if case["periodic_model"] == "infinite_physical"
        else "legacy_not_decomposed"
    )
    lower_boundary = (
        "e_bottom_zero"
        if case["periodic_model"] == "infinite_physical"
        else "legacy_implicit"
    )
    lines = [
        "checkpoint_schema_version=11",
        "mesh_fingerprint=1111111111111111",
        "field_backend=fmm",
        "field_normalization=si",
        "field_source_model=triangle_p0",
        "field_kernel_id=triangle_p0_exact_p2m_near",
        f"mesh_nelem={mesh_nelem}",
        f"mesh_count={mesh_count}",
        "mpi_world_size=6",
        "processed_particles=10",
        "absorbed=6",
        "escaped=4",
        f"batches={batches}",
        "escaped_boundary=3",
        "survived_max_step=1",
        "last_rel_change=1.0e-3",
        f"periodic2_nonzero_mode_backend={backend}",
        f"periodic2_zero_mode_policy={zero_mode}",
        f"periodic2_lower_boundary_model={lower_boundary}",
        "periodic2_generation_tolerance=1.0e-8",
        f"electrostatic_split_periodic_active={'T' if backend == 'cached_kneq0' else 'F'}",
        f"electrostatic_status={'cached_kneq0' if backend == 'cached_kneq0' else 'not_applicable'}",
        "gauss_residual_C=0.0",
        "charge_ledger_nspecies=1",
        f"charge_ledger_batch_count={batches}",
        "charge_ledger_residual_C=0.0",
    ]
    build_origin = manifest.get("build_origin")
    if build_origin is not None:
        lines[1:1] = [
            f"{key}={build_origin[key]}"
            for key in (
                "build_info_schema_version",
                "build_version",
                "build_version_mode",
                "build_source_commit",
                "build_id",
            )
        ]
    if cache_hit is not None:
        cache_path = validation_root / "cache/periodic2/operator.bin"
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_bytes(b"deterministic-cached-operator\n")
        lines.extend(
            [
                f"periodic2_cache_hit={'T' if cache_hit else 'F'}",
                f"periodic2_operator_build_count={build_count}",
                f"periodic2_cache_fingerprint={fingerprint}",
                f"periodic2_cache_path={cache_path}",
            ]
        )
    (out / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    charge_rows = ["elem_idx,charge_C"]
    charge_rows.extend(f"{index},1.0e-15" for index in range(1, mesh_nelem + 1))
    (out / "charges.csv").write_text("\n".join(charge_rows) + "\n", encoding="utf-8")
    potential_rows = ["elem_idx,potential_V"]
    potential_rows.extend(f"{index},2.0" for index in range(1, mesh_nelem + 1))
    (out / "mesh_potential.csv").write_text(
        "\n".join(potential_rows) + "\n", encoding="utf-8"
    )
    source_rows = [
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count"
    ]
    source_rows.extend(
        ",".join(
            str(source[key])
            for key in (
                "mesh_id",
                "source_kind",
                "template_kind",
                "surface_model",
                "epsilon_r",
                "elem_count",
            )
        )
        for source in mesh_sources
    )
    (out / "mesh_sources.csv").write_text(
        "\n".join(source_rows) + "\n", encoding="utf-8"
    )
    triangle_rows = [
        "elem_idx,v0_x_m,v0_y_m,v0_z_m,v1_x_m,v1_y_m,v1_z_m,"
        "v2_x_m,v2_y_m,v2_z_m,charge_C,mesh_id"
    ]
    triangle_rows.extend(
        f"{index},{2 * (index - 1)},0,0,{2 * (index - 1) + 1},0,0,"
        f"{2 * (index - 1)},1,0,1.0e-15,{mesh_id}"
        for index, mesh_id in enumerate(element_mesh_ids, start=1)
    )
    (out / "mesh_triangles.csv").write_text(
        "\n".join(triangle_rows) + "\n", encoding="utf-8"
    )
    (out / "charge_ledger.csv").write_text(
        "batch,species_idx,injected_from_remote_C,emitted_from_surface_C,"
        "absorbed_on_surface_C,escaped_to_infinity_C,discarded_unresolved_C,"
        "injected_count,"
        "emitted_count,absorbed_count,escaped_count,discarded_unresolved_count\n"
        f"{batches},1,0,0,0,0,0,10,0,6,3,1\n",
        encoding="utf-8",
    )
    # Schema 11 checkpoint (src/runtime/bem_checkpoint_contract.f90): global
    # macro residuals and a completion manifest, no RNG state files.
    (out / "macro_residuals.csv").write_text(
        "species_idx,face,residual\n1,0,0.0\n",
        encoding="utf-8",
    )
    (out / "checkpoint_complete.txt").write_text(
        "schema_version=1\nstate=complete\n"
        f"batches={batches}\nmpi_world_size=6\n"
        "macro_residuals_present=T\ncharge_ledger_present=T\n",
        encoding="utf-8",
    )
    stride = int(case["history_stride"])
    start = 1
    if case["restart_from"] is not None:
        previous = _load_toml(Path(case["config_path"]))["output"]["restart_from"]
        previous_case = next(
            value
            for value in manifest["cases"].values()
            if value["output_dir"] == previous
        )
        start = int(previous_case["batch_count"]) + 1
    history_batches = [
        batch for batch in range(start, batches + 1) if (batch - 1) % stride == 0
    ]
    history = ["batch,processed_particles,rel_change,elem_idx,charge_C"]
    history.extend(
        f"{batch},10,1.0e-3,{element},1.0e-15"
        for batch in history_batches
        for element in range(1, mesh_nelem + 1)
    )
    (out / "charge_history.csv").write_text("\n".join(history) + "\n", encoding="utf-8")
    return out


def _write_submission_provenance(
    validation_root: Path,
    manifest: dict[str, object],
) -> dict[str, str]:
    jobs = _TEST_SUBMITTED_JOBS
    job_ids = {
        "submission_complete": True,
        **{name: job_id for name, (job_id, _job, _resources) in jobs.items()},
    }
    (validation_root / "submit/job_ids.json").write_text(
        json.dumps(job_ids), encoding="utf-8"
    )
    source = manifest["source"]
    source_snapshot = manifest["source_snapshot"]
    binary = manifest["binary"]
    library = manifest["analysis_library"]
    assert isinstance(source, dict)
    assert isinstance(source_snapshot, dict)
    assert isinstance(binary, dict)
    assert isinstance(library, dict)
    cases = manifest["cases"]
    assert isinstance(cases, dict)
    for name, (job_id, job_name, resources) in jobs.items():
        token = f"{job_id}.{job_name}"
        module_path = validation_root / "provenance/modules" / f"{token}.txt"
        hash_path = validation_root / "provenance/hashes" / f"{token}.sha256"
        module_path.parent.mkdir(parents=True, exist_ok=True)
        hash_path.parent.mkdir(parents=True, exist_ok=True)
        module_path.write_text(
            "SysA/2022\nintel/2023.2\nintelmpi/2023.2\n",
            encoding="utf-8",
        )
        artifact = (
            library["staged_path"] if name == "analysis" else binary["staged_path"]
        )
        config_lines = [
            f"{cases[case_name]['config_path']}: OK"
            for case_name, role in _TEST_CASE_PRODUCER_ROLES.items()
            if role == name
        ]
        hash_path.write_text(
            "\n".join(
                (
                    f"{source_snapshot['hash_file']}: OK",
                    f"{artifact}: OK",
                    *config_lines,
                )
            )
            + "\n",
            encoding="utf-8",
        )
        if name == "analysis":
            continue
        status_path = validation_root / "provenance/jobs" / f"{token}.status"
        status_path.parent.mkdir(parents=True, exist_ok=True)
        status_path.write_text(
            "\n".join(
                (
                    f"job_id={job_id}",
                    f"job_name={job_name}",
                    f"source_commit={source['commit']}",
                    f"resources={resources}",
                    "exit_code=0",
                )
            )
            + "\n",
            encoding="utf-8",
        )
    return {name: job_id for name, (job_id, _job, _resources) in jobs.items()}


def _fake_periodic_kernel_oracle() -> dict[str, object]:
    expected_decay_ratio = float(np.exp(-np.pi * 0.25))
    cache_identities = {
        group: {
            "hit": True,
            "build_count": 0,
            "fingerprint": f"fake-triangle_p0-{group}",
            "path": f"__CACHE_PATH_{group}__",
            "sha256": f"__CACHE_SHA256_{group}__",
        }
        for group in ("4", "8")
    }
    group_for_label = {
        label: group
        for group, labels in {
            "4": ("uniform_4", "cosine_4"),
            "8": ("cosine_8",),
        }.items()
        for label in labels
    }
    return {
        "field_source_model": "triangle_p0",
        "field_kernel_id": "triangle_p0_exact_p2m_near",
        "effective_field_contract": {
            "requested_periodic_model": "infinite_physical",
            "configured_far_correction": "cached_kneq0",
            "effective_far_correction": "cached_kneq0",
            "nonzero_mode_backend": "cached_kneq0",
            "zero_mode_policy": "exclude_k0",
            "lower_boundary_model": "e_bottom_zero",
        },
        "cache_diagnostics": copy.deepcopy(cache_identities["4"]),
        "cache_identities": cache_identities,
        "cache_evaluations": [
            {
                "label": label,
                "hit": True,
                "build_count": 0,
                "fingerprint": cache_identities[group_for_label[label]]["fingerprint"],
                "path": cache_identities[group_for_label[label]]["path"],
                "sha256": cache_identities[group_for_label[label]]["sha256"],
            }
            for label in ("uniform_4", "cosine_4", "cosine_8")
        ],
        "uniform_plane": {
            "closure": "e_bottom_zero",
            "cells_per_axis": 4,
            "below_absolute_error_V_m": 0.0,
            "nonzero_relative_error": 0.01,
            "transverse_absolute_error_V_m": 0.0,
            "potential_gauge_absolute_error_V": 0.0,
            "potential_nonzero_relative_error": 0.01,
            "force_relative_errors": [0.01, 0.01],
            "force_transverse_relative_errors": [0.0, 0.0],
            "torque_relative_errors": [0.0, 0.0],
            "quadrature_relative_difference": 0.0,
            "quadrature_relative_tolerance": 0.01,
            "relative_tolerance": 0.12,
            "field_cells_per_axis": 4,
            "target_integration": "gauss_duffy_order_3_and_7",
            "quadrature_order": [3, 7],
            "interpretation": (
                "Maxwell traction under E_bottom=0; not universal free-space self force"
            ),
            "wrench_refinement": [
                {
                    "cells_per_axis": 4,
                    "force_relative_error": 0.01,
                    "transverse_relative_error": 0.0,
                    "torque_relative_error": 0.0,
                    "other_objects_normalized_absolute": 0.0,
                    "external_uniform_normalized_absolute": 0.0,
                    "total_minus_images_normalized_absolute": 0.0,
                    "primary_free_subtraction_normalized_absolute": 1.0e-12,
                    "component_consistency_relative_error": 1.0e-12,
                },
                {
                    "cells_per_axis": 4,
                    "force_relative_error": 0.01,
                    "transverse_relative_error": 0.0,
                    "torque_relative_error": 0.0,
                    "other_objects_normalized_absolute": 0.0,
                    "external_uniform_normalized_absolute": 0.0,
                    "total_minus_images_normalized_absolute": 0.0,
                    "primary_free_subtraction_normalized_absolute": 1.0e-12,
                    "component_consistency_relative_error": 1.0e-12,
                },
            ],
        },
        "neutral_cosine_plane": {
            "errors": [
                {
                    "cells_per_axis": 4,
                    "field_relative_error": 0.12,
                    "potential_relative_error": 0.11,
                    "field_decay_ratio": expected_decay_ratio * 1.10,
                    "potential_decay_ratio": expected_decay_ratio * 0.92,
                    "field_decay_ratio_relative_error": 0.10,
                    "potential_decay_ratio_relative_error": 0.08,
                    "charge_neutrality_ratio": 1.0e-16,
                },
                {
                    "cells_per_axis": 8,
                    "field_relative_error": 0.04,
                    "potential_relative_error": 0.03,
                    "field_decay_ratio": expected_decay_ratio * 1.03,
                    "potential_decay_ratio": expected_decay_ratio * 0.98,
                    "field_decay_ratio_relative_error": 0.03,
                    "potential_decay_ratio_relative_error": 0.02,
                    "charge_neutrality_ratio": 1.0e-16,
                },
            ],
            "sample_abs_z_m": [0.25, 0.5],
            "expected_decay_ratio": expected_decay_ratio,
            "decay_ratio_relative_tolerance": 0.18,
            "fine_relative_tolerance": 0.08,
            "expected_decay": "exp(-k*abs(z-z0))",
        },
    }


def _fake_periodic_kernel_oracles() -> dict[str, dict[str, object]]:
    return {"triangle_p0": _fake_periodic_kernel_oracle()}


def _write_periodic_oracle_receipt(
    tool,
    validation_root: Path,
    library: Path,
    *,
    overrides: dict[str, object] | None = None,
    kernel_oracles: dict[str, dict[str, object]] | None = None,
    shared_operator: bool = False,
) -> Path:
    cache_dir = validation_root / "cache/oracles"
    cache_dir.mkdir(parents=True, exist_ok=True)
    if kernel_oracles is None:
        kernel_oracles = _fake_periodic_kernel_oracles()
    cache_artifacts = {
        label: {
            group: cache_dir / f"operator-{label}-{group}.bin" for group in ("4", "8")
        }
        for label in kernel_oracles
    }
    for label, grouped_paths in cache_artifacts.items():
        for group, path in grouped_paths.items():
            path.write_bytes(f"deterministic-oracle-cache:{label}:{group}\n".encode())
    if not cache_artifacts:
        (cache_dir / "operator-missing-model-fixture.bin").write_bytes(
            b"deterministic-oracle-cache:missing-model\n"
        )
    config_path = validation_root / "provenance/oracles/periodic_plane.toml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    tool._oracle_panel_config(config_path)
    for label, kernel_oracle in kernel_oracles.items():
        identities = kernel_oracle["cache_identities"]
        for group, identity in identities.items():
            default_cache_path = cache_artifacts[label][group]
            if str(identity["path"]).startswith("__CACHE_PATH_"):
                identity["path"] = str(default_cache_path.resolve())
            if str(identity["sha256"]).startswith("__CACHE_SHA256_"):
                identity["sha256"] = tool._sha256(default_cache_path)
        if shared_operator:
            identities["8"] = copy.deepcopy(identities["4"])
        kernel_oracle["cache_diagnostics"] = copy.deepcopy(identities["4"])
        group_for_evaluation = {
            evaluation_label: group
            for group, evaluation_labels in tool.ORACLE_CACHE_EVALUATION_GROUPS.items()
            for evaluation_label in evaluation_labels
        }
        for evaluation in kernel_oracle.get("cache_evaluations", []):
            identity = identities[group_for_evaluation[evaluation["label"]]]
            if shared_operator:
                evaluation["fingerprint"] = identity["fingerprint"]
            if str(evaluation["path"]).startswith("__CACHE_PATH_"):
                evaluation["path"] = identity["path"]
            if str(evaluation["sha256"]).startswith("__CACHE_SHA256_"):
                evaluation["sha256"] = identity["sha256"]
    report: dict[str, object] = {
        "receipt_schema_version": 1,
        "oracle_schema_version": 3,
        "status": "qualified",
        "manifest_sha256": tool._sha256(validation_root / "manifest.json"),
        "library": str(library.resolve()),
        "library_sha256": tool._sha256(library),
        "config": str(config_path.resolve()),
        "config_sha256": tool._sha256(config_path),
        "kernel_configs": {
            "triangle_p0": {
                "path": str(config_path.resolve()),
                "sha256": tool._sha256(config_path),
            },
        },
        "cache_dir": str(cache_dir.resolve()),
        "cache_files": tool._output_inventory(cache_dir),
        "execution_job_id": "9006",
        "kernel_oracles": kernel_oracles,
        "verified_at": "2026-07-12T00:00:00+00:00",
    }
    manifest = json.loads(
        (validation_root / "manifest.json").read_text(encoding="utf-8")
    )
    if manifest.get("build_origin") is not None:
        report["library_build_origin"] = manifest["build_origin"]
    if overrides:
        report.update(overrides)
    receipt = validation_root / "provenance/oracles/periodic_plane.json"
    tool._publish_execution_receipt(receipt, report)
    return receipt


def _wrench_component(
    force: np.ndarray,
    torque: np.ndarray,
    energy: float,
) -> SimpleNamespace:
    return SimpleNamespace(
        force_N=np.asarray(force, dtype=float),
        torque_Nm=np.asarray(torque, dtype=float),
        potential_energy_J=energy,
    )
