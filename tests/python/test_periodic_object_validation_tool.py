from __future__ import annotations

import copy
import csv
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import numpy as np
import tomli_w


ROOT = Path(__file__).resolve().parents[2]
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


def _legacy_current_rows() -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    wrenches: list[dict[str, object]] = []
    curves: list[dict[str, object]] = []
    paths: list[dict[str, object]] = []
    radius = 3.5e-5
    for batch, mesh_id in LEGACY_NATIVE_KEYS:
        charge = (mesh_id - 5) * 1.0e-15
        wrenches.append(
            {
                "case": "archived_v1_3",
                "periodic_model": "configured",
                "resolved_batch": batch,
                "mesh_id": mesh_id,
                "component": "total_external",
                "total_charge_C": charge,
                "status": "available",
            }
        )
        for displacement in (0.0, 2.0 * radius):
            for component, scale in (
                ("other_objects_all_images", 0.5),
                ("target_periodic_images", 0.25),
                ("total_external", 0.75),
            ):
                curves.append(
                    {
                        "case": "archived_v1_3",
                        "periodic_model": "configured",
                        "resolved_batch": batch,
                        "mesh_id": mesh_id,
                        "component": component,
                        "displacement_m": displacement,
                        "force_z_N": scale
                        * (mesh_id * 1.0e-12 + batch * 1.0e-18),
                        "status": "converged",
                    }
                )
        paths.append(
            {
                "case": "archived_v1_3",
                "periodic_model": "configured",
                "resolved_batch": batch,
                "mesh_id": mesh_id,
                "radius_m": radius,
                "endpoint_work_J": 4.0e-16,
                "status": "available",
            }
        )
    return wrenches, curves, paths


def _load_tool():
    spec = importlib.util.spec_from_file_location(
        "beach_periodic_object_validation", TOOL_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_tool_import_does_not_require_tomli_w_at_runtime() -> None:
    script = f"""
import importlib.abc
import importlib.util
import sys

class BlockTomliWriter(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path, target=None):
        if fullname == "tomli_w":
            raise ModuleNotFoundError("blocked runtime-only import test")
        return None

sys.modules.pop("tomli_w", None)
sys.meta_path.insert(0, BlockTomliWriter())
spec = importlib.util.spec_from_file_location(
    "beach_periodic_object_validation_without_tomli_w", {str(TOOL_PATH)!r}
)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.fixture
def archive_run(tmp_path: Path) -> Path:
    run = tmp_path / "archive" / "R20260625-0002"
    (run / "input").mkdir(parents=True)
    (run / "work" / "latest").mkdir(parents=True)
    (run / "input" / "beach.toml").write_text(
        """
[sim]
dt = 1.979898987322333e-12
batch_duration = 6.060915267313266
batch_count = 280000
max_step = 100000
tol_rel = 1.0e-8
softening = 9.899494936611663e-11
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [9.899494936611664e-05, 9.899494936611664e-05, 0.0009899494936611664]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_far_correction = "none"

[particles]
[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 5000
inject_face = "z_high"
pos_low = [0.0, 0.0, 0.0009899494936611664]
pos_high = [9.899494936611664e-05, 9.899494936611664e-05, 0.0009899494936611664]
drift_velocity = [0.0, 0.0, -4.0e5]

[mesh]
mode = "template"
[[mesh.templates]]
kind = "sphere"
enabled = true
radius = 3.5e-05
center = [2.474873734152916e-05, 2.474873734152916e-05, 0.00013597484835343895]
n_lon = 24
n_lat = 10

[output]
write_files = true
dir = "work/latest"
history_stride = 1000
write_mesh_potential = true
resume = true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (run / "manifest.toml").write_text(
        """
[simulator_source]
executable = "/nonexistent/archived/beach"
exe_hash = "sha256:b82ac96cdcb6d0d5b200d4816908714a18432f52ad5d220242193ecfed2f7c63"

[job]
processes = 6
threads = 112
cores = 112
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (run / "work" / "SIMULATOR_VERSION.txt").write_text(
        "executable: /nonexistent/archived/beach\n1.3.0-v1.3.0\n",
        encoding="utf-8",
    )
    (run / "work" / "latest" / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,sphere,insulator,1.0,1\n",
        encoding="utf-8",
    )
    _write_release_mechanics_fixture(run)
    _write_legacy_estimator_fixture(run)
    return run


@pytest.fixture
def binary(tmp_path: Path) -> Path:
    path = tmp_path / "build" / "beach"
    path.parent.mkdir()
    path.write_bytes(b"deterministic-beach-test-binary\n")
    path.chmod(0o755)
    return path


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
        return tomllib.load(stream)


def _write_toml(path: Path, data: dict[str, object]) -> None:
    path.write_text(tomli_w.dumps(data), encoding="utf-8")


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


def test_stage_creates_fresh_canonical_and_immutable_restart_inputs(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"

    manifest = tool.stage_validation(archive_run, validation_root, binary)

    finite = _load_toml(validation_root / "input/full/finite_configured.toml")
    infinite = _load_toml(validation_root / "input/full/infinite_physical.toml")
    finite_140 = _load_toml(
        validation_root
        / "input/full/segments/finite_configured_140000.toml"
    )
    finite_280 = _load_toml(
        validation_root
        / "input/full/segments/finite_configured_280000.toml"
    )
    infinite_280 = _load_toml(
        validation_root
        / "input/full/segments/infinite_physical_280000.toml"
    )

    assert finite["sim"]["field_periodic_far_correction"] == "none"
    assert finite["sim"]["rng_seed"] == 12345
    assert finite["output"]["resume"] is False
    assert finite["sim"]["batch_count"] == 280000
    assert infinite["sim"]["field_periodic_far_correction"] == "cached_kneq0"
    assert infinite["sim"]["field_periodic_ewald_layers"] == 4
    assert infinite["sim"]["field_periodic_generation_tolerance"] == 1.0e-8
    assert infinite["sim"]["rng_seed"] == 12345
    assert infinite["output"]["resume"] is False

    assert finite_140["sim"]["batch_count"] == 140000
    assert finite_140["output"]["resume"] is False
    assert "restart_from" not in finite_140["output"]
    assert finite_280["sim"]["batch_count"] == 280000
    assert finite_280["output"]["resume"] is True
    assert finite_280["output"]["restart_from"].endswith(
        "/run/full/finite_configured/140000"
    )
    assert finite_280["output"]["dir"].endswith(
        "/run/full/finite_configured/280000"
    )
    assert infinite_280["output"]["restart_from"].endswith(
        "/run/full/infinite_physical/140000"
    )

    assert manifest["schema_version"] == 1
    assert manifest["resources"] == {
        "partition": "gr20001a",
        "mpi_processes": 6,
        "openmp_threads": 112,
        "cores_per_process": 112,
    }
    assert manifest["archive"]["exact_executable_available"] is False
    assert {
        relative
        for relative in manifest["archive"]["analysis_inputs"]
        if relative.startswith("analysis/local_release/")
    } == {
        "analysis/local_release/release_model_summary.json",
        "analysis/local_release/moving_sphere_model_summary.json",
        "analysis/local_release/moving_sphere_force_curves.csv",
        "analysis/local_release/moving_sphere_release_summary.csv",
        "analysis/local_release/force_timeseries.csv",
    }
    assert manifest["execution"]["fresh_resume"] is False
    assert manifest["execution"]["segment_batches"] == [140000, 280000]
    assert "output.resume" in manifest["allowed_differences"]
    assert Path(manifest["binary"]["staged_path"]).read_bytes() == binary.read_bytes()
    source_snapshot = manifest["source_snapshot"]
    assert Path(source_snapshot["tool"]).is_file()
    assert Path(source_snapshot["root"], "beach/__init__.py").is_file()
    assert Path(source_snapshot["root"], "src/physics").is_dir()
    assert Path(source_snapshot["root"], "app").is_dir()
    assert Path(source_snapshot["root"], "fpm.toml").is_file()
    assert len(source_snapshot["hash_file_sha256"]) == 64

    smoke = (validation_root / "submit/smoke_sysa.sh").read_text(encoding="utf-8")
    submit_chain = (validation_root / "submit/submit_chain.sh").read_text(
        encoding="utf-8"
    )
    assert "#SBATCH --rsc p=6:t=112:c=112" in smoke
    assert "OMP_NUM_THREADS=112" in smoke
    assert 'module switch "${current_sys}" SysA/2022' in smoke
    assert "grep -Eq '^Sys(B|C|CL|G)/'" in smoke
    assert str(source_snapshot["root"]) in smoke
    assert "sha256sum --check" in smoke
    assert "BINARY_SHA256=" in smoke
    assert "input hash mismatch before execution" in smoke
    assert "/usr/bin/time -v" in smoke
    assert "srun \"${BINARY}\"" in smoke
    assert "mpiexec" not in smoke
    assert "mpirun" not in smoke
    assert "--ntasks" not in smoke
    assert "--cpus-per-task" not in smoke
    assert "export PYTHONPATH=" in smoke
    assert "SOURCE_COMMIT=" in smoke
    assert "date -Iseconds" in smoke
    assert "trap record_status EXIT" in smoke
    assert "exit_code=" in smoke
    assert "p=6:t=112:c=112" in smoke
    assert '--producer-job-role "smoke"' in smoke
    assert "--dependency=afterok:" in submit_chain
    assert "job_ids.json" in submit_chain
    assert "verify-inputs" in submit_chain
    assert "submission_complete" in submit_chain
    assert submit_chain.count("write_job_ids") >= 3
    finite_restart = (
        validation_root / "submit/full_finite_280000_sysa.sh"
    ).read_text(encoding="utf-8")
    assert "full/finite_configured/140000" in finite_restart
    assert '--expected-batches "${restart_batches}"' in finite_restart
    assert '"140000"' in finite_restart
    assert '--producer-job-role "finite_280000"' in finite_restart
    parent_verify = next(
        line for line in finite_restart.splitlines() if "--require-existing-receipt" in line
    )
    assert "--producer-job-role" not in parent_verify

    report = tool.verify_inputs(validation_root)
    assert report["status"] == "ok"
    assert report["case_count"] == 9

    staged_tool = Path(source_snapshot["tool"])
    staged_tool.chmod(0o644)
    staged_tool.write_text("# changed\n", encoding="utf-8")
    with pytest.raises(tool.ValidationError, match="source snapshot hash mismatch"):
        tool.verify_inputs(validation_root)


def test_stage_freezes_archive_mesh_source_contract(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    archive_sources = archive_run / "work/latest/mesh_sources.csv"
    archive_sources.write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,plane,insulator,  1.00000000E+00,800\n"
        "2,template,sphere,insulator,  1.00000000E+00,432\n"
        "3,template,sphere,insulator,  1.00000000E+00,432\n"
        "4,template,sphere,insulator,  1.00000000E+00,432\n"
        "5,template,sphere,insulator,  1.00000000E+00,432\n"
        "6,template,sphere,insulator,  1.00000000E+00,432\n"
        "7,template,sphere,insulator,  1.00000000E+00,432\n",
        encoding="utf-8",
    )

    manifest = tool.stage_validation(archive_run, validation_root, binary)

    assert manifest["archive"]["mesh_source_contract"] == {
        "path": str(archive_sources.resolve()),
        "sha256": tool._sha256(archive_sources),
        "ordered_mesh_ids": [1, 2, 3, 4, 5, 6, 7],
        "by_mesh_id": {
            "1": {
                "source_kind": "template",
                "template_kind": "plane",
                "surface_model": "insulator",
                "epsilon_r": 1.0,
                "elem_count": 800,
            },
            **{
                str(mesh_id): {
                    "source_kind": "template",
                    "template_kind": "sphere",
                    "surface_model": "insulator",
                    "epsilon_r": 1.0,
                    "elem_count": 432,
                }
                for mesh_id in range(2, 8)
            },
        },
    }


def test_stage_writes_explicit_finite_periodic_image_layer_to_every_case(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"

    manifest = tool.stage_validation(archive_run, validation_root, binary)

    for case in manifest["cases"].values():
        config = _load_toml(Path(case["config_path"]))
        assert config["sim"]["field_periodic_image_layers"] == 1


@pytest.mark.parametrize(
    ("relative_path", "mutate"),
    [
        (
            "input/full/finite_configured.toml",
            lambda data: data["sim"].__setitem__("dt", 9.0),
        ),
        (
            "input/full/infinite_physical.toml",
            lambda data: data["sim"].__setitem__("rng_seed", 7),
        ),
        (
            "input/full/infinite_physical.toml",
            lambda data: data["sim"].__setitem__(
                "field_periodic_ewald_layers", 5
            ),
        ),
        (
            "input/full/finite_configured.toml",
            lambda data: data["mesh"]["templates"][0].__setitem__(
                "radius", 9.0e-5
            ),
        ),
    ],
)
def test_verify_inputs_rejects_undeclared_physics_changes(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    relative_path: str,
    mutate,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    path = validation_root / relative_path
    data = _load_toml(path)
    mutate(data)
    _write_toml(path, data)

    with pytest.raises(tool.ValidationError, match="input mismatch"):
        tool.verify_inputs(validation_root)


def test_verify_inputs_rejects_mpi_metadata_change(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["resources"]["mpi_processes"] = 5
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="MPI resource metadata"):
        tool.verify_inputs(validation_root)


def test_verify_inputs_rejects_case_graph_change(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    del manifest["cases"]["smoke_infinite_physical"]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="case graph"):
        tool.verify_inputs(validation_root)


@pytest.mark.parametrize(
    "mutation",
    [
        "history_stride",
        "cache_expectation",
        "role",
        "restart_parent",
        "output_path",
        "config_path",
    ],
)
def test_verify_inputs_rejects_case_graph_contract_mutations(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    case = manifest["cases"]["full_finite_configured_280000"]
    if mutation == "history_stride":
        case["history_stride"] = 2000
    elif mutation == "cache_expectation":
        manifest["cases"]["smoke_infinite_physical"]["cache_expectation"] = "miss"
    elif mutation == "role":
        case["role"] = "canonical"
    elif mutation == "restart_parent":
        case["restart_from"] = manifest["cases"]["smoke_finite_configured"][
            "output_dir"
        ]
    elif mutation == "output_path":
        case["output_dir"] = manifest["cases"]["full_infinite_physical_280000"][
            "output_dir"
        ]
    elif mutation == "config_path":
        case["config_path"] = manifest["cases"]["full_infinite_physical_280000"][
            "config_path"
        ]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="case graph"):
        tool.verify_inputs(validation_root)


def test_verify_inputs_can_recheck_static_contract_after_outputs_exist(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )

    assert tool.verify_inputs(
        validation_root, require_empty_outputs=False
    )["status"] == "ok"
    with pytest.raises(tool.ValidationError, match="fresh output is not empty"):
        tool.verify_inputs(validation_root)


def test_verify_inputs_rejects_extra_source_snapshot_file(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    manifest = tool.stage_validation(archive_run, validation_root, binary)
    snapshot_root = Path(manifest["source_snapshot"]["root"])
    snapshot_root.chmod(0o755)
    (snapshot_root / "sitecustomize.py").write_text(
        "raise RuntimeError('unexpected')\n", encoding="utf-8"
    )

    with pytest.raises(tool.ValidationError, match="source snapshot inventory"):
        tool.verify_inputs(validation_root)


def test_verify_inputs_requires_exact_generated_script_set(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["scripts"] = {}
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="script set"):
        tool.verify_inputs(validation_root)


def test_stage_rejects_archive_mpi_resource_mismatch(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    manifest_path = archive_run / "manifest.toml"
    manifest = _load_toml(manifest_path)
    manifest["job"]["processes"] = 5
    _write_toml(manifest_path, manifest)

    with pytest.raises(tool.ValidationError, match="archived MPI resources"):
        tool.stage_validation(archive_run, tmp_path / "validation", binary)


def test_stage_rejects_nonempty_validation_root(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    (validation_root / "cache/periodic2").mkdir(parents=True)
    (validation_root / "cache/periodic2/operator.bin").write_bytes(b"stale")

    with pytest.raises(tool.ValidationError, match="must be empty"):
        tool.stage_validation(archive_run, validation_root, binary)


@pytest.mark.parametrize("location", ["repository", "archive"])
def test_stage_rejects_validation_root_inside_source_or_archive(
    archive_run: Path,
    binary: Path,
    location: str,
) -> None:
    tool = _load_tool()
    validation_root = ROOT if location == "repository" else archive_run

    with pytest.raises(tool.ValidationError, match="outside the repository and archive"):
        tool.stage_validation(archive_run, validation_root, binary)


@pytest.mark.parametrize(
    "mutation",
    [
        "root_metadata",
        "config_dotdot",
        "config_double_slash",
        "output_trailing_slash",
        "provenance_ancestor_symlink",
    ],
)
def test_verify_inputs_rejects_noncanonical_or_symlinked_validation_paths(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if mutation == "root_metadata":
        manifest["validation_root"] = str(validation_root / ".." / "validation")
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    elif mutation == "config_dotdot":
        manifest["cases"]["smoke_finite_configured"]["config_path"] = str(
            validation_root
            / "input/smoke"
            / ".."
            / "smoke/finite_configured.toml"
        )
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    elif mutation == "config_double_slash":
        case = manifest["cases"]["smoke_finite_configured"]
        case["config_path"] = str(case["config_path"]).replace(
            "/input/", "//input/", 1
        )
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    elif mutation == "output_trailing_slash":
        case = manifest["cases"]["smoke_finite_configured"]
        case["output_dir"] = str(case["output_dir"]) + "/"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        provenance = validation_root / "provenance"
        real_provenance = validation_root / "provenance-real"
        provenance.rename(real_provenance)
        provenance.symlink_to(real_provenance, target_is_directory=True)

    with pytest.raises(tool.ValidationError, match="validation root|canonical|symlink"):
        tool.verify_inputs(validation_root)


def test_stage_rejects_validation_root_symlink_before_writing(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    real_root = tmp_path / "real-validation"
    real_root.mkdir()
    validation_root = tmp_path / "validation-link"
    validation_root.symlink_to(real_root, target_is_directory=True)

    with pytest.raises(tool.ValidationError, match="validation root.*symlink"):
        tool.stage_validation(archive_run, validation_root, binary)

    assert not any(real_root.iterdir())


@pytest.mark.parametrize(
    "mutation",
    ["input_path", "output_path", "input_ancestor_symlink"],
)
def test_verify_inputs_rejects_archive_metadata_path_substitution_and_symlinks(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if mutation == "input_path":
        relative = "input/release_kernel_base.toml"
        replacement = tmp_path / "replacement.toml"
        replacement.write_bytes((archive_run / relative).read_bytes())
        manifest["archive"]["analysis_inputs"][relative]["path"] = str(replacement)
    elif mutation == "output_path":
        relative = "work/latest/mesh_sources.csv"
        replacement = tmp_path / "replacement.csv"
        replacement.write_bytes((archive_run / relative).read_bytes())
        manifest["archive"]["analysis_outputs"][relative]["path"] = str(replacement)
    else:
        local_release = archive_run / "analysis/local_release"
        replacement = archive_run / "analysis/local_release-real"
        local_release.rename(replacement)
        local_release.symlink_to(replacement, target_is_directory=True)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="archive analysis.*canonical|symlink"):
        tool.verify_inputs(validation_root)


def test_stage_with_library_snapshots_kernel_and_adds_dependent_analysis_job(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    library = tmp_path / "build/libbeach_field_kernel.so"
    library.parent.mkdir(exist_ok=True)
    library.write_bytes(b"deterministic-field-kernel\n")

    manifest = tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=library,
    )

    staged = Path(manifest["analysis_library"]["staged_path"])
    assert staged.read_bytes() == library.read_bytes()
    assert manifest["analysis_library"]["sha256"] == tool._sha256(staged)
    analysis_job = validation_root / "submit/analysis_sysa.sh"
    analysis_text = analysis_job.read_text(encoding="utf-8")
    assert "#SBATCH --rsc p=1:t=28:c=28" in analysis_text
    assert "--require-complete" in analysis_text
    assert str(staged) in analysis_text
    assert str(archive_run.resolve()) in analysis_text
    assert "sha256sum --check" in analysis_text
    oracle_command = 'srun python3.11 "${TOOL}" probe-periodic-oracles'
    analyze_command = 'srun python3.11 "${TOOL}" analyze'
    assert analysis_text.count(oracle_command) == 1
    assert analysis_text.count(analyze_command) == 1
    assert analysis_text.index(oracle_command) < analysis_text.index(analyze_command)
    assert oracle_command + " \\\n  --validation-root \"${VALIDATION_ROOT}\" \\\n  --library \"${LIBRARY}\"" in analysis_text
    smoke_text = (validation_root / "submit/smoke_sysa.sh").read_text(
        encoding="utf-8"
    )
    assert "probe-library" in smoke_text
    assert str(staged) in smoke_text
    subprocess.run(["bash", "-n", analysis_job], check=True)
    chain = (validation_root / "submit/submit_chain.sh").read_text(
        encoding="utf-8"
    )
    assert "afterok:${finite_280}:${infinite_280}" in chain
    assert '"analysis": "%s"' in chain
    for script_name in manifest["scripts"]:
        script_text = (validation_root / "submit" / script_name).read_text(
            encoding="utf-8"
        )
        assert script_text.count("unset PYTHONHOME") == 1
        assert script_text.count("export PYTHONNOUSERSITE=1") == 1
        assert script_text.count('export PYTHONPATH="${SOURCE_ROOT}"') == 1
        assert "${PYTHONPATH:+" not in script_text
    assert tool.verify_inputs(validation_root)["status"] == "ok"


@pytest.mark.parametrize(
    "unsafe_name",
    ["validation root", "validation@root", "validation$root", "validation'root", "validation\nroot"],
)
def test_stage_rejects_unsafe_path_characters_before_filesystem_use(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    unsafe_name: str,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / unsafe_name
    validation_root.mkdir()
    (validation_root / "sentinel").write_text("do not stage\n", encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="unsafe validation root path"):
        tool.stage_validation(archive_run, validation_root, binary)


def test_verify_inputs_rejects_analysis_script_with_oracle_after_analysis(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    library = tmp_path / "build/libbeach_field_kernel.so"
    library.parent.mkdir(exist_ok=True)
    library.write_bytes(b"deterministic-field-kernel\n")
    tool.stage_validation(
        archive_run,
        validation_root,
        binary,
        library=library,
    )
    script = validation_root / "submit/analysis_sysa.sh"
    text = script.read_text(encoding="utf-8")
    oracle_token = '"${TOOL}" probe-periodic-oracles'
    analyze_token = '"${TOOL}" analyze'
    assert text.count(oracle_token) == 1
    assert text.count(analyze_token) == 1
    text = text.replace(oracle_token, '"${TOOL}" __temporary__', 1)
    text = text.replace(analyze_token, oracle_token, 1)
    text = text.replace('"${TOOL}" __temporary__', analyze_token, 1)
    script.write_text(text, encoding="utf-8")
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["scripts"]["analysis_sysa.sh"]["sha256"] = tool._sha256(script)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="plane oracles"):
        tool.verify_inputs(validation_root)


def test_submit_chain_preserves_partial_ids_and_refuses_duplicate_submission(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()
    scripts = {
        "module": "#!/bin/bash\nif [[ \" $* \" == *\" list \"* ]]; then echo SysA/2022; fi\n",
        "spartition": "#!/bin/bash\nprintf 'Partition State\\ngr20001a UP\\n'\n",
        "qgroup": "#!/bin/bash\nexit 0\n",
        "python3.11": "#!/bin/bash\nexit 0\n",
        "sbatch": (
            "#!/bin/bash\n"
            "value=0\n"
            "[ ! -f \"${SBATCH_COUNTER}\" ] || value=$(cat \"${SBATCH_COUNTER}\")\n"
            "value=$((value + 1))\n"
            "printf '%s\\n' \"${value}\" > \"${SBATCH_COUNTER}\"\n"
            "if [ \"${value}\" = \"${SBATCH_FAIL_AT:-0}\" ]; then exit 1; fi\n"
            "printf '900%s\\n' \"${value}\"\n"
        ),
    }
    for name, content in scripts.items():
        path = fake_bin / name
        path.write_text(content, encoding="utf-8")
        path.chmod(0o755)
    counter = tmp_path / "sbatch-count"
    environment = {
        **os.environ,
        "PATH": f"{fake_bin}:{os.environ['PATH']}",
        "SBATCH_COUNTER": str(counter),
        "SBATCH_FAIL_AT": "3",
    }
    chain = validation_root / "submit/submit_chain.sh"

    first = subprocess.run(
        ["bash", chain],
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert first.returncode != 0
    assert counter.read_text(encoding="utf-8").strip() == "3"
    job_ids = (validation_root / "submit/job_ids.json").read_bytes()
    job_report = json.loads(job_ids)
    assert job_report["submission_complete"] is False
    assert job_report["smoke"] == "9001"
    assert job_report["finite_140000"] == "9002"
    assert job_report["infinite_140000"] == ""
    second = subprocess.run(
        ["bash", chain],
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert second.returncode == 2
    assert "refusing to resubmit" in second.stderr
    assert counter.read_text(encoding="utf-8").strip() == "3"
    assert (validation_root / "submit/job_ids.json").read_bytes() == job_ids


def test_stage_cli_accepts_analysis_library() -> None:
    tool = _load_tool()
    args = tool._parser().parse_args(
        [
            "stage",
            "--archive-run",
            "/archive",
            "--validation-root",
            "/validation",
            "--binary",
            "/build/beach",
            "--library",
            "/build/libbeach_field_kernel.so",
        ]
    )

    assert args.library == Path("/build/libbeach_field_kernel.so")


def test_probe_periodic_oracles_cli_dispatches_validation_root_and_library(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    calls: list[tuple[Path, Path]] = []

    def fake_probe(validation_root: Path, library: Path) -> dict[str, str]:
        calls.append((validation_root, library))
        return {"status": "qualified"}

    monkeypatch.setattr(tool, "probe_periodic_oracles", fake_probe)
    validation_root = tmp_path / "validation"
    library = tmp_path / "libbeach_field_kernel.so"

    status = tool.main(
        [
            "probe-periodic-oracles",
            "--validation-root",
            str(validation_root),
            "--library",
            str(library),
        ]
    )

    assert status == 0
    assert calls == [(validation_root, library)]


def test_verify_run_cli_accepts_producer_job_role() -> None:
    tool = _load_tool()

    args = tool._parser().parse_args(
        [
            "verify-run",
            "--case-dir",
            "/validation/run/smoke/finite_configured",
            "--expected-batches",
            "100",
            "--producer-job-role",
            "smoke",
        ]
    )

    assert args.producer_job_role == "smoke"


def test_verify_inputs_rejects_wrong_generated_producer_role(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    script = validation_root / "submit/smoke_sysa.sh"
    text = script.read_text(encoding="utf-8")
    self_verify = next(
        line
        for line in text.splitlines()
        if "verify-run" in line and '--expected-batches "${batches}"' in line
    )
    if "--producer-job-role" in self_verify:
        invalid_verify = re.sub(
            r'--producer-job-role "[^"]+"',
            '--producer-job-role "finite_140000"',
            self_verify,
        )
    else:
        invalid_verify = self_verify + ' --producer-job-role "finite_140000"'
    script.write_text(text.replace(self_verify, invalid_verify), encoding="utf-8")
    manifest_path = validation_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["scripts"]["smoke_sysa.sh"]["sha256"] = tool._sha256(script)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="producer job role"):
        tool.verify_inputs(validation_root)


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
    not os.environ.get("BEACH_VALIDATION_BINARY")
    or not os.environ.get("BEACH_FIELD_KERNEL_LIB"),
    reason=(
        "native build-origin integration requires BEACH_VALIDATION_BINARY and "
        "BEACH_FIELD_KERNEL_LIB"
    ),
)
def test_native_binary_and_library_match_clean_current_source() -> None:
    tool = _load_tool()
    binary_value = os.environ.get("BEACH_VALIDATION_BINARY")
    library_value = os.environ.get("BEACH_FIELD_KERNEL_LIB")
    assert binary_value is not None
    assert library_value is not None
    binary_path = Path(binary_value).expanduser().resolve()
    library_path = Path(library_value).expanduser().resolve()

    assert binary_path.is_file()
    assert library_path.is_file()
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


def test_probe_library_executes_basic_field_and_potential(
    binary: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    calls: list[Path] = []

    class FakeKernel:
        def __init__(self, positions, charges, *, library_path):
            assert np.asarray(positions).shape == (1, 3)
            assert np.asarray(charges).shape == (1,)
            calls.append(Path(library_path))

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def eval_e(self, points):
            assert np.asarray(points).shape == (1, 3)
            return np.array([[1.0, 0.0, 0.0]])

        def eval_phi(self, points):
            return np.array([2.0])

    monkeypatch.setattr(tool, "FieldKernel", FakeKernel, raising=False)
    report = tool.probe_library(binary)

    assert report["status"] == "ok"
    assert report["library_sha256"] == tool._sha256(binary)
    assert calls == [binary.resolve()]


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
    manifest = json.loads((validation_root / "manifest.json").read_text(encoding="utf-8"))
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
    model_fingerprint = (
        "ABCDEF0123456789"
        if case["periodic_model"] == "infinite_physical"
        else "0123456789ABCDEF"
    )
    lines = [
        "checkpoint_schema_version=3",
        f"model_fingerprint={model_fingerprint}",
        "mesh_fingerprint=1111111111111111",
        "species_fingerprint=2222222222222222",
        "field_backend=fmm",
        "field_normalization=si",
        "field_source_model=point",
        "field_kernel_id=softened_point",
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
        "electrostatic_status=applicable",
        "interface_potential_V=0.0",
        "interface_normal_field_V_m=0.0",
        "last_outer_update_batch=0",
        "outer_applicability_status=0",
        "outer_nonlinear_iterations=0",
        "outer_nonlinear_residual=0.0",
        "outer_infinity_potential_V=0.0",
        "outer_debye_length_m=0.0",
        "outer_linearity_ratio=0.0",
        "outer_max_linearity_ratio=0.0",
        "outer_integrated_charge_per_area_C_m2=0.0",
        "outer_electron_current_density_A_m2=0.0",
        "outer_ion_current_density_A_m2=0.0",
        "outer_photoelectron_current_density_A_m2=0.0",
        "outer_total_current_density_A_m2=0.0",
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
    charge_rows.extend(
        f"{index},1.0e-15" for index in range(1, mesh_nelem + 1)
    )
    (out / "charges.csv").write_text(
        "\n".join(charge_rows) + "\n", encoding="utf-8"
    )
    potential_rows = ["elem_idx,potential_V"]
    potential_rows.extend(
        f"{index},2.0" for index in range(1, mesh_nelem + 1)
    )
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
        "interface_outward_gross_C,interface_returned_gross_C,injected_count,"
        "emitted_count,absorbed_count,escaped_count,discarded_unresolved_count\n"
        f"{batches},1,0,0,0,0,0,0,0,10,0,6,3,1\n",
        encoding="utf-8",
    )
    (out / "macro_residuals.csv").write_text(
        "species_idx,residual_macro_particles\n1,0.0\n",
        encoding="utf-8",
    )
    for rank in range(6):
        (out / f"rng_state_rank{rank:05d}.txt").write_text(
            "1\n12345\n", encoding="utf-8"
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
        batch
        for batch in range(start, batches + 1)
        if (batch - 1) % stride == 0
    ]
    history = ["batch,processed_particles,rel_change,elem_idx,charge_C"]
    history.extend(
        f"{batch},10,1.0e-3,{element},1.0e-15"
        for batch in history_batches
        for element in range(1, mesh_nelem + 1)
    )
    (out / "charge_history.csv").write_text(
        "\n".join(history) + "\n", encoding="utf-8"
    )
    return out


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
        ("field_source_model", "triangle_p0"),
        ("field_kernel_id", None),
        ("field_kernel_id", "triangle_p0_exact_p2m_near"),
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
        validation_root
        / "provenance/hashes"
        / f"{selected_job_id}.{job_name}.sha256"
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
        "config_sha256": manifest["cases"]["smoke_finite_configured"][
            "config_sha256"
        ],
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
            validation_root
            / "provenance/hashes/9001.beach-periodic-smoke.sha256"
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
            validation_root
            / "provenance/verified/smoke_finite_configured.json"
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
        validation_root
        / "provenance/modules/9006.beach-periodic-analysis.txt"
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
        validation_root
        / "provenance/modules/9006.beach-periodic-analysis.txt"
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


def _fake_periodic_kernel_oracle(
    field_source_model: str,
    field_kernel_id: str,
) -> dict[str, object]:
    expected_decay_ratio = float(np.exp(-np.pi * 0.25))
    cache_labels = (
        ["uniform_4", "uniform_8", "point_primary", "cosine_4", "cosine_8"]
        if field_source_model == "point"
        else ["uniform_4", "cosine_4", "cosine_8"]
    )
    cache_groups = (
        {
            "4": ("uniform_4", "cosine_4"),
            "8": ("uniform_8", "cosine_8", "point_primary"),
        }
        if field_source_model == "point"
        else {
            "4": ("uniform_4", "cosine_4"),
            "8": ("cosine_8",),
        }
    )
    cache_identities = {
        group: {
            "hit": True,
            "build_count": 0,
            "fingerprint": f"fake-{field_source_model}-{group}",
            "path": f"__CACHE_PATH_{group}__",
            "sha256": f"__CACHE_SHA256_{group}__",
        }
        for group in ("4", "8")
    }
    group_for_label = {
        label: group for group, labels in cache_groups.items() for label in labels
    }
    result: dict[str, object] = {
        "field_source_model": field_source_model,
        "field_kernel_id": field_kernel_id,
        "effective_field_contract": {
            "requested_periodic_model": "infinite_physical",
            "configured_far_correction": "none",
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
                "fingerprint": cache_identities[group_for_label[label]][
                    "fingerprint"
                ],
                "path": cache_identities[group_for_label[label]]["path"],
                "sha256": cache_identities[group_for_label[label]]["sha256"],
            }
            for label in cache_labels
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
            "field_cells_per_axis": (
                8 if field_source_model == "point" else 4
            ),
            "interpretation": (
                "Maxwell traction under E_bottom=0; not universal free-space "
                "self force"
            ),
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
                    "field_parity_relative_error": 1.0e-12,
                    "potential_parity_relative_error": 1.0e-12,
                },
                {
                    "cells_per_axis": 8,
                    "field_relative_error": 0.04,
                    "potential_relative_error": 0.03,
                    "field_decay_ratio": expected_decay_ratio * 1.03,
                    "potential_decay_ratio": expected_decay_ratio * 0.98,
                    "field_decay_ratio_relative_error": 0.03,
                    "potential_decay_ratio_relative_error": 0.02,
                    "field_parity_relative_error": 1.0e-12,
                    "potential_parity_relative_error": 1.0e-12,
                },
            ],
            "sample_abs_z_m": [0.25, 0.5],
            "expected_decay_ratio": expected_decay_ratio,
            "decay_ratio_relative_tolerance": 0.18,
            "fine_relative_tolerance": 0.08,
            "expected_decay": "exp(-k*abs(z-z0))",
        },
    }
    if field_source_model == "point":
        result["uniform_plane"]["force_relative_errors"] = [0.005, 0.001]
        result["uniform_plane"][
            "force_transverse_relative_errors"
        ] = [0.004, 0.001]
        result["uniform_plane"]["torque_relative_errors"] = [0.003, 0.001]
        result["uniform_plane"]["quadrature_relative_difference"] = 0.005
        result["uniform_plane"]["target_integration"] = "point_centroid"
        result["uniform_plane"]["quadrature_order"] = None
        result["uniform_plane"][
            "source_refinement_force_relative_difference"
        ] = 0.005
        result["uniform_plane"]["wrench_refinement"] = [
            {
                "cells_per_axis": 4,
                "force_relative_error": 0.005,
                "transverse_relative_error": 0.004,
                "torque_relative_error": 0.003,
                "other_objects_normalized_absolute": 0.0,
                "external_uniform_normalized_absolute": 0.0,
                "total_minus_images_normalized_absolute": 0.0,
                "primary_free_subtraction_normalized_absolute": 0.036,
                "component_consistency_relative_error": 0.036,
            },
            {
                "cells_per_axis": 8,
                "force_relative_error": 0.001,
                "transverse_relative_error": 0.001,
                "torque_relative_error": 0.001,
                "other_objects_normalized_absolute": 0.0,
                "external_uniform_normalized_absolute": 0.0,
                "total_minus_images_normalized_absolute": 0.0,
                "primary_free_subtraction_normalized_absolute": 0.001,
                "component_consistency_relative_error": 0.001,
            },
        ]
        result["uniform_plane"]["primary_self_subtraction"] = {
            "force_normalized_absolute": 0.0,
            "potential_relative_error": 1.0e-13,
            "relative_tolerance": 1.0e-11,
        }
        result["uniform_plane"]["point_source_refinement"] = [
            {
                "cells_per_axis": 4,
                "off_surface_modulation_rms_V_m": 0.03,
            },
            {
                "cells_per_axis": 8,
                "off_surface_modulation_rms_V_m": 0.002,
            },
        ]
        for row in result["neutral_cosine_plane"]["errors"]:
            row["charge_neutrality_ratio"] = 1.0e-16
            row["field_parity_relative_error"] = 1.0e-12
            row["potential_parity_relative_error"] = 1.0e-12
        result["softened_point_micro_oracle"] = {
            "sample_count": 4,
            "softening_to_period_ratio": 1.0e-6,
            "target_distance_over_softening": 1.0,
            "field_relative_error": 1.0e-13,
            "potential_relative_error": 1.0e-13,
            "ordinary_field_relative_mismatch": 1.0e-13,
            "ordinary_potential_relative_mismatch": 1.0e-13,
            "self_field_normalized_absolute": 0.0,
            "relative_tolerance": 1.0e-11,
        }
    else:
        result["uniform_plane"][
            "target_integration"
        ] = "gauss_duffy_order_3_and_7"
        result["uniform_plane"]["quadrature_order"] = [3, 7]
        result["uniform_plane"][
            "source_refinement_force_relative_difference"
        ] = None
        result["uniform_plane"]["wrench_refinement"] = [
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
        ]
    return result


def _fake_periodic_kernel_oracles() -> dict[str, dict[str, object]]:
    return {
        "triangle_p0": _fake_periodic_kernel_oracle(
            "triangle_p0", "triangle_p0_exact_p2m_near"
        ),
        "production_point": _fake_periodic_kernel_oracle(
            "point", "softened_point"
        ),
    }


@pytest.mark.parametrize("mutation", ["refinement", "parity", "component"])
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
    kernel_oracles = _fake_periodic_kernel_oracles()
    production = kernel_oracles["production_point"]
    if mutation == "refinement":
        rows = production["uniform_plane"]["point_source_refinement"]
        rows[1]["off_surface_modulation_rms_V_m"] = rows[0][
            "off_surface_modulation_rms_V_m"
        ]
    elif mutation == "parity":
        production["neutral_cosine_plane"]["errors"][1][
            "field_parity_relative_error"
        ] = 0.09
    else:
        row = production["uniform_plane"]["wrench_refinement"][0]
        row["total_minus_images_normalized_absolute"] = 1.0e-5
        row["component_consistency_relative_error"] = 1.0e-5

    def fake_run(*, cache_dir: Path, field_source_model: str, **_kwargs):
        label = "production_point" if field_source_model == "point" else "triangle_p0"
        result = copy.deepcopy(kernel_oracles[label])
        identities = result["cache_identities"]
        for group, identity in identities.items():
            cache_path = cache_dir / f"operator-{label}-{group}.bin"
            cache_path.write_bytes(f"cache:{label}:{group}\n".encode())
            identity["path"] = str(cache_path.resolve())
            identity["sha256"] = tool._sha256(cache_path)
        result["cache_diagnostics"] = copy.deepcopy(identities["4"])
        groups = tool.ORACLE_CACHE_EVALUATION_GROUPS[label]
        group_for_evaluation = {
            evaluation_label: group
            for group, evaluation_labels in groups.items()
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
        match="production point-plane|production_point",
    ):
        tool.probe_periodic_oracles(validation_root, library)

    assert not receipt.exists()


def _write_periodic_oracle_receipt(
    tool,
    validation_root: Path,
    library: Path,
    *,
    overrides: dict[str, object] | None = None,
    kernel_oracles: dict[str, dict[str, object]] | None = None,
) -> Path:
    cache_dir = validation_root / "cache/oracles"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_artifacts = {
        label: {
            group: cache_dir / f"operator-{label}-{group}.bin"
            for group in ("4", "8")
        }
        for label in ("triangle_p0", "production_point")
    }
    for label, grouped_paths in cache_artifacts.items():
        for group, path in grouped_paths.items():
            path.write_bytes(
                f"deterministic-oracle-cache:{label}:{group}\n".encode()
            )
    config_path = validation_root / "provenance/oracles/periodic_plane.toml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    tool._oracle_panel_config(config_path)
    point_config_path = (
        validation_root / "provenance/oracles/periodic_plane_point.toml"
    )
    tool._oracle_point_config(point_config_path)
    if kernel_oracles is None:
        kernel_oracles = _fake_periodic_kernel_oracles()
    for label, kernel_oracle in kernel_oracles.items():
        identities = kernel_oracle["cache_identities"]
        for group, identity in identities.items():
            default_cache_path = cache_artifacts[label][group]
            if str(identity["path"]).startswith("__CACHE_PATH_"):
                identity["path"] = str(default_cache_path.resolve())
            if str(identity["sha256"]).startswith("__CACHE_SHA256_"):
                identity["sha256"] = tool._sha256(default_cache_path)
        kernel_oracle["cache_diagnostics"] = copy.deepcopy(identities["4"])
        groups = tool.ORACLE_CACHE_EVALUATION_GROUPS[label]
        group_for_evaluation = {
            evaluation_label: group
            for group, evaluation_labels in groups.items()
            for evaluation_label in evaluation_labels
        }
        for evaluation in kernel_oracle.get("cache_evaluations", []):
            identity = identities[group_for_evaluation[evaluation["label"]]]
            if str(evaluation["path"]).startswith("__CACHE_PATH_"):
                evaluation["path"] = identity["path"]
            if str(evaluation["sha256"]).startswith("__CACHE_SHA256_"):
                evaluation["sha256"] = identity["sha256"]
    legacy_triangle = kernel_oracles.get("triangle_p0") or _fake_periodic_kernel_oracle(
        "triangle_p0", "triangle_p0_exact_p2m_near"
    )
    report: dict[str, object] = {
        "receipt_schema_version": 1,
        "oracle_schema_version": 2,
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
            "production_point": {
                "path": str(point_config_path.resolve()),
                "sha256": tool._sha256(point_config_path),
            },
        },
        "cache_dir": str(cache_dir.resolve()),
        "cache_files": tool._output_inventory(cache_dir),
        "execution_job_id": "9006",
        "kernel_oracles": kernel_oracles,
        # Temporary aliases keep old receipts readable while schema 2 is rolled out.
        "uniform_plane": legacy_triangle["uniform_plane"],
        "neutral_cosine_plane": legacy_triangle["neutral_cosine_plane"],
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


def test_periodic_oracle_receipt_is_bound_to_library_and_cache_inventory(
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

    report = tool._verify_periodic_oracle_receipt(validation_root, library)

    assert report["status"] == "qualified"
    assert report["receipt_sha256"] == tool._sha256(receipt)
    assert set(report["kernel_oracles"]) == {"triangle_p0", "production_point"}
    assert report["kernel_oracles"]["triangle_p0"]["field_source_model"] == (
        "triangle_p0"
    )
    assert report["kernel_oracles"]["production_point"]["field_source_model"] == (
        "point"
    )
    assert report["kernel_oracles"]["production_point"]["field_kernel_id"] == (
        "softened_point"
    )
    assert [
        row["label"]
        for row in report["kernel_oracles"]["triangle_p0"]["cache_evaluations"]
    ] == ["uniform_4", "cosine_4", "cosine_8"]
    assert [
        row["label"]
        for row in report["kernel_oracles"]["production_point"][
            "cache_evaluations"
        ]
    ] == [
        "uniform_4",
        "uniform_8",
        "point_primary",
        "cosine_4",
        "cosine_8",
    ]
    for oracle in report["kernel_oracles"].values():
        identities = oracle["cache_identities"]
        assert oracle["cache_diagnostics"] == identities["4"]
        assert identities["4"]["fingerprint"] != identities["8"]["fingerprint"]
        assert identities["4"]["path"] != identities["8"]["path"]


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
    cosine = kernel_oracles["production_point"]["neutral_cosine_plane"]
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
    else:
        cosine["decay_ratio_relative_tolerance"] = 0.19
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(
        tool.ValidationError,
        match="cosine|decay|sample|threshold|contract",
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

    assert tool._verify_periodic_oracle_receipt(validation_root, library)[
        "status"
    ] == "qualified"
    current["build_id"] = f"{'a' * 40}:dirty"
    with pytest.raises(tool.ValidationError, match="build origin"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize("alias_name", ["uniform_plane", "neutral_cosine_plane"])
def test_periodic_oracle_receipt_rejects_legacy_alias_mismatch(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    alias_name: str,
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
    nested_triangle = _fake_periodic_kernel_oracle(
        "triangle_p0", "triangle_p0_exact_p2m_near"
    )
    alias = json.loads(json.dumps(nested_triangle[alias_name]))
    if alias_name == "uniform_plane":
        alias["nonzero_relative_error"] = 0.02
    else:
        alias["errors"][1]["field_relative_error"] = 0.05
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        overrides={alias_name: alias},
    )

    with pytest.raises(
        tool.ValidationError,
        match=r"legacy.*alias|alias.*triangle_p0",
    ):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    ("field_name", "tampered_value"),
    [
        ("requested_periodic_model", "finite_images"),
        ("configured_far_correction", "cached_kneq0"),
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
    production = kernel_oracles["production_point"]
    if field_name in production["effective_field_contract"]:
        production["effective_field_contract"][field_name] = tampered_value
    else:
        key = field_name.removeprefix("cache_")
        production["cache_identities"]["4"][key] = tampered_value
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


@pytest.mark.parametrize("missing_model", ["triangle_p0", "production_point"])
def test_periodic_oracle_receipt_requires_both_kernel_models(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    missing_model: str,
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
    del kernel_oracles[missing_model]
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match=missing_model):
        tool._verify_periodic_oracle_receipt(validation_root, library)


def test_periodic_oracle_receipt_rejects_point_config_tampering(
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
    point_config = validation_root / "provenance/oracles/periodic_plane_point.toml"
    point_config.write_text(
        point_config.read_text(encoding="utf-8") + "\n# tampered\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="production_point.*config"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    ("field_name", "tampered_value"),
    [
        ("field_source_model", "triangle_p0"),
        ("field_kernel_id", "triangle_p0_exact_p2m_near"),
    ],
)
def test_periodic_oracle_receipt_rejects_nonproduction_point_identity(
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
    kernel_oracles["production_point"][field_name] = tampered_value
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
        triangle["target_integration"] = "point_centroid"
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
    row = kernel_oracles["triangle_p0"]["uniform_plane"][
        "wrench_refinement"
    ][1]
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
def test_periodic_oracle_receipt_rejects_production_point_threshold_excess(
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
    production = kernel_oracles["production_point"]
    if oracle_kind == "uniform":
        production["uniform_plane"]["nonzero_relative_error"] = 0.13
    else:
        production["neutral_cosine_plane"]["errors"][1][
            "field_relative_error"
        ] = 0.09
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="threshold"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    ("summary_key", "row_key"),
    [
        ("force_relative_errors", "force_relative_error"),
        (
            "force_transverse_relative_errors",
            "transverse_relative_error",
        ),
        ("torque_relative_errors", "torque_relative_error"),
    ],
)
def test_periodic_oracle_receipt_rejects_worse_fine_point_wrench_error(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    summary_key: str,
    row_key: str,
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
    uniform = kernel_oracles["production_point"]["uniform_plane"]
    coarse_error = float(uniform[summary_key][0])
    worse_fine_error = coarse_error + 1.0e-3
    uniform[summary_key][1] = worse_fine_error
    uniform["wrench_refinement"][1][row_key] = worse_fine_error
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(
        tool.ValidationError,
        match="production point-plane|production_point",
    ):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize("collision_kind", ["fingerprint", "path"])
def test_periodic_oracle_receipt_rejects_cross_model_cache_collision(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    collision_kind: str,
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
    triangle_cache = kernel_oracles["triangle_p0"]["cache_identities"]["4"]
    point_oracle = kernel_oracles["production_point"]
    point_cache = point_oracle["cache_identities"]["4"]
    if collision_kind == "fingerprint":
        point_cache["fingerprint"] = triangle_cache["fingerprint"]
        for evaluation in point_oracle["cache_evaluations"]:
            if evaluation["label"] in {"uniform_4", "cosine_4"}:
                evaluation["fingerprint"] = triangle_cache["fingerprint"]
    else:
        point_cache["path"] = str(
            (
                validation_root
                / "cache/oracles/operator-triangle_p0-4.bin"
            ).resolve()
        )
        for evaluation in point_oracle["cache_evaluations"]:
            if evaluation["label"] in {"uniform_4", "cosine_4"}:
                evaluation["path"] = point_cache["path"]
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="cache"):
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
        "group_collision",
    ],
)
def test_periodic_oracle_receipt_binds_every_model_cache_evaluation(
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
    evaluations = kernel_oracles["production_point"]["cache_evaluations"]
    if mutation == "missing_label":
        evaluations.pop()
    elif mutation == "duplicate_label":
        evaluations[-1]["label"] = evaluations[0]["label"]
    elif mutation == "fingerprint_drift":
        evaluations[-1]["fingerprint"] = "drifted-point-cache"
    elif mutation == "path_drift":
        evaluations[-1]["path"] = str(
            (
                validation_root
                / "cache/oracles/operator-triangle-p0.bin"
            ).resolve()
        )
    elif mutation == "sha256_drift":
        evaluations[-1]["sha256"] = "0" * 64
    elif mutation == "cache_miss":
        evaluations[-1]["hit"] = False
    elif mutation == "cache_build":
        evaluations[-1]["build_count"] = 1
    else:
        production = kernel_oracles["production_point"]
        production["cache_identities"]["8"] = copy.deepcopy(
            production["cache_identities"]["4"]
        )
        coarse = production["cache_identities"]["4"]
        for evaluation in evaluations:
            if evaluation["label"] in {"uniform_8", "cosine_8", "point_primary"}:
                for key in ("fingerprint", "path", "sha256"):
                    evaluation[key] = coarse[key]
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="cache evaluation|4/8 cache"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    "mutation", ["missing", "field_error", "scale", "ordinary"]
)
def test_periodic_oracle_receipt_requires_softened_point_micro_oracle(
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
    production = kernel_oracles["production_point"]
    if mutation == "missing":
        del production["softened_point_micro_oracle"]
    elif mutation == "field_error":
        production["softened_point_micro_oracle"]["field_relative_error"] = 2.0e-11
    elif mutation == "scale":
        production["softened_point_micro_oracle"]["softening_to_period_ratio"] = 2.0e-6
    else:
        production["softened_point_micro_oracle"][
            "ordinary_field_relative_mismatch"
        ] = 2.0e-11
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(tool.ValidationError, match="softened-point"):
        tool._verify_periodic_oracle_receipt(validation_root, library)


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_refinement",
        "nondecreasing_refinement",
        "charge",
        "parity",
        "integration",
        "potential",
        "field_cells",
        "component",
        "component_identity",
        "primary",
    ],
)
def test_periodic_oracle_receipt_gates_point_plane_structure(
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
    production = kernel_oracles["production_point"]
    if mutation == "missing_refinement":
        del production["uniform_plane"]["point_source_refinement"]
    elif mutation == "nondecreasing_refinement":
        production["uniform_plane"]["point_source_refinement"][1][
            "off_surface_modulation_rms_V_m"
        ] = 0.04
    elif mutation == "charge":
        production["neutral_cosine_plane"]["errors"][1][
            "charge_neutrality_ratio"
        ] = 1.0e-8
    elif mutation == "parity":
        production["neutral_cosine_plane"]["errors"][1][
            "field_parity_relative_error"
        ] = 0.09
    elif mutation == "integration":
        production["uniform_plane"]["target_integration"] = "gauss_duffy"
    elif mutation == "potential":
        production["uniform_plane"]["potential_nonzero_relative_error"] = 0.13
    elif mutation == "field_cells":
        production["uniform_plane"]["field_cells_per_axis"] = 4
    elif mutation == "component":
        production["uniform_plane"]["wrench_refinement"][1][
            "component_consistency_relative_error"
        ] = 0.02
    elif mutation == "component_identity":
        production["uniform_plane"]["wrench_refinement"][1][
            "other_objects_normalized_absolute"
        ] = 1.0e-5
    else:
        production["uniform_plane"]["primary_self_subtraction"][
            "potential_relative_error"
        ] = 2.0e-11
    _write_periodic_oracle_receipt(
        tool,
        validation_root,
        library,
        kernel_oracles=kernel_oracles,
    )

    with pytest.raises(
        tool.ValidationError,
        match="production point-plane|production_point",
    ):
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
            "model_fingerprint=0123456789ABCDEF",
            "model_fingerprint=FEDCBA9876543210",
        ),
        encoding="utf-8",
    )
    with pytest.raises(tool.ValidationError, match="restart parent"):
        tool.verify_run(second, 280000)


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


def test_verify_run_rejects_cache_and_particle_ledger_mismatch(
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
        cache_hit=True,
        build_count=0,
    )

    with pytest.raises(tool.ValidationError, match="cache miss"):
        tool.verify_run(prime, 1)

    bad = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    summary = (bad / "summary.txt").read_text(encoding="utf-8")
    (bad / "summary.txt").write_text(
        summary.replace("processed_particles=10", "processed_particles=11"),
        encoding="utf-8",
    )
    with pytest.raises(tool.ValidationError, match="processed_particles"):
        tool.verify_run(bad, 100)

    bad = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    ledger_path = bad / "charge_ledger.csv"
    ledger = ledger_path.read_text(encoding="utf-8")
    ledger_path.write_text(ledger.replace(",6,3,1\n", ",6,4,0\n"), encoding="utf-8")
    with pytest.raises(tool.ValidationError, match="outcome counts"):
        tool.verify_run(bad, 100)


def test_verify_run_rejects_nonfinite_mesh_and_history_gap(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
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
    triangles = (out / "mesh_triangles.csv").read_text(encoding="utf-8")
    (out / "mesh_triangles.csv").write_text(
        triangles.replace("1,0,0,0", "1,nan,0,0"), encoding="utf-8"
    )
    with pytest.raises(tool.ValidationError, match="mesh_triangles.csv"):
        tool.verify_run(out, 100)

    (out / "mesh_triangles.csv").write_text(triangles, encoding="utf-8")
    history = (out / "charge_history.csv").read_text(encoding="utf-8").splitlines()
    (out / "charge_history.csv").write_text(
        "\n".join([history[0], *history[2:]]) + "\n", encoding="utf-8"
    )
    with pytest.raises(tool.ValidationError, match="history batch sequence"):
        tool.verify_run(out, 100)


def test_verify_run_requires_configured_mesh_potential_output(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
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
    (out / "mesh_potential.csv").unlink()

    with pytest.raises(tool.ValidationError, match="mesh_potential.csv"):
        tool.verify_run(out, 100)

    out = _write_run_output(
        validation_root,
        "smoke_finite_configured",
        batches=100,
        cache_hit=None,
        build_count=None,
    )
    (out / "mesh_sources.csv").unlink()
    with pytest.raises(tool.ValidationError, match="mesh_sources.csv"):
        tool.verify_run(out, 100)


def test_verify_run_requires_complete_six_rank_checkpoint(
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
    (out / "rng_state_rank00005.txt").unlink()

    with pytest.raises(tool.ValidationError, match="checkpoint"):
        tool.verify_run(out, 1)


def test_verify_run_requires_complete_infinite_electrostatic_restart_state(
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
            if not line.startswith("outer_total_current_density_A_m2=")
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="electrostatic restart state"):
        tool.verify_run(out, 1)


def test_verify_run_requires_conditional_outer_checkpoint_files(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    archive_input = archive_run / "input/beach.toml"
    config = _load_toml(archive_input)
    config["outer_plasma"] = {
        "model": "kinetic_1d",
        "photoelectron_histogram_enabled": True,
    }
    _write_toml(archive_input, config)
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    out = _write_run_output(
        validation_root,
        "cache_prime",
        batches=1,
        cache_hit=False,
        build_count=1,
    )

    with pytest.raises(tool.ValidationError, match="outer_plasma_profile.csv"):
        tool.verify_run(out, 1)
    (out / "outer_plasma_profile.csv").write_text(
        "point,z_m,potential_V,field_V_m,charge_density_C_m3\n"
        "1,0.0,0.0,0.0,0.0\n",
        encoding="utf-8",
    )
    with pytest.raises(tool.ValidationError, match="photoelectron_histogram.csv"):
        tool.verify_run(out, 1)


def test_legacy_estimator_audit_compares_exact_native_archive_keys(
    archive_run: Path,
) -> None:
    tool = _load_tool()
    wrenches, curves, paths = _legacy_current_rows()

    rows, summary = tool._legacy_estimator_audit(
        archive_run,
        wrenches,
        curves,
        paths,
        strict=True,
    )

    assert summary["status"] == "complete"
    assert summary["covered_native_keys"] == [
        [149001, 7],
        [180001, 6],
        [279001, 6],
        [279001, 7],
    ]
    assert len(rows) == 36
    assert {int(row["batch"]) for row in rows} == {149001, 180001, 279001}
    assert 280000 not in {int(row["batch"]) for row in rows}
    assert {
        row["estimator"] for row in rows
    } == {
        "direct_object_z",
        "local_pair_proxy",
        "moving_top_mesh_pairwise_coulomb",
    }
    assert {
        row["component"]
        for row in rows
        if row["comparison_kind"] == "legacy_moving_sphere_force"
    } == {
        "other_objects_all_images",
        "target_periodic_images",
        "total_external",
    }
    assert all(row["status"] == "computed" for row in rows)


@pytest.mark.parametrize(
    "mutation",
    [
        "missing",
        "schema",
        "duplicate",
        "nonfinite",
        "charge",
        "radius",
        "curve_coverage",
        "endpoint",
        "work",
    ],
)
def test_legacy_estimator_audit_strictly_validates_inputs_and_mapping(
    archive_run: Path,
    mutation: str,
) -> None:
    tool = _load_tool()
    analysis = archive_run / "analysis/local_release"
    wrenches, curves, paths = _legacy_current_rows()
    if mutation == "missing":
        (analysis / "force_timeseries.csv").unlink()
    elif mutation == "schema":
        path = analysis / "moving_sphere_force_curves.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text(
            lines[0].replace(",net_force_z_n", "")
            + "\n"
            + "\n".join(lines[1:])
            + "\n",
            encoding="utf-8",
        )
    elif mutation == "duplicate":
        path = analysis / "moving_sphere_force_curves.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("\n".join([*lines, lines[1]]) + "\n", encoding="utf-8")
    elif mutation == "nonfinite":
        path = analysis / "force_timeseries.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        fields = lines[1].split(",")
        fields[10] = "nan"
        lines[1] = ",".join(fields)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    elif mutation == "charge":
        wrenches[0]["total_charge_C"] = 9.0e-15
    elif mutation == "radius":
        paths[0]["radius_m"] = 9.0e-5
    elif mutation == "curve_coverage":
        curves.pop()
    else:
        path = analysis / "moving_sphere_release_summary.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        fields = lines[1].split(",")
        fields[12 if mutation == "endpoint" else 19] = "9.0e-6"
        lines[1] = ",".join(fields)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(tool.ValidationError, match="legacy estimator"):
        tool._legacy_estimator_audit(
            archive_run,
            wrenches,
            curves,
            paths,
            strict=True,
        )


def test_legacy_estimator_audit_non_strict_reports_unavailable(
    archive_run: Path,
) -> None:
    tool = _load_tool()
    (archive_run / "analysis/local_release/force_timeseries.csv").unlink()

    rows, summary = tool._legacy_estimator_audit(
        archive_run,
        [],
        [],
        [],
        strict=False,
    )

    assert summary["status"] == "unavailable"
    assert rows[0]["status"] == "unavailable"


def test_analyze_archive_only_writes_stable_explicit_missing_artifacts(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    (archive_run / "work/latest/summary.txt").write_text(
        "mesh_nelem=3392\nprocessed_particles=100\nabsorbed=60\n"
        "escaped=40\nbatches=280000\nlast_rel_change=1.0e-3\n",
        encoding="utf-8",
    )
    (archive_run / "work/latest/charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-15\n",
        encoding="utf-8",
    )

    class FakeBeach:
        def __init__(self, output_dir, *, config_path=None):
            assert Path(output_dir) == archive_run / "work/latest"
            assert Path(config_path) == archive_run / "input/beach.toml"
            self.result = SimpleNamespace(
                mesh_nelem=3392,
                processed_particles=100,
                absorbed=60,
                escaped=40,
                batches=280000,
                escaped_boundary=39,
                survived_max_step=1,
                last_rel_change=1.0e-3,
                charges=[1.0e-15, -2.0e-15],
                mesh_ids=None,
                mesh_potential_v=None,
                charge_ledger=None,
                charge_ledger_residual_c=None,
                model_fingerprint=None,
                mesh_fingerprint=None,
                species_fingerprint=None,
            )

    monkeypatch.setattr(tool, "Beach", FakeBeach)
    report = tool.analyze_validation(
        archive_run,
        validation_root,
        library=binary,
    )

    analysis = validation_root / "analysis"
    expected = {
        "run_summary.csv",
        "charge_history_pair.csv",
        "particle_ledger_pair.csv",
        "mesh_potential_pair.csv",
        "snapshot_manifest.csv",
        "object_wrench.csv",
        "object_path_curves.csv",
        "object_path_summary.csv",
        "finite_shell_convergence.csv",
        "comparison_matrix.csv",
        "legacy_estimator_comparison.csv",
        "analysis_summary.json",
        "review_ja.md",
        "artifacts.json",
    }
    assert expected.issubset({path.name for path in analysis.iterdir()})
    assert report["cases"]["archived_v1_3"]["status"] == "available"
    assert report["cases"]["new_finite_configured"]["status"] == "missing"
    assert report["cases"]["new_infinite_physical"]["status"] == "missing"
    assert report["schema_version"] == 2
    qualification = report["numerical_qualification_for_local_frozen_model"]
    assert qualification["status"] == "not_qualified"
    assert qualification["status_semantics"] == (
        "path_work_shell_on_fixed_saved_discretization"
    )
    assert qualification["saved_sphere_mesh_refinement_status"] == "not_evaluated"
    assert qualification["source_discretization_refinement_status"] == "not_evaluated"
    assert qualification["plane_oracle_used_as_sphere_error_bar"] is False
    review = (analysis / "review_ja.md").read_text(encoding="utf-8")
    assert "旧結果は二つの実行バイナリ" in review
    assert "新 finite: 未実行" in review
    assert "fixed-discretization path/work/shell gate" in review
    assert "保存済み sphere mesh/source refinement: not_evaluated" in review
    assert "平面 oracle の誤差を sphere error bar に流用しません" in review
    with (analysis / "comparison_matrix.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    assert rows[0]["status"] == "missing_comparison"


def test_analyze_available_archive_uses_public_wrench_path_and_shell_apis(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    tool.stage_validation(archive_run, validation_root, binary)
    cache_path = validation_root / "cache/periodic2/operator.bin"
    cache_path.write_bytes(b"analysis-cache\n")
    (archive_run / "work/latest/summary.txt").write_text(
        "mesh_nelem=2\nprocessed_particles=100\nabsorbed=60\n"
        "escaped=40\nbatches=280000\nlast_rel_change=1.0e-3\n",
        encoding="utf-8",
    )
    (archive_run / "work/latest/charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-15\n2,2.0e-15\n",
        encoding="utf-8",
    )
    calls: dict[str, list[object]] = {
        "models": [],
        "steps": [],
        "shell_meshes": [],
        "path_grids": [],
    }

    class FakeRelease:
        def __init__(self, grid: np.ndarray):
            self.gravity_work_J = grid * 0.1
            self.adhesion_work_J = grid * 0.2
            self.available_energy_J = grid * 0.3
            self.electrostatic_only_speed_m_s = grid + 1.0
            self.gravity_corrected_speed_m_s = grid + 0.5
            self.speed_m_s = grid + 0.25
            self.barrier_free_from_rest = True
            self.endpoint_positive = True
            self.endpoint_reachable_from_rest = True
            self.endpoint_speed_m_s = 2.0
            self.maximum_reachable_speed_m_s = 3.0
            self.instantaneous_force_margin_N = 4.0
            self.minimum_available_energy_J = -5.0e-19
            self.first_inaccessible_displacement_m = None
            self.endpoint_available_energy_J = 6.0e-18

    class FakePath:
        def __init__(self, grid: np.ndarray, *, qualified: bool):
            self.displacement_m = grid
            self.force_N = np.column_stack((grid * 0.0, grid * 0.0, grid + 3.0))
            self.torque_Nm = np.column_stack((grid + 1.0, grid * 0.0, grid * 0.0))
            self.electrostatic_work_J = grid * 3.0
            self.potential_difference_work_J = grid * 3.0 if qualified else None
            self.component_force_N = {
                "other_objects_all_images": self.force_N * 0.25,
                "target_periodic_images": self.force_N * 0.5,
                "external_uniform": self.force_N * 0.25,
                "total_external": self.force_N,
            }
            self.component_torque_Nm = {
                "other_objects_all_images": self.torque_Nm * 0.25,
                "target_periodic_images": self.torque_Nm * 0.5,
                "external_uniform": self.torque_Nm * 0.25,
                "total_external": self.torque_Nm,
            }
            self.status = "converged" if qualified else "not_converged"
            self.work_relative_mismatch = 0.0
            torque_origins = np.broadcast_to(
                np.array([1.0, 2.0, 3.0]), (grid.size, 3)
            ).copy()
            torque_origins[:, 2] += grid
            self.numerical_metadata = {
                "target_geometry_representation": "periodic2_mesh_connected",
                "vertex_bounding_radius_m": 3.5e-5,
                "torque_origin_policy": "moving_geometric_area_centroid",
                "torque_origin_m": torque_origins,
            }

        def evaluate_release(self, **_kwargs):
            return FakeRelease(self.displacement_m)

    class FakeProbe:
        def __init__(self, mesh_id: int, periodic_model: str):
            self.mesh_id = mesh_id
            self.periodic_model = periodic_model
            self.vertex_bounding_radius_m = 3.5e-5
            self.target_geometry_representation = "periodic2_mesh_connected"
            self.geometric_area_centroid_m = np.array([1.0, 2.0, 3.0])

        def wrench(self, **_kwargs):
            total = SimpleNamespace(
                force_N=np.array([1.0, 2.0, 3.0]),
                torque_Nm=np.array([4.0, 5.0, 6.0]),
                potential_energy_J=7.0,
            )
            other = _wrench_component(
                total.force_N * 0.25,
                total.torque_Nm * 0.25,
                total.potential_energy_J * 0.25,
            )
            images = _wrench_component(
                total.force_N * 0.5,
                total.torque_Nm * 0.5,
                total.potential_energy_J * 0.5,
            )
            uniform = _wrench_component(
                total.force_N * 0.25,
                total.torque_Nm * 0.25,
                total.potential_energy_J * 0.25,
            )
            return SimpleNamespace(
                total_charge_C=8.0e-15,
                force_N=total.force_N,
                torque_Nm=total.torque_Nm,
                torque_origin_m=np.array([1.0, 2.0, 3.0]),
                components={
                    "other_objects_all_images": other,
                    "target_periodic_images": images,
                    "external_uniform": uniform,
                    "total_external": total,
                },
                numerical_metadata={
                    "effective_far_correction": (
                        "none"
                        if self.periodic_model == "configured"
                        else "cached_kneq0"
                    ),
                    "target_geometry_representation": (
                        "periodic2_mesh_connected"
                    ),
                    "vertex_bounding_radius_m": self.vertex_bounding_radius_m,
                    "torque_origin_policy": "moving_geometric_area_centroid",
                    "torque_origin_m": np.array([1.0, 2.0, 3.0]),
                    "periodic_cache": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "hit": True,
                            "build_count": 0,
                            "fingerprint": "analysis-cache-fingerprint",
                            "path": cache_path,
                        }
                    ),
                    "periodic_kneq0": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "force_N": np.array([0.1, 0.2, 0.3]),
                            "torque_Nm": np.array([0.4, 0.5, 0.6]),
                            "potential_energy_J": 0.7,
                        }
                    ),
                    "physical_k0": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "force_N": np.array([0.0, 0.0, 0.8]),
                            "torque_Nm": np.zeros(3),
                            "potential_energy_J": 0.9,
                        }
                    ),
                    "cached_kneq0_trace_correction": (
                        None
                        if self.periodic_model == "configured"
                        else {
                            "force_N": np.array([9.0, 8.0, 7.0]),
                            "torque_Nm": np.array([6.0, 5.0, 4.0]),
                            "potential_energy_J": 3.0,
                        }
                    ),
                    "primary_free_subtraction": {
                        "force_N": np.array([-0.1, -0.2, -0.3]),
                        "torque_Nm": np.array([-0.4, -0.5, -0.6]),
                        "potential_energy_J": -0.7,
                    },
                },
            )

        def vertical_path(self, displacement_m, **_kwargs):
            calls["path_grids"].append(np.asarray(displacement_m, dtype=float))
            return FakePath(
                np.asarray(displacement_m, dtype=float),
                qualified=self.periodic_model == "configured",
            )

    class FakeSnapshot:
        def __init__(self, model: str):
            self.model = model

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def object_probe(self, mesh_id: int):
            return FakeProbe(mesh_id, self.model)

    triangles = np.array(
        [
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]],
        ]
    )

    class FakeBeach:
        def __init__(self, output_dir, *, config_path=None):
            assert Path(output_dir) == archive_run / "work/latest"
            assert Path(config_path) == archive_run / "input/beach.toml"
            self.result = SimpleNamespace(
                mesh_nelem=2,
                processed_particles=100,
                absorbed=60,
                escaped=40,
                batches=280000,
                escaped_boundary=39,
                survived_max_step=1,
                last_rel_change=1.0e-3,
                charges=np.array([1.0e-15, 2.0e-15]),
                mesh_ids=np.array([6, 7]),
                triangles=triangles,
                mesh_potential_v=None,
                charge_ledger=None,
                charge_ledger_residual_c=None,
                model_fingerprint=None,
                mesh_fingerprint=None,
                species_fingerprint=None,
                history=SimpleNamespace(
                    has_data=True,
                    batch_indices=np.array([149001, 180001, 279001]),
                ),
            )

        def object_interaction_snapshot(self, *, periodic_model, step, **_kwargs):
            calls["models"].append(periodic_model)
            calls["steps"].append(step)
            return FakeSnapshot(periodic_model)

    def fake_shell(_snapshot, probe, displacement_m, **_kwargs):
        calls["shell_meshes"].append(probe.mesh_id)
        assert np.asarray(displacement_m)[0] == 0.0
        return SimpleNamespace(
            image_layers=np.array([0, 1, 2, 3]),
            force_increment_error_N=np.array([1.0, 0.1, 0.01]),
            work_increment_error_J=np.array([2.0, 0.2, 0.02]),
            force_tail_proxy_N=np.array([10.0, 1.0, 0.1]),
            work_tail_proxy_J=np.array([20.0, 2.0, 0.2]),
            increment_converged=np.array([False, True, True]),
            reference_model="infinite_physical",
            reference_force_error_N=np.array([3.0, 2.0, 1.0, 0.1]),
            reference_work_error_J=np.array([6.0, 4.0, 2.0, 0.2]),
            reference_converged=np.array([False, False, True, True]),
            status="converged",
            selected_image_layers=3,
        )

    monkeypatch.setattr(tool, "Beach", FakeBeach)
    monkeypatch.setattr(tool, "finite_shell_convergence", fake_shell)
    report = tool.analyze_validation(
        archive_run,
        validation_root,
        library=binary,
    )

    assert calls["models"] == [
        "configured",
        "infinite_physical",
        "configured",
        "infinite_physical",
        "configured",
        "infinite_physical",
        "configured",
        "infinite_physical",
    ]
    assert calls["steps"] == [
        149001,
        149001,
        180001,
        180001,
        279001,
        279001,
        None,
        None,
    ]
    assert calls["shell_meshes"] == [6, 7]
    assert calls["path_grids"]
    assert all(grid[0] == 0.0 for grid in calls["path_grids"])
    assert all(grid[-1] == pytest.approx(7.0e-5) for grid in calls["path_grids"])
    with (validation_root / "analysis/object_wrench.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        wrench_rows = list(csv.DictReader(stream))
    assert any(
        row["mesh_id"] == "6"
        and row["resolved_batch"] == "280000"
            and row["effective_far_correction"] == "cached_kneq0"
            and row["zero_mode_policy"] == "exclude_k0"
        and row["lower_boundary_model"] == "e_bottom_zero"
            and row["periodic_cache_hit"] == "True"
            and row["periodic_cache_build_count"] == "0"
            and row["periodic_cache_fingerprint"]
            == "analysis-cache-fingerprint"
            and row["periodic_cache_path"] == str(cache_path)
            and row["periodic_cache_path_sha256"] == tool._sha256(cache_path)
        and row["periodic_model"] == "infinite_physical"
        and row["component"] == "numerical:physical_k0"
        and row["component_kind"] == "numerical_decomposition"
        and row["status"] == "available"
        for row in wrench_rows
    )
    assert all(
        row["torque_origin_policy"] == "moving_geometric_area_centroid"
        and float(row["torque_origin_x_m"]) == pytest.approx(1.0)
        and float(row["torque_origin_y_m"]) == pytest.approx(2.0)
        and float(row["torque_origin_z_m"]) == pytest.approx(3.0)
        for row in wrench_rows
        if row["status"] == "available"
    )
    assert all(
        float(row["model_radius_m"]) == pytest.approx(3.5e-5)
        and row["radius_source"]
        == "release_model_summary.assumptions.radius_m"
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["radius_relative_mismatch"]) == pytest.approx(0.0)
        and float(row["radius_relative_tolerance"]) == pytest.approx(5.0e-3)
        and row["target_geometry_representation"]
        == "periodic2_mesh_connected"
        for row in wrench_rows
        if row["status"] == "available"
    )
    assert any(
        row["component"] == "total_external"
        and row["component_kind"] == "physical_total"
        for row in wrench_rows
    )
    assert any(
        row["component"] == "target_periodic_images"
        and row["component_kind"] == "physical_additive"
        for row in wrench_rows
    )
    assert any(
        row["component"] == "numerical:cached_kneq0_trace_correction"
        and row["component_kind"] == "numerical_diagnostic_included"
        for row in wrench_rows
    )
    wrench_groups: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in wrench_rows:
        if row["status"] != "available":
            continue
        key = (
            row["case"],
            row["periodic_model"],
            row["resolved_batch"],
            row["mesh_id"],
        )
        wrench_groups.setdefault(key, []).append(row)
    for rows in wrench_groups.values():
        physical = {
            row["component"]
            for row in rows
            if row["component_kind"].startswith("physical_")
        }
        assert physical == set(tool.PHYSICAL_OBJECT_COMPONENTS)
        numerical = {
            row["component"].removeprefix("numerical:")
            for row in rows
            if row["component"].startswith("numerical:")
        }
        expected_numerical = (
            set(tool.CACHED_NUMERICAL_COMPONENTS)
            if rows[0]["effective_far_correction"] == "cached_kneq0"
            else {"primary_free_subtraction"}
        )
        assert numerical == expected_numerical
    with (validation_root / "analysis/object_path_curves.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        curve_rows = list(csv.DictReader(stream))
    curve_keys = [
        (
            row["case"],
            row["periodic_model"],
            row["resolved_batch"],
            row["mesh_id"],
            row["component"],
            row["displacement_m"],
        )
        for row in curve_rows
    ]
    assert len(curve_keys) == len(set(curve_keys))
    assert all(
        row["torque_origin_policy"] == "moving_geometric_area_centroid"
        and float(row["torque_origin_x_m"]) == pytest.approx(1.0)
        and float(row["torque_origin_y_m"]) == pytest.approx(2.0)
        and float(row["torque_origin_z_m"])
        == pytest.approx(3.0 + float(row["displacement_m"]))
        for row in curve_rows
    )
    curve_groups: dict[tuple[str, ...], dict[str, list[tuple[str, ...]]]] = {}
    for row in curve_rows:
        key = (
            row["case"],
            row["periodic_model"],
            row["resolved_batch"],
            row["mesh_id"],
        )
        signature = (
            row["displacement_m"],
            row["torque_origin_policy"],
            row["torque_origin_x_m"],
            row["torque_origin_y_m"],
            row["torque_origin_z_m"],
        )
        curve_groups.setdefault(key, {}).setdefault(row["component"], []).append(
            signature
        )
    for components in curve_groups.values():
        assert set(components) == set(tool.PHYSICAL_OBJECT_COMPONENTS)
        reference_signature = components["total_external"]
        assert all(values == reference_signature for values in components.values())
    assert all(
        float(row["model_radius_m"]) == pytest.approx(3.5e-5)
        and row["radius_source"]
        == "release_model_summary.assumptions.radius_m"
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["radius_relative_mismatch"]) == pytest.approx(0.0)
        and float(row["radius_relative_tolerance"]) == pytest.approx(5.0e-3)
        and row["target_geometry_representation"]
        == "periodic2_mesh_connected"
        for row in curve_rows
    )
    with (validation_root / "analysis/object_path_summary.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        path_rows = list(csv.DictReader(stream))
    assert any(
        row["mesh_id"] == "7"
        and row["resolved_batch"] == "280000"
        and row["periodic_model"] == "configured"
        and row["status"] == "available"
        and row["potential_work_available"] == "True"
        and row["numerically_qualified"] == "True"
            and row["barrier_free_from_rest"] == "True"
            and row["minimum_available_energy_J"] == "-5e-19"
            and row["endpoint_available_energy_J"] == "6e-18"
        for row in path_rows
    )
    declared_rows = [row for row in path_rows if row["status"] == "available"]
    assert declared_rows
    expected_mass = 4.0 * np.pi * 3000.0 * (3.5e-5) ** 3 / 3.0
    assert all(float(row["radius_m"]) == pytest.approx(3.5e-5) for row in declared_rows)
    assert all(float(row["mass_kg"]) == pytest.approx(expected_mass) for row in declared_rows)
    assert all(
        row["radius_source"] == "release_model_summary.assumptions.radius_m"
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["radius_relative_mismatch"]) == pytest.approx(0.0)
        and float(row["radius_relative_tolerance"]) == pytest.approx(5.0e-3)
        and row["target_geometry_representation"] == "periodic2_mesh_connected"
        and row["torque_origin_policy"] == "moving_geometric_area_centroid"
        and float(row["initial_torque_origin_x_m"]) == pytest.approx(1.0)
        and float(row["initial_torque_origin_y_m"]) == pytest.approx(2.0)
        and float(row["initial_torque_origin_z_m"]) == pytest.approx(3.0)
        and float(row["path_start_m"]) == pytest.approx(0.0)
        and float(row["path_end_m"]) == pytest.approx(7.0e-5)
        for row in declared_rows
    )
    assert any(
        row["mesh_id"] == "7"
        and row["periodic_model"] == "infinite_physical"
        and row["status"] == "available"
        and row["path_status"] == "not_converged"
        and row["potential_work_available"] == "False"
        and row["numerically_qualified"] == "False"
        and row["speed_status"] == "numerically_unqualified"
        for row in path_rows
    )
    with (validation_root / "analysis/finite_shell_convergence.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        shell_rows = list(csv.DictReader(stream))
    assert shell_rows[-1]["selected_image_layers"] == "3"
    assert shell_rows[-1]["force_tail_proxy_N"] == "0.1"
    assert shell_rows[-1]["work_tail_proxy_J"] == "0.2"
    assert shell_rows[-1]["reference_model"] == "infinite_physical"
    assert shell_rows[-1]["reference_force_error_N"] == "0.1"
    assert shell_rows[-1]["reference_work_error_J"] == "0.2"
    assert shell_rows[-1]["reference_converged"] == "True"
    assert all(
        float(row["model_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["geometry_radius_m"]) == pytest.approx(3.5e-5)
        and float(row["path_start_m"]) == pytest.approx(0.0)
        and float(row["path_end_m"]) == pytest.approx(7.0e-5)
        for row in shell_rows
    )
    assert report["physics_evaluation"]["status"] == "available"
    review = (validation_root / "analysis/review_ja.md").read_text(encoding="utf-8")
    assert "mesh 6" in review
    assert "batch 280000" in review
    assert "fixed-discretization path/work/shell gate: path_work=True" in review
    assert "barrier on fixed saved discretization=True" in review
    assert "平面 oracle の誤差を sphere error bar に流用しません" in review
    assert "physical_k0 Fz" in review
    with (validation_root / "analysis/comparison_matrix.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        comparison_rows = list(csv.DictReader(stream))
    assert any(
        row["comparison_kind"] == "frozen_field_override"
        and row["metric"] == "force_z_N"
        and row["resolved_batch"] == "280000"
        and row["mesh_id"] == "6"
        for row in comparison_rows
    )
    assert any(
        row["comparison_kind"] == "frozen_field_override"
        and row["metric"] == "barrier_free_from_rest"
        for row in comparison_rows
    )
    assert any(
        row["comparison_kind"] == "frozen_field_override"
        and row["metric"] == "minimum_available_energy_J"
        for row in comparison_rows
    )

    strict_run = FakeBeach(
        archive_run / "work/latest",
        config_path=archive_run / "input/beach.toml",
    )
    _wrench, _curves, _paths, _shells, strict_physics = (
        tool._evaluate_object_physics_at_step(
            archive_run=archive_run,
            validation_root=validation_root,
            library=binary,
            run_rows=[{"case": "archived_v1_3", "status": "available"}],
            runs={"archived_v1_3": strict_run},
            step=None,
            evaluation_target_ids=(6,),
            run_shell=False,
            allow_geometry_radius_fallback=False,
            require_complete_contract=True,
        )
    )
    assert strict_physics["status"] == "available"
    assert strict_physics["successful_target_models"] == 2
    assert strict_physics["failures"] == []


def test_release_parameters_reject_invalid_vdw_config(tmp_path: Path) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "input/release_kernel_base.toml").write_text(
        '[adhesion]\nmodel = "vdw_work"\nhamaker_constant = 1.0e-19\n',
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="invalid vdw adhesion"):
        tool._release_parameters(archive, 3.5e-5)


def test_release_parameters_rejects_declared_geometry_radius_mismatch(
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    _write_release_mechanics_fixture(archive)

    within_tolerance = tool._release_parameters(
        archive,
        3.51e-5,
        allow_geometry_radius_fallback=False,
    )
    assert within_tolerance["radius_m"] == pytest.approx(3.5e-5)
    assert within_tolerance["geometry_radius_m"] == pytest.approx(3.51e-5)
    assert within_tolerance["mass_kg"] == pytest.approx(
        4.0 * np.pi * 3000.0 * (3.5e-5) ** 3 / 3.0
    )
    declared_coefficient = 3.0 * 0.1 * 1.0e-19 * (3.5e-5 / 2.0) / 6.0
    assert within_tolerance["adhesion_force_N"] == pytest.approx(
        declared_coefficient / (0.4e-9) ** 2
    )
    assert within_tolerance["adhesion_work_J"] == pytest.approx(
        0.5 * declared_coefficient * (1.0 / 0.4e-9 - 1.0 / 10.0e-9)
    )

    with pytest.raises(tool.ValidationError, match="radius.*mismatch"):
        tool._release_parameters(
            archive,
            3.6e-5,
            allow_geometry_radius_fallback=False,
        )


def test_release_parameters_geometry_fallback_is_nonstrict_and_missing_only(
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "input/release_kernel_base.toml").write_text(
        '[adhesion]\nmodel = "none"\n',
        encoding="utf-8",
    )

    mechanics = tool._release_parameters(
        archive,
        3.5e-5,
        allow_geometry_radius_fallback=True,
    )

    assert mechanics["radius_m"] == pytest.approx(3.5e-5)
    assert mechanics["geometry_radius_m"] == pytest.approx(3.5e-5)
    assert mechanics["radius_source"] == "probe.vertex_bounding_radius_m_fallback"
    with pytest.raises(tool.ValidationError, match="declared dust radius is missing"):
        tool._release_parameters(
            archive,
            3.5e-5,
            allow_geometry_radius_fallback=False,
        )

    summary = archive / "analysis/local_release/release_model_summary.json"
    summary.parent.mkdir(parents=True)
    summary.write_text(
        json.dumps({"assumptions": {"radius_m": None}}),
        encoding="utf-8",
    )
    with pytest.raises(tool.ValidationError, match="declared dust radius"):
        tool._release_parameters(
            archive,
            3.5e-5,
            allow_geometry_radius_fallback=True,
        )


@pytest.mark.parametrize("invalid_radius", [True, np.array([3.5e-5])])
def test_probe_geometry_contract_rejects_non_scalar_radius(invalid_radius: object) -> None:
    tool = _load_tool()
    probe = SimpleNamespace(
        target_geometry_representation="periodic2_mesh_connected",
        vertex_bounding_radius_m=invalid_radius,
    )

    with pytest.raises(tool.ValidationError, match="bounding radius"):
        tool._probe_geometry_contract(probe)


def test_torque_origin_provenance_wraps_array_conversion_errors() -> None:
    tool = _load_tool()

    with pytest.raises(tool.ValidationError, match="torque origin metadata"):
        tool._torque_origin_provenance(
            {
                "torque_origin_policy": "moving_geometric_area_centroid",
                "torque_origin_m": object(),
            }
        )


def test_production_torque_origin_contract_tracks_final_displacement_grid() -> None:
    tool = _load_tool()
    probe = SimpleNamespace(geometric_area_centroid_m=np.array([1.0, 2.0, 3.0]))
    displacement = np.array([0.0, 0.1, 0.25])
    origins = np.broadcast_to(np.array([1.0, 2.0, 3.0]), (3, 3)).copy()
    origins[:, 2] += displacement

    tool._validate_production_torque_origin_contract(
        probe=probe,
        wrench_policy="moving_geometric_area_centroid",
        wrench_origin_m=origins[0],
        path_policy="moving_geometric_area_centroid",
        path_origin_m=origins,
        displacement_m=displacement,
    )
    origins[1, 0] += 1.0e-3
    with pytest.raises(tool.ValidationError, match="moving torque origins"):
        tool._validate_production_torque_origin_contract(
            probe=probe,
            wrench_policy="moving_geometric_area_centroid",
            wrench_origin_m=np.array([1.0, 2.0, 3.0]),
            path_policy="moving_geometric_area_centroid",
            path_origin_m=origins,
            displacement_m=displacement,
        )


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


def test_strict_wrench_component_contract_requires_complete_additive_sets() -> None:
    tool = _load_tool()
    first = _wrench_component(np.array([1.0, 0.0, 0.0]), np.zeros(3), 1.0)
    second = _wrench_component(np.array([0.0, 2.0, 0.0]), np.zeros(3), 2.0)
    third = _wrench_component(np.array([0.0, 0.0, 3.0]), np.zeros(3), 3.0)
    total = _wrench_component(np.array([1.0, 2.0, 3.0]), np.zeros(3), 6.0)
    components = {
        "other_objects_all_images": first,
        "target_periodic_images": second,
        "external_uniform": third,
        "total_external": total,
    }
    diagnostic = {
        "force_N": np.zeros(3),
        "torque_Nm": np.zeros(3),
        "potential_energy_J": 0.0,
    }
    cached_metadata = {
        name: diagnostic
        for name in (
            "periodic_kneq0",
            "physical_k0",
            "primary_free_subtraction",
            "cached_kneq0_trace_correction",
        )
    }

    def validate(
        component_values: dict[str, SimpleNamespace],
        metadata_values: dict[str, object],
        *,
        effective_far_correction: str,
        wrench_force_N: np.ndarray = total.force_N,
        wrench_torque_Nm: np.ndarray = total.torque_Nm,
    ) -> None:
        tool._validate_strict_wrench_component_contract(
            component_values,
            metadata_values,
            effective_far_correction=effective_far_correction,
            wrench_force_N=wrench_force_N,
            wrench_torque_Nm=wrench_torque_Nm,
        )

    validate(
        components,
        cached_metadata,
        effective_far_correction="cached_kneq0",
    )
    with pytest.raises(tool.ValidationError, match="physical wrench component set"):
        validate(
            {name: value for name, value in components.items() if name != "external_uniform"},
            cached_metadata,
            effective_far_correction="cached_kneq0",
        )
    with pytest.raises(tool.ValidationError, match="cached numerical wrench component set"):
        validate(
            components,
            {
                name: value
                for name, value in cached_metadata.items()
                if name != "physical_k0"
            },
            effective_far_correction="cached_kneq0",
        )
    nonfinite_metadata = dict(cached_metadata)
    nonfinite_metadata["physical_k0"] = {
        "force_N": np.zeros(3),
        "torque_Nm": np.zeros(3),
        "potential_energy_J": np.nan,
    }
    with pytest.raises(tool.ValidationError, match="numerical.*potential energy"):
        validate(
            components,
            nonfinite_metadata,
            effective_far_correction="cached_kneq0",
        )
    with pytest.raises(tool.ValidationError, match="additive physical wrench"):
        inconsistent = dict(components)
        inconsistent["total_external"] = _wrench_component(
            np.array([1.0, 2.0, 4.0]), np.zeros(3), 6.0
        )
        validate(
            inconsistent,
            cached_metadata,
            effective_far_correction="cached_kneq0",
        )
    with pytest.raises(tool.ValidationError, match="finite numerical wrench component set"):
        validate(
            components,
            cached_metadata,
            effective_far_correction="none",
        )
    validate(
        components,
        {"primary_free_subtraction": diagnostic},
        effective_far_correction="none",
    )
    with pytest.raises(tool.ValidationError, match="primary_free_subtraction"):
        validate(
            components,
            {},
            effective_far_correction="none",
        )
    with pytest.raises(tool.ValidationError, match="wrench total"):
        validate(
            components,
            cached_metadata,
            effective_far_correction="cached_kneq0",
            wrench_force_N=np.array([1.0, 2.0, 4.0]),
        )
    with pytest.raises(tool.ValidationError, match="wrench total"):
        validate(
            components,
            cached_metadata,
            effective_far_correction="cached_kneq0",
            wrench_torque_Nm=np.array([0.0, 0.0, 1.0]),
        )


def test_additive_identity_tolerance_scales_with_cancelling_components() -> None:
    tool = _load_tool()

    tool._require_additive_identity(
        np.array([1.0]),
        [np.array([1.0e16]), np.array([1.0]), np.array([-1.0e16])],
        label="cancelling test",
    )
    tool._require_additive_identity(
        np.array([2.0]),
        [np.array([1.0e308]), np.array([-1.0e308]), np.array([2.0])],
        label="large cancelling test",
    )


@pytest.mark.parametrize(
    "components",
    [
        [np.array([1.0e16]), np.array([1.0]), np.array([-1.0e16])],
        [np.array([1.0e308]), np.array([-1.0e308]), np.array([2.0])],
    ],
)
def test_additive_identity_rejects_false_total_under_cancellation(
    components: list[np.ndarray],
) -> None:
    tool = _load_tool()

    with pytest.raises(tool.ValidationError, match="additive components"):
        tool._require_additive_identity(
            np.array([0.0]),
            components,
            label="cancelling test",
        )


@pytest.mark.parametrize(
    ("total_charge", "potential_energy"),
    [(np.nan, 1.0), (1.0, np.inf)],
)
def test_wrench_row_rejects_nonfinite_scalar_fields(
    total_charge: float,
    potential_energy: float,
) -> None:
    tool = _load_tool()
    mechanics = {
        "radius_m": 3.5e-5,
        "radius_source": "release_model_summary.assumptions.radius_m",
        "geometry_radius_m": 3.5e-5,
        "radius_relative_mismatch": 0.0,
        "radius_relative_tolerance": 5.0e-3,
    }

    with pytest.raises(tool.ValidationError, match="finite"):
        tool._wrench_row(
            case="case",
            periodic_model="configured",
            effective_far_correction="none",
            zero_mode_policy="legacy_not_decomposed",
            lower_boundary_model="legacy_implicit",
            periodic_cache=None,
            step_selector="final",
            resolved_batch=1,
            mesh_id=6,
            mechanics=mechanics,
            target_geometry_representation="periodic2_mesh_connected",
            component="total_external",
            component_kind="physical_total",
            force=np.zeros(3),
            torque=np.zeros(3),
            torque_origin_policy="moving_geometric_area_centroid",
            torque_origin_m=np.zeros(3),
            total_charge=total_charge,
            potential_energy=potential_energy,
        )


def test_strict_path_component_contract_requires_additive_components() -> None:
    tool = _load_tool()
    first = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    second = np.array([[0.0, 2.0, 0.0], [0.0, 3.0, 0.0]])
    third = np.array([[0.0, 0.0, 3.0], [0.0, 0.0, 4.0]])
    total = first + second + third
    zeros = np.zeros_like(total)
    path = SimpleNamespace(
        displacement_m=np.array([0.0, 1.0]),
        force_N=total,
        torque_Nm=zeros,
        component_force_N={
            "other_objects_all_images": first,
            "target_periodic_images": second,
            "external_uniform": third,
            "total_external": total,
        },
        component_torque_Nm={
            "other_objects_all_images": zeros,
            "target_periodic_images": zeros,
            "external_uniform": zeros,
            "total_external": zeros,
        },
    )

    tool._validate_strict_path_component_contract(
        path,
        expected_start_m=0.0,
        expected_end_m=1.0,
    )
    path.component_force_N = dict(path.component_force_N)
    path.component_force_N.pop("external_uniform")
    with pytest.raises(tool.ValidationError, match="path force component set"):
        tool._validate_strict_path_component_contract(
            path,
            expected_start_m=0.0,
            expected_end_m=1.0,
        )


@pytest.mark.parametrize(
    "displacement",
    [
        np.array([1.0e-6, 1.0]),
        np.array([0.0, 0.9]),
        np.array([0.0, 0.75, 0.5, 1.0]),
    ],
)
def test_strict_path_component_contract_rejects_invalid_final_grid(
    displacement: np.ndarray,
) -> None:
    tool = _load_tool()
    zeros = np.zeros((displacement.size, 3))
    path = SimpleNamespace(
        displacement_m=displacement,
        force_N=zeros,
        torque_Nm=zeros,
        component_force_N={name: zeros for name in tool.PHYSICAL_OBJECT_COMPONENTS},
        component_torque_Nm={name: zeros for name in tool.PHYSICAL_OBJECT_COMPONENTS},
    )

    with pytest.raises(tool.ValidationError, match="path displacement"):
        tool._validate_strict_path_component_contract(
            path,
            expected_start_m=0.0,
            expected_end_m=1.0,
        )


def test_shell_reference_contract_requires_selected_physical_reference() -> None:
    tool = _load_tool()
    shell = SimpleNamespace(
        image_layers=np.array([0, 1, 2, 3]),
        increment_converged=np.array([False, True, True]),
        status="converged",
        selected_image_layers=3,
        reference_model="infinite_physical",
        reference_converged=np.array([False, False, True, True]),
        reference_force_error_N=np.array([3.0, 2.0, 1.0, 0.1]),
        reference_work_error_J=np.array([6.0, 4.0, 2.0, 0.2]),
        force_increment_error_N=np.array([1.0, 0.1, 0.01]),
        work_increment_error_J=np.array([2.0, 0.2, 0.02]),
        force_tail_proxy_N=np.array([10.0, 1.0, 0.1]),
        work_tail_proxy_J=np.array([20.0, 2.0, 0.2]),
    )

    tool._validate_shell_reference_contract(shell)
    shell.image_layers = np.array([1, 2, 3, 4])
    with pytest.raises(tool.ValidationError, match="start at zero"):
        tool._validate_shell_reference_contract(shell)
    shell.image_layers = np.array([0, 1, 3, 3])
    with pytest.raises(tool.ValidationError, match="consecutive"):
        tool._validate_shell_reference_contract(shell)
    shell.image_layers = np.array([0, 1, 2, 3])
    shell.status = "invalid"
    with pytest.raises(tool.ValidationError, match="status"):
        tool._validate_shell_reference_contract(shell)
    shell.status = "converged"
    shell.increment_converged = np.array([False, False, True])
    with pytest.raises(tool.ValidationError, match="successive"):
        tool._validate_shell_reference_contract(shell)
    shell.increment_converged = np.array([False, True, True])
    shell.reference_converged = np.array([False, False, False, True])
    with pytest.raises(tool.ValidationError, match="successive"):
        tool._validate_shell_reference_contract(shell)
    shell.reference_converged = np.array([False, False, True, False])
    with pytest.raises(tool.ValidationError, match="physical reference"):
        tool._validate_shell_reference_contract(shell)
    shell.reference_converged = np.array([False, False, True, True])
    shell.reference_force_error_N = np.array([3.0, 2.0, -1.0, 0.1])
    with pytest.raises(tool.ValidationError, match="nonnegative"):
        tool._validate_shell_reference_contract(shell)
    shell.reference_force_error_N = np.array([3.0, 2.0, 1.0, 0.1])
    shell.force_tail_proxy_N = np.array([10.0, -1.0, 0.1])
    with pytest.raises(tool.ValidationError, match="increment/tail"):
        tool._validate_shell_reference_contract(shell)
    shell.force_tail_proxy_N = np.array([10.0, 1.0])
    with pytest.raises(tool.ValidationError, match="increment/tail"):
        tool._validate_shell_reference_contract(shell)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("increment_converged", np.array([0, 1, 1])),
        ("increment_converged", np.array([False, np.nan, True])),
        ("reference_converged", np.array([0, 0, 1, 1])),
        ("reference_converged", np.array([False, False, np.inf, True])),
    ],
)
def test_shell_reference_contract_rejects_nonboolean_gate_arrays(
    field: str,
    value: np.ndarray,
) -> None:
    tool = _load_tool()
    shell = SimpleNamespace(
        image_layers=np.array([0, 1, 2, 3]),
        increment_converged=np.array([False, True, True]),
        status="converged",
        selected_image_layers=3,
        reference_model="infinite_physical",
        reference_converged=np.array([False, False, True, True]),
        reference_force_error_N=np.array([3.0, 2.0, 1.0, 0.1]),
        reference_work_error_J=np.array([6.0, 4.0, 2.0, 0.2]),
        force_increment_error_N=np.array([1.0, 0.1, 0.01]),
        work_increment_error_J=np.array([2.0, 0.2, 0.02]),
        force_tail_proxy_N=np.array([10.0, 1.0, 0.1]),
        work_tail_proxy_J=np.array([20.0, 2.0, 0.2]),
    )
    setattr(shell, field, value)

    with pytest.raises(tool.ValidationError, match="boolean"):
        tool._validate_shell_reference_contract(shell)


def test_release_parameters_rejects_malformed_assumption_files(
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "analysis/local_release").mkdir(parents=True)
    (archive / "input/release_kernel_base.toml").write_text(
        "[adhesion\nmodel = 'none'\n", encoding="utf-8"
    )

    with pytest.raises(tool.ValidationError, match="failed to read TOML"):
        tool._release_parameters(archive, 3.5e-5)

    (archive / "input/release_kernel_base.toml").unlink()
    (archive / "analysis/local_release/release_model_summary.json").write_text(
        "{invalid", encoding="utf-8"
    )
    with pytest.raises(tool.ValidationError, match="release model summary"):
        tool._release_parameters(archive, 3.5e-5)


def test_object_physics_requires_both_target_meshes(tmp_path: Path) -> None:
    tool = _load_tool()

    class IncompleteRun:
        result = SimpleNamespace(mesh_ids=np.array([6]))

        def object_interaction_snapshot(self, **_kwargs):
            raise AssertionError("snapshot must not be built for incomplete targets")

    _wrench, _curves, _paths, _shells, physics = tool._evaluate_object_physics(
        archive_run=tmp_path / "archive",
        validation_root=tmp_path / "validation",
        library=tmp_path / "kernel.so",
        run_rows=[{"case": "archived_v1_3", "status": "available"}],
        runs={"archived_v1_3": IncompleteRun()},
    )

    assert physics["status"] == "not_evaluated"
    assert physics["successful_target_models"] == 0
    assert "required target mesh ids" in physics["failures"][0]["message"]


@pytest.mark.parametrize(
    ("metadata_override", "expected_message"),
    [
        ({"effective_far_correction": "free"}, "effective far correction mismatch"),
        (
            {"target_geometry_representation": "raw_saved"},
            "target geometry representation mismatch",
        ),
        ({"vertex_bounding_radius_m": 4.0e-5}, "metadata radius"),
        (
            {"torque_origin_m": np.array([1.0, 0.0, 0.0])},
            "metadata torque origin does not match",
        ),
    ],
)
def test_object_physics_rejects_invalid_wrench_contract_metadata(
    tmp_path: Path,
    metadata_override: dict[str, object],
    expected_message: str,
) -> None:
    tool = _load_tool()
    triangles = np.array(
        [
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]],
        ]
    )

    class FakeProbe:
        def __init__(self, periodic_model: str):
            self.periodic_model = periodic_model
            self.vertex_bounding_radius_m = 3.5e-5
            self.target_geometry_representation = "periodic2_mesh_connected"

        def wrench(self, **_kwargs):
            metadata = {
                "effective_far_correction": (
                    "none"
                    if self.periodic_model == "configured"
                    else "cached_kneq0"
                ),
                "target_geometry_representation": "periodic2_mesh_connected",
                "vertex_bounding_radius_m": 3.5e-5,
                "torque_origin_policy": "moving_geometric_area_centroid",
                "torque_origin_m": np.zeros(3),
            }
            metadata.update(metadata_override)
            return SimpleNamespace(
                total_charge_C=1.0e-15,
                force_N=np.zeros(3),
                torque_Nm=np.zeros(3),
                torque_origin_m=np.zeros(3),
                components={},
                numerical_metadata=metadata,
            )

        def vertical_path(self, *_args, **_kwargs):
            raise AssertionError("mismatched effective model must stop before path")

    class FakeSnapshot:
        def __init__(self, periodic_model: str):
            self.periodic_model = periodic_model

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def object_probe(self, _mesh_id: int):
            return FakeProbe(self.periodic_model)

    class FakeRun:
        result = SimpleNamespace(
            batches=280000,
            mesh_ids=np.array([6, 7]),
            triangles=triangles,
        )

        def object_interaction_snapshot(self, **kwargs):
            return FakeSnapshot(kwargs["periodic_model"])

    _wrench, _curves, _paths, _shells, physics = (
        tool._evaluate_object_physics_at_step(
            archive_run=tmp_path / "archive",
            validation_root=tmp_path / "validation",
            library=tmp_path / "kernel.so",
            run_rows=[{"case": "archived_v1_3", "status": "available"}],
            runs={"archived_v1_3": FakeRun()},
            step=None,
            evaluation_target_ids=(6,),
            run_shell=False,
        )
    )

    assert physics["status"] == "not_evaluated"
    assert physics["successful_target_models"] == 0
    assert all(
        expected_message in failure["message"]
        for failure in physics["failures"]
    )


def test_new_infinite_final_charge_gets_infinite_shell_reference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    archive_run = tmp_path / "archive"
    (archive_run / "input").mkdir(parents=True)
    _write_release_mechanics_fixture(archive_run)
    calls: list[tuple[object, str, int]] = []
    triangles = np.array(
        [
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]],
        ]
    )

    class FakeSnapshot:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def object_probe(self, mesh_id: int):
            return SimpleNamespace(
                mesh_id=mesh_id,
                vertex_bounding_radius_m=3.5e-5,
                target_geometry_representation="periodic2_mesh_connected",
                geometric_area_centroid_m=np.array([1.0, 2.0, 3.0]),
            )

    class FakeRun:
        result = SimpleNamespace(
            batches=280000,
            mesh_ids=np.array([6, 7]),
            triangles=triangles,
        )

        def object_interaction_snapshot(self, *, step, periodic_model, **_kwargs):
            calls.append((step, periodic_model, -1))
            return FakeSnapshot()

    def fake_shell(_snapshot, probe, displacement_m, **_kwargs):
        calls.append((None, "shell", probe.mesh_id))
        assert len(displacement_m) == 17
        assert displacement_m[-1] == pytest.approx(7.0e-5)
        return SimpleNamespace(
            image_layers=np.array([0, 1, 2, 3]),
            force_increment_error_N=np.array([1.0, 0.1, 0.01]),
            work_increment_error_J=np.array([2.0, 0.2, 0.02]),
            force_tail_proxy_N=np.array([10.0, 1.0, 0.1]),
            work_tail_proxy_J=np.array([20.0, 2.0, 0.2]),
            increment_converged=np.array([False, True, True]),
            reference_model="infinite_physical",
            reference_force_error_N=np.array([3.0, 2.0, 1.0, 0.1]),
            reference_work_error_J=np.array([6.0, 4.0, 2.0, 0.2]),
            reference_converged=np.array([False, False, True, True]),
            status="converged",
            selected_image_layers=3,
        )

    monkeypatch.setattr(tool, "finite_shell_convergence", fake_shell)
    rows, failures = tool._evaluate_new_infinite_shell_reference(
        archive_run=archive_run,
        validation_root=tmp_path / "validation",
        library=tmp_path / "kernel.so",
        run=FakeRun(),
    )

    assert calls == [
        (None, "infinite_physical", -1),
        (None, "shell", 6),
        (None, "shell", 7),
    ]
    assert failures == []
    assert len(rows) == 8
    assert {(row["mesh_id"], row["resolved_batch"]) for row in rows} == {
        (6, 280000),
        (7, 280000),
    }
    assert all(row["periodic_model"] == "infinite_physical" for row in rows)
    assert all(row["effective_far_correction"] == "cached_kneq0" for row in rows)
    assert all(row["zero_mode_policy"] == "exclude_k0" for row in rows)
    assert all(row["lower_boundary_model"] == "e_bottom_zero" for row in rows)
    assert all(float(row["model_radius_m"]) == pytest.approx(3.5e-5) for row in rows)
    assert all(float(row["geometry_radius_m"]) == pytest.approx(3.5e-5) for row in rows)
    assert all(float(row["path_start_m"]) == pytest.approx(0.0) for row in rows)
    assert all(float(row["path_end_m"]) == pytest.approx(7.0e-5) for row in rows)

    def fake_shell_without_reference(*args, **kwargs):
        shell = fake_shell(*args, **kwargs)
        shell.reference_model = None
        return shell

    calls.clear()
    monkeypatch.setattr(tool, "finite_shell_convergence", fake_shell_without_reference)
    invalid_rows, invalid_failures = tool._evaluate_new_infinite_shell_reference(
        archive_run=archive_run,
        validation_root=tmp_path / "validation",
        library=tmp_path / "kernel.so",
        run=FakeRun(),
    )
    assert invalid_rows == []
    assert len(invalid_failures) == 2
    assert all("reference_model" in failure["message"] for failure in invalid_failures)


def test_charge_vector_comparison_preserves_spatial_difference() -> None:
    tool = _load_tool()
    left = np.array([1.0, 0.0])
    right = np.array([0.0, 1.0])

    rows = tool._charge_vector_comparison_rows(
        comparison="left_to_right",
        comparison_kind="boundary_history_response",
        sample_kind="history",
        resolved_batch=1001,
        scope="object",
        mesh_id=6,
        left_snapshot_id="left:1001",
        right_snapshot_id="right:1001",
        left=left,
        right=right,
        interpretation="test",
    )
    by_metric = {row["metric"]: row for row in rows}

    assert by_metric["charge_sum_C"]["signed_difference"] == 0.0
    assert by_metric["charge_sum_C"]["absolute_difference"] == 0.0
    assert by_metric["element_charge_l1_C"]["left_value"] == 1.0
    assert by_metric["element_charge_l1_C"]["right_value"] == 1.0
    assert by_metric["element_charge_l1_C"]["absolute_difference"] == 2.0
    assert by_metric["element_charge_l2_C"]["absolute_difference"] == pytest.approx(
        np.sqrt(2.0)
    )
    assert by_metric["element_charge_linf_C"]["absolute_difference"] == 1.0


def test_geometry_comparison_allows_declared_roundoff_not_order_change() -> None:
    tool = _load_tool()
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0e-3, 0.0, 0.0], [0.0, 1.0e-3, 0.0]],
            [[0.0, 0.0, 1.0e-3], [1.0e-3, 0.0, 1.0e-3], [0.0, 1.0e-3, 1.0e-3]],
        ]
    )
    left = SimpleNamespace(mesh_ids=np.array([6, 7]), triangles=triangles)
    preliminary = tool._geometry_comparison(left, left)
    perturbed = triangles.copy()
    perturbed[0, 0, 0] += 0.5 * preliminary["coordinate_tolerance_m"]
    within = tool._geometry_comparison(
        left,
        SimpleNamespace(mesh_ids=np.array([6, 7]), triangles=perturbed),
    )
    reordered = tool._geometry_comparison(
        left,
        SimpleNamespace(mesh_ids=np.array([7, 6]), triangles=triangles),
    )

    assert within["status"] == "match"
    assert within["max_coordinate_difference_m"] > 0.0
    assert reordered["status"] == "mesh_id_order_mismatch"


def _complete_local_qualification_fixture(tool):
    paths: list[dict[str, object]] = []
    wrenches: list[dict[str, object]] = []
    for case, model, batch, mesh_id in sorted(
        tool._expected_object_evaluation_keys()
    ):
        common = {
            "case": case,
            "periodic_model": model,
            "resolved_batch": batch,
            "mesh_id": mesh_id,
            "status": "available",
        }
        paths.append({**common, "numerically_qualified": True})
        wrenches.append({**common, "component": "total_external"})
    shell_rows = [
        {
            "case": case,
            "periodic_model": "infinite_physical",
            "resolved_batch": 280000,
            "mesh_id": mesh_id,
            "status": "converged",
        }
        for case in (
            "archived_v1_3",
            "new_finite_configured",
            "new_infinite_physical",
        )
        for mesh_id in (6, 7)
    ]
    cases = {
        case: {"status": "available"}
        for case in (
            "archived_v1_3",
            "new_finite_configured",
            "new_infinite_physical",
        )
    }
    return paths, wrenches, shell_rows, cases


def test_local_qualification_reports_fixed_discretization_scope_when_qualified() -> None:
    tool = _load_tool()
    paths, wrenches, shell_rows, cases = _complete_local_qualification_fixture(tool)

    qualification = tool._local_model_numerical_qualification(
        {"status": "available"},
        wrenches,
        paths,
        shell_rows,
        cases,
    )

    assert qualification["status"] == "qualified"
    assert qualification["status_semantics"] == (
        "path_work_shell_on_fixed_saved_discretization"
    )
    assert qualification["verified_numerical_axes"] == [
        "fixed_saved_discretization_coverage",
        "path_integration",
        "work_potential_consistency",
        "force_and_barrier_decision_resolution",
        "finite_shell_to_infinite_reference",
    ]
    assert qualification["unverified_numerical_axes"] == [
        "saved_sphere_mesh_refinement",
        "source_discretization_refinement",
        "sphere_absolute_force_error",
        "sphere_absolute_torque_error",
    ]
    assert qualification["saved_sphere_mesh_refinement_status"] == "not_evaluated"
    assert qualification["source_discretization_refinement_status"] == "not_evaluated"
    assert qualification["sphere_absolute_force_error_status"] == "not_evaluated"
    assert qualification["sphere_absolute_torque_error_status"] == "not_evaluated"
    assert qualification["plane_oracle_used_as_sphere_error_bar"] is False
    assert qualification["claim_scope"] == (
        "fixed_saved_mesh_and_source_discretization; "
        "local_frozen_field_0_to_2R_path_work_shell_only; "
        "not full_discretization_or_escape_to_infinity"
    )


def test_local_qualification_rejects_duplicate_keys_masking_missing_state() -> None:
    tool = _load_tool()
    paths, wrenches, shell_rows, cases = _complete_local_qualification_fixture(tool)
    paths[-1] = dict(paths[0])
    wrenches[-1] = dict(wrenches[0])

    qualification = tool._local_model_numerical_qualification(
        {"status": "available"},
        wrenches,
        paths,
        shell_rows,
        cases,
    )

    assert qualification["status"] == "not_qualified"
    assert qualification["duplicate_path_keys"]
    assert qualification["missing_path_keys"]
    assert qualification["duplicate_wrench_keys"]
    assert qualification["missing_wrench_keys"]


def test_comparison_contract_rejects_wrong_frozen_snapshot_semantics() -> None:
    tool = _load_tool()
    snapshots = [
        {"snapshot_id": "case:history:1", "status": "available"},
        {"snapshot_id": "other:history:1", "status": "available"},
    ]
    effective = {
        "archive_version_drift": ("none", "none"),
        "frozen_field_override": ("none", "cached_kneq0"),
        "boundary_history_response_common_evaluator": (
            "cached_kneq0",
            "cached_kneq0",
        ),
        "actual_end_to_end": ("none", "cached_kneq0"),
    }
    rows = [
        {
            "comparison": kind,
            "comparison_kind": kind,
            "metric": metric,
            "left_snapshot_id": "case:history:1",
            "right_snapshot_id": "case:history:1",
            "left_effective_far_correction": pair[0],
            "right_effective_far_correction": pair[1],
            "status": "computed",
        }
        for kind, pair in effective.items()
        for metric in (
            "force_z_N",
            "endpoint_work_J",
            "minimum_available_energy_J",
            "barrier_free_from_rest",
            "endpoint_reachable_from_rest",
        )
    ]
    assert tool._comparison_artifact_contract(rows, snapshots)["status"] == "complete"
    frozen = next(
        row
        for row in rows
        if row["comparison_kind"] == "frozen_field_override"
    )
    frozen["right_snapshot_id"] = "other:history:1"

    contract = tool._comparison_artifact_contract(rows, snapshots)
    assert contract["status"] == "incomplete"
    assert "different charge snapshots" in contract["semantic_errors"][0]


def test_run_comparisons_use_canonical_comparison_kinds() -> None:
    tool = _load_tool()
    common = {
        "status": "available",
        "processed_particles": 100,
        "absorbed": 60,
        "escaped": 40,
        "total_charge_C": 1.0,
    }
    rows = tool._comparison_rows(
        [
            {"case": "archived_v1_3", **common},
            {"case": "new_finite_configured", **common},
            {"case": "new_infinite_physical", **common},
        ]
    )

    assert {row["comparison_kind"] for row in rows} == {
        "archive_version_drift",
        "actual_end_to_end",
    }


@pytest.mark.parametrize(
    ("model", "contact", "cutoff", "message"),
    [
        ("unknown", 1.0e-9, 2.0e-9, "unsupported adhesion model"),
        ("vdw_work", 2.0e-9, 1.0e-9, "cutoff_distance"),
    ],
)
def test_release_parameters_rejects_unknown_or_nonphysical_adhesion(
    tmp_path: Path,
    model: str,
    contact: float,
    cutoff: float,
    message: str,
) -> None:
    tool = _load_tool()
    archive = tmp_path / "archive"
    (archive / "input").mkdir(parents=True)
    (archive / "input/release_kernel_base.toml").write_text(
        "\n".join(
            [
                "[adhesion]",
                f'model = "{model}"',
                "hamaker_constant = 1.0e-19",
                f"contact_distance = {contact}",
                f"cutoff_distance = {cutoff}",
                "roughness_factor = 0.5",
                "contact_count = 2.0",
                "peel_factor = 1.0",
                'contact_geometry = "sphere_sphere"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match=message):
        tool._release_parameters(archive, 3.5e-5)


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


def test_complete_run_contract_rejects_finite_infinite_species_mismatch(
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
    definitions = (
        ("cache_prime", 1, False, 1),
        ("smoke_finite_configured", 100, None, None),
        ("smoke_infinite_physical", 100, True, 0),
        ("full_finite_configured_140000", 140000, None, None),
        ("full_finite_configured_280000", 280000, None, None),
        ("full_infinite_physical_140000", 140000, True, 0),
        ("full_infinite_physical_280000", 280000, True, 0),
    )
    infinite_cases = {
        "cache_prime",
        "smoke_infinite_physical",
        "full_infinite_physical_140000",
        "full_infinite_physical_280000",
    }
    for name, batches, hit, count in definitions:
        output = _write_run_output(
            validation_root,
            name,
            batches=batches,
            cache_hit=hit,
            build_count=count,
        )
        if name in infinite_cases:
            summary = (output / "summary.txt").read_text(encoding="utf-8")
            (output / "summary.txt").write_text(
                summary.replace(
                    "species_fingerprint=2222222222222222",
                    "species_fingerprint=3333333333333333",
                ),
                encoding="utf-8",
            )
        tool.verify_run(output, batches)

    with pytest.raises(tool.ValidationError, match="pair species_fingerprint"):
        tool._verify_complete_runs(validation_root, manifest)


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
                mesh_nelem=1,
                processed_particles=10,
                absorbed=6,
                escaped=4,
                batches=280000,
                escaped_boundary=3,
                survived_max_step=1,
                last_rel_change=1.0e-3,
                charges=np.array([1.0e-15]),
                mesh_ids=np.array([1]),
                triangles=np.array(
                    [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
                ),
                mesh_potential_v=None,
                charge_ledger=None,
                charge_ledger_residual_c=0.0,
                model_fingerprint="0123456789ABCDEF",
                mesh_fingerprint="1111111111111111",
                species_fingerprint="2222222222222222",
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
    assert (
        validation_root / "analysis/legacy_estimator_comparison.csv"
    ).is_file()
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


def test_analyze_cli_returns_nonzero_for_partial_physics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    monkeypatch.setattr(
        tool,
        "analyze_validation",
        lambda *_args, **_kwargs: {"physics_evaluation": {"status": "partial"}},
    )

    status = tool.main(
        [
            "analyze",
            "--archive-run",
            str(tmp_path / "archive"),
            "--validation-root",
            str(tmp_path / "validation"),
            "--library",
            str(tmp_path / "kernel.so"),
        ]
    )

    assert status == 2


def test_analyze_cli_require_complete_returns_nonzero_for_not_evaluated(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    monkeypatch.setattr(
        tool,
        "analyze_validation",
        lambda *_args, **_kwargs: {
            "cases": {
                "archived_v1_3": {"status": "available"},
                "new_finite_configured": {"status": "available"},
                "new_infinite_physical": {"status": "available"},
            },
            "physics_evaluation": {"status": "not_evaluated"},
        },
    )

    status = tool.main(
        [
            "analyze",
            "--archive-run",
            str(tmp_path / "archive"),
            "--validation-root",
            str(tmp_path / "validation"),
            "--library",
            str(tmp_path / "kernel.so"),
            "--require-complete",
        ]
    )

    assert status == 2


def test_analyze_cli_require_complete_rejects_unqualified_local_model(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tool = _load_tool()
    monkeypatch.setattr(
        tool,
        "analyze_validation",
        lambda *_args, **_kwargs: {
            "physics_evaluation": {"status": "available"},
            "numerical_qualification_for_local_frozen_model": {
                "status": "not_qualified"
            },
        },
    )

    status = tool.main(
        [
            "analyze",
            "--archive-run",
            str(tmp_path / "archive"),
            "--validation-root",
            str(tmp_path / "validation"),
            "--library",
            str(tmp_path / "kernel.so"),
            "--require-complete",
        ]
    )

    assert status == 2
