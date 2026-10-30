"""Stage, verify-inputs, and submit contracts of the periodic validation tool."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import tomli_w
from validation_support import (
    ROOT,
    TOOL_PATH,
    _load_toml,
    _load_tool,
    _write_run_output,
)

pytestmark = pytest.mark.workflow


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


@pytest.mark.smoke
def test_stage_creates_fresh_canonical_and_immutable_restart_inputs(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"

    manifest = tool.stage_validation(archive_run, validation_root, binary)

    import tomllib
    grouped_input = tomllib.loads((validation_root / "input/full/finite_configured.toml").read_text())
    assert 'sim' not in grouped_input
    assert grouped_input['fields']['periodic']['backend'] == 'finite_images'
    finite = _load_toml(validation_root / "input/full/finite_configured.toml")
    infinite = _load_toml(validation_root / "input/full/infinite_physical.toml")
    finite_140 = _load_toml(
        validation_root / "input/full/segments/finite_configured_140000.toml"
    )
    finite_280 = _load_toml(
        validation_root / "input/full/segments/finite_configured_280000.toml"
    )
    infinite_280 = _load_toml(
        validation_root / "input/full/segments/infinite_physical_280000.toml"
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
    assert finite_280["output"]["dir"].endswith("/run/full/finite_configured/280000")
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
    assert 'srun "${BINARY}"' in smoke
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
    finite_restart = (validation_root / "submit/full_finite_280000_sysa.sh").read_text(
        encoding="utf-8"
    )
    assert "full/finite_configured/140000" in finite_restart
    assert '--expected-batches "${restart_batches}"' in finite_restart
    assert '"140000"' in finite_restart
    assert '--producer-job-role "finite_280000"' in finite_restart
    parent_verify = next(
        line
        for line in finite_restart.splitlines()
        if "--require-existing-receipt" in line
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


def _write_toml(path: Path, data: dict[str, object]) -> None:
    from beach.config import to_grouped_layout
    if 'sim' in data and 'domain' in data:
        data = to_grouped_layout(data)
    path.write_text(tomli_w.dumps(data), encoding="utf-8")


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
            lambda data: data["sim"].__setitem__("field_periodic_ewald_layers", 5),
        ),
        (
            "input/full/finite_configured.toml",
            lambda data: data["mesh"]["parts"][0].__setitem__("radius", 9.0e-5),
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


@pytest.mark.smoke
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

    assert (
        tool.verify_inputs(validation_root, require_empty_outputs=False)["status"]
        == "ok"
    )
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


@pytest.mark.parametrize(
    "removed_setting",
    [
        "softening = 1.0e-10",
        'field_source_model = "point"',
        'field_kernel_id = "softened_point"',
    ],
)
def test_stage_rejects_removed_point_source_configuration(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    removed_setting: str,
) -> None:
    tool = _load_tool()
    config_path = archive_run / "input/beach.toml"
    config = config_path.read_text(encoding="utf-8")
    config_path.write_text(
        config.replace("[sim]\n", f"[sim]\n{removed_setting}\n", 1),
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match="removed point-source"):
        tool.stage_validation(archive_run, tmp_path / "validation", binary)


def test_stage_rejects_removed_point_element_kernel_table(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
) -> None:
    tool = _load_tool()
    config_path = archive_run / "input/beach.toml"
    config_path.write_text(
        config_path.read_text(encoding="utf-8")
        + '\n[field]\nelement_kernel = "point"\n',
        encoding="utf-8",
    )

    with pytest.raises(tool.ValidationError, match=r"removed point-source.*\[field\]"):
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
    assert (validation_root / "cache/periodic2/operator.bin").read_bytes() == b"stale"


@pytest.mark.parametrize("root_exists", [False, True])
def test_failed_stage_attempt_leaves_root_ready_for_retry(
    archive_run: Path,
    binary: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    root_exists: bool,
) -> None:
    tool = _load_tool()
    validation_root = tmp_path / "validation"
    if root_exists:
        validation_root.mkdir()

    # Job scripts are written after the configs and the read-only source snapshot.
    def fail_after_partial_staging(*_args, **_kwargs):
        raise tool.ValidationError("injected job-script failure")

    with monkeypatch.context() as patch:
        patch.setattr(tool, "_write_jobs", fail_after_partial_staging)
        with pytest.raises(tool.ValidationError, match="injected"):
            tool.stage_validation(archive_run, validation_root, binary)

    assert validation_root.exists() is root_exists
    assert not root_exists or not any(validation_root.iterdir())
    manifest = tool.stage_validation(archive_run, validation_root, binary)
    assert (validation_root / "manifest.json").is_file()
    assert manifest["validation_root"] == str(validation_root.resolve())


@pytest.mark.parametrize("location", ["repository", "archive"])
def test_stage_rejects_validation_root_inside_source_or_archive(
    archive_run: Path,
    binary: Path,
    location: str,
) -> None:
    tool = _load_tool()
    validation_root = ROOT if location == "repository" else archive_run

    with pytest.raises(
        tool.ValidationError, match="outside the repository and archive"
    ):
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
            validation_root / "input/smoke" / ".." / "smoke/finite_configured.toml"
        )
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    elif mutation == "config_double_slash":
        case = manifest["cases"]["smoke_finite_configured"]
        case["config_path"] = str(case["config_path"]).replace("/input/", "//input/", 1)
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

    with pytest.raises(
        tool.ValidationError, match="archive analysis.*canonical|symlink"
    ):
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
    assert (
        oracle_command
        + ' \\\n  --validation-root "${VALIDATION_ROOT}" \\\n  --library "${LIBRARY}"'
        in analysis_text
    )
    smoke_text = (validation_root / "submit/smoke_sysa.sh").read_text(encoding="utf-8")
    assert "probe-library" in smoke_text
    assert str(staged) in smoke_text
    subprocess.run(["bash", "-n", analysis_job], check=True)
    chain = (validation_root / "submit/submit_chain.sh").read_text(encoding="utf-8")
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
    [
        "validation root",
        "validation@root",
        "validation$root",
        "validation'root",
        "validation\nroot",
    ],
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


@pytest.mark.smoke
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
        "module": '#!/bin/bash\nif [[ " $* " == *" list "* ]]; then echo SysA/2022; fi\n',
        "spartition": "#!/bin/bash\nprintf 'Partition State\\ngr20001a UP\\n'\n",
        "qgroup": "#!/bin/bash\nexit 0\n",
        "python3.11": "#!/bin/bash\nexit 0\n",
        "sbatch": (
            "#!/bin/bash\n"
            "value=0\n"
            '[ ! -f "${SBATCH_COUNTER}" ] || value=$(cat "${SBATCH_COUNTER}")\n'
            "value=$((value + 1))\n"
            'printf \'%s\\n\' "${value}" > "${SBATCH_COUNTER}"\n'
            'if [ "${value}" = "${SBATCH_FAIL_AT:-0}" ]; then exit 1; fi\n'
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
