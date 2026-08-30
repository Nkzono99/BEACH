from __future__ import annotations

import os
import subprocess
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools" / "run_physics_release_gate.sh"


def release_gate_env(**overrides: str) -> dict[str, str]:
    env = os.environ.copy()
    for name in (
        "SLURM_JOB_ID", "BEACH_RELEASE_MPI_RUNNER",
        "BEACH_RELEASE_MAX_RSS_KB", "BEACH_RELEASE_TIME_COMMAND",
    ):
        env.pop(name, None)
    env.update({
        "BEACH_RELEASE_GATE_HOSTNAME": "test-compute-node",
        "FPM_FC": "/bin/true", "MPI_FC": "/bin/true",
    })
    env.update(overrides)
    return env


def read_manifest(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def run_fake_gate(
    tmp_path: Path, make_body: str, *, budget_kb: int | None = None,
) -> subprocess.CompletedProcess[str]:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_make = fake_bin / "make"
    fake_make.write_text(f"#!/usr/bin/env bash\n{make_body}", encoding="utf-8")
    fake_make.chmod(0o755)
    env = release_gate_env(BEACH_FAKE_MAKE_LOG=str(tmp_path / "make.log"))
    if budget_kb is not None:
        env["BEACH_RELEASE_MAX_RSS_KB"] = str(budget_kb)
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    return subprocess.run(
        [str(SCRIPT), "--manifest", str(tmp_path / "manifest.txt")],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def test_release_gate_refuses_kudpc_login_node(tmp_path: Path) -> None:
    result = subprocess.run(
        [str(SCRIPT), "--dry-run", "--manifest", str(tmp_path / "manifest.txt")],
        cwd=ROOT,
        env=release_gate_env(BEACH_RELEASE_GATE_HOSTNAME="laurel31"),
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode != 0
    assert "compute node" in result.stderr


def test_release_gate_dry_run_writes_reproducible_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.txt"
    fake_fc = tmp_path / "fake-fpm-fc"
    fake_fc.write_text("#!/usr/bin/env bash\necho FPM-FC-TEST-VERSION\n", encoding="utf-8")
    fake_fc.chmod(0o755)
    result = subprocess.run(
        [str(SCRIPT), "--dry-run", "--manifest", str(manifest)],
        cwd=ROOT,
        env=release_gate_env(FC="/bin/false", FPM_FC=str(fake_fc)),
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    values = read_manifest(manifest)
    expected = {
        "schema_version": "2",
        "status": "planned",
        "fortran_compiler": "FPM-FC-TEST-VERSION",
        "test_l3.command": "make test-fortran-release-correctness",
        "far_correction.command": "make test-fortran-far-correction",
        "mpi.command": "make test-mpi",
        "mpi_cache.command": "make test-mpi-periodic-cache",
        "convergence.required_categories": "boris_dt,panel_fmm_order,rough_panel_mesh",
        "test_l3.max_rss_kb": "0",
    }
    assert {key: values[key] for key in expected} == expected
    assert int(values["budget.max_rss_kb"]) > 0
    assert Path(values["convergence_csv"]) == tmp_path / "convergence.csv"
    assert Path(values["test_l3.target_timings_csv"]) == tmp_path / "test_l3-target-timings.csv"
    assert Path(values["far_correction.target_timings_csv"]) == tmp_path / "far_correction-target-timings.csv"
    assert values["git_commit"]


def test_release_gate_removes_stale_artifacts_and_enforces_memory_budget(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.txt"
    stale_artifacts = [
        tmp_path / "convergence.csv",
        tmp_path / "test_l3-target-timings.csv",
        tmp_path / "far_correction-target-timings.csv",
    ]
    for artifact in stale_artifacts:
        artifact.write_text("stale-success\n", encoding="utf-8")
    result = run_fake_gate(tmp_path, "exit 0\n", budget_kb=1)

    assert result.returncode != 0
    assert all(not artifact.exists() for artifact in stale_artifacts)
    values = read_manifest(manifest)
    assert values["test_l3.status"] == "failed_memory_budget"
    assert int(values["test_l3.max_rss_kb"]) > 1
    assert values["status"] == "failed"


def test_release_gate_propagates_failed_make_to_manifest(tmp_path: Path) -> None:
    result = run_fake_gate(tmp_path, "exit 7\n", budget_kb=8_388_608)

    assert result.returncode != 0
    values = read_manifest(tmp_path / "manifest.txt")
    assert values["test_l3.status"] == "failed"
    assert values["status"] == "failed"


def test_release_gate_collects_prefixed_convergence_and_runs_each_gate(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.txt"
    make_log = tmp_path / "make.log"
    result = run_fake_gate(
        tmp_path,
        """
printf '%s\n' "$*" >>"$BEACH_FAKE_MAKE_LOG"
if [[ "$1" == "test-fortran-release-correctness" ]]; then
  printf '[1/3] dt ... BEACH_CONVERGENCE,boris_dt,fixture,1,2,pass\n'
  printf '[2/3] order ... BEACH_CONVERGENCE,panel_fmm_order,fixture,1,2,pass\n'
  printf '[3/3] mesh ... BEACH_CONVERGENCE,rough_panel_mesh,fixture,1,2,pass\n'
fi
""",
    )

    assert result.returncode == 0, result.stderr
    assert make_log.read_text(encoding="utf-8").splitlines() == [
        "test-fortran-release-correctness",
        "test-fortran-far-correction",
        "test-mpi",
        "test-mpi-periodic-cache",
    ]
    assert (tmp_path / "convergence.csv").read_text(encoding="utf-8").splitlines() == [
        "category,configuration,metric_1,metric_2,metric_3,acceptance",
        "boris_dt,fixture,1,2,pass",
        "panel_fmm_order,fixture,1,2,pass",
        "rough_panel_mesh,fixture,1,2,pass",
    ]
    values = read_manifest(manifest)
    assert values["status"] == "passed"
    assert values["convergence.status"] == "passed"
    assert values["convergence.rows"] == "3"
    assert all(values[f"{gate}.status"] == "passed" for gate in ("test_l3", "far_correction", "mpi", "mpi_cache"))


def test_portable_physics_contract_workflow_runs_l2() -> None:
    workflow = (ROOT / ".github" / "workflows" / "physics-contracts.yml").read_text(encoding="utf-8")
    assert "uses: fortran-lang/setup-fpm@" in workflow
    assert "MPI_FC: mpifort" in workflow
    assert "libopenmpi-dev openmpi-bin" in workflow
    assert "make test-l2" in workflow


def test_release_mpi_build_defaults_to_ifx() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    release_gate = SCRIPT.read_text(encoding="utf-8")
    camphor_profile = (ROOT / "env" / "camphor.env").read_text(encoding="utf-8")

    assert "MPI_FC ?= mpiifx" in makefile.splitlines()
    assert 'mpi_fc_command="${MPI_FC:-mpiifx}"' in release_gate
    assert ': "${FC:=mpiifx}"' in camphor_profile


def test_makefile_separates_release_correctness_from_benchmark() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    with (ROOT / "fpm.toml").open("rb") as handle:
        manifest = tomllib.load(handle)
    benchmark = "benchmark_periodic2_runtime"
    l3_targets = makefile.split("FORTRAN_L3_TARGETS ?=", 1)[1].split(
        "FORTRAN_RELEASE_CONVERGENCE_TARGETS ?=", 1
    )[0]
    far_targets = makefile.split("FORTRAN_FAR_CORRECTION_TARGETS ?=", 1)[1].split(
        "FORTRAN_BENCHMARK_TARGETS ?=", 1
    )[0]

    assert (
        "$(call run_fortran_targets,$(FORTRAN_RELEASE_CONVERGENCE_TARGETS) "
        "$(FORTRAN_L3_TARGETS),debug,test)"
    ) in makefile
    assert "$(call run_fortran_targets,$(FORTRAN_FAR_CORRECTION_TARGETS),debug,test)" in makefile
    assert "$(call run_fortran_targets,$(FORTRAN_BENCHMARK_TARGETS),release,run,example)" in makefile
    assert benchmark not in l3_targets
    assert benchmark not in far_targets
    assert benchmark in {item["name"] for item in manifest["example"]}
    assert benchmark not in {item["name"] for item in manifest["test"]}
    assert benchmark not in {item["name"] for item in manifest["executable"]}
    assert 'tools/run_physics_release_gate.sh --manifest "$(PHYSICS_RELEASE_MANIFEST)"' in makefile
