from __future__ import annotations

import os
import re
import subprocess
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools" / "run_physics_release_gate.sh"


def test_release_gate_refuses_kudpc_login_node(tmp_path: Path) -> None:
    env = os.environ.copy()
    env.pop("SLURM_JOB_ID", None)
    env["BEACH_RELEASE_GATE_HOSTNAME"] = "laurel31"
    result = subprocess.run(
        [str(SCRIPT), "--dry-run", "--manifest", str(tmp_path / "manifest.txt")],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode != 0
    assert "compute node" in result.stderr


def test_release_gate_dry_run_writes_reproducible_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.txt"
    env = os.environ.copy()
    env["BEACH_RELEASE_GATE_HOSTNAME"] = "test-compute-node"
    result = subprocess.run(
        [str(SCRIPT), "--dry-run", "--manifest", str(manifest)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    text = manifest.read_text(encoding="utf-8")
    assert "schema_version=2" in text
    assert "status=planned" in text
    assert "test_l3.command=make test-fortran-release-correctness" in text
    assert "far_correction.command=make test-fortran-far-correction" in text
    assert "mpi.command=make test-mpi" in text
    assert "mpi_cache.command=make test-mpi-periodic-cache" in text
    assert "budget.max_rss_kb=" in text
    assert "convergence_csv=" in text
    assert "test_l3.target_timings_csv=" in text
    assert "far_correction.target_timings_csv=" in text
    assert "test_l3.max_rss_kb=0" in text
    assert "git_commit=" in text


def test_release_gate_removes_stale_convergence_before_failed_run(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.txt"
    convergence = tmp_path / "convergence.csv"
    convergence.write_text("stale-success\n", encoding="utf-8")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_make = fake_bin / "make"
    fake_make.write_text("#!/usr/bin/env bash\nexit 7\n", encoding="utf-8")
    fake_make.chmod(0o755)
    env = os.environ.copy()
    env.update(
        {
            "BEACH_RELEASE_GATE_HOSTNAME": "test-compute-node",
            "FC": "/bin/true",
            "MPI_FC": "/bin/true",
            "PATH": f"{fake_bin}:{env['PATH']}",
        }
    )

    result = subprocess.run(
        [str(SCRIPT), "--manifest", str(manifest)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert not convergence.exists()


def test_release_manifest_reports_fpm_fortran_compiler(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.txt"
    fake_fc = tmp_path / "fake-fpm-fc"
    fake_fc.write_text("#!/usr/bin/env bash\necho FPM-FC-TEST-VERSION\n", encoding="utf-8")
    fake_fc.chmod(0o755)
    env = os.environ.copy()
    env.update(
        {
            "BEACH_RELEASE_GATE_HOSTNAME": "test-compute-node",
            "FC": "/bin/false",
            "FPM_FC": str(fake_fc),
            "MPI_FC": "/bin/true",
        }
    )

    result = subprocess.run(
        [str(SCRIPT), "--dry-run", "--manifest", str(manifest)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "fortran_compiler=FPM-FC-TEST-VERSION" in manifest.read_text(encoding="utf-8")


def test_portable_physics_contract_workflow_runs_l2() -> None:
    workflow = (ROOT / ".github" / "workflows" / "physics-contracts.yml").read_text(encoding="utf-8")
    assert "fortran-lang/setup-fpm@v10" in workflow
    assert "MPI_FC: mpifort" in workflow
    assert "libopenmpi-dev openmpi-bin" in workflow
    assert "make test-l2" in workflow
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    assert "test-physics-release:" in makefile
    assert "tools/run_physics_release_gate.sh" in makefile


def test_release_mpi_build_defaults_to_ifx() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    release_gate = SCRIPT.read_text(encoding="utf-8")
    camphor_profile = (ROOT / "env" / "camphor.env").read_text(encoding="utf-8")

    assert re.search(r"^MPI_FC \?= mpiifx$", makefile, flags=re.MULTILINE)
    assert 'mpi_fc_command="${MPI_FC:-mpiifx}"' in release_gate
    assert ': "${FC:=mpiifx}"' in camphor_profile


def test_release_gate_uses_minimal_correctness_subset() -> None:
    script = SCRIPT.read_text(encoding="utf-8")
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")

    assert "run_gate test_l3 make test-l3" not in script
    assert "test-fortran-release-correctness" in script
    assert "FORTRAN_RELEASE_CONVERGENCE_TARGETS" in makefile
    assert "test-fortran-release-correctness:" in makefile


def test_far_correction_diagnostic_and_benchmark_are_opt_in() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    match = re.search(
        r"FORTRAN_FAR_CORRECTION_TARGETS \?=(.*?)(?=^FORTRAN_[A-Z_]+_TARGETS \?=)",
        makefile,
        flags=re.MULTILINE | re.DOTALL,
    )

    assert match is not None
    assert "test_periodic2_flat_oracle_diag" not in match.group(1)
    assert "FORTRAN_FAR_CORRECTION_DIAGNOSTIC_TARGETS" in makefile
    assert "test-fortran-far-correction-diagnostics:" in makefile
    assert "FORTRAN_BENCHMARK_TARGETS" in makefile
    assert "test-fortran-benchmark:" in makefile
    assert "$(call run_fortran_targets,$(FORTRAN_BENCHMARK_TARGETS),release,run,example)" in makefile


def test_runtime_loop_is_not_part_of_debug_correctness_test() -> None:
    correctness = (ROOT / "tests" / "fortran" / "test_periodic2_infinite_operator.f90").read_text(encoding="utf-8")
    benchmark = ROOT / "benchmarks" / "fortran" / "benchmark_periodic2_runtime.f90"
    with (ROOT / "fpm.toml").open("rb") as handle:
        manifest = tomllib.load(handle)

    assert "cached_runtime_cost_vs_explicit_shell" not in correctness
    assert benchmark.exists()
    assert "benchmark_periodic2_runtime" in {item["name"] for item in manifest["example"]}
    assert "benchmark_periodic2_runtime" not in {item["name"] for item in manifest["executable"]}
    assert "benchmark_periodic2_runtime" not in {item["name"] for item in manifest["test"]}


def test_operator_cache_does_not_build_a_probe_plan() -> None:
    source = (ROOT / "tests" / "fortran" / "test_periodic2_operator_cache.f90").read_text(encoding="utf-8")

    assert "probe_plan" not in source
    assert "system_clock" in source


def test_release_gate_accepts_prefixed_convergence_markers() -> None:
    script = SCRIPT.read_text(encoding="utf-8")
    assert "s/^.*BEACH_CONVERGENCE,//p" in script
    assert "grep -h '^BEACH_CONVERGENCE,'" not in script
