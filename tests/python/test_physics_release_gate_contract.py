from __future__ import annotations

import os
import subprocess
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
    assert "test_l3.command=make test-l3" in text
    assert "far_correction.command=make test-fortran-far-correction" in text
    assert "mpi.command=make test-mpi" in text
    assert "mpi_cache.command=make test-mpi-periodic-cache" in text
    assert "budget.max_rss_kb=" in text
    assert "convergence_csv=" in text
    assert "test_l3.max_rss_kb=0" in text
    assert "git_commit=" in text


def test_portable_physics_contract_workflow_runs_l2() -> None:
    workflow = (ROOT / ".github" / "workflows" / "physics-contracts.yml").read_text(encoding="utf-8")
    assert "fortran-lang/setup-fpm@v10" in workflow
    assert "make test-l2" in workflow
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    assert "test-physics-release:" in makefile
    assert "tools/run_physics_release_gate.sh" in makefile


def test_release_gate_accepts_prefixed_convergence_markers() -> None:
    script = SCRIPT.read_text(encoding="utf-8")
    assert "s/^.*BEACH_CONVERGENCE,//p" in script
    assert "grep -h '^BEACH_CONVERGENCE,'" not in script
