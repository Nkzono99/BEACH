from __future__ import annotations

import csv
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "tools" / "run_fortran_targets.sh"


def _write_fake_build(path: Path) -> None:
    path.write_text(
        "#!/usr/bin/env bash\n"
        "set -eu\n"
        "target=\"${!#}\"\n"
        "printf '%s\\n' \"$target\" >>\"$FAKE_BUILD_LOG\"\n"
        "[[ \"$target\" != \"$FAKE_FAIL_TARGET\" ]]\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _run(
    tmp_path: Path,
    *targets: str,
    fail_target: str = "never",
    action: str = "test",
    target_kind: str = "target",
) -> subprocess.CompletedProcess[str]:
    fake_build = tmp_path / "fake-build.sh"
    _write_fake_build(fake_build)
    env = os.environ.copy()
    env.update(
        {
            "BUILD_SH": str(fake_build),
            "BEACH_FORTRAN_TIMING_CSV": str(tmp_path / "timings.csv"),
            "FAKE_BUILD_LOG": str(tmp_path / "order.log"),
            "FAKE_FAIL_TARGET": fail_target,
            "FPM_PROFILE": "debug",
            "FPM_ACTION": action,
            "FPM_TARGET_KIND": target_kind,
        }
    )
    return subprocess.run(
        [str(RUNNER), *targets],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def test_runner_records_targets_in_sequential_order(tmp_path: Path) -> None:
    result = _run(tmp_path, "alpha", "beta")

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "order.log").read_text(encoding="utf-8").splitlines() == ["alpha", "beta"]
    with (tmp_path / "timings.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["target"] for row in rows] == ["alpha", "beta"]
    assert [row["profile"] for row in rows] == ["debug", "debug"]
    assert [row["status"] for row in rows] == ["passed", "passed"]
    assert all(re.fullmatch(r"\d+\.\d{3}", row["elapsed_seconds"]) for row in rows)
    assert all(float(row["elapsed_seconds"]) >= 0.0 for row in rows)


def test_runner_records_failure_and_stops(tmp_path: Path) -> None:
    result = _run(tmp_path, "alpha", "broken", "not-run", fail_target="broken")

    assert result.returncode != 0
    assert (tmp_path / "order.log").read_text(encoding="utf-8").splitlines() == ["alpha", "broken"]
    with (tmp_path / "timings.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["status"] for row in rows] == ["passed", "failed"]


def test_runner_reports_requested_fpm_action(tmp_path: Path) -> None:
    result = _run(tmp_path, "benchmark", action="run", target_kind="example")

    assert result.returncode == 0, result.stderr
    assert "fpm run --example --target benchmark" in result.stdout
