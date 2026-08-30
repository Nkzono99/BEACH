from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "tools" / "run_fortran_targets.sh"


def _write_fake_build(path: Path) -> None:
    path.write_text(
        "#!/usr/bin/env bash\n"
        "set -eu\n"
        "target=\"${!#}\"\n"
        "printf '%s\\t%s\\t%s\\n' \"$FPM_ACTION\" \"$FPM_PROFILE\" \"$*\" "
        ">>\"$FAKE_BUILD_LOG\"\n"
        "if [[ \"$target\" == \"$FAKE_FAIL_TARGET\" ]]; then exit 7; fi\n"
        "exit 0\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _run(
    tmp_path: Path,
    *targets: str,
    fail_target: str = "never",
    action: str | None = "test",
    target_kind: str = "target",
    profile: str | None = "debug",
    write_timings: bool = True,
) -> subprocess.CompletedProcess[str]:
    fake_build = tmp_path / "fake-build.sh"
    _write_fake_build(fake_build)
    env = os.environ.copy()
    env.pop("BEACH_FORTRAN_TIMING_CSV", None)
    env.pop("FPM_ACTION", None)
    env.pop("FPM_PROFILE", None)
    env.update(
        {
            "BUILD_SH": str(fake_build),
            "FAKE_BUILD_LOG": str(tmp_path / "build.log"),
            "FAKE_FAIL_TARGET": fail_target,
            "FPM_TARGET_KIND": target_kind,
        }
    )
    if action is not None:
        env["FPM_ACTION"] = action
    if profile is not None:
        env["FPM_PROFILE"] = profile
    if write_timings:
        env["BEACH_FORTRAN_TIMING_CSV"] = str(tmp_path / "timings.csv")
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
    assert (tmp_path / "build.log").read_text(encoding="utf-8").splitlines() == [
        "test\tdebug\t--target alpha",
        "test\tdebug\t--target beta",
    ]
    with (tmp_path / "timings.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["target"] for row in rows] == ["alpha", "beta"]
    assert [row["profile"] for row in rows] == ["debug", "debug"]
    assert [row["status"] for row in rows] == ["passed", "passed"]
    assert all(float(row["elapsed_seconds"]) >= 0.0 for row in rows)


def test_runner_passes_its_defaults_to_the_build(tmp_path: Path) -> None:
    result = _run(tmp_path, "alpha", action=None, profile=None)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "build.log").read_text(encoding="utf-8").splitlines() == [
        "test\tdebug\t--target alpha"
    ]


def test_runner_records_failure_and_stops(tmp_path: Path) -> None:
    result = _run(tmp_path, "alpha", "broken", "not-run", fail_target="broken")

    assert result.returncode == 7
    assert (tmp_path / "build.log").read_text(encoding="utf-8").splitlines() == [
        "test\tdebug\t--target alpha",
        "test\tdebug\t--target broken",
    ]
    with (tmp_path / "timings.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [(row["target"], row["status"]) for row in rows] == [
        ("alpha", "passed"),
        ("broken", "failed"),
    ]


def test_runner_forwards_run_example_without_requiring_timings(tmp_path: Path) -> None:
    result = _run(
        tmp_path,
        "benchmark",
        action="run",
        target_kind="example",
        profile="release",
        write_timings=False,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "build.log").read_text(encoding="utf-8").splitlines() == [
        "run\trelease\t--example --target benchmark"
    ]
    assert not (tmp_path / "timings.csv").exists()


def test_runner_rejects_example_for_test_action_before_build(tmp_path: Path) -> None:
    result = _run(
        tmp_path,
        "alpha",
        action="test",
        target_kind="example",
    )

    assert result.returncode == 2
    assert "example targets require FPM_ACTION=run" in result.stderr
    assert not (tmp_path / "build.log").exists()
