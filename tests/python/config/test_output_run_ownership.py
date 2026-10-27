"""Run the freshly built ``beach`` twice in one output directory (issues #79, #80).

A fresh run must not leave the previous run's checkpoints or histories behind, and a
resume from an older checkpoint must rewrite every history so the strict readers accept
it. ``make test-config-contract`` supplies the executable; each run tracks one electron
for a few batches.
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess

import numpy as np
import tomli_w

from beach.config._toml import load_toml_file
from beach.io import load_run

ROOT = Path(__file__).resolve().parents[3]


def _run(beach: Path, tmp_path: Path, output: Path, **settings) -> str:
    config = load_toml_file(ROOT / "examples/smoke_single_electron.toml")
    config["run"]["batch_count"] = settings["batch_count"]
    if settings.get("resume"):
        config["run"]["restart"] = {}
    config["output"]["dir"] = str(output)
    config["output"]["history"] = {"stride_batches": settings["history_stride"]}
    config["output"]["checkpoint"] = {"stride_batches": settings["checkpoint_stride"]}
    path = tmp_path / "beach.toml"
    path.write_text(tomli_w.dumps(config), encoding="utf-8")
    env = dict(os.environ, OMP_NUM_THREADS="1", BEACH_CACHE_DIR=str(tmp_path / "cache"))
    result = subprocess.run(
        [str(beach), str(path)], cwd=tmp_path, env=env,
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def _complete(directory: Path) -> bool:
    manifest = directory / "checkpoint_complete.txt"
    return manifest.is_file() and "state=complete" in manifest.read_text(encoding="utf-8")


def test_fresh_run_owns_a_reused_output_directory(native_checker: Path, tmp_path: Path) -> None:
    output = tmp_path / "out"
    _run(native_checker, tmp_path, output, batch_count=3, history_stride=1, checkpoint_stride=1)
    assert (output / "charge_history.csv").is_file()
    assert _complete(output / "checkpoints" / "slot0")

    _run(native_checker, tmp_path, output, batch_count=2, history_stride=0, checkpoint_stride=0)
    result = load_run(output)
    assert result.summary.batches == 2
    assert result.history.charge is None, "run A's charge history must not survive run B"
    assert not (output / "checkpoint_latest.txt").exists()
    for slot in ("slot0", "slot1"):
        assert not _complete(output / "checkpoints" / slot)

    stdout = _run(native_checker, tmp_path, output, batch_count=3, history_stride=0,
                  checkpoint_stride=0, resume=True)
    assert "resuming_from_batches=2" in stdout, "resume must continue run B, not run A's slot"


def test_resume_from_an_older_checkpoint_rewrites_histories(
    native_checker: Path, tmp_path: Path,
) -> None:
    output = tmp_path / "out"
    _run(native_checker, tmp_path, output, batch_count=4, history_stride=1, checkpoint_stride=2)
    reference = load_run(output).history.charge
    np.testing.assert_array_equal(reference.batch_indices, [1, 2, 3, 4])
    final_charge = reference[4]

    # Crash after batch 4: only the batch-2 slot is complete and the history ends mid-row.
    for directory in (output, output / "checkpoints" / "slot1"):
        (directory / "checkpoint_complete.txt").unlink()
    with (output / "charge_history.csv").open("a", encoding="utf-8") as stream:
        stream.write("5,6,0.0,1")

    stdout = _run(native_checker, tmp_path, output, batch_count=4, history_stride=1,
                  checkpoint_stride=2, resume=True)
    assert "resuming_from_batches=2" in stdout
    history = load_run(output).history.charge
    np.testing.assert_array_equal(history.batch_indices, [1, 2, 3, 4])
    np.testing.assert_array_equal(history[4], final_charge)
