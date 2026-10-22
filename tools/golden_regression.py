#!/usr/bin/env python3
"""Run reduced example cases and compare their outputs before and after a refactor.

The cases are small copies of examples/ that run in seconds to a minute with one
OpenMP thread and a fixed seed, so a behaviour-preserving change reproduces every
output file bit for bit.

    tools/golden_regression.py run OUT_DIR          # build release in this checkout, run all cases
    tools/golden_regression.py compare BASE NEW     # bitwise, or --rtol for a relative tolerance

Run on a compute node; ``run`` calls ``./build.sh`` with ``FPM_PROFILE=release``.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib
import tomli_w

ROOT = Path(__file__).resolve().parents[1]

# Lines carrying build identity or timing differ between builds and are not compared.
VOLATILE_KEYS = ("build_", "version", "elapsed", "time_s", "wall", "timestamp", "date", "host", "git_", "perf")


def _output(cfg: dict[str, Any], out_dir: Path, case: str, potential: bool = False) -> None:
    history: dict[str, Any] = {"stride_batches": 1}
    if potential:
        history["potential"] = True
    cfg["output"] = {"enabled": True, "dir": str(out_dir / case), "history": history}


def _shrink_sampling(cfg: dict[str, Any], macros: int) -> None:
    cfg["particles"]["tracking"]["max_steps_per_particle"] = 20000
    for species in cfg["particles"]["species"]:
        sampling = species.get("sampling", {})
        if sampling.get("target_macro_particles_per_batch", 0) > 0:
            sampling["target_macro_particles_per_batch"] = macros
        if "rays_per_batch" in sampling:
            sampling["rays_per_batch"] = macros


def case_configs(out_dir: Path) -> dict[str, dict[str, Any]]:
    """Return the reduced cases: free-space direct, periodic2 FMM with photoelectrons, zero-current refresh."""

    def load(name: str) -> dict[str, Any]:
        return tomllib.loads((ROOT / "examples" / name).read_text())

    cases: dict[str, dict[str, Any]] = {}

    cfg = load("panel_direct.toml")
    cfg["run"]["batch_count"] = 3
    cfg["run"]["rng_seed"] = 4242
    _output(cfg, out_dir, "g1_free_direct")
    cases["g1_free_direct"] = cfg

    cfg = load("closed_photoelectron.toml")
    cfg["run"]["batch_count"] = 3
    cfg["mesh"]["parts"][0].update(nx=6, ny=6)
    cfg["mesh"]["parts"][1].update(n_lon=8, n_lat=4)
    _shrink_sampling(cfg, 60)
    _output(cfg, out_dir, "g2_periodic_photo", potential=True)
    cases["g2_periodic_photo"] = cfg

    cfg = load("zero_current_refresh.toml")
    cfg["run"]["batch_count"] = 4
    cfg["mesh"]["parts"][0].update(nx=4, ny=4)
    _shrink_sampling(cfg, 80)
    _output(cfg, out_dir, "g3_zero_current_refresh")
    cases["g3_zero_current_refresh"] = cfg
    return cases


def run_cases(out_dir: Path) -> int:
    out_dir = out_dir.resolve()
    config_dir = out_dir / "cfg"
    config_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ | {"OMP_NUM_THREADS": "1", "FPM_PROFILE": "release"}
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout
    print(f"golden run head={head.strip()} dirty_files={len(dirty.splitlines())} out={out_dir}")
    with (out_dir / "build.log").open("w") as log:
        built = subprocess.run(["./build.sh"], cwd=ROOT, env=env | {"FPM_ACTION": "build"}, stdout=log, stderr=log)
    if built.returncode != 0:
        print(f"BUILD FAILED, see {out_dir / 'build.log'}")
        return 1
    status = 0
    for name, cfg in case_configs(out_dir).items():
        path = config_dir / f"{name}.toml"
        path.write_text(tomli_w.dumps(cfg))
        start = time.monotonic()
        with (out_dir / f"{name}.log").open("w") as log:
            result = subprocess.run(
                ["./build.sh", "--target", "beach", "--", str(path)],
                cwd=ROOT, env=env | {"FPM_ACTION": "run"}, stdout=log, stderr=log,
            )
        state = "PASS" if result.returncode == 0 else "FAIL"
        status |= result.returncode != 0
        print(f"{state} {name} {time.monotonic() - start:.0f}s")
    return status


def _tokens(line: str) -> list[float | str]:
    values: list[float | str] = []
    for token in line.replace(",", " ").replace("=", " ").split():
        try:
            values.append(float(token.replace("D", "E")))
        except ValueError:
            values.append(token)
    return values


def _stable_lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if not any(key in line.lower() for key in VOLATILE_KEYS)]


def compare_trees(base: Path, new: Path, rtol: float) -> int:
    """Print one line per case; return nonzero when a file differs beyond ``rtol`` (0 means bitwise)."""

    status = 0
    worst_overall = 0.0
    for case in sorted(p for p in base.iterdir() if p.is_dir() and p.name.startswith("g")):
        files = sorted(f for f in case.rglob("*") if f.is_file() and f.suffix in (".csv", ".txt"))
        notes: list[str] = []
        for f in files:
            other = new / case.name / f.relative_to(case)
            if not other.exists():
                notes.append(f"MISSING {f.name}")
                status = 1
                continue
            lines_a, lines_b = _stable_lines(f), _stable_lines(other)
            if lines_a == lines_b:
                continue
            if len(lines_a) != len(lines_b):
                notes.append(f"{f.name}: line count {len(lines_a)} vs {len(lines_b)}")
                status = 1
                continue
            worst, text_diffs = 0.0, 0
            for line_a, line_b in zip(lines_a, lines_b):
                if line_a == line_b:
                    continue
                values_a, values_b = _tokens(line_a), _tokens(line_b)
                if len(values_a) != len(values_b):
                    text_diffs += 1
                    continue
                for a, b in zip(values_a, values_b):
                    if isinstance(a, float) and isinstance(b, float):
                        if a != b:
                            worst = max(worst, abs(a - b) / max(abs(a), abs(b), 1.0e-300))
                    elif a != b:
                        text_diffs += 1
            worst_overall = max(worst_overall, worst)
            if text_diffs or worst > rtol:
                status = 1
            notes.append(f"{f.name}: max_rel={worst:.3e} text_diffs={text_diffs}")
        print(f"[{case.name}] files={len(files)} " + ("BITWISE IDENTICAL" if not notes else "; ".join(notes)))
    print(f"overall max relative difference: {worst_overall:.3e} ({'OK' if status == 0 else 'DIFFERS'} at rtol={rtol:g})")
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="build this checkout and run every golden case into OUT_DIR")
    run.add_argument("out_dir", type=Path)
    compare = sub.add_parser("compare", help="compare two OUT_DIR trees")
    compare.add_argument("base", type=Path)
    compare.add_argument("new", type=Path)
    compare.add_argument("--rtol", type=float, default=0.0, help="relative tolerance (default 0: bitwise)")
    args = parser.parse_args(argv)
    if args.command == "run":
        return run_cases(args.out_dir)
    return compare_trees(args.base, args.new, args.rtol)


if __name__ == "__main__":
    sys.exit(main())
