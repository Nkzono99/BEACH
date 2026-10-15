"""runops (>= 2.0) adapter for BEACH, published as the ``beach`` entry point.

runops discovers :class:`Adapter` through the ``runops.adapters`` entry point
group declared in ``pyproject.toml``. The adapter is only called from inside
runops; runops helpers are imported lazily.

Run directory conventions (the job's working directory is the run directory):

* ``input/beach.toml``: the case's ``beach.toml`` (or ``beach_template.toml``)
  with run parameters applied and ``[output].dir`` set to ``work/latest``
* ``work/latest/``: BEACH output directory (``summary.txt`` marks completion)

Agent skills and references come from ``plugins/beach-context``; wheels ship
them as the data package ``beach.agent_context``.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path
from typing import Any

INPUT = "input/beach.toml"
OUTPUT = "work/latest"
ERROR_WORDS = ("error", "fatal", "killed", "oom")
BATCH_RE = re.compile(r"\bbatch\s+(\d+)\s*/\s*(\d+)\b", re.IGNORECASE)
TAIL_BYTES = 256 * 1024

_HERE = Path(__file__).resolve().parent


def agent_context_path() -> Path | None:
    """Directory with ``skills/`` and ``references/`` for runops to deploy."""
    for candidate in (
        _HERE / "agent_context",
        _HERE.parent / "plugins" / "beach-context",
    ):
        if (candidate / "skills").is_dir():
            return candidate
    return None


class Adapter:
    """BEACH adapter for runops (ported from the runops v1 ``beach`` adapter)."""

    name = "beach"

    # ---- required ----------------------------------------------------
    def render_inputs(self, params: dict[str, Any], run_dir: Path) -> list[str]:
        """Apply dotted-key ``params`` to ``input/beach.toml`` and set the output dir."""
        from runops.adapters import toml_override

        path = run_dir / INPUT
        template = run_dir / "input" / "beach_template.toml"
        if not path.is_file() and template.is_file():
            shutil.move(str(template), str(path))
        if not path.is_file():
            raise ValueError("case must contain beach.toml")
        config = toml_override(path, {**params, "output.dir": OUTPUT})
        written = [INPUT]
        obj = str(config.get("mesh", {}).get("obj_path", ""))
        if obj and Path(obj).is_absolute() and Path(obj).is_file():
            dest = run_dir / "input" / Path(obj).name
            shutil.copy2(obj, dest)
            toml_override(path, {"mesh.obj_path": f"input/{dest.name}"})
            written.append(f"input/{dest.name}")
        (run_dir / OUTPUT).mkdir(parents=True, exist_ok=True)
        return written

    def build_command(self, runtime: Any, run_dir: Path) -> list[str]:
        executable = str(runtime.executable)
        return [executable if executable not in ("", ".") else "beach", INPUT]

    def detect_status(self, run_dir: Path) -> str:
        """``failed`` on error words in the job's stderr or a missing ``summary.txt``."""
        work = run_dir / "work"
        if not work.is_dir():
            return "unknown"
        job_id = _job_id(run_dir)
        logs = [work / f"stderr.{job_id}.log", work / f"{job_id}.err"] if job_id else []
        logs = [p for p in logs if p.is_file()] or sorted(
            [*work.glob("stderr.*.log"), *work.glob("*.err")],
            key=lambda p: p.stat().st_mtime,
        )[-1:]
        if any(any(word in _tail(p).lower() for word in ERROR_WORDS) for p in logs):
            return "failed"
        if _summary_file(run_dir) is not None:
            return "completed"
        return "failed" if any(work.iterdir()) else "unknown"

    # ---- optional ----------------------------------------------------
    def case_template(self) -> dict[str, str]:
        return {
            "beach.toml": (
                "# Base BEACH configuration. [output].dir is set by the adapter.\n"
                "# Override values with dotted keys, e.g. --set sim.batch_count=200\n"
                "# or --set particles.species.0.pos_low=[0.0,0.0,1.0e-3].\n"
                "# See examples/beach.toml for a complete configuration.\n\n"
                "[sim]\n"
                "dt = 1.0e-6\n"
                "batch_count = 100\n"
            )
        }

    def pip_packages(self) -> list[str]:
        return ["beach-bem"]

    def keep_patterns(self) -> list[str]:
        return ["summary.txt", "charges.csv", "stdout.*", "stderr.*", "*.out", "*.err"]

    def summarize(self, run_dir: Path) -> dict[str, Any]:
        summary: dict[str, Any] = {"status": self.detect_status(run_dir)}
        path = _summary_file(run_dir)
        if path is not None:
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                key, sep, value = line.partition("=")
                if sep:
                    summary[key.strip()] = _number(value.strip())
        progress = _progress(run_dir)
        if progress is not None:
            summary["batch"], summary["batch_count"] = progress
        return summary

    def agent_context_dir(self) -> Path | None:
        return agent_context_path()


def _number(value: str) -> Any:
    for cast in (int, float):
        try:
            return cast(value)
        except ValueError:
            pass
    return value


def _job_id(run_dir: Path) -> str:
    if sys.version_info >= (3, 11):
        import tomllib
    else:  # pragma: no cover
        import tomli as tomllib
    try:
        data = tomllib.loads((run_dir / "manifest.toml").read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError):
        return ""
    return str(data.get("job", {}).get("job_id", ""))


def _summary_file(run_dir: Path) -> Path | None:
    work = run_dir / "work"
    for directory in (
        work / "latest",
        work / "outputs" / "latest",
        work / "outputs",
        work,
    ):
        if (directory / "summary.txt").is_file():
            return directory / "summary.txt"
    return None


def _progress(run_dir: Path) -> tuple[int, int] | None:
    work = run_dir / "work"
    logs = sorted(
        [*work.glob("stdout.*.log"), *work.glob("*.out")],
        key=lambda p: p.stat().st_mtime,
    )
    if not logs:
        return None
    matches = BATCH_RE.findall(_tail(logs[-1]))
    return (int(matches[-1][0]), int(matches[-1][1])) if matches else None


def _tail(path: Path) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, 2)
            handle.seek(max(0, handle.tell() - TAIL_BYTES))
            return handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""
