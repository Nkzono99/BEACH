"""runops (>= 2.0) adapter for BEACH, published as the ``beach`` entry point.

runops discovers :class:`Adapter` through the ``runops.adapters`` entry point
group declared in ``pyproject.toml``. The adapter is only called from inside
runops; runops helpers are imported lazily.

Run directory conventions (the job's working directory is the run directory):

* ``input/beach.toml``: the case's ``beach.toml`` (or ``beach_template.toml``)
  in the grouped layout (``run``, ``domain``, ``mesh``, ``particles``, ``fields``,
  ``sheath``, ``output``), with run parameters applied as dotted keys (e.g.
  ``run.batch_count``) and ``[output].dir`` set to ``work/latest``
* OBJ meshes in ``[[mesh.parts]]`` are referenced from the run directory: a
  relative path to a file shipped with the case becomes ``input/<path>``, and
  an absolute path is copied into ``input/``
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
        from beach.config._layout import is_grouped

        if not is_grouped(config):
            raise ValueError(
                "input/beach.toml uses the pre-3.0 flat layout ([sim] ...); BEACH 3 accepts "
                "only the grouped tables run, domain, mesh, particles, fields, sheath and "
                "output. Start from `beachx config init` or the files in examples/."
            )
        written = [INPUT]
        moves = {}
        parts = config.get("mesh", {}).get("parts", [])
        for index, part in enumerate(parts if isinstance(parts, list) else []):
            if (
                not isinstance(part, dict)
                or part.get("kind") != "obj"
                or not part.get("path")
            ):
                continue
            new_path = _localize_obj(Path(str(part["path"])), run_dir, index)
            if new_path is not None:
                moves[f"mesh.parts.{index}.path"] = new_path
                if new_path not in written:
                    written.append(new_path)
        if moves:
            toml_override(path, moves)
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
        """The same starting config as ``beachx config init`` (grouped layout)."""
        from beach.config import default_config, dump_beach_toml

        lines = dump_beach_toml(default_config()).splitlines(keepends=True)
        note = [
            "# runops sets [output].dir for each run. Vary values per run with dotted keys,\n",
            "# e.g. `runo new CASE --set run.batch_count=200`. More starting points: examples/.\n",
        ]
        return {"beach.toml": "".join(lines[:2] + note + lines[2:])}

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


def _localize_obj(path: Path, run_dir: Path, index: int) -> str | None:
    """Return a run-directory-relative OBJ path to write back, or None to keep it."""
    if path.is_absolute():
        if not path.is_file():
            return None
        dest = run_dir / "input" / path.name
        if dest.exists() and dest.read_bytes() != path.read_bytes():
            dest = run_dir / "input" / f"part{index}-{path.name}"
        shutil.copy2(path, dest)
        return dest.relative_to(run_dir).as_posix()
    if (run_dir / path).is_file():
        return None  # already valid from the run directory
    if (run_dir / "input" / path).is_file():
        return (Path("input") / path).as_posix()
    return None


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
