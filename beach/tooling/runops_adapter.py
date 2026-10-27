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
* ``work/latest/``: BEACH output directory. A run is complete when its
  ``summary.txt`` reports ``batches`` equal to ``run.batch_count``
* ``runo new --from RUN --continue``: ``work/restart`` links to the source run's
  output directory and ``[run.restart].from`` points at it; ``input/restart_source.txt``
  records where the checkpoint came from

Agent skills and references come from ``plugins/beach-context``; wheels ship
them as the data package ``beach.agent_context``.
"""

from __future__ import annotations

import filecmp
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any

from beach.io.summary import parse_summary_text

INPUT = "input/beach.toml"
OUTPUT = "work/latest"
RESTART_LINK = "work/restart"
RESTART_NOTE = "input/restart_source.txt"
# Files of the final checkpoint in work/latest/, kept by `runo tidy --purge-work`
# so a purged run can still be continued.
CHECKPOINT_FILES = (
    "rng_state*.txt",
    "macro_residuals.csv",
    "charge_ledger.csv",
    "checkpoint_complete.txt",
)
ERROR_WORDS = ("error", "fatal", "killed", "oom")
BATCH_RE = re.compile(r"\bbatch\s+(\d+)\s*/\s*(\d+)\b", re.IGNORECASE)
RESUME_RE = re.compile(r"^resuming_from_batches\s*=\s*(\d+)", re.MULTILINE)
TAIL_BYTES = 256 * 1024

_PACKAGE_DIR = Path(__file__).resolve().parents[1]


def agent_context_path() -> Path | None:
    """Directory with ``skills/`` and ``references/`` for runops to deploy."""
    for candidate in (
        _PACKAGE_DIR / "agent_context",
        _PACKAGE_DIR.parent / "plugins" / "beach-context",
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
        """``failed`` on error words in the job's stderr, or when ``summary.txt`` is
        missing or reports fewer ``batches`` than ``run.batch_count``."""
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
        path = _summary_file(run_dir)
        if path is not None:
            batches = _summary_batches(path)
            target = _batch_count(run_dir)
            if batches is not None and target is not None and batches < target:
                return "failed"
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
        return [
            "summary.txt",
            "charges.csv",
            *(f"latest/{name}" for name in CHECKPOINT_FILES),
            "stdout.*",
            "stderr.*",
            "*.out",
            "*.err",
        ]

    def validate_params(self, params: dict[str, Any]) -> list[str]:
        """Warn about dotted keys that the BEACH schema does not define."""
        try:
            from beach.config.schema import load_schema

            schema, _ = load_schema()
        except (ImportError, OSError, ValueError):
            return []
        warnings = []
        for key in params:
            if not _schema_has_key(schema, key):
                warnings.append(
                    f"{key!r} is not a BEACH 3 input key; check the spelling and the table "
                    "(see beach-references Parameters / beach.schema.json)"
                )
        return warnings

    def extend_inputs(
        self, src_run: Path, run_dir: Path, params: dict[str, Any]
    ) -> list[str]:
        """Resume from the source run's checkpoint (``runo new --from RUN --continue``).

        ``run.batch_count`` is cumulative: it must exceed the source's ``batches``.
        BEACH itself checks the mesh fingerprint and the MPI world size at start.
        """
        from runops.adapters import toml_override

        source = _checkpoint_dir(src_run)
        if source is None:
            raise ValueError(
                f"{src_run} has no BEACH checkpoint (summary.txt) under work/ to continue from"
            )
        summaries = [
            source / "summary.txt",
            *source.glob("checkpoints/slot*/summary.txt"),
        ]
        counts = [
            b
            for p in summaries
            if p.is_file() and (b := _summary_batches(p)) is not None
        ]
        saved = max(counts, default=None)
        target = _batch_count(run_dir)
        if saved is not None and (target is None or target <= saved):
            raise ValueError(
                f"run.batch_count is cumulative; the checkpoint already has {saved} batches, "
                f"so continue with --set run.batch_count=<more than {saved}>"
            )
        link = run_dir / RESTART_LINK
        link.parent.mkdir(parents=True, exist_ok=True)
        link.symlink_to(source.resolve(), target_is_directory=True)
        toml_override(
            run_dir / INPUT,
            {"run.restart.from": RESTART_LINK, "output.enabled": True},
        )
        (run_dir / RESTART_NOTE).write_text(
            f"source = {source.resolve()}\nbatches = {saved if saved is not None else ''}\n",
            encoding="utf-8",
        )
        return [INPUT, RESTART_NOTE, RESTART_LINK]

    def summarize(self, run_dir: Path) -> dict[str, Any]:
        summary: dict[str, Any] = {"status": self.detect_status(run_dir)}
        path = _summary_file(run_dir)
        if path is not None:
            values = parse_summary_text(
                path.read_text(encoding="utf-8", errors="replace"), source=str(path)
            )
            summary.update({key: _number(value) for key, value in values.items()})
        progress = _progress(run_dir)
        if progress is not None:
            summary["batch"], summary["batch_count"] = progress
        return summary

    def agent_context_dir(self) -> Path | None:
        return agent_context_path()

    def progress(self, run_dir: Path) -> dict[str, Any] | None:
        """Batches of the latest job for ``runo sync`` (runops >= 2.5).

        ``start`` is the batch a resumed run continued from (``resuming_from_batches=``
        at the top of stdout), so runops measures speed from the resume point.
        """
        log = _latest_stdout(run_dir)
        if log is None:
            return None
        matches = BATCH_RE.findall(_tail(log))
        if not matches:
            return None
        done, total = int(matches[-1][0]), int(matches[-1][1])
        resumed = RESUME_RE.search(_head(log))
        return {
            "done": done,
            "total": total,
            "unit": "batch",
            "start": int(resumed.group(1)) if resumed else 0,
        }


def _localize_obj(path: Path, run_dir: Path, index: int) -> str | None:
    """Return a run-directory-relative OBJ path to write back, or None to keep it."""
    if path.is_absolute():
        if not path.is_file():
            return None
        dest = run_dir / "input" / path.name
        if dest.exists() and not (
            os.path.samefile(path, dest) or filecmp.cmp(path, dest, shallow=False)
        ):
            dest = run_dir / "input" / f"part{index}-{path.name}"
        # The path may already name the destination itself, directly or by a symlink.
        if not (dest.exists() and os.path.samefile(path, dest)):
            shutil.copy2(path, dest)
        return dest.relative_to(run_dir).as_posix()
    if (run_dir / path).is_file():
        return None  # already valid from the run directory
    if (run_dir / "input" / path).is_file():
        return (Path("input") / path).as_posix()
    return None


def _checkpoint_dir(run_dir: Path) -> Path | None:
    """The source run's output directory that holds a final or periodic checkpoint."""
    work = run_dir / "work"
    for directory in (
        work / "latest",
        work / "outputs" / "latest",
        work / "outputs",
        work,
    ):
        if (directory / "summary.txt").is_file() or any(
            directory.glob("checkpoints/slot*/summary.txt")
        ):
            return directory
    return None


def _summary_batches(path: Path) -> int | None:
    for line in _tail(path).splitlines():
        key, sep, value = line.partition("=")
        if sep and key.strip() == "batches":
            try:
                return int(value.strip())
            except ValueError:
                return None
    return None


def _batch_count(run_dir: Path) -> int | None:
    if sys.version_info >= (3, 11):
        import tomllib
    else:  # pragma: no cover
        import tomli as tomllib
    try:
        config = tomllib.loads((run_dir / INPUT).read_text(encoding="utf-8"))
        return int(config["run"]["batch_count"])
    except (OSError, tomllib.TOMLDecodeError, KeyError, TypeError, ValueError):
        return None


_INDEX_RE = re.compile(r"\[(\d+)\]")


def _schema_has_key(schema: dict[str, Any], key: str) -> bool:
    """Whether the dotted ``key`` (``a.b``, ``a.0.b``, ``a[0].b``) names a schema property."""
    steps: list[str | int] = []
    for part in key.split("."):
        name = _INDEX_RE.sub("", part)
        if name:
            steps.append(int(name) if name.isdigit() else name)
        steps += [int(i) for i in _INDEX_RE.findall(part)]
    rules = [schema]
    for step in steps:
        rules = _expand(schema, rules)
        if isinstance(step, int):
            rules = [r["items"] for r in rules if isinstance(r.get("items"), dict)]
        else:
            found = [
                r["properties"][name]
                for r in rules
                for name in (step, step.lower())
                if name in r.get("properties", {})
            ]
            if not found:
                extra = [
                    r["additionalProperties"]
                    for r in rules
                    if "additionalProperties" in r
                ]
                found = [
                    x for x in extra if isinstance(x, dict)
                ]  # a map of named entries
                if not found:
                    return not any(x is False for x in extra)
            rules = found
        if not rules:
            return False
    return True


def _expand(schema: dict[str, Any], rules: list[Any]) -> list[dict[str, Any]]:
    """Follow ``$ref`` and flatten ``allOf``/``anyOf``/``oneOf``/``then``/``else`` branches."""
    out: list[dict[str, Any]] = []
    stack = list(rules)
    seen = 0
    while stack and seen < 10000:
        seen += 1
        rule = stack.pop()
        if not isinstance(rule, dict):
            continue
        ref = rule.get("$ref")
        if isinstance(ref, str) and ref.startswith("#/"):
            node: Any = schema
            for part in ref[2:].split("/"):
                node = node.get(part.replace("~1", "/").replace("~0", "~"), {})
            stack.append(node)
        for name in ("allOf", "anyOf", "oneOf"):
            stack.extend(rule.get(name, []))
        for name in ("then", "else"):
            if isinstance(rule.get(name), dict):
                stack.append(rule[name])
        out.append(rule)
    return out


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


def _latest_stdout(run_dir: Path) -> Path | None:
    """stdout of the manifest's job if present, else the newest stdout in ``work/``."""
    work = run_dir / "work"
    if not work.is_dir():
        return None
    job_id = _job_id(run_dir)
    for name in (f"stdout.{job_id}.log", f"{job_id}.out") if job_id else ():
        if (work / name).is_file():
            return work / name
    logs = sorted(
        [*work.glob("stdout.*.log"), *work.glob("*.out")],
        key=lambda p: p.stat().st_mtime,
    )
    return logs[-1] if logs else None


def _progress(run_dir: Path) -> tuple[int, int] | None:
    log = _latest_stdout(run_dir)
    if log is None:
        return None
    matches = BATCH_RE.findall(_tail(log))
    return (int(matches[-1][0]), int(matches[-1][1])) if matches else None


def _head(path: Path, size: int = 16 * 1024) -> str:
    try:
        with path.open("rb") as handle:
            return handle.read(size).decode("utf-8", errors="replace")
    except OSError:
        return ""


def _tail(path: Path) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, 2)
            handle.seek(max(0, handle.tell() - TAIL_BYTES))
            return handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""
