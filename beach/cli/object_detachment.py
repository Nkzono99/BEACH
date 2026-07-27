"""CLI for frozen-source periodic object detachment diagnostics."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from beach import (
    AdhesionProfile,
    FieldKernelError,
    ObjectInteractionSnapshot,
    load_fortran_result,
)

from ._shared import (
    configure_entry_parser,
    require_finite,
    require_nonnegative_finite,
    require_positive_finite,
)

COMMAND_NAME = "object-detachment"
_INSTANTANEOUS_FIELDS = (
    "component",
    "component_kind",
    "force_x_N",
    "force_y_N",
    "force_z_N",
    "torque_x_Nm",
    "torque_y_Nm",
    "torque_z_Nm",
    "potential_energy_J",
)
_PATH_COMPONENTS = (
    "other_objects_all_images",
    "target_periodic_images",
    "external_uniform",
    "total_external",
)
_PATH_BASE_FIELDS = (
    "displacement_m",
    "force_x_N",
    "force_y_N",
    "force_z_N",
    "torque_x_Nm",
    "torque_y_Nm",
    "torque_z_Nm",
    "electrostatic_work_J",
    "potential_energy_J",
    "potential_difference_work_J",
    "gravity_work_J",
    "adhesion_work_J",
    "available_energy_J",
    "electrostatic_only_speed_m_s",
    "gravity_corrected_speed_m_s",
    "speed_m_s",
)
_PATH_FIELDS = _PATH_BASE_FIELDS + tuple(
    f"{component}_{kind}_{axis}_{unit}"
    for component in _PATH_COMPONENTS
    for kind, unit in (("force", "N"), ("torque", "Nm"))
    for axis in ("x", "y", "z")
)
_ARTIFACT_NAMES = (
    "instantaneous_wrench.csv",
    "path.csv",
    "summary.json",
    "report.md",
)


def _configure_parser(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("output_dir", help="Fortran output directory")
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="beach.toml used by the simulation",
    )
    parser.add_argument("--target-mesh-id", type=int, required=True)
    parser.add_argument(
        "--periodic-model",
        choices=("configured", "infinite-physical"),
        required=True,
        help="use the saved finite-image policy or cached infinite-periodic field",
    )
    parser.add_argument("--z-max-m", type=float, required=True)
    parser.add_argument("--z-points", type=int, required=True)
    parser.add_argument("--mass-kg", type=float, required=True)
    parser.add_argument(
        "--gravity-m-s2",
        type=float,
        default=1.62,
        help="gravity used for release mechanics (default: lunar 1.62 m/s^2)",
    )
    parser.add_argument("--adhesion-force-n", type=float, default=None)
    parser.add_argument("--adhesion-range-m", type=float, default=None)
    parser.add_argument("--eta-translation", type=float, default=1.0)
    parser.add_argument(
        "--torque-origin",
        default="geometric-area-centroid",
        help="geometric-area-centroid, origin, or X,Y,Z in metres",
    )
    parser.add_argument(
        "--step",
        type=_parse_step,
        default=None,
        metavar="final|BATCH",
        help=(
            "charge snapshot: final uses charges.csv (default); an integer selects "
            "that stored history batch, and -1 selects the last stored history row"
        ),
    )
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--library", type=Path, default=None)
    parser.add_argument(
        "--generation-tolerance",
        type=float,
        default=None,
        help="cached periodic-operator generation tolerance",
    )
    parser.add_argument(
        "--fixed-grid",
        action="store_true",
        help="disable adaptive midpoint refinement of the initial z grid",
    )
    parser.add_argument("--relative-tolerance", type=float, default=5.0e-3)
    parser.add_argument(
        "--force-absolute-tolerance-n",
        type=float,
        default=1.0e-12,
    )
    parser.add_argument(
        "--work-absolute-tolerance-j",
        type=float,
        default=1.0e-18,
    )
    parser.add_argument("--max-refinement", type=int, default=8)
    parser.add_argument(
        "--output-dir",
        dest="artifact_dir",
        type=Path,
        default=None,
        help="artifact directory (default: OUTPUT_DIR/object_detachment_mesh_<id>)",
    )


def build_parser(*, prog: str | None = COMMAND_NAME) -> argparse.ArgumentParser:
    """Build the parser for the object-detachment CLI."""

    parser = argparse.ArgumentParser(prog=prog)
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def add_subparser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command under the unified ``beachx`` CLI."""

    parser = subparsers.add_parser(
        COMMAND_NAME,
        help="evaluate periodic object force, work, and release speed",
    )
    _configure_parser(parser)
    return configure_entry_parser(parser, run)


def _parse_step(value: str) -> int | None:
    normalized = value.strip().lower()
    if normalized == "final":
        return None
    try:
        return int(normalized)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "step must be final or an integer batch."
        ) from exc


def run(args: argparse.Namespace) -> None:
    """Evaluate one object and write deterministic review artifacts."""

    parser = args._parser
    require_positive_finite(parser, args.mass_kg, "--mass-kg")
    require_positive_finite(parser, args.z_max_m, "--z-max-m")
    require_nonnegative_finite(parser, args.gravity_m_s2, "--gravity-m-s2")
    require_finite(parser, args.eta_translation, "--eta-translation")
    if not 0.0 <= args.eta_translation <= 1.0:
        parser.error("--eta-translation must be between 0 and 1.")
    if args.z_points < 2:
        parser.error("--z-points must be >= 2.")
    if args.max_refinement < 0:
        parser.error("--max-refinement must be >= 0.")
    require_nonnegative_finite(
        parser,
        args.relative_tolerance,
        "--relative-tolerance",
    )
    require_nonnegative_finite(
        parser,
        args.force_absolute_tolerance_n,
        "--force-absolute-tolerance-n",
    )
    require_nonnegative_finite(
        parser,
        args.work_absolute_tolerance_j,
        "--work-absolute-tolerance-j",
    )
    require_positive_finite(
        parser,
        args.generation_tolerance,
        "--generation-tolerance",
    )
    adhesion = _adhesion_from_args(parser, args)
    torque_origin = _parse_torque_origin(parser, args.torque_origin)

    run_dir = Path(args.output_dir)
    artifact_dir = args.artifact_dir or (
        run_dir / f"object_detachment_mesh_{args.target_mesh_id}"
    )
    requested_model = args.periodic_model
    api_model = requested_model.replace("-", "_")
    displacement = np.linspace(0.0, args.z_max_m, args.z_points)

    try:
        result = load_fortran_result(run_dir)
        charge_source, resolved_charge_batch = _resolve_charge_provenance(
            result,
            args.step,
        )
        with ObjectInteractionSnapshot.from_result(
            result,
            step=args.step,
            config_path=args.config,
            periodic_model=api_model,
            cache_dir=args.cache_dir,
            generation_tolerance=args.generation_tolerance,
            library_path=args.library,
        ) as snapshot:
            probe = snapshot.object_probe(args.target_mesh_id)
            wrench = probe.wrench(torque_origin=torque_origin, components=True)
            path = probe.vertical_path(
                displacement,
                adaptive=not args.fixed_grid,
                relative_tolerance=args.relative_tolerance,
                force_absolute_tolerance_N=args.force_absolute_tolerance_n,
                work_absolute_tolerance_J=args.work_absolute_tolerance_j,
                max_refinement=args.max_refinement,
                torque_origin=torque_origin,
                components=True,
            )
            resolved_snapshot_step = snapshot.step
        release = path.evaluate_release(
            mass_kg=args.mass_kg,
            gravity_m_s2=args.gravity_m_s2,
            adhesion=adhesion,
            eta_translation=args.eta_translation,
        )
    except (FileNotFoundError, ValueError, FieldKernelError) as exc:
        raise SystemExit(str(exc)) from exc

    summary = _build_summary(
        args=args,
        run_dir=run_dir,
        artifact_dir=artifact_dir,
        requested_model=requested_model,
        resolved_snapshot_step=resolved_snapshot_step,
        charge_source=charge_source,
        resolved_charge_batch=resolved_charge_batch,
        adhesion=adhesion,
        wrench=wrench,
        path=path,
        release=release,
    )
    artifacts = {
        "instantaneous_wrench.csv": _instantaneous_csv_text(wrench),
        "path.csv": _path_csv_text(path, release),
        "summary.json": _json_text(summary),
        "report.md": _render_report(summary),
    }
    artifact_dir.mkdir(parents=True, exist_ok=True)
    _replace_artifact_set(artifact_dir, artifacts)
    print(f"saved_dir={artifact_dir}")
    print(f"force_z_N={wrench.force_N[2]:.6e}")
    print(f"barrier_free_from_rest={release.barrier_free_from_rest}")
    release_summary = summary["release"]
    assert isinstance(release_summary, Mapping)
    print(f"numerically_qualified={release_summary['numerically_qualified']}")


def _resolve_charge_provenance(
    result: Any,
    step: int | None,
) -> tuple[str, int]:
    if step is None:
        return "charges.csv", int(result.batches)
    if step != -1:
        return "charge_history.csv", int(step)
    history = result.history
    if history is not None and history.has_data:
        return "charge_history.csv", int(history.batch_indices[-1])
    return "charges.csv", int(result.batches)


def _adhesion_from_args(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> AdhesionProfile:
    force = args.adhesion_force_n
    extent = args.adhesion_range_m
    if (force is None) != (extent is None):
        parser.error(
            "--adhesion-force-n and --adhesion-range-m must be supplied together."
        )
    if force is None:
        return AdhesionProfile.none()
    require_positive_finite(parser, force, "--adhesion-force-n")
    require_positive_finite(parser, extent, "--adhesion-range-m")
    return AdhesionProfile.finite_range_constant(force, extent)


def _parse_torque_origin(
    parser: argparse.ArgumentParser,
    value: str,
) -> str | np.ndarray:
    normalized = value.strip().lower()
    if normalized == "geometric-area-centroid":
        return "geometric_area_centroid"
    if normalized == "origin":
        return "origin"
    try:
        result = np.asarray(
            [float(part.strip()) for part in value.split(",")],
            dtype=np.float64,
        )
    except ValueError:
        result = np.empty(0)
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        parser.error(
            "--torque-origin must be geometric-area-centroid, origin, or X,Y,Z."
        )
    return result


def _instantaneous_csv_text(wrench: Any) -> str:
    rows: list[dict[str, object]] = []
    for name, component in wrench.components.items():
        kind = "physical_total" if name == "total_external" else "physical_additive"
        rows.append(_component_row(name, kind, component))

    metadata = wrench.numerical_metadata
    for name in (
        "periodic_kneq0",
        "physical_k0",
        "primary_free_subtraction",
        "cached_kneq0_trace_correction",
    ):
        component = metadata.get(name)
        if not isinstance(component, Mapping):
            continue
        kind = (
            "numerical_diagnostic_included"
            if name == "cached_kneq0_trace_correction"
            else "numerical_decomposition"
        )
        rows.append(_metadata_component_row(name, kind, component))
    return _csv_text(_INSTANTANEOUS_FIELDS, rows)


def _component_row(name: str, kind: str, component: Any) -> dict[str, object]:
    force = np.asarray(component.force_N)
    torque = np.asarray(component.torque_Nm)
    return {
        "component": name,
        "component_kind": kind,
        "force_x_N": force[0],
        "force_y_N": force[1],
        "force_z_N": force[2],
        "torque_x_Nm": torque[0],
        "torque_y_Nm": torque[1],
        "torque_z_Nm": torque[2],
        "potential_energy_J": component.potential_energy_J,
    }


def _metadata_component_row(
    name: str,
    kind: str,
    component: Mapping[str, object],
) -> dict[str, object]:
    force = np.asarray(component["force_N"])
    torque = np.asarray(component["torque_Nm"])
    return {
        "component": name,
        "component_kind": kind,
        "force_x_N": force[0],
        "force_y_N": force[1],
        "force_z_N": force[2],
        "torque_x_Nm": torque[0],
        "torque_y_Nm": torque[1],
        "torque_z_Nm": torque[2],
        "potential_energy_J": component.get("potential_energy_J"),
    }


def _path_csv_text(path: Any, release: Any) -> str:
    potential = path.potential_energy_J
    potential_work = path.potential_difference_work_J
    rows = []
    for index, displacement in enumerate(path.displacement_m):
        row = {
            "displacement_m": displacement,
            "force_x_N": path.force_N[index, 0],
            "force_y_N": path.force_N[index, 1],
            "force_z_N": path.force_N[index, 2],
            "torque_x_Nm": path.torque_Nm[index, 0],
            "torque_y_Nm": path.torque_Nm[index, 1],
            "torque_z_Nm": path.torque_Nm[index, 2],
            "electrostatic_work_J": path.electrostatic_work_J[index],
            "potential_energy_J": None if potential is None else potential[index],
            "potential_difference_work_J": (
                None if potential_work is None else potential_work[index]
            ),
            "gravity_work_J": release.gravity_work_J[index],
            "adhesion_work_J": release.adhesion_work_J[index],
            "available_energy_J": release.available_energy_J[index],
            "electrostatic_only_speed_m_s": (
                release.electrostatic_only_speed_m_s[index]
            ),
            "gravity_corrected_speed_m_s": (
                release.gravity_corrected_speed_m_s[index]
            ),
            "speed_m_s": release.speed_m_s[index],
        }
        for component in _PATH_COMPONENTS:
            force_values = path.component_force_N.get(component)
            torque_values = path.component_torque_Nm.get(component)
            for axis_index, axis in enumerate(("x", "y", "z")):
                row[f"{component}_force_{axis}_N"] = (
                    None if force_values is None else force_values[index, axis_index]
                )
                row[f"{component}_torque_{axis}_Nm"] = (
                    None if torque_values is None else torque_values[index, axis_index]
                )
        rows.append(row)
    return _csv_text(_PATH_FIELDS, rows)


def _build_summary(
    *,
    args: argparse.Namespace,
    run_dir: Path,
    artifact_dir: Path,
    requested_model: str,
    resolved_snapshot_step: int | None,
    charge_source: str,
    resolved_charge_batch: int,
    adhesion: AdhesionProfile,
    wrench: Any,
    path: Any,
    release: Any,
) -> dict[str, object]:
    potential_available = path.potential_difference_work_J is not None
    unavailable_reason = None
    if not potential_available:
        unavailable_reason = (
            "The selected field evaluator did not provide a potential-energy path."
        )
    numerical = wrench.numerical_metadata
    numerically_qualified = bool(
        path.status == "converged"
        and potential_available
        and release.numerically_qualified
    )
    return {
        "schema_version": 1,
        "inputs": {
            "run_dir": str(run_dir.resolve()),
            "config": str(args.config.resolve()),
            "artifact_dir": str(artifact_dir.resolve()),
            "target_mesh_id": args.target_mesh_id,
            "requested_step": "final" if args.step is None else args.step,
            "resolved_snapshot_step": resolved_snapshot_step,
            "charge_source": charge_source,
            "resolved_charge_batch": resolved_charge_batch,
            "mass_kg": args.mass_kg,
            "gravity_m_s2": args.gravity_m_s2,
            "eta_translation": args.eta_translation,
            "initial_z_points": args.z_points,
            "z_max_m": args.z_max_m,
        },
        "physics_policy": {
            "requested_periodic_model": requested_model,
            "effective_far_correction": numerical.get("effective_far_correction"),
            "self_policy": "exclude_primary_keep_images",
            "surface_trace": "principal_value",
            "source_geometry_policy": "frozen",
            "target_motion": "vertical_translation",
            "source_model": "triangle_p0",
            "target_integration": numerical.get("target_integration"),
            "quadrature_order": numerical.get("quadrature_order"),
            "uniform_potential_has_no_half_factor": True,
            "cached_trace_correction_semantics": (
                "diagnostic_already_included_in_periodic_kneq0_not_additive"
            ),
        },
        "configuration": {
            "saved_config": str(args.config.resolve()),
            "periodic_override": (
                None
                if requested_model == "configured"
                else {"field_periodic_far_correction": "cached_kneq0"}
            ),
            "cache_dir_override": (
                None if args.cache_dir is None else str(args.cache_dir.resolve())
            ),
            "generation_tolerance_override": args.generation_tolerance,
        },
        "adhesion": {
            "kind": adhesion.kind,
            "force_N": (
                adhesion.constant_force_N
                if adhesion.kind == "finite_range_constant"
                else None
            ),
            "range_m": (
                adhesion.range_m
                if adhesion.kind == "finite_range_constant"
                else None
            ),
        },
        "instantaneous_wrench": {
            "mesh_id": wrench.mesh_id,
            "step": wrench.step,
            "total_charge_C": wrench.total_charge_C,
            "force_N": wrench.force_N,
            "torque_Nm": wrench.torque_Nm,
            "torque_origin_m": wrench.torque_origin_m,
            "components": wrench.components,
            "numerical_metadata": numerical,
        },
        "path": {
            "status": path.status,
            "status_reason": path.numerical_metadata.get("status_reason"),
            "point_count": len(path.displacement_m),
            "refinement_count": path.refinement_count,
            "work_relative_mismatch": path.work_relative_mismatch,
            "work_absolute_mismatch_J": path.work_absolute_mismatch_J,
            "potential_difference_available": potential_available,
            "potential_difference_unavailable_reason": unavailable_reason,
            "numerical_metadata": path.numerical_metadata,
        },
        "release": {
            "numerically_qualified": numerically_qualified,
            "mechanics_record_numerically_qualified": (
                release.numerically_qualified
            ),
            "source_path_status": release.source_path_status,
            "instantaneous_force_margin_N": release.instantaneous_force_margin_N,
            "minimum_available_energy_J": release.minimum_available_energy_J,
            "barrier_free_from_rest": release.barrier_free_from_rest,
            "first_inaccessible_displacement_m": (
                release.first_inaccessible_displacement_m
            ),
            "endpoint_available_energy_J": release.endpoint_available_energy_J,
            "endpoint_positive": release.endpoint_positive,
            "endpoint_reachable_from_rest": release.endpoint_reachable_from_rest,
            "endpoint_speed_m_s": release.endpoint_speed_m_s,
            "maximum_reachable_speed_m_s": release.maximum_reachable_speed_m_s,
        },
        "assumptions": [
            "The source geometry and charges are frozen while the target moves.",
            "Only the target primary cell is removed; its periodic images remain.",
            "Adhesion is an explicit user model and defaults to none.",
            "A positive endpoint energy alone does not prove from-rest accessibility.",
            "Execution success alone does not establish physical validity.",
        ],
    }


def _render_report(summary: Mapping[str, object]) -> str:
    wrench = summary["instantaneous_wrench"]
    path = summary["path"]
    release = summary["release"]
    assert isinstance(wrench, Mapping)
    assert isinstance(path, Mapping)
    assert isinstance(release, Mapping)
    force = np.asarray(wrench["force_N"])
    return "\n".join(
        (
            "# BEACH object detachment diagnostic",
            "",
            f"- Target mesh: `{wrench['mesh_id']}`",
            f"- Instantaneous Fz: `{force[2]:.9e} N`",
            f"- Path status: `{path['status']}`",
            f"- Numerically qualified: `{release['numerically_qualified']}`",
            f"- From-rest barrier free: `{release['barrier_free_from_rest']}`",
            f"- Endpoint speed: `{release['endpoint_speed_m_s']:.9e} m/s`",
            "",
            "The result is conditional on the frozen-source, gravity, adhesion, and",
            "periodic-boundary assumptions recorded in `summary.json`. A successful",
            "run is not by itself evidence of numerical convergence or physical validity.",
            "",
        )
    )


def _csv_text(
    fieldnames: Sequence[str],
    rows: Sequence[Mapping[str, object]],
) -> str:
    validated = [
        {name: _strict_csv_value(row.get(name)) for name in fieldnames}
        for row in rows
    ]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(validated)
    return stream.getvalue()


def _strict_csv_value(value: object) -> object:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("CSV output contains a non-finite number.")
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"unsupported CSV value: {type(value).__name__}")


def _json_text(value: Mapping[str, object]) -> str:
    serializable = _json_value(value)
    return (
        json.dumps(
            serializable,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    )


def _replace_artifact_set(
    directory: Path,
    artifacts: Mapping[str, str],
) -> None:
    if tuple(artifacts) != _ARTIFACT_NAMES:
        raise ValueError("object-detachment artifact set is incomplete or unordered.")
    staged: dict[str, Path] = {}
    backups: dict[str, Path] = {}
    installed: list[str] = []
    try:
        for name in _ARTIFACT_NAMES:
            staged[name] = _stage_text(directory, name, artifacts[name])
        for name in _ARTIFACT_NAMES:
            target = directory / name
            if target.exists():
                backups[name] = _backup_file(directory, name, target)
        for name in _ARTIFACT_NAMES:
            os.replace(staged[name], directory / name)
            installed.append(name)
    except BaseException as exc:
        rollback_errors: list[OSError] = []
        for name in reversed(installed):
            target = directory / name
            try:
                backup = backups.pop(name, None)
                if backup is None:
                    target.unlink(missing_ok=True)
                else:
                    os.replace(backup, target)
            except OSError as rollback_error:
                rollback_errors.append(rollback_error)
        if rollback_errors:
            raise RuntimeError(
                "artifact replacement failed and rollback was incomplete: "
                + "; ".join(str(error) for error in rollback_errors)
            ) from exc
        raise
    finally:
        for path in (*staged.values(), *backups.values()):
            path.unlink(missing_ok=True)


def _stage_text(directory: Path, name: str, content: str) -> Path:
    descriptor, raw_path = tempfile.mkstemp(
        prefix=f".{name}.",
        suffix=".tmp",
        dir=directory,
        text=True,
    )
    path = Path(raw_path)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        path.unlink(missing_ok=True)
        raise
    return path


def _backup_file(directory: Path, name: str, source: Path) -> Path:
    descriptor, raw_path = tempfile.mkstemp(
        prefix=f".{name}.",
        suffix=".backup",
        dir=directory,
    )
    os.close(descriptor)
    path = Path(raw_path)
    try:
        shutil.copy2(source, path)
    except BaseException:
        path.unlink(missing_ok=True)
        raise
    return path


def _json_value(value: object) -> object:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("JSON output contains a non-finite number.")
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if hasattr(value, "force_N") and hasattr(value, "torque_Nm"):
        return {
            "force_N": _json_value(value.force_N),
            "torque_Nm": _json_value(value.torque_Nm),
            "potential_energy_J": _json_value(value.potential_energy_J),
        }
    raise TypeError(f"unsupported JSON value: {type(value).__name__}")


def main(argv: Sequence[str] | None = None) -> None:
    """Run the object-detachment CLI entry point."""

    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
