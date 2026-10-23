"""CLI adapter for BEACH particle-workload estimation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from beach.tooling.workload import (
    completed_batches_from_resume_config,
    estimate_workload,
    load_toml,
    read_macro_residuals,
)

from . import _options
from ._output import emit


def _default_threads() -> int:
    value = os.environ.get("OMP_NUM_THREADS", "").strip()
    if value:
        try:
            parsed = int(value)
            if parsed > 0:
                return parsed
        except ValueError:
            pass
    return 1


def _default_mpi_ranks() -> int:
    for key in ("PMI_SIZE", "OMPI_COMM_WORLD_SIZE", "SLURM_NTASKS"):
        value = os.environ.get(key, "").strip()
        if not value:
            continue
        try:
            parsed = int(value)
            if parsed > 0:
                return parsed
        except ValueError:
            pass
    return 1


def _default_mpi_rank() -> int:
    for key in ("PMI_RANK", "OMPI_COMM_WORLD_RANK", "SLURM_PROCID"):
        value = os.environ.get(key, "").strip()
        if not value:
            continue
        try:
            parsed = int(value)
            if parsed >= 0:
                return parsed
        except ValueError:
            pass
    return 0


def configure(parser: argparse.ArgumentParser) -> None:
    """Add the ``workload`` options."""

    parser.add_argument("config", type=Path, help="path to beach.toml")
    parser.add_argument(
        "--threads",
        type=int,
        default=_default_threads(),
        help="OpenMP thread count for per-thread estimate (default: OMP_NUM_THREADS or 1)",
    )
    parser.add_argument(
        "--macro-residuals",
        type=Path,
        default=None,
        help="optional macro_residuals.csv to start from resume state",
    )
    parser.add_argument(
        "--mpi-ranks",
        type=int,
        default=_default_mpi_ranks(),
        help="MPI world size used for local(rank) workload estimate (default: MPI env or 1)",
    )
    parser.add_argument(
        "--mpi-rank",
        type=int,
        default=_default_mpi_rank(),
        help="rank index used for local(rank) workload estimate (default: MPI env or 0)",
    )
    parser.add_argument(
        "--show-batches",
        type=int,
        default=10,
        help="number of head batches to print in detail (default: 10)",
    )
    _options.add_json(parser)


def run(args: argparse.Namespace) -> None:
    """Execute the workload-estimation command."""

    if not args.config.exists():
        raise SystemExit(f"config file not found: {args.config}")

    config = load_toml(args.config)
    species_raw = config.get("particles", {}).get("species", [])
    initial_residuals = read_macro_residuals(args.macro_residuals, len(species_raw))
    completed_batches = completed_batches_from_resume_config(config)

    result = estimate_workload(
        config=config,
        threads=args.threads,
        initial_residuals=initial_residuals,
        mpi_ranks=args.mpi_ranks,
        mpi_rank=args.mpi_rank,
        completed_batches=completed_batches,
    )

    batch_totals = result["batch_totals"]
    batch_thread_min = result["batch_thread_min"]
    batch_thread_max = result["batch_thread_max"]
    total_particles = result["total_particles"]
    batch_count = result["batch_count"]
    threads = result["threads"]

    report: dict[str, object] = {
        "config": args.config,
        "threads": threads,
        "mpi_ranks": result["mpi_ranks"],
        "mpi_rank": result["mpi_rank"],
        "estimate_scope": "local_rank",
    }
    if result["secondary_emission_excluded"]:
        report["secondary_emission_scope"] = "excluded_collision_generated_particles"
    report.update(
        {
            "target_batch_count": result["target_batch_count"],
            "completed_batches": result["completed_batches"],
            "remaining_batch_count": batch_count,
            "resolved_batch_duration": result["resolved_batch_duration"],
            "total_particles": total_particles,
            "global_total_particles": result["global_total_particles"],
            "local_reservoir_particles": result["local_reservoir_particles"],
            "global_reservoir_particles": result["global_reservoir_particles"],
            "particles_per_batch_min": min(batch_totals, default=0),
            "particles_per_batch_max": max(batch_totals, default=0),
            "particles_per_batch_avg": (
                total_particles / batch_count if batch_count > 0 else 0.0
            ),
            "per_thread_particles_min": min(batch_thread_min, default=0),
            "per_thread_particles_max": max(batch_thread_max, default=0),
            "per_thread_particles_avg": (
                total_particles / (batch_count * threads) if batch_count > 0 else 0.0
            ),
            "final_macro_residuals": result["final_residuals"],
            "batch_details": [
                {
                    "batch": batch_idx + 1,
                    "total": batch_totals[batch_idx],
                    "global_total": result["global_batch_totals"][batch_idx],
                    "per_thread_min": batch_thread_min[batch_idx],
                    "per_thread_max": batch_thread_max[batch_idx],
                    "species": result["species_per_batch"][batch_idx],
                    "global_species": result["global_species_per_batch"][batch_idx],
                }
                for batch_idx in range(max(0, min(args.show_batches, batch_count)))
            ],
        }
    )
    emit(args, report)
