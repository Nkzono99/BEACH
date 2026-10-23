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
    global_total_particles = result["global_total_particles"]
    batch_count = result["batch_count"]
    target_batch_count = result["target_batch_count"]
    completed_batches = result["completed_batches"]
    threads = result["threads"]
    mpi_ranks = result["mpi_ranks"]
    mpi_rank = result["mpi_rank"]

    print(f"config={args.config}")
    print(f"threads={threads}")
    print(f"mpi_ranks={mpi_ranks}")
    print(f"mpi_rank={mpi_rank}")
    print("estimate_scope=local_rank")
    if result["secondary_emission_excluded"]:
        print("secondary_emission_scope=excluded_collision_generated_particles")
    print(f"target_batch_count={target_batch_count}")
    print(f"completed_batches={completed_batches}")
    print(f"remaining_batch_count={batch_count}")
    print(f"resolved_batch_duration={result['resolved_batch_duration']}")
    print(f"total_particles={total_particles}")
    print(f"global_total_particles={global_total_particles}")
    print(f"local_reservoir_particles={result['local_reservoir_particles']}")
    print(f"global_reservoir_particles={result['global_reservoir_particles']}")
    if batch_count > 0:
        print(f"particles_per_batch_min={min(batch_totals)}")
        print(f"particles_per_batch_max={max(batch_totals)}")
        print(f"particles_per_batch_avg={total_particles / batch_count:.3f}")
        print(f"per_thread_particles_min={min(batch_thread_min)}")
        print(f"per_thread_particles_max={max(batch_thread_max)}")
        print(
            "per_thread_particles_avg="
            f"{total_particles / (batch_count * threads):.3f}"
        )
    else:
        print("particles_per_batch_min=0")
        print("particles_per_batch_max=0")
        print("particles_per_batch_avg=0.000")
        print("per_thread_particles_min=0")
        print("per_thread_particles_max=0")
        print("per_thread_particles_avg=0.000")
    print(f"final_macro_residuals={result['final_residuals']}")

    n_detail = max(0, min(args.show_batches, batch_count))
    if n_detail > 0:
        print("batch_details=")
        for batch_idx in range(n_detail):
            species_counts = result["species_per_batch"][batch_idx]
            global_species_counts = result["global_species_per_batch"][batch_idx]
            print(
                "  "
                f"batch={batch_idx + 1} "
                f"total={batch_totals[batch_idx]} "
                f"global_total={result['global_batch_totals'][batch_idx]} "
                f"per_thread=[{batch_thread_min[batch_idx]},{batch_thread_max[batch_idx]}] "
                f"species={species_counts} "
                f"global_species={global_species_counts}"
            )

