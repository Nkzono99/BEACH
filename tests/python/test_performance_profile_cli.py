from __future__ import annotations

from pathlib import Path

import pytest

from beach.cli_plot_performance_profile import main as plot_performance_main


def _write_profile_fixture(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# BEACH performance profile",
                "# mpi_world_size=4",
                "# omp_max_threads=8",
                "# detail_enabled=T",
                "# use rank_max_s of simulation_total for scaling comparisons",
                "region,calls_sum,calls_mean,rank_min_s,rank_mean_s,rank_max_s,imbalance_ratio",
                "program_total,1,0.25,1.2,1.4,1.8,1.286",
                "simulation_total,1,0.25,1.0,1.2,1.6,1.333",
                "field_refresh,4,1.0,0.3,0.4,0.6,1.500",
                "particle_batch,4,1.0,0.4,0.5,0.7,1.400",
                "particle_field_eval,1000,250.0,0.9,1.1,1.5,1.364",
            ]
        ),
        encoding="utf-8",
    )


def test_plot_performance_profile_saves_png(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    profile_path = run_dir / "performance_profile.csv"
    save_path = tmp_path / "performance_profile.png"
    _write_profile_fixture(profile_path)

    plot_performance_main([str(run_dir), "--save", str(save_path), "--top", "4"])

    assert save_path.exists()
