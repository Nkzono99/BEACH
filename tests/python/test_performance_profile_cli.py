from __future__ import annotations

from pathlib import Path

import pytest

from beach.cli.main import main as beachx_main


def _write_profile_fixture(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# BEACH performance profile",
                "# mpi_world_size=4",
                "# omp_max_threads=8",
                "# use rank_max_s of simulation_total for scaling comparisons",
                "region,calls_sum,calls_mean,rank_min_s,rank_mean_s,rank_max_s,imbalance_ratio",
                "program_total,1,0.25,1.2,1.4,1.8,1.286",
                "simulation_total,1,0.25,1.0,1.2,1.6,1.333",
                "field_refresh,4,1.0,0.3,0.4,0.6,1.500",
                "particle_batch,4,1.0,0.4,0.5,0.7,1.400",
                "commit_charge,4,1.0,0.2,0.3,0.5,1.667",
            ]
        ),
        encoding="utf-8",
    )


def test_beachx_profile_saves_png_and_reports_top_regions(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    profile_path = run_dir / "performance_profile.csv"
    save_path = tmp_path / "performance_profile.png"
    _write_profile_fixture(profile_path)

    beachx_main(
        ["profile", str(run_dir), "--save", str(save_path), "--top", "4"]
    )

    streams = capsys.readouterr()
    assert save_path.is_file()
    assert save_path.stat().st_size > 0
    assert f"saved={save_path}" in streams.out
    assert "top_regions_by_rank_max_s=" in streams.out
    assert streams.out.index("program_total:") < streams.out.index("simulation_total:")
    assert streams.out.index("simulation_total:") < streams.out.index("particle_batch:")
    assert streams.out.index("particle_batch:") < streams.out.index("field_refresh:")
    assert "commit_charge:" not in streams.out


def test_beachx_profile_rejects_missing_required_column(tmp_path: Path) -> None:
    profile_path = tmp_path / "performance_profile.csv"
    profile_path.write_text(
        "region,calls_sum,calls_mean,rank_min_s,rank_mean_s,imbalance_ratio\n"
        "simulation_total,1,1,1,1,1\n",
        encoding="utf-8",
    )
    expected = (
        f'Failed to parse "{profile_path}": Missing required column: rank_max_s'
    )

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["profile", str(profile_path)])

    assert excinfo.value.code == expected


def test_beachx_profile_reports_empty_plot_as_cli_error(tmp_path: Path) -> None:
    profile_path = tmp_path / "performance_profile.csv"
    profile_path.write_text(
        "region,calls_sum,calls_mean,rank_min_s,rank_mean_s,rank_max_s,"
        "imbalance_ratio\n"
        "simulation_total,0,0,0,0,0,0\n",
        encoding="utf-8",
    )

    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["profile", str(profile_path)])

    assert excinfo.value.code == (
        f'Failed to plot "{profile_path}": No non-zero rows are available to plot.'
    )
