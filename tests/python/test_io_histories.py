"""Public contracts for the shared batch history and auxiliary CSV readers."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from beach import BatchHistory, load_run
from beach.io import (
    load_fixed_current_history,
    load_matching_plane_spectrum_history,
    load_performance_profile,
    load_top_reference_history,
    summary_fields_markdown,
)
from result_fixtures import _write_minimal_result_fixture


@pytest.mark.parametrize("quantity", ["charge", "potential"])
def test_batch_history_reads_strided_batches_by_element_id(
    tmp_path: Path, quantity: str
) -> None:
    path = tmp_path / f"{quantity}_history.csv"
    column = "charge_C" if quantity == "charge" else "potential_V"
    metadata = ",processed_particles,rel_change" if quantity == "charge" else ""
    prefix = ",10,0.5" if quantity == "charge" else ""
    path.write_text(
        f"batch{metadata},elem_idx,{column}\n"
        f"2{prefix},2,-2.0\n2{prefix},1,1.0\n"
        f"7{prefix},1,3.0\n7{prefix},2,-4.0\n",
        encoding="utf-8",
    )
    history = BatchHistory(path, mesh_nelem=2, value_column=column)

    np.testing.assert_array_equal(history.batch_indices, [2, 7])
    np.testing.assert_array_equal(history.get_step(-1), [3.0, -4.0])
    np.testing.assert_array_equal(history.as_array(), [[1.0, 3.0], [-2.0, -4.0]])
    if quantity == "potential":
        assert history.processed_particles_by_batch is None
        assert history.rel_change_by_batch is None
    with pytest.raises(ValueError, match="step=3 is not found"):
        history.get_step(3)


def test_load_run_exposes_saved_potential_and_auxiliary_histories(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run"
    out.mkdir()
    _write_minimal_result_fixture(out)
    (out / "potential_history.csv").write_text(
        "batch,elem_idx,potential_V\n5,1,7.5\n5,2,-2.5\n", encoding="utf-8"
    )
    (out / "top_reference_history.csv").write_text(
        "batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,potential_min_V,potential_max_V\n"
        "5,0.1,2.0,16,3.0,0.2,2.0,4.0\n",
        encoding="utf-8",
    )
    (out / "matching_plane_spectrum_history.csv").write_text(
        "batch,simulated_time_s,bin,energy_low_eV,energy_high_eV,photoelectron_outward_flux_m2_s\n"
        "5,0.1,2,0.001,0.002,4e6\n",
        encoding="utf-8",
    )
    (out / "fixed_current_history.csv").write_text(
        "batch,simulated_time_s,window_batches,species_idx,channel,tracked_charge_C,target_charge_C,target_over_tracked\n"
        "5,0.1,5,1,absorbed,0,-1e-8,nan\n",
        encoding="utf-8",
    )
    result = load_run(out)

    assert isinstance(result.history.potential, BatchHistory)
    np.testing.assert_array_equal(result.history.potential[-1], [7.5, -2.5])
    assert result.top_reference_history[0].sample_n == 16
    assert result.top_reference_history[0].switch_phi_rms_mismatch_v == 0.0
    assert result.matching_plane_spectrum_history[0].energy_high_ev == 0.002
    fixed = result.fixed_current_history[0]
    assert fixed.target_charge_c == -1e-8
    assert math.isnan(fixed.target_over_tracked)


@pytest.mark.parametrize(
    "reader,header",
    [
        (
            load_top_reference_history,
            "batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,potential_min_V,potential_max_V",
        ),
        (
            load_matching_plane_spectrum_history,
            "batch,simulated_time_s,bin,energy_low_eV,energy_high_eV,photoelectron_outward_flux_m2_s",
        ),
        (
            load_fixed_current_history,
            "batch,simulated_time_s,window_batches,species_idx,channel,tracked_charge_C,target_charge_C,target_over_tracked",
        ),
    ],
)
def test_auxiliary_history_distinguishes_absent_from_header_only(
    tmp_path: Path, reader, header: str
) -> None:
    path = tmp_path / "history.csv"
    assert reader(path) is None
    path.write_text(header + "\n", encoding="utf-8")
    assert reader(path) == ()
    path.write_text("unknown\n0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing required columns"):
        reader(path)


def test_profile_reader_is_available_without_the_cli(tmp_path: Path) -> None:
    path = tmp_path / "performance_profile.csv"
    path.write_text(
        "# mpi_world_size=4\n"
        "region,calls_sum,calls_mean,rank_min_s,rank_mean_s,rank_max_s,imbalance_ratio\n"
        "simulation_total,4,1,2,3,4,1.333\n",
        encoding="utf-8",
    )
    metadata, rows = load_performance_profile(path)
    assert metadata == {"mpi_world_size": "4"}
    assert rows[0]["region"] == "simulation_total"
    assert rows[0]["rank_max_s"] == 4.0


def test_generated_summary_reference_matches_the_reader_contract() -> None:
    path = Path(__file__).resolve().parents[2] / "docs/reference/OutputReference.md"
    table = path.read_text(encoding="utf-8").split(
        "<!-- BEGIN GENERATED SUMMARY FIELDS -->", 1
    )[1]
    table = table.split("<!-- END GENERATED SUMMARY FIELDS -->", 1)[0]
    assert table.strip() == summary_fields_markdown()
