from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from beach.outer_kinetic import (
    ELEMENTARY_CHARGE,
    RAW_HEADER,
    CertificationConfig,
    KineticQuery,
    OuterKineticConfig,
    ReservoirConfig,
    VelocityGridConfig,
    _certify_time_history,
    _finalize_classification,
    _make_photoelectron_state,
    _make_reservoir_state,
    advect_v_finite_volume,
    advect_z_finite_volume,
    compute_field_profile,
    convert_kinetic_table,
    run_kinetic_query,
    write_kinetic_atlas,
)


def _balanced_config() -> OuterKineticConfig:
    grid = VelocityGridConfig(nv=16, vmin_mps=-4.0, vmax_mps=4.0)
    common = {
        "number_density_m3": 1.0e6,
        "temperature_ev": 1.0,
        "drift_velocity_mps": 0.0,
        "mass_kg": ELEMENTARY_CHARGE,
        "grid": grid,
    }
    return OuterKineticConfig(
        matching_plane_z_m=2.0,
        z_length_m=1.0,
        nz=16,
        cfl=0.4,
        max_time_s=1.0,
        initial_condition="reservoir",
        electron=ReservoirConfig(charge_c=-ELEMENTARY_CHARGE, **common),
        ion=ReservoirConfig(charge_c=ELEMENTARY_CHARGE, **common),
        photoelectron_grid=grid,
        certification=CertificationConfig(
            warmup_time_s=0.2,
            averaging_window_s=0.4,
            sample_interval_s=0.05,
            stationarity_rtol=1.0e-10,
            drift_rtol=1.0e-10,
            sem_rtol=1.0e-10,
            steady_fluctuation_rtol=1.0e-10,
            autocorrelation_windows=2.0,
            far_field_abs_v_m=1.0e-10,
            far_charge_rel=1.0e-10,
            gauss_abs_c_m2=1.0e-20,
            charge_budget_rtol=1.0e-10,
            velocity_loss_rtol=1.0e-10,
        ),
    )


def _zero_query() -> KineticQuery:
    return KineticQuery(
        displacement_c_m2=0.0,
        photoelectron_outward_number_flux_m2_s=0.0,
        photoelectron_outward_mean_normal_energy_ev=0.0,
    )


def test_field_integration_matches_zero_and_uniform_charge() -> None:
    zero = compute_field_profile(
        np.zeros(4), displacement_c_m2=3.0 * 8.8541878128e-12, dz_m=0.5
    )
    np.testing.assert_allclose(zero.electric_faces_v_m, 3.0)
    np.testing.assert_allclose(zero.potential_faces_v, [6.0, 4.5, 3.0, 1.5, 0.0])
    assert zero.gauss_residual_c_m2 == pytest.approx(0.0, abs=1.0e-26)

    rho = np.full(4, 2.0e-12)
    charged = compute_field_profile(rho, displacement_c_m2=0.0, dz_m=0.25)
    expected = np.arange(5) * 2.0e-12 * 0.25 / 8.8541878128e-12
    np.testing.assert_allclose(charged.electric_faces_v_m, expected)
    assert charged.gauss_residual_c_m2 == pytest.approx(0.0, abs=1.0e-26)


def test_finite_volume_transport_shifts_one_cell_and_preserves_positivity() -> None:
    distribution = np.array(
        [[1.0, 4.0], [2.0, 3.0], [3.0, 2.0], [4.0, 1.0]]
    )
    shifted_z, z_inventory_change = advect_z_finite_volume(
        distribution,
        np.array([1.0, -1.0]),
        dt_s=1.0,
        dz_m=1.0,
        dv_mps=1.0,
        left_inflow=np.zeros(2),
        right_inflow=np.zeros(2),
    )
    np.testing.assert_allclose(
        shifted_z,
        [[0.0, 3.0], [1.0, 2.0], [2.0, 1.0], [3.0, 0.0]],
    )
    assert np.sum(shifted_z) - np.sum(distribution) == pytest.approx(
        z_inventory_change
    )

    velocity_distribution = np.array([[1.0, 2.0, 3.0, 4.0]])
    shifted_v, v_inventory_change = advect_v_finite_volume(
        velocity_distribution,
        np.array([1.0]),
        dt_s=1.0,
        dv_mps=1.0,
        dz_m=1.0,
    )
    np.testing.assert_allclose(shifted_v, [[0.0, 1.0, 2.0, 3.0]])
    assert np.sum(shifted_v) - np.sum(velocity_distribution) == pytest.approx(
        v_inventory_change
    )
    assert np.min(shifted_z) >= 0.0
    assert np.min(shifted_v) >= 0.0


def test_discrete_reservoir_density_and_photoelectron_flux_are_exact() -> None:
    config = _balanced_config()
    state = _make_reservoir_state(
        "electron",
        config.electron,
        nz=config.nz,
        dz_m=config.z_length_m / config.nz,
        initial_condition="reservoir",
    )
    assert np.sum(state.right_inflow) * state.dv_mps == pytest.approx(
        config.electron.number_density_m3, rel=2.0e-15
    )

    query = replace(
        _zero_query(),
        photoelectron_outward_number_flux_m2_s=2.5e12,
        photoelectron_outward_mean_normal_energy_ev=2.2,
    )
    physical_grid = replace(
        config,
        photoelectron_grid=VelocityGridConfig(
            nv=32, vmin_mps=-5.0e6, vmax_mps=5.0e6
        ),
    )
    photoelectron = _make_photoelectron_state(
        physical_grid, query, dz_m=physical_grid.z_length_m / physical_grid.nz
    )
    positive = photoelectron.velocity_mps > 0.0
    represented_flux = np.sum(
        photoelectron.velocity_mps[positive]
        * photoelectron.left_inflow[positive]
    ) * photoelectron.dv_mps
    assert represented_flux == pytest.approx(
        query.photoelectron_outward_number_flux_m2_s, rel=2.0e-15
    )


def test_stationarity_classifier_distinguishes_steady_and_secular() -> None:
    config = _balanced_config()
    times = np.arange(0.0, 1.0001, 0.05)
    history = {"time_s": times}
    for name in (
        "matching_potential_v",
        "electron_inward_number_flux_m2_s",
        "ion_inward_number_flux_m2_s",
        "photoelectron_return_number_flux_m2_s",
        "photoelectron_escape_number_flux_m2_s",
        "total_outer_charge_c_m2",
    ):
        history[name] = np.ones_like(times)
    classification, *_ = _certify_time_history(config, _zero_query(), history)
    assert classification == "steady"

    history["matching_potential_v"] = 100.0 * times
    classification, *_ = _certify_time_history(config, _zero_query(), history)
    assert classification == "secular"


def test_far_boundary_gate_does_not_hide_secular_time_history() -> None:
    certification = _balanced_config().certification
    classification, reason = _finalize_classification(
        "secular",
        numerical_failure=None,
        budget_residual=0.0,
        velocity_loss=0.0,
        gauss_residual=0.0,
        far_field=10.0 * certification.far_field_abs_v_m,
        far_charge=10.0 * certification.far_charge_rel,
        certification=certification,
    )
    assert classification == "secular"
    assert reason == "time history has secular drift"


def test_balanced_zero_field_query_is_certified_and_conservative() -> None:
    result = run_kinetic_query(_balanced_config(), _zero_query())
    assert result.classification == "steady"
    assert result.response[0] == pytest.approx(0.0, abs=1.0e-13)
    assert result.response[1] == pytest.approx(result.response[2], rel=1.0e-14)
    assert result.far_field_v_m == pytest.approx(0.0, abs=1.0e-13)
    assert result.far_charge_imbalance == pytest.approx(0.0, abs=1.0e-13)
    assert result.charge_budget_residual < 1.0e-12
    assert result.max_velocity_boundary_loss_fraction == 0.0


def _rewrite_manifest_hash(manifest_path: Path, raw_path: Path, rows: int) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["raw_csv_sha256"] = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    manifest["query_count"] = rows
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def test_atlas_manifest_and_table_conversion_fail_closed(tmp_path: Path) -> None:
    atlas = tmp_path / "atlas"
    results = write_kinetic_atlas(_balanced_config(), [_zero_query()], atlas)
    assert results[0].classification == "steady"
    raw_path = atlas / "kinetic_response_raw.csv"
    manifest_path = atlas / "kinetic_response_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for key in (
        "git_commit",
        "solver_source_sha256",
        "config_fingerprint_sha256",
        "nz",
        "nv_e",
        "nv_i",
        "z_max_m",
        "velocity_ranges_m_s",
        "stationarity_tolerances",
        "far_boundary_tolerances",
        "raw_csv_sha256",
    ):
        assert key in manifest

    table_path = tmp_path / "response.csv"
    assert convert_kinetic_table(raw_path, manifest_path, table_path) == 1
    assert table_path.read_text(encoding="utf-8").startswith(
        "# matching_plane_z_m=2"
    )

    raw_path.write_text(raw_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="does not match its manifest"):
        convert_kinetic_table(raw_path, manifest_path, tmp_path / "tampered.csv")


def test_table_conversion_rejects_unresolved_and_incomplete_grid(
    tmp_path: Path,
) -> None:
    atlas = tmp_path / "atlas"
    write_kinetic_atlas(_balanced_config(), [_zero_query()], atlas)
    raw_path = atlas / "kinetic_response_raw.csv"
    manifest_path = atlas / "kinetic_response_manifest.json"
    with raw_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    rows[0]["classification"] = "unresolved_transient"
    with raw_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RAW_HEADER)
        writer.writeheader()
        writer.writerows(rows)
    _rewrite_manifest_hash(manifest_path, raw_path, 1)
    with pytest.raises(ValueError, match="uncertified point"):
        convert_kinetic_table(raw_path, manifest_path, tmp_path / "unresolved.csv")

    base = rows[0]
    base["classification"] = "steady"
    incomplete = []
    for displacement, pe_flux in ((0.0, 0.0), (0.0, 1.0), (1.0, 0.0)):
        row = dict(base)
        row["displacement_c_m2"] = str(displacement)
        row["photoelectron_outward_number_flux_m2_s"] = str(pe_flux)
        incomplete.append(row)
    with raw_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RAW_HEADER)
        writer.writeheader()
        writer.writerows(incomplete)
    _rewrite_manifest_hash(manifest_path, raw_path, 3)
    with pytest.raises(ValueError, match="not a complete Cartesian product"):
        convert_kinetic_table(raw_path, manifest_path, tmp_path / "incomplete.csv")
