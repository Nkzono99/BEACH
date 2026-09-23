"""Keep the independent failure diagnostic honest on captured physical endpoints."""

import json
from pathlib import Path
import subprocess
import sys

import pytest


pytest.importorskip("scipy")

REPO_ROOT = Path(__file__).resolve().parents[2]

# The captured runs' physical inputs; geometry is irrelevant to the outer model.
CONFIG = """
[surface_current_model]
response_backend = "zhao_online"
photoelectron_closure = "energy_spectrum"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"

[[particles.species]]
species_key = "electron"
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
temperature_ev = 10.0
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
species_key = "ion"
q_particle = 1.602176634e-19
m_particle = 1.67262192369e-27
number_density_cm3 = 5.0
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
species_key = "photoelectron"
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
"""


@pytest.mark.parametrize(
    ("spectrum_name", "displacement_before", "expected_phi", "native_tolerance"),
    [
        pytest.param(
            "matching_plane_pe_captured_spectrum.csv",
            0.0,
            [5.829194512642798, 5.945909126703211, 6.055444426495559],
            None,
            id="three-physical-endpoints",
        ),
        pytest.param(
            "matching_plane_pe_batch14_spectrum.csv",
            1.9174651506096518e-11,
            [6.432956538216116],
            None,
            id="batch14-bisection-tolerance-miss",
        ),
        pytest.param(
            "matching_plane_pe_batch154_spectrum.csv",
            1.9063843552418728e-11,
            [6.433217419077733539385019749],
            2.629e-19,
            id="physical-profile-with-unresolved-BE-residual",
        ),
        pytest.param(None, 0.0, [], None, id="zero-source-absence"),
    ],
)
def test_failure_diagnostic_classifies_physical_endpoints(
    tmp_path: Path,
    spectrum_name: str | None,
    displacement_before: float,
    expected_phi: list[float],
    native_tolerance: float | None,
) -> None:
    config = tmp_path / "beach.toml"
    config.write_text(CONFIG)
    if spectrum_name is None:
        spectrum = tmp_path / "zero-source.csv"
        spectrum.write_text("energy_low_ev,energy_high_ev,flux_m2_s\n0.0,1.0,0.0\n")
    else:
        spectrum = REPO_ROOT / "tests" / "fixtures" / spectrum_name
    if native_tolerance is not None:
        # Preserve the actual failure receipt's requested tolerance. A physical
        # profile alone must not certify the returned double-precision BE state.
        spectrum_text = spectrum.read_text()
        spectrum = tmp_path / "failed-run.err"
        spectrum.write_text(
            f"search_unresolved; F_tol= {native_tolerance:.3E}C/m2\n"
            "matching-plane failed endpoint: D_before, D_seed [C/m2], duration [s]="
            f" {displacement_before:.17E} {displacement_before:.17E} 2.0\n"
            "matching-plane failed feedback: PE flux, PE energy, electron flux, ion flux="
            " 1.7485673005233697E+13 2.4900318074921848E+00 0.0 0.0\n"
            "matching-plane failed spectrum begin\n"
            + spectrum_text
            + "matching-plane failed spectrum end\n"
        )
    output = tmp_path / "diagnosis"
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "tools" / "diagnose_matching_plane_failure.py"),
            str(spectrum),
            "--config", str(config),
            "--duration", "2",
            "--displacement-before", str(displacement_before),
            "--output", str(output),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    diagnosis = json.loads(output.with_suffix(".json").read_text())
    physical_roots = [
        candidate for candidate in diagnosis["candidates"]
        if candidate["kind"] == "backward_euler" and candidate["physical"]
    ]
    assert diagnosis["physical_BE_candidates"] == len(expected_phi)
    assert sorted(root["phi_H"] for root in physical_roots) == pytest.approx(
        expected_phi, rel=0.0, abs=1e-6,
    )
    assert diagnosis["absence_certificate"]["complete"] is (not expected_phi)
    confirmation = diagnosis["BE_residual_confirmation"]
    assert confirmation["requested_tolerance_C_m2"] == native_tolerance
    if native_tolerance is not None:
        assert diagnosis["diagnosis"] == "physical_profile_found_BE_residual_unresolved"
        assert confirmation["physical_candidates_within_tolerance"] == 0
        assert all(abs(root["BE_residual_C_m2"]) > native_tolerance for root in physical_roots)
        assert all(root["independent_BE_residual_within_native_tolerance"] is False for root in physical_roots)
    elif expected_phi:
        assert diagnosis["diagnosis"] == "physical_profile_found_BE_tolerance_unverified"
        assert confirmation["physical_candidates_within_tolerance"] is None
        assert all(root["independent_BE_residual_within_native_tolerance"] is None for root in physical_roots)
    else:
        assert diagnosis["diagnosis"] == "type_B_absence_certified_in_model_with_roundoff_margin"
        assert diagnosis["absence_certificate"]["flat_state_excluded"] is True
